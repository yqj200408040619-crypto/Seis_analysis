import json
import numpy as np
import pandas as pd
import pytest
from ndstudio.examples import inversion_example
from ndstudio.inversion import validate_cases,invert_tables,SelectionError
from ndstudio.project import new_project
from ndstudio.candidates import import_curve_csv,auto_selection

@pytest.fixture
def known(tmp_path):
    root=inversion_example(tmp_path/'known');return pd.read_csv(root/'curves.csv'),json.loads((root/'selection.json').read_text(encoding='utf-8-sig')),root

def test_reject_duplicate_path_overlap(known):
    frame,s,_=known;s['cases'][0]['bands']*=2
    with pytest.raises(SelectionError,match='重复计权'):validate_cases(frame,s)

@pytest.mark.parametrize('column,value',[('instrument_response_corrected',False),('background_discriminated',False),('is_model_prediction',True)])
def test_quality_or_prediction_cannot_be_bypassed(known,column,value):
    frame,s,_=known;frame[column]=value
    with pytest.raises(SelectionError):validate_cases(frame,s)

def test_modal_assumption_is_required(known):
    frame,s,_=known;del s['cases'][0]['bands'][0]['mode']
    with pytest.raises(SelectionError,match='模态'):validate_cases(frame,s)

def test_auto_groups_preserve_user_search_bounds(known,tmp_path):
    frame,s,r=known;p=new_project(tmp_path/'project');import_curve_csv(p,r/'curves.csv')
    (p/'inversion_settings.json').write_text(json.dumps(s['settings']))
    generated=auto_selection(p,max_cases=1,minimum_wavelengths=0)
    assert generated['settings']['thickness_bounds_m']==[10.,150.]

def test_wavelength_floor_is_explicit_and_does_not_relax_phase_quality(known,tmp_path):
    frame,s,r=known;frame['distance_over_wavelength']=2.5
    source=tmp_path/'curves.csv';frame.to_csv(source,index=False)
    p=new_project(tmp_path/'project');import_curve_csv(p,source)
    assert not auto_selection(p,minimum_wavelengths=3)['cases']
    selected=auto_selection(p,minimum_wavelengths=2)
    assert selected['cases'] and selected['selection_notes']['minimum_distance_over_wavelength']==2
    for bad in (-1,float('nan')):
        with pytest.raises(ValueError):auto_selection(p,minimum_wavelengths=bad)

def test_known_layered_structure_recovers_and_roots_agree(known):
    frame,s,r=known;models=invert_tables(frame,s,r/'output')
    assert models,'Known layered data yielded no compatible model'
    best=min(models,key=lambda m:m['phase_weighted_rms_fraction'])
    assert best['phase_weighted_rms_fraction']<.002
    np.testing.assert_allclose(best['thickness_m'],[35.],rtol=.04)
    np.testing.assert_allclose(best['vs_m_s'],[350.,800.],rtol=.03)
    if best['independent_solver_available']:assert best['independent_solver_pass']
    pred=pd.read_csv(r/'output/model_predictions.csv');assert not pred.is_observation.any()
    obs=pd.read_csv(r/'output/selected_observed_dispersion.csv');assert obs.is_observation.all()
    assert not best['unique_structure_identified']
    layers=pd.read_csv(r/'output/structure_layers.csv')
    q=layers[layers.model_id==best['model_id']]
    assert q.top_depth_m.iloc[0]==0
    assert q.bottom_depth_m.iloc[0]==pytest.approx(sum(best['thickness_m']))
    assert q.layer_type.iloc[-1]=='halfspace' and pd.isna(q.bottom_depth_m.iloc[-1])

def test_missing_independent_solver_stops_when_required(known,monkeypatch):
    import sys
    frame,s,r=known;s['settings']['require_independent_solver']=True;monkeypatch.setitem(sys.modules,'pysurf96',None)
    with pytest.raises(ImportError,match='pysurf96'):invert_tables(frame,s,r/'output')

def test_missing_modal_roots_reject_attempt_without_aborting(known,monkeypatch):
    import ndstudio.inversion as module
    frame,s,r=known;s['settings'].update(starts=1,max_evaluations=2)
    def absent(*args,**kwargs):raise ValueError('Missing modal root')
    monkeypatch.setattr(module,'forward',absent)
    assert invert_tables(frame,s,r/'output')==[]
    assert pd.read_csv(r/'output/inversion_attempts.csv').reason.iloc[0]=='missing_modal_root'
    assert json.loads((r/'output/inversion_complete.json').read_text(encoding='utf-8-sig'))['compatible_models']==0


def test_independent_solver_missing_root_is_recorded_without_aborting(known,monkeypatch):
    wrapper=pytest.importorskip('pysurf96.wrapper')
    import ndstudio.inversion as module
    frame,s,r=known;s['settings'].update(starts=1,require_independent_solver=True)
    def absent(*args,**kwargs):raise wrapper.Surf96Error('No valid independent modal root')
    monkeypatch.setattr(module,'independent_phase',absent)
    assert invert_tables(frame,s,r/'output')==[]
    attempts=pd.read_csv(r/'output/inversion_attempts.csv')
    assert attempts.reason.eq('independent_solver_mismatch_or_missing_root').any()
    assert not attempts.accepted.any()

def test_automatic_joint_groups_require_pairwise_compatibility(known,tmp_path):
    frame,s,r=known;pieces=[]
    for i,cid in enumerate(['A','B','C']):
        part=frame.copy();part['curve_id']=cid;part['pair_index']=i;part['pair']=cid+'--X';pieces.append(part)
    combined=pd.concat(pieces);path=tmp_path/'group_input.csv';combined.to_csv(path,index=False)
    root=new_project(tmp_path/'clique_project');import_curve_csv(root,path)
    pd.DataFrame({'curve_a':['A','B'],'curve_b':['B','C']}).to_csv(root/'output/cross_pair_compatible_curves.csv',index=False)
    result=auto_selection(root,max_cases=3,minimum_wavelengths=0)
    assert result['cases']
    for case in result['cases']:assert len(case['bands'])<=2,'A–B–C chain was incorrectly accepted as a compatible joint group'

def test_selection_change_invalidates_completed_structure(known,tmp_path):
    from ndstudio.inversion import selection_digest,code_digest,inversion_current
    frame,s,r=known;root=new_project(tmp_path/'freshness');(root/'selection.json').write_text(json.dumps(s))
    marker={'selection_sha256':selection_digest(s),'studio_code_sha256':code_digest()}
    (root/'output/inversion_complete.json').parent.mkdir(exist_ok=True);(root/'output/inversion_complete.json').write_text(json.dumps(marker))
    assert inversion_current(root)
    s['cases'][0]['bands'][0]['mode']=1;(root/'selection.json').write_text(json.dumps(s));assert not inversion_current(root)

def test_string_false_quality_flag_is_not_truthy(known):
    frame,s,r=known;frame['instrument_response_corrected']='no'
    with pytest.raises(SelectionError,match='响应'):validate_cases(frame,s)
