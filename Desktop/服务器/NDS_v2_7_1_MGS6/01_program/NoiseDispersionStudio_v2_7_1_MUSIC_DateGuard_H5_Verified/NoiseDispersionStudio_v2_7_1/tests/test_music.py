"""Physical checks for MUSIC, changing azimuth, two sources and both handoffs."""
import json
import numpy as np
import pandas as pd
import pytest
from ndstudio.direction import settings,estimate,pair_gate
from ndstudio.music_front import frequency_reference,correct_candidates

XY=np.array([[0,0],[600,50],[150,800],[900,500],[-300,1100],[1200,-450]],float)


def snapshots(f,azimuth,velocity,windows=5,count=48,noise=.005):
    rng=np.random.default_rng(483);axis=np.array([np.sin(np.radians(azimuth)),np.cos(np.radians(azimuth))])
    a=np.exp(-2j*np.pi*f*(XY@axis)/velocity)
    source=rng.normal(size=(windows,count))+1j*rng.normal(size=(windows,count))
    return a[:,None,None]*source[None]+noise*(rng.normal(size=(len(XY),windows,count))+1j*rng.normal(size=(len(XY),windows,count)))


def test_covariance_music_recovers_frequency_varying_direction_above_2000():
    c=settings(dict(velocity_min_m_s=1500.,velocity_max_m_s=4500.,velocity_points=101))
    for f,angle,speed in [(1.,10.,3200.),(1.2,35.,3000.),(1.5,70.,2800.)]:
        r,_=estimate(snapshots(f,angle,speed)[None],XY,[f],c)
        assert r['quality'].all() and (r['music_source_count']==1).all()
        np.testing.assert_allclose((r['azimuth_deg']-angle+180)%360-180,0,atol=.5)
        np.testing.assert_allclose(r['velocity_m_s'],speed,rtol=.01)
        assert r['music_snapshot_count'][0]==48


def test_two_wavefronts_have_two_signal_eigenvectors_and_are_not_combined():
    rng=np.random.default_rng(78);f=1.;count=100
    axes=np.array([[.5,np.sqrt(3)/2],[-.6,.8]])
    A=np.exp(-2j*np.pi*f*(XY@axes.T)/2000.)
    source=rng.normal(size=(2,3,count))+1j*rng.normal(size=(2,3,count))
    z=np.einsum('sq,qwk->swk',A,source)+.01*(rng.normal(size=(6,3,count))+1j*rng.normal(size=(6,3,count)))
    r,_=estimate(z[None],XY,[f],settings(dict(velocity_min_m_s=1200.,velocity_max_m_s=3500.,velocity_points=101)))
    assert (r['music_source_count']==2).all()
    assert not r['quality'].any() and (r['reason']=='multiple_sources').all()


def test_pure_noise_and_line_array_never_make_a_true_velocity():
    rng=np.random.default_rng(91);z=rng.normal(size=(1,6,4,80))+1j*rng.normal(size=(1,6,4,80))
    r,_=estimate(z,XY,[1.],settings());assert not r['quality'].any()
    line=XY.copy();line[:,1]=0
    r,_=estimate(z,line,[1.],settings());assert (r['reason']=='geometry').all()
    with pytest.raises(ValueError,match='仅支持MUSIC'):settings({'method':'fk'})


def test_front_has_no_45_degree_gate_and_uses_each_frequency_projection():
    c=settings();f=np.arange(1.,1.51,.1);delta=np.array([0.,1000.]);refs={};rows=[]
    for i,k in enumerate(f):
        angle=20.+8*i;speed=700.-20*i;axis=np.array([np.sin(np.radians(angle)),np.cos(np.radians(angle))])
        r=dict(azimuth_deg=np.full(10,angle),velocity_m_s=np.full(10,speed),slow_e_s_m=np.full(10,axis[0]/speed),slow_n_s_m=np.full(10,axis[1]/speed),
            direction_uncertainty_deg=np.full(10,2.5),quality=np.ones(10,bool),reason=np.full(10,'accepted',dtype='U32'),direction_group=np.full(10,round(angle/15),int))
        ref=frequency_reference(r,np.arange(10)*100,100,delta,c,'positive',False)
        assert ref is not None
        if angle>45:assert frequency_reference(r,np.arange(10)*100,100,delta,c,'positive',True) is None
        refs[('A--B','positive',round(k,9))]=ref
        rows.append(dict(curve_id='C',pair='A--B',branch='positive',pair_index=0,segment_id='S',cohort='all',frequency_hz=k,
            phase_velocity_m_s=speed/np.cos(np.radians(angle)),group_velocity_m_s=800.,source_phase_rad=0.,constant_cycle=0,background_discriminated=True))
    q,a,cat=correct_candidates(pd.DataFrame(rows),pd.DataFrame({'distance_m':[1000.]},index=[0]),refs,dict(step_hz=.1,min_points=3,min_band_hz=.2))
    assert len(q)==6 and len(cat)==1 and q.projection_applied.all()
    np.testing.assert_allclose(q.phase_velocity_m_s,700.-20*np.arange(6))
    assert q.projected_distance_m.nunique()==6
    assert not np.allclose(q.group_velocity_m_s,pd.DataFrame(rows).group_velocity_m_s*np.cos(np.radians(q.propagation_azimuth_deg)))


def test_front_native_stage_and_goertzel_export_do_not_double_project(tmp_path):
    from test_direction import fixture
    from ndstudio.music_front import run_music_front,current
    from ndstudio.phase_stack import export_targets
    from dispersion.common import write_json
    root,h=fixture(tmp_path/'p')
    write_json(root/'phase_settings.json',dict(direction=dict(velocity_min_m_s=200.,velocity_max_m_s=1000.,velocity_points=41)))
    r=run_music_front(root,aligned_h5=h)
    assert r['corrected_points']>0 and current(root)
    q=pd.read_csv(root/'output/music_front/dispersion.csv');t,_=export_targets(root)
    np.testing.assert_allclose(t.phase_velocity_m_s,q.phase_velocity_m_s)
    np.testing.assert_allclose(t.prediction_reconstruction_error_rad,0,atol=1e-6)
    assert t.projection_applied.all() and t.projected_distance_m.notna().all()


def test_large_station_distance_adapts_above_40_without_capping_velocity(tmp_path):
    from ndstudio.project import new_project,adapt_geometry_settings
    from dispersion.common import write_json
    root=new_project(tmp_path/'p')
    pd.DataFrame([dict(station='S'+str(i),file='unused.npz',seed_id='XX.S'+str(i)+'..HHZ',x_m=x,y_m=y,precorrected=True,units='m/s',polarity=1) for i,(x,y) in enumerate([(0,0),(110000,0),(50000,20000)])]).to_csv(root/'stations.csv',index=False)
    adapt_geometry_settings(root);c=json.loads((root/'config.json').read_text())
    assert c['scan']['arrival_max_s']>=61
    assert c['phase']['velocity_bounds_m_s'][1]>2000


def test_empty_front_can_only_supply_screening_priors_and_cannot_be_inverted(tmp_path,monkeypatch):
    from test_direction import fixture
    from ndstudio import music_front
    from ndstudio.phase_stack import export_targets
    from ndstudio.candidates import auto_selection
    from ndstudio.inversion import invert_project,SelectionError
    root,_=fixture(tmp_path/'p');out=root/'output';front=out/'music_front';front.mkdir();(front/'music_complete.json').write_text('{}')
    raw=pd.read_csv(out/'conditional_dispersion.csv')
    distances=pd.read_csv(out/'pairs.csv').set_index('pair_index').distance_m
    raw['distance_over_wavelength']=raw.pair_index.map(distances)*raw.frequency_hz/raw.phase_velocity_m_s
    raw.to_csv(out/'conditional_dispersion.csv',index=False);raw.iloc[:0].to_csv(front/'dispersion.csv',index=False)
    monkeypatch.setattr(music_front,'current',lambda root:True)
    with pytest.raises(ValueError,match='No conditional'):export_targets(root)
    targets,_=export_targets(root,screening_priors_only=True)
    assert len(targets)==len(raw) and targets.target_role.eq('apparent_screening_prior_only').all()
    assert targets.is_model_prediction.all() and targets.velocity_type.eq('conditional_pair_velocity').all()
    cat=raw.groupby('curve_id').agg(segment_id=('segment_id','first'),cohort=('cohort','first'),pair_index=('pair_index','first'),branch=('branch','first'),constant_cycle=('constant_cycle','first'),source_phase_rad=('source_phase_rad','first'),fmin_hz=('frequency_hz','min'),fmax_hz=('frequency_hz','max'),points=('frequency_hz','size')).reset_index()
    cat['background_discriminated']=True;cat['compatible_other_pairs']=0;cat.to_csv(out/'curve_catalog.csv',index=False)
    selection=auto_selection(root,max_cases=1,minimum_wavelengths=0,screening_priors_only=True)
    assert selection['cases'] and selection['selection_notes']['screening_priors_only']
    monkeypatch.setattr('ndstudio.candidates.read_curves',lambda project:(raw,cat))
    with pytest.raises(SelectionError,match='先验不能直接反演'):invert_project(root)
