"""Known signed waves, geometry degeneracy, angle limits and pipeline handoff."""
import json, os
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
from ndstudio.direction import settings,estimate,pair_gate,independent_indices,joint_estimate
from ndstudio.phase_bands import assign_bands

XY=np.array([[0,0],[67,12],[21,94],[133,48],[-72,173],[219,159],[87,-123]],float)


def wave(xy,azimuth=30.,velocity=500.,f=1.5,windows=10):
    s=np.array([np.sin(np.radians(azimuth)),np.cos(np.radians(azimuth))])/velocity
    return np.exp(-2j*np.pi*f*(xy@s))[:,None]*np.exp(1j*np.arange(windows))[None,:]


@pytest.mark.parametrize('angle',[0.,30.,90.,210.,355.])
def test_known_propagation_sign_rotation_and_velocity(angle):
    c=settings(dict(velocity_min_m_s=200.,velocity_max_m_s=1000.,velocity_points=41))
    r,_=estimate(wave(XY,angle)[None],XY,[1.5],c)
    assert r['quality'].all()
    np.testing.assert_allclose((r['azimuth_deg']-angle+180)%360-180,0,atol=1e-5)
    np.testing.assert_allclose(r['velocity_m_s'],500,rtol=1e-7)


def test_angle_45_limit_uncertainty_and_reverse_projection():
    c=settings();r=dict(azimuth_deg=np.array([0.,0.,0.,180.]),direction_uncertainty_deg=np.array([0.,0.,5.,0.]),
        slow_e_s_m=np.zeros(4),slow_n_s_m=np.array([1.,1.,1.,-1.])/500,velocity_m_s=np.full(4,500.),
        quality=np.ones(4,bool),reason=np.full(4,'accepted',dtype='U24'))
    g=pair_gate(r,100/np.sqrt(2),100/np.sqrt(2),c)
    assert g['direction_valid'].tolist()==[True,True,False,True]
    assert g['direction_reason'][2]=='angle_boundary'
    assert g['projected_distance_m'][0]>0 and g['projected_distance_m'][3]<0
    assert not pair_gate(r,100,1,c)['direction_valid'].any()


def test_line_array_and_missing_subarray_do_not_invent_direction():
    xy=np.array([[0,0],[100,0],[200,0],[300,0]],float);z=wave(xy)
    z[:2,1]=np.nan
    r,_=estimate(z[None],xy,[1.5],settings())
    assert not r['quality'].any();assert r['reason'][0]=='geometry';assert r['reason'][1]=='few_stations'


def test_neighbor_frequency_power_average_preserves_opposed_bin_phases():
    f=np.array([1.49,1.5,1.51]);c=settings(dict(velocity_min_m_s=200.,velocity_max_m_s=1000.,velocity_points=41))
    z=np.array([wave(XY,f=k)*np.exp(1j*phase) for k,phase in zip(f,[0,np.pi,np.pi/2])])
    r,_=estimate(z,XY,f,c)
    assert r['quality'].all();np.testing.assert_allclose(r['velocity_m_s'],500,rtol=1e-6)


def test_joint_phase_cycles_use_signed_projection_and_independent_windows():
    c=settings(dict(velocity_min_m_s=200.,velocity_max_m_s=1000.,velocity_points=41))
    z=wave(XY,windows=100);r,_=estimate(z[None],XY,[1.5],c)
    rows,pairs=joint_estimate(z,XY,np.arange(100)*10,100,1.5,r,c,False)
    assert rows and rows[0]['accepted'] and rows[0]['independent_windows']==10
    assert rows[0]['phase_velocity_m_s']==pytest.approx(500,rel=1e-7)
    assert pairs and max(abs(p['final_phase_deg']) for p in pairs)>0
    assert len(independent_indices(np.arange(100)*10,100,np.ones(100,bool)))==10


def test_continuous_bands_split_rejected_frequency_and_direction_change():
    f=np.arange(1.,2.01,.1);q=pd.DataFrame(dict(curve_id='A',frequency_hz=f,accepted=True,selected_direction_group=1))
    q.loc[4,'accepted']=False;q.loc[8:,'selected_direction_group']=2
    result=assign_bands(q,f,dict(step_hz=.1,min_points=3,min_band_hz=.2))
    assert result.loc[3,'band_id']!=result.loc[5,'band_id']!=result.loc[8,'band_id']
    assert not result.loc[4,'accepted'];assert result.accepted.sum()==10


def fixture(root):
    import h5py
    from ndstudio.project import new_project
    from dispersion.common import write_json,sha
    root=new_project(root);out=root/'output';out.mkdir(exist_ok=True)
    fs=20.;duration=600.;t=np.arange(round(duration*fs))/fs;f=np.arange(1.3,1.81,.1)
    names=['S'+str(i) for i in range(len(XY))];az=30.;v=500.;slow=np.array([.5,np.sqrt(3)/2])/v
    # Broadband delayed stochastic wave; real-domain generation has independent spectral phases.
    rng=np.random.default_rng(510);freq=np.fft.rfftfreq(len(t),1/fs);spectrum=np.fft.rfft(rng.normal(size=len(t)))
    spectrum[(freq<.8)|(freq>2.5)]=0
    data=np.array([np.fft.irfft(spectrum*np.exp(-2j*np.pi*freq*(xy@slow)),n=len(t)) for xy in XY])
    stations=pd.DataFrame(dict(station=names,file=[s+'.npz' for s in names],seed_id=['XX.'+s+'..HHZ' for s in names],x_m=XY[:,0],y_m=XY[:,1],precorrected=True,units='m/s',polarity=1))
    stations.to_csv(root/'stations.csv',index=False)
    c=json.loads((root/'config.json').read_text());c['scan'].update(min_hz=1.,max_hz=2.,step_hz=.1,min_band_hz=.25,min_points=3);write_json(root/'config.json',c)
    curves=[];measured=[];pairs=[]
    for k,(a,b) in enumerate([(0,2),(0,1),(1,3)]):
        d=XY[b]-XY[a];distance=np.linalg.norm(d);pair=names[a]+'--'+names[b];pairs.append(dict(pair_index=k,pair=pair,distance_m=distance))
        for frequency in f:
            phase=-2*np.pi*frequency*(d@slow);cycle=1;source=np.pi/4
            conditional=2*np.pi*frequency*distance/(source-phase+2*np.pi*cycle)
            curves.append(dict(curve_id=pair,segment_id=pair,pair_index=k,pair=pair,branch='positive',frequency_hz=frequency,
                phase_velocity_m_s=conditional,group_velocity_m_s=500.,source_phase_rad=source,constant_cycle=cycle,
                cohort='synthetic',instrument_response_corrected=True,background_discriminated=True,source_phase_verified=False,is_model_prediction=False))
            measured.append(dict(segment_id=pair,frequency_hz=frequency,measured_unwrapped_phase_rad=phase,relative_phase_q025_rad=-.4,relative_phase_q975_rad=.4))
    pd.DataFrame(curves).to_csv(out/'conditional_dispersion.csv',index=False);pd.DataFrame(measured).to_csv(out/'measured_phase_segments.csv',index=False);pd.DataFrame(pairs).to_csv(out/'pairs.csv',index=False)
    path=root/'aligned.h5'
    with h5py.File(path,'w') as h:
        h.attrs.update(fs_hz=fs,start_utc='2020-01-01T00:00:00.000000Z',end_utc_exclusive='2020-01-01T00:10:00.000000Z',units='m/s',response_corrected=True,preprocessing_config_sha256=sha(root/'config.json'))
        for name,x in zip(names,data):
            g=h.create_group('stations/'+name);g.create_dataset('data',data=x);g.create_dataset('valid',data=np.ones(len(x),np.uint8))
    return root,path


def test_real_worker_chain_shared_windows_angle_masks_reports_and_band_handoff(tmp_path,monkeypatch):
    from ndstudio.phase_stack import run_phase_stack
    from ndstudio.phase_audit import audit
    root,h=fixture(tmp_path/'p');direction=dict(velocity_min_m_s=200.,velocity_max_m_s=1000.,velocity_points=41)
    r=run_phase_stack(root,aligned_h5=h,stabilization_s=0,repeat_cycles=40,slide_cycles=2,taper='hann',frequency_method='exact',direction=direction,amplitude_rejection=False,tolerance_deg=20)
    assert r['shared_utc_windows'] and r['accepted_window_evaluations']>0 and r['inversion_handoff_points']>0
    out=root/'output/phase_guided';assert audit(out)['pass_audit']
    trace=json.loads((out/'phase_trace.json').read_text());assert trace['complete']
    assert (out/'direction/report.html').exists() and (out/'inversion_curves.csv').exists()
    q=pd.read_csv(out/'phase_gate_summary.csv');assert (q.direction_eligible_windows<=q.valid_windows).all()
    assert q[q.pair=='S0--S1'].accepted_windows.sum()==0
    assert 'href="velocity_comparison/report.html"' in (out/'phase_report.html').read_text()
    starts={}
    for row in q.itertuples():
        with np.load(out/row.windows_file) as z:
            if row.frequency_hz in starts:np.testing.assert_array_equal(starts[row.frequency_hz],z['start_sample'])
            starts[row.frequency_hz]=z['start_sample']
            if z['accepted'].any():assert np.max(z['pair_angle_upper_deg'][z['accepted']])<=45+1e-8
    joint=pd.read_csv(out/'direction/joint_dispersion.csv')
    assert joint.accepted.any();assert abs(joint.phase_velocity_m_s.median()-500)<30
    from ndstudio.candidates import read_curves,review,auto_selection
    (root/'output/families_complete.json').write_text('{}')
    inherited,catalog=read_curves(root)
    assert len(inherited)==r['inversion_handoff_points'] and inherited.parent_curve_id.nunique()>0
    assert inherited.constant_cycle.eq(1).all() and inherited.source_phase_rad.eq(np.pi/4).all()
    assert len(review(root))==len(catalog)
    # Changing only the angle must reuse spectra/FK, then change downstream masks.
    import ndstudio.direction as module
    original=module.estimate
    def bounded_only(z,*args):
        assert z.shape[2]<=1,'A full FK scan was unnecessarily recomputed'
        return original(z,*args)
    monkeypatch.setattr(module,'estimate',bounded_only)
    narrowed=run_phase_stack(root,aligned_h5=h,stabilization_s=0,repeat_cycles=40,slide_cycles=2,taper='hann',frequency_method='exact',
        direction=direction|{'max_pair_angle_deg':10.},amplitude_rejection=False,tolerance_deg=20)
    assert narrowed['accepted_window_evaluations']<r['accepted_window_evaluations']
    trace=json.loads((out/'phase_trace.json').read_text())
    stages=[x for x in trace['stages'] if x['stage'].startswith('方向估计')]
    assert stages and all(x['results']['cache_reused'] for x in stages)
    inherited,catalog=read_curves(root)
    assert inherited.empty and catalog.empty
    assert not auto_selection(root,minimum_wavelengths=0)['cases']


def test_direction_json_edits_and_ui_edits_merge_without_old_values(tmp_path):
    from PySide6.QtWidgets import QApplication
    from ndstudio.gui import MainWindow
    from ndstudio.project import new_project
    app=QApplication.instance() or QApplication([]);w=MainWindow(new_project(tmp_path/'p'))
    c=json.loads(w.directioneditor.toPlainText());c['max_pair_angle_deg']=35.;c['min_coherence']=.8;w.directioneditor.setPlainText(json.dumps(c))
    w.directionfields['direction_step_deg'].setValue(10)
    used=w.phase_settings()['direction']
    assert used['max_pair_angle_deg']==35 and used['min_coherence']==.8 and used['direction_step_deg']==10
    assert w.directionfields['max_pair_angle_deg'].value()==35
    w.close()
