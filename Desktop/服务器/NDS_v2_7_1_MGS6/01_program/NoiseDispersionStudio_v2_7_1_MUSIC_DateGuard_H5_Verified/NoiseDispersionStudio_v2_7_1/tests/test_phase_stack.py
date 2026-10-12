import numpy as np
import pandas as pd
import pytest
from ndstudio.phase_stack import goertzel_windows,phase_interval_mask,wrap,stack_stats,export_targets,robust_flags

@pytest.mark.parametrize('f,n,taper',[ (2.37,137,'none'),(4.13,101,'hann'),(1.1,97,'none')])
def test_exact_goertzel_matches_dft(f,n,taper):
    rng=np.random.default_rng(73);x=rng.normal(size=700);starts=np.array([0,7,83,251],dtype=np.int64)
    w=np.ones(n) if taper=='none' else np.hanning(n)
    actual=goertzel_windows(x,starts,n,f,20.,w)
    expected=np.array([np.sum((x[a:a+n]-x[a:a+n].mean())*w*np.exp(-2j*np.pi*f*np.arange(n)/20)) for a in starts])
    np.testing.assert_allclose(actual,expected,rtol=1e-11,atol=1e-10)

def test_circle_wrap_and_large_arc():
    phases=np.radians([179,-179,150,0])
    assert phase_interval_mask(phases,np.pi,np.radians(-5),np.radians(5)).tolist()==[True,True,False,False]
    assert phase_interval_mask(phases,0,-np.pi,np.pi).all()
    assert phase_interval_mask(np.array([np.nan]),0,-np.pi,np.pi).tolist()==[False]
    with pytest.raises(ValueError):phase_interval_mask(phases,0,1,-1)

def test_delayed_wave_cross_sign():
    fs=20.;f=2.5;t=np.arange(4000)/fs;delay=.073
    a=np.cos(2*np.pi*f*t);b=np.cos(2*np.pi*f*(t-delay));starts=np.array([0,83,251])
    ga=goertzel_windows(a,starts,160,f,fs,np.ones(160));gb=goertzel_windows(b,starts,160,f,fs,np.ones(160))
    np.testing.assert_allclose(wrap(np.angle(gb*np.conj(ga))+2*np.pi*f*delay),0,atol=1e-12)

def test_missing_and_zero():
    x=np.ones(120);x[25]=np.nan;z=goertzel_windows(x,np.array([0,30]),20,2.3,20.,np.ones(20))
    assert np.isfinite(z[0]);assert np.isfinite(z[1]);assert abs(z[1])==0
    z=goertzel_windows(x,np.array([10]),20,2.3,20.,np.ones(20));assert np.isnan(z[0])
    assert stack_stats(np.array([0j]),np.array([False]))['unit_phase_deg'] is None

def test_changed_export_invalidates_old_completion(tmp_path):
    o=tmp_path/'output';o.mkdir();f=2.;d=100.;theta=-.2;phi=0.;cycle=1
    pd.DataFrame([dict(curve_id='N1',segment_id='s',pair_index=0,pair='a--b',branch='positive',frequency_hz=f,phase_velocity_m_s=2*np.pi*f*d/(phi-theta+2*np.pi*cycle),source_phase_rad=phi,constant_cycle=cycle)]).to_csv(o/'conditional_dispersion.csv',index=False)
    pd.DataFrame([dict(pair_index=0,distance_m=d)]).to_csv(o/'pairs.csv',index=False)
    pd.DataFrame([dict(segment_id='s',frequency_hz=f,measured_unwrapped_phase_rad=theta,relative_phase_q025_rad=-.1,relative_phase_q975_rad=.1)]).to_csv(o/'measured_phase_segments.csv',index=False)
    _,dest=export_targets(tmp_path,tolerance_deg=15)
    (dest/'phase_stack_complete.json').write_text('{}');(dest/'phase_report.html').write_text('old')
    export_targets(tmp_path,tolerance_deg=5)
    assert not (dest/'phase_stack_complete.json').exists();assert not (dest/'phase_report.html').exists()

def test_no_amplitude_censor():
    z=np.exp(1j*np.radians([1,2,3,4,5]));z[-1]*=1e9
    base=abs(z)>0;accepted=base&phase_interval_mask(np.angle(z),0,-.2,.2)
    large,_=robust_flags(abs(z));assert large[-1] and accepted[-1]
    result=stack_stats(z,accepted);assert result['count']==5 and result['unit_resultant_length']>.99

@pytest.mark.parametrize('branch,sign',[('positive',1),('negative',-1)])
def test_export_source_cycle_and_branch(tmp_path,branch,sign):
    o=tmp_path/'output';o.mkdir();f=2.5;d=1000.;theta=-1.3;phi=np.pi/4
    curves=[]
    for cycle in [1,2]:
        v=2*np.pi*f*d/(phi-theta+2*np.pi*cycle)
        curves.append(dict(curve_id=f'N{cycle}',segment_id='s',pair_index=0,pair='a--b',branch=branch,frequency_hz=f,phase_velocity_m_s=v,source_phase_rad=phi,constant_cycle=cycle))
    pd.DataFrame(curves).to_csv(o/'conditional_dispersion.csv',index=False)
    pd.DataFrame([dict(pair_index=0,distance_m=d)]).to_csv(o/'pairs.csv',index=False)
    pd.DataFrame([dict(segment_id='s',frequency_hz=f,measured_unwrapped_phase_rad=theta,relative_phase_q025_rad=-.1,relative_phase_q975_rad=.2)]).to_csv(o/'measured_phase_segments.csv',index=False)
    q,_=export_targets(tmp_path,tolerance_deg=0,range_source='measured')
    np.testing.assert_allclose(q.predicted_phase_deg,np.degrees(sign*theta),atol=1e-10)
    np.testing.assert_allclose(q.phase_offset_low_deg,np.degrees(-.1 if sign>0 else -.2))
    np.testing.assert_allclose(q.phase_offset_high_deg,np.degrees(.2 if sign>0 else .1))
    np.testing.assert_allclose(q.predicted_phase_deg.iloc[0],q.predicted_phase_deg.iloc[1],atol=1e-10)
    with pytest.raises(ValueError):export_targets(tmp_path,tolerance_deg=180)

@pytest.mark.parametrize('native_fs,start_fraction',[(20.,0.),(100.,.013)])
def test_alignment_preserves_utc_gap_and_native_rate(tmp_path,monkeypatch,native_fs,start_fraction):
    import json,h5py
    from obspy import Stream,Trace,UTCDateTime
    import ndstudio.phase_stack as mod
    from dispersion.common import DEFAULTS,merge
    c=merge(DEFAULTS,{'response':{'output':'VEL'}})
    (tmp_path/'config.json').write_text(json.dumps(c))
    monkeypatch.setattr(mod,'load_config',lambda path:c)
    row=dict(station='A',seed_id='XX.A..HHZ',precorrected=True,units='m/s',polarity=1)
    monkeypatch.setattr(mod,'station_table',lambda c:(pd.DataFrame([row]),None,None))
    origin=UTCDateTime('2020-01-01T00:00:00');start=origin-100+start_fraction
    t=np.arange(round(800*native_fs))/native_fs
    x=np.cos(2*np.pi*2*(t+float(start-origin)));x[(t+float(start-origin)>=200)&(t+float(start-origin)<210)]=np.nan
    tr=Trace(x,dict(sampling_rate=native_fs,starttime=start));stream=Stream([tr])
    monkeypatch.setattr(mod,'read_station',lambda row,c:(stream.copy(),[]))
    monkeypatch.setattr(mod,'response_table',lambda row,st,c:(None,{'response_status':'physical_input_user_declared'}))
    arr,fs,utc,meta=mod.align_project(tmp_path,{'A'},str(origin),str(origin+600),tmp_path/'a.h5')
    assert fs==20. and utc==origin and len(arr['A'])==12000
    assert np.isnan(arr['A'][4000:4200]).all()
    target=np.cos(2*np.pi*2*np.arange(12000)/20.)
    assert np.corrcoef(arr['A'][1000:2000],target[1000:2000])[0,1]>.999
    with h5py.File(tmp_path/'a.h5') as h:
        assert h.attrs['response_corrected'];assert h.attrs['units']=='m/s'
        assert h['stations/A/valid'][4000:4200].sum()==0
