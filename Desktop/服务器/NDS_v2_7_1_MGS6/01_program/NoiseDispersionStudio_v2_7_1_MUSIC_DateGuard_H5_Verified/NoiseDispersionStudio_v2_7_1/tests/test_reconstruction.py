"""Focused checks of restored low-frequency projection and inversion handoff."""
import json,os
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
import numpy as np
import pandas as pd
import pytest
from ndstudio.reconstruction import settings,low_targets,align_cross,phase_speed,shifted_null,comparison_tables

@pytest.mark.parametrize('sign',[1,-1])
def test_signed_projection(sign):
    f=.7;v=600.;d=1000.;projection=sign*d*np.cos(np.radians([15.,25.,35.]));reference=sign*d*np.cos(np.radians(25.))
    raw=np.exp(-2j*np.pi*f*projection/v)*[1,3,2];aligned,_=align_cross(raw,f,projection,reference,v);anchor=-2*np.pi*f*reference/v
    phase=np.angle(aligned.sum());phase+=2*np.pi*round((anchor-phase)/(2*np.pi))
    assert phase_speed(f,phase,reference)==pytest.approx(v)
    assert phase_speed(f,phase,sign*d)!=pytest.approx(v)
    np.testing.assert_allclose(np.angle(aligned*np.exp(-1j*anchor)),0,atol=1e-14)

def test_prediction_and_background_null():
    q=pd.DataFrame(dict(curve_id=['C']*3,frequency_hz=[1,1.1,1.2],phase_velocity_m_s=[600,570,550],background_discriminated=True,is_model_prediction=False))
    spec=dict(bandpass_hz=[.1,4],direction_half_width=0,velocity_bounds=[150,3000]);opt=settings(dict(enabled=True,lower_frequency_hz=.3,step_hz=.1))
    r=low_targets(q,opt,spec);low=r[r.frequency_hz<1]
    assert low.frequency_hz.min()==pytest.approx(.3) and low.is_model_prediction.all() and not low.background_discriminated.any()
    with pytest.raises(ValueError):low_targets(q.assign(phase_velocity_m_s=[500,550,600]),opt,spec)
    rng=np.random.default_rng(541);a=np.exp(1j*rng.uniform(-np.pi,np.pi,100));center=-2*np.pi/5
    args=(np.arange(100)*40,40,np.ones(100,bool),np.full(100,100.),100.,1.,500.,center,-.2,.2,opt,False)
    assert shifted_null(a,a*np.exp(1j*center),*args)[0]<=.05
    assert shifted_null(a,np.exp(1j*rng.uniform(-np.pi,np.pi,100)),*args)[0]>.05

def test_common_and_offset_grids():
    b=pd.DataFrame(dict(curve_id=['C']*4,frequency_hz=[1.,1.1,1.2,1.3],phase_velocity_m_s=[600,570,540,510],group_velocity_m_s=500))
    a=pd.DataFrame(dict(curve_id=['C__B1']*7,parent_curve_id=['C']*7,frequency_hz=[.7,.8,.9,1.,1.1,1.2,1.3],phase_velocity_m_s=[690,660,630,600,570,540,510],group_velocity_m_s=500))
    scan=dict(step_hz=.1,min_points=3,min_band_hz=.25);before,after=comparison_tables(b,a,scan)
    assert len(before)==len(after)==4
    pd.testing.assert_frame_equal(before[['curve_id','frequency_hz']],after[['curve_id','frequency_hz']])
    bb=b.copy();bb['curve_id']='D';bb['frequency_hz']+=.05
    aa=a.copy();aa['curve_id']='D__B1';aa['parent_curve_id']='D';aa['frequency_hz']+=.05
    before,after=comparison_tables(pd.concat([b,bb]),pd.concat([a,aa]),scan);assert len(before)==len(after)==8

@pytest.fixture(scope='module')
def native(tmp_path_factory):
    import h5py
    from test_direction import fixture,XY
    from ndstudio.phase_stack import run_phase_stack
    root,hpath=fixture(tmp_path_factory.mktemp('projected_native'));fs=20.;count=28000;f=np.fft.rfftfreq(count,1/fs);axis=np.array([.5,np.sqrt(3)/2]);v=500*(np.maximum(f,.1)/1.5)**(-.25)
    rng=np.random.default_rng(510);z=np.fft.rfft(rng.normal(size=count));z[(f<.6)|(f>2.5)]=0
    data=[np.fft.irfft(z*np.exp(-2j*np.pi*f*(xy@axis)/v),n=count) for xy in XY]
    with h5py.File(hpath,'r+') as h:
        h.attrs['end_utc_exclusive']='2020-01-01T00:23:20.000000Z'
        for i,x in enumerate(data):
            g=h['stations/S'+str(i)];del g['data'];del g['valid'];g.create_dataset('data',data=x);g.create_dataset('valid',data=np.ones(count,np.uint8))
    out=root/'output';q=pd.read_csv(out/'conditional_dispersion.csv');q=q[q.pair=='S0--S2'].copy();q['frequency_hz']=q.frequency_hz.round(9)
    d=np.linalg.norm(XY[2]-XY[0]);projection=(XY[2]-XY[0])@axis;speed=500*(q.frequency_hz/1.5)**(-.25)
    q['source_phase_rad']=0.;q['constant_cycle']=0;q['phase_velocity_m_s']=speed*d/projection;q['group_velocity_m_s']=speed/1.25;q['distance_over_wavelength']=d*q.frequency_hz/q.phase_velocity_m_s;q.to_csv(out/'conditional_dispersion.csv',index=False)
    measured=pd.read_csv(out/'measured_phase_segments.csv');measured=measured[measured.segment_id=='S0--S2'].copy();measured['frequency_hz']=measured.frequency_hz.round(9);measured['measured_unwrapped_phase_rad']=-2*np.pi*measured.frequency_hz*projection/(500*(measured.frequency_hz/1.5)**(-.25));measured.to_csv(out/'measured_phase_segments.csv',index=False)
    options=dict(stabilization_s=0,repeat_cycles=20,slide_cycles=20,taper='hann',frequency_method='exact',direction=dict(velocity_min_m_s=200,velocity_max_m_s=1000,velocity_points=41),amplitude_rejection=False,reconstruction=dict(enabled=True,lower_frequency_hz=.8,step_hz=.1))
    marker=run_phase_stack(root,aligned_h5=hpath,**options);(out/'families_complete.json').write_text('{}');return root,hpath,options,marker

def test_native_measurements_and_audit(native):
    from ndstudio.phase_audit import audit
    from ndstudio.candidates import read_curves,auto_selection
    from ndstudio.inversion import validate_cases
    root,_,_,m=native;q,_=read_curves(root)
    assert m['low_frequency_measured_points']==5 and len(q)==11
    assert not q.is_model_prediction.any() and not q.is_prior_prediction.any() and q.projection_applied.all()
    assert (q.background_null_p<=.05).all() and (np.diff(q.phase_velocity_m_s)<=0).all()
    assert abs(q.phase_velocity_m_s.median()-500*(q.frequency_hz.median()/1.5)**(-.25))<15
    assert audit(root/'output/phase_guided')['pass_audit']
    selection=auto_selection(root,max_cases=1);assert len(selection['cases'])==len(validate_cases(q,selection,3,.25)[1])==1

def test_native_paired_inversion(native):
    from ndstudio.candidates import auto_selection
    from ndstudio.reconstruction import paired_compare
    root,_,_,_=native;selection=auto_selection(root,max_cases=1);selection['settings'].update(layers=[1],starts=1,max_evaluations=12)
    dest=root/'paired_solver_smoke';r=paired_compare(root,selection,dest)
    b=pd.read_csv(dest/'before_common_curves.csv');a=pd.read_csv(dest/'after_common_curves.csv');extended=pd.read_csv(dest/'after_extended_curves.csv')
    pd.testing.assert_frame_equal(b[['curve_id','frequency_hz']],a[['curve_id','frequency_hz']]);assert extended.frequency_hz.min()<a.frequency_hz.min()
    for name in ('before_common','after_common','after_extended'):assert name in r and (dest/name/'structure_layers.csv').exists()
    r=paired_compare(root,dict(settings={},cases=[]),dest)
    assert r['common_points']==r['extended_points']==0
    assert all(not (dest/n).exists() for n in ('before_common','after_common','after_extended'))

def test_direction_unavailable(native,monkeypatch):
    from ndstudio.phase_stack import run_phase_stack
    from ndstudio.direction import DirectionSession
    root,h,options,_=native;original=DirectionSession.get
    def unavailable(self,*args):
        d=original(self,*args);d['quality'][:]=False;return d
    monkeypatch.setattr(DirectionSession,'get',unavailable)
    r=run_phase_stack(root,output=root/'rejected_phase',aligned_h5=h,**options)
    assert r['inversion_handoff_points']==r['low_frequency_measured_points']==0

def test_ui_and_mode_requirements(tmp_path):
    from PySide6.QtWidgets import QApplication
    from ndstudio.gui import MainWindow
    from ndstudio.project import new_project
    from ndstudio.phase_stack import run_phase_stack
    from dispersion.common import write_json
    app=QApplication.instance() or QApplication([]);root=new_project(tmp_path/'ui');w=MainWindow(root)
    assert not w.phasereconstruct.isChecked();w.phasereconstruct.setChecked(True);w.phaselow.setValue(.25);w.phaselowstep.setValue(.025)
    opts=w.phase_settings();assert opts['frequency_method']=='exact' and opts['direction']['mode']=='filter'
    write_json(root/'phase_settings.json',opts);w.load(root);assert w.phase_settings()['reconstruction']==opts['reconstruction'];w.close()
    with pytest.raises(ValueError,match='方向筛窗'):run_phase_stack(root,reconstruction={'enabled':True})
