"""Full time axes and an independent translation of the user's MATLAB core."""
import json,time,os
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
from ndstudio.phase_stack import goertzel_windows,run_phase_stack
from ndstudio.phase_screening import screen_cross,accumulation_series,geometric_fraction
from ndstudio.phase_records import available_phase_extent,phase_window_geometry,full_window_starts
from ndstudio.project import new_project
from dispersion.common import write_json,parse_utc


def matlab_core(x,y,f,fs,n,step):
    """Literal scalar recurrence, MAD delta[2:end], sign split and +/-90 gate."""
    k=int(np.floor(f*n/fs+.5));w=2*np.pi*k/n;cross=[]
    for first in range(0,len(x)-n+1,step):
        gs=[]
        for data in (x[first:first+n],y[first:first+n]):
            data=data-data.mean();s1=s2=0.
            for v in data:s0=v+2*np.cos(w)*s1-s2;s2=s1;s1=s0
            gs.append(s1-np.exp(-1j*w)*s2)
        cross.append(gs[1]*np.conj(gs[0]))
    c=np.array(cross);amp=abs(c);median=np.median(amp);scale=1.4826*np.median(abs(amp-median))
    ath=median+5*scale if scale>0 else 5*median
    delta=np.r_[0.,abs(np.diff(amp))];median=np.median(delta[1:]);scale=1.4826*np.median(abs(delta[1:]-median))
    dth=median+5*scale if scale>0 else 5*median
    good=(amp<=ath)&(delta<=dth);phase=np.angle(c);masks=[]
    for sign in (phase>=0,phase<0):
        base=good&sign;total=c[base].sum();center=np.angle(total) if abs(total)>0 else np.nan
        difference=np.angle(c*np.exp(-1j*center));masks.append(base&(abs(difference)<=np.pi/2))
    return c,masks


def test_matlab_scalar_core_matches_goertzel_and_branch_screening():
    rng=np.random.default_rng(122);fs=20.;f=2.37;n=int(np.floor(20*fs/f+.5));step=9
    x=rng.normal(size=7000);y=.4*x+rng.normal(size=len(x));y[4000:4050]*=90
    cross,(positive,negative)=matlab_core(x,y,f,fs,n,step)
    starts=full_window_starts(len(x),n,step);fk=int(np.floor(f*n/fs+.5))*fs/n
    ga=goertzel_windows(x,starts,n,fk,fs,np.ones(n));gb=goertzel_windows(y,starts,n,fk,fs,np.ones(n))
    actual=gb*np.conj(ga);np.testing.assert_allclose(actual,cross,rtol=1e-10,atol=1e-9)
    result=screen_cross(actual,.1,-.3,.3,mode='matlab')
    np.testing.assert_array_equal(result['accepted_positive'],positive);np.testing.assert_array_equal(result['accepted_negative'],negative)
    for key,mask in [('accepted_positive',positive),('accepted_negative',negative)]:
        np.testing.assert_allclose(actual[result[key]].sum(),cross[mask].sum(),rtol=1e-11,atol=1e-7)


def test_curve_gate_keeps_original_complex_sign_and_separate_branches():
    z=np.exp(1j*np.radians([178,179,-179,-178,0,4,120]));r=screen_cross(z,np.pi,-.1,.1,amplitude_rejection=False)
    assert r['accepted_positive'].tolist()==[True,True,False,False,False,False,False]
    assert r['accepted_negative'].tolist()==[False,False,True,True,False,False,False]
    assert np.angle(z[r['accepted_positive']].sum())>0 and np.angle(z[r['accepted_negative']].sum())<0
    assert geometric_fraction(np.pi,-.1,.1,'nonnegative')==pytest.approx(.1/(2*np.pi))
    assert r['accepted'].sum()==2


def test_cumulative_display_includes_all_accepted_windows():
    rng=np.random.default_rng(11);z=rng.normal(size=200123)+1j*rng.normal(size=200123);mask=rng.random(len(z))>.2
    r=accumulation_series(z,mask,limit=31);all_z=z[mask]
    assert len(r['count'])==31 and r['count'][-1]==len(all_z)
    np.testing.assert_allclose(r['raw'],np.cumsum(all_z)[r['count']-1],rtol=1e-10,atol=1e-9)
    np.testing.assert_allclose(r['unit'],np.cumsum(all_z/abs(all_z))[r['count']-1],rtol=1e-10,atol=1e-9)


def raw_fixture(directory):
    root=new_project(directory);origin=parse_utc('2020-01-01');data=root/'canonical';data.mkdir();rows=[]
    # A--B must start 350 s earlier than A--C, independent of the other pair.
    for station,start,duration,delay in [('A',0.,1800.,0.),('B',50.,1750.,.06),('C',400.,1300.,-.06)]:
        t=np.arange(round(duration*20))/20+start
        x=np.cos(2*np.pi*.2*(t-delay))+.7*np.cos(2*np.pi*.4*(t-delay))
        if station=='B':x[(t>=1100)&(t<1120)]=np.nan
        np.savez_compressed(data/(station+'.npz'),data=x,fs_hz=20.,start_utc=str(origin+start))
        rows.append(dict(station=station,file=station+'.npz',seed_id='XX.'+station+'..BHZ',x_m=len(rows)*100,y_m=0,units='m/s',precorrected=True,polarity=1))
    pd.DataFrame(rows).to_csv(root/'stations.csv',index=False)
    c=json.loads((root/'config.json').read_text());c['start_utc']=str(origin+900);c['end_utc']=str(origin+1000);c['bandpass_hz']=[.05,7];c['response']['pre_filt_hz']=[.005,.02,7.5,9];write_json(root/'config.json',c)
    out=root/'output';out.mkdir();curves=[];measured=[]
    for index,pair,delay in [(0,'A--B',.06),(1,'A--C',-.06)]:
        for f in (.2,.4):
            theta=-2*np.pi*f*delay;d=100*(index+1);velocity=2*np.pi*f*d/(np.pi/4-theta+2*np.pi)
            curves.append(dict(curve_id=pair,segment_id=pair,pair_index=index,pair=pair,branch='positive',frequency_hz=f,
                phase_velocity_m_s=velocity,velocity_min_m_s=velocity*.9,velocity_max_m_s=velocity*1.1,source_phase_rad=np.pi/4,constant_cycle=1))
            measured.append(dict(segment_id=pair,frequency_hz=f,measured_unwrapped_phase_rad=theta,relative_phase_q025_rad=-.2,relative_phase_q975_rad=.2))
    pd.DataFrame(curves).to_csv(out/'conditional_dispersion.csv',index=False);pd.DataFrame(measured).to_csv(out/'measured_phase_segments.csv',index=False)
    pd.DataFrame([dict(pair_index=i,distance_m=100*(i+1)) for i in (0,1)]).to_csv(out/'pairs.csv',index=False)
    return root


def test_full_raw_extent_pair_specific_starts_gaps_and_actual_phase(tmp_path):
    from ndstudio.phase_audit import audit
    root=raw_fixture(tmp_path/'project');origin,end,bounds=available_phase_extent(root,{'A','B','C'})
    assert end-origin==1800 and len(bounds)==3
    r=run_phase_stack(root,stabilization_s=600,amplitude_rejection=False)
    assert r['time_mode']=='full' and r['record_duration_s']==1800
    intervals={q['pair']:q for q in r['pair_processing_intervals']}
    assert parse_utc(intervals['A--B']['processing_start_utc'])-origin==650
    assert parse_utc(intervals['A--C']['processing_start_utc'])-origin==1000
    gates=pd.read_csv(root/'output/phase_guided/phase_gate_summary.csv')
    for g in gates.itertuples():
        with np.load(root/'output/phase_guided'/g.windows_file) as z:
            start=round((parse_utc(intervals[g.pair]['processing_start_utc'])-origin)*20)
            assert int(z['start_sample'][0])==start
            assert np.all(np.diff(z['start_sample'])==int(z['slide_samples']))
            assert int(z['start_sample'][-1])+int(z['window_samples'])<=int(z['processing_end_sample_exclusive'])
            if g.pair=='A--B':assert not z['valid'][(z['start_sample']<1120*20)&(z['start_sample']+int(z['window_samples'])>1100*20)].any()
        assert abs(g.after_raw_phase_deg-g.predicted_phase_deg)<.1
    assert audit(root/'output/phase_guided')['pass_audit']
    assert not list((root/'output/phase_guided').glob('.phase_work_*'))


def test_old_ui_times_do_not_limit_full_mode_and_primary_button_routes_phase(tmp_path,monkeypatch):
    from PySide6.QtWidgets import QApplication
    from ndstudio.gui import MainWindow
    app=QApplication.instance() or QApplication([]);root=new_project(tmp_path/'project')
    write_json(root/'phase_settings.json',dict(start='2020-01-01',end='2020-01-02',repeat_cycles=10))
    w=MainWindow(root);w.show();app.processEvents();w.tabs.setCurrentWidget(w.phasepage)
    assert w.phase_settings()['time_mode']=='full' and w.phase_settings()['start'] is None
    assert not w.phaseinterval.isVisible();assert not w.phaseextras.isChecked()
    assert w.phasestabilization.value()==600 and w.phaseamp.isChecked()
    calls=[];monkeypatch.setattr(w,'start',lambda *args,**kw:calls.append(args));w.primary_analysis();assert calls==[('phase-stack',)]
    w.phaseextras.setChecked(True);w.phasetimemode.setCurrentIndex(1);app.processEvents();assert w.phaseinterval.isVisible()
    assert w.phase_settings()['start']=='2020-01-01'
    w.close()


def test_actual_gui_worker_processes_full_records_without_utc_entry(tmp_path,monkeypatch):
    from PySide6.QtWidgets import QApplication,QMessageBox
    from PySide6.QtCore import QProcess
    from ndstudio.gui import MainWindow
    app=QApplication.instance() or QApplication([]);errors=[];monkeypatch.setattr(QMessageBox,'warning',lambda *a:errors.append(a[-1]))
    root=raw_fixture(tmp_path/'project');w=MainWindow(root);w.tabs.setCurrentWidget(w.phasepage);w.phaseamp.setChecked(False)
    w.primary_analysis();assert w.current_job=='phase-stack';deadline=time.monotonic()+35
    while w.process.state()!=QProcess.ProcessState.NotRunning and time.monotonic()<deadline:app.processEvents();time.sleep(.01)
    app.processEvents();assert w.process.state()==QProcess.ProcessState.NotRunning and w.process.exitCode()==0,w.buffer
    r=json.loads((root/'output/phase_guided/phase_stack_complete.json').read_text())
    assert r['time_mode']=='full' and r['record_duration_s']==1800 and '600' in w.phaselast.text()
    saved=json.loads((root/'phase_settings.json').read_text());assert saved['start'] is None and saved['end'] is None and saved['stabilization_s']==600
    assert not errors;w.close()


def test_blocked_native_reader_keeps_phase_across_six_hour_boundary(tmp_path,monkeypatch):
    import ndstudio.phase_stack as mod
    from dispersion.common import DEFAULTS,merge
    c=merge(DEFAULTS,{'data_root':str(tmp_path),'response':{'output':'VEL'}});write_json(tmp_path/'config.json',c)
    origin=parse_utc('2020-01-01');n=round(21780*20);t=np.arange(n)/20
    np.savez_compressed(tmp_path/'A.npz',data=np.cos(2*np.pi*2*t),fs_hz=20.,start_utc=str(origin))
    row=dict(station='A',seed_id='XX.A..HHZ',file='A.npz',precorrected=True,units='m/s',polarity=1,clock_offset_s=0.)
    monkeypatch.setattr(mod,'load_config',lambda _:c);monkeypatch.setattr(mod,'station_table',lambda _:(pd.DataFrame([row]),None,None))
    work=tmp_path/'work';work.mkdir();resources=[]
    arrays,fs,utc,meta=mod.align_project(tmp_path,{'A'},str(origin),str(origin+21780),tmp_path/'aligned.h5',work,resources)
    seam=21600*20;x=arrays['A'][seam-400:seam+400]
    assert np.isfinite(x).all();assert np.corrcoef(x,np.cos(2*np.pi*2*t[seam-400:seam+400]))[0,1]>.999999
    assert len(list(work.glob('*_data.npy')))==1
    for resource in reversed(resources):
        if isinstance(resource,np.memmap):resource._mmap.close()
        else:resource.close()


def test_unavailable_windows_are_reported_without_inventing_a_phase(tmp_path):
    root=raw_fixture(tmp_path/'project')
    r=run_phase_stack(root,stabilization_s=2000.)
    assert r['accepted_window_evaluations']==0
    gates=pd.read_csv(root/'output/phase_guided/phase_gate_summary.csv')
    assert gates.after_raw_phase_deg.isna().all() and gates.unavailable_reason.str.contains('不足').all()
    assert not (root/'output/phase_guided/velocity_comparison/report.html').exists()


def test_native_failure_closes_disk_mappings_and_hides_stale_report(tmp_path,monkeypatch):
    import ndstudio.phase_stack as mod
    root=raw_fixture(tmp_path/'project');dest=root/'output/phase_guided';dest.mkdir()
    (dest/'phase_report.html').write_text('STALE')
    def unavailable(*args):raise ValueError('Calibration unavailable')
    monkeypatch.setattr(mod,'response_table',unavailable)
    with pytest.raises(ValueError,match='Calibration unavailable'):mod.run_phase_stack(root)
    assert not list(dest.glob('.phase_work_*')) and not (dest/'phase_report.html').exists()
    assert not (dest/'phase_stack_complete.json').exists()
