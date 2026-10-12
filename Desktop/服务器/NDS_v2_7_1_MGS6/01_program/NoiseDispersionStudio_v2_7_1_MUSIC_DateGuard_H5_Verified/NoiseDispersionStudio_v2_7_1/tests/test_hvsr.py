"""Scientific acceptance checks: units, timing, gaps, rotation and non-uniqueness."""
from pathlib import Path
from io import BytesIO
import json
import numpy as np
import pandas as pd
import h5py
import pytest
from obspy import Trace,Stream,UTCDateTime
from obspy.core.inventory.response import Response
from ndstudio.hvsr import (discover,run_hvsr,settings,spectral_ratio,orientation_matrix,
                           peak_table,_read_interval,_independent_indices)
from ndstudio.hvsr_example import create_example
from ndstudio.hvsr_structure import analyze_structure,ellipticity_forward,ellipticity_search


def simple(folder,starts=(0,0,0),locations=('','',''),channels=('HHN','HHE','HHZ'),units='m/s',n=12000,fs=20):
    folder.mkdir(parents=True,exist_ok=True);p=folder/'wave.h5';rng=np.random.default_rng(8);base=rng.normal(size=n)
    with h5py.File(p,'w') as h:
        for j,ch in enumerate(channels):
            d=h.create_dataset(str(j),data=base*(2 if j<2 else 1))
            d.attrs.update(role='waveform',seed_id=f'XX.A.{locations[j]}.{ch}',sampling_rate_hz=fs,
                           start_utc=str(UTCDateTime('2026-01-01')+starts[j]),units=units,precorrected=units!='counts',
                           azimuth_deg=90. if j==1 else 0.,dip_deg=-90. if j==2 else 0.)
    return p


def test_full_spectral_example_recovers_peak_rejects_injected_windows_and_preserves_sources(tmp_path):
    root=create_example(tmp_path/'demo');out=root/'output/hvsr';s=json.loads((out/'hvsr_summary.json').read_text())
    assert s['accepted_windows']==28 and s['independent_windows']==28
    assert abs(s['selected_peak']['f0_hz']-3)<.1
    assert s['selected_peak']['clear_peak'] and s['selected_peak']['clarity_pass_count']==6
    qc=pd.read_csv(out/'window_qc.csv');assert not qc.iloc[4].accepted and not qc.iloc[12].accepted
    assert s['structure']['depth_at_f0_range_m'][0]<20<s['structure']['depth_at_f0_range_m'][1]
    assert (out/'hvsr_complete.json').exists() and (out/'hvsr_report.html').exists()
    import hashlib
    for p,digest in s['source_sha256'].items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==digest


@pytest.mark.parametrize('locations,channels',[(('','','01'),('HHN','HHE','HHZ')),(('','',''),('HHN','HHE','BHZ'))])
def test_never_mix_instrument_location_or_family(tmp_path,locations,channels):
    simple(tmp_path/'raw',locations=locations,channels=channels)
    assert not any(r['ready'] for r in discover(tmp_path/'raw')['instruments'])


def test_missing_vertical_has_no_complete_instrument(tmp_path):
    p=simple(tmp_path/'raw')
    with h5py.File(p,'a') as h:del h['2']
    with pytest.raises(ValueError,match='三分量'):run_hvsr(tmp_path/'raw',tmp_path/'out',config=dict(fmax_hz=7))


def test_fractional_sample_clock_error_is_not_silently_interpolated(tmp_path):
    simple(tmp_path/'raw',starts=(0,.0125,0))
    with pytest.raises(ValueError,match='有效时窗'):run_hvsr(tmp_path/'raw',tmp_path/'out',config=dict(fmin_hz=.3,fmax_hz=7))
    qc=pd.read_csv(tmp_path/'out/window_qc.csv');assert qc.reason.str.contains('时间网格不齐').all()


@pytest.mark.parametrize('method',['rms','total'])
def test_rms_and_total_are_invariant_to_horizontal_rotation(method):
    x=np.random.default_rng(7).normal(size=(3,6000));c=settings(dict(fmin_hz=.3,fmax_hz=15,horizontal=method))
    angle=.61;r=np.array([[np.cos(angle),np.sin(angle),0],[-np.sin(angle),np.cos(angle),0],[0,0,1]])
    np.testing.assert_allclose(spectral_ratio(x,50,c)[1],spectral_ratio(r@x,50,c)[1],rtol=2e-14)


def test_oriented_1_2_recover_neu_and_coplanar_components_fail():
    rows=[dict(azimuth_deg=30.,dip_deg=0.),dict(azimuth_deg=120.,dip_deg=0.),dict(azimuth_deg=0.,dip_deg=-90.)]
    m=orientation_matrix(rows);x=np.random.default_rng(2).normal(size=(3,100))
    np.testing.assert_allclose(np.linalg.solve(m,m@x),x,atol=1e-15)
    with pytest.raises(ValueError,match='退化'):orientation_matrix([rows[0],rows[0],rows[2]])


def test_same_samples_in_overlap_allowed_conflicting_overlap_fails(tmp_path):
    p=simple(tmp_path/'raw');records=discover(tmp_path/'raw')['instruments'][0]['channels']['HHN'];r=records[0]
    y,_=_read_interval(records+records,r['start_ns'],100,20,{})
    with h5py.File(p,'a') as h:
        h.copy(h['0'],'conflict');h['conflict'][50]=9999
    bad=dict(r,dataset='conflict')
    with pytest.raises(ValueError,match='重叠波形冲突'):_read_interval(records+[bad],r['start_ns'],100,20,{})


def test_overlapping_windows_do_not_inflate_statistical_window_count():
    assert _independent_indices(np.arange(0,121,15),60).tolist()==[0,4,8]


def test_flat_curve_and_edge_peak_do_not_generate_structure(tmp_path):
    f=np.geomspace(.3,15,120);curves=np.ones((10,len(f)))
    assert peak_table(f,np.ones(len(f)),np.ones(len(f)),curves,60,10)==[]
    r=analyze_structure(f,np.ones(len(f)),np.ones(len(f)),None,settings(),tmp_path)
    assert r['status']=='no_reliable_peak' and not (tmp_path/'depth_vs_tradeoff.csv').exists()


def test_counts_need_complete_response_and_individual_gains_are_removed(tmp_path):
    from test_download_handoff import inventory,DATE
    from h5data_handoff import export_h5
    inv=inventory(stations=('A',),channels=('HHN','HHE','HHZ'),fs=20)
    gains=[1000,500,2000];physical=[2,4,1];st=Stream();base=np.round(np.random.default_rng(9).normal(size=12000)*100).astype(np.int32)
    for j,ch in enumerate(inv[0][0]):
        ch.azimuth=90 if j==1 else 0
        ch.response=Response.from_paz([],[],stage_gain=gains[j],normalization_factor=1.,input_units='M/S',output_units='COUNTS')
        st+=Trace(base*gains[j]*physical[j],dict(network='XX',station='A',location='',channel=ch.code,sampling_rate=20,starttime=DATE))
    raw=tmp_path/'raw';raw.mkdir();mseed=raw/'native.mseed';st.write(str(mseed),format='MSEED');buf=BytesIO();inv.write(buf,format='STATIONXML')
    export_h5(mseed,mseed.with_suffix('.h5'),stationxml=buf.getvalue(),strict=True)
    run_hvsr(raw,tmp_path/'out','XX.A..HH',dict(fmin_hz=.3,fmax_hz=7,window_s=60))
    d=pd.read_csv(tmp_path/'out/hvsr_curve.csv')
    np.testing.assert_allclose(d.hv,np.sqrt(10),rtol=.003)
    np.testing.assert_allclose(d.n_over_v,2,rtol=.003);np.testing.assert_allclose(d.e_over_v,4,rtol=.003)
    p=simple(tmp_path/'incomplete',units='counts')
    assert discover(p)['errors'] and not discover(p)['instruments']


def test_band_above_native_nyquist_is_rejected(tmp_path):
    simple(tmp_path/'raw')
    with pytest.raises(ValueError,match='0.45'):run_hvsr(tmp_path/'raw',tmp_path/'out',config=dict(fmax_hz=12))


def test_quarter_wave_retains_vs_depth_nonuniqueness(tmp_path):
    c=settings(dict(vs_min_m_s=200,vs_max_m_s=400));f=np.geomspace(.3,15,100)
    p=dict(clear_peak=True,f0_hz=2,sigma_f_hz=.1,left_drop=True,right_drop=True)
    s=analyze_structure(f,np.ones(100),np.ones(100),p,c,tmp_path)
    assert s['depth_at_f0_range_m']==[25,50]
    assert 'vs30' not in s


def test_disba_conditional_ellipticity_search_recovers_an_injected_model_trial(tmp_path):
    c=settings(dict(ellipticity_trials=20,fmin_hz=.8,fmax_hz=10,vs_min_m_s=200,vs_max_m_s=350,thickness_min_m=5,thickness_max_m=40))
    rng=np.random.default_rng(c['ellipticity_seed']);v=np.exp(rng.uniform(np.log(200),np.log(350)));h=np.exp(rng.uniform(np.log(5),np.log(40)))
    f=np.geomspace(.8,10,48);pred=ellipticity_forward(f,h,v,1000,c)
    r=ellipticity_search(f,pred,np.ones(48)*1.1,c,tmp_path)
    assert r['best_model']['log_rms']<1e-6 and r['acceptable_fit']
    assert r['best_model']['thickness_m']==pytest.approx(h)
    assert r['best_model']['vs_m_s']==pytest.approx(v)


@pytest.mark.parametrize('left,right',[(False,True),(True,False),(False,False)])
def test_sesame_clear_but_one_sided_peak_does_not_generate_structural_depth(tmp_path,left,right):
    f=np.geomspace(.05,20,100)
    peak=dict(clear_peak=True,clarity_pass_count=5,f0_hz=.424,left_drop=left,right_drop=right)
    result=analyze_structure(f,np.ones(100),np.ones(100),peak,settings(),tmp_path)
    assert result['status']=='incomplete_peak_shape'
    assert 'depth_at_f0_range_m' not in result
    assert not (tmp_path/'depth_vs_tradeoff.csv').exists()


def test_response_padding_loss_is_distinguished_from_missing_core_waveform(tmp_path):
    simple(tmp_path/'raw');rows=discover(tmp_path/'raw')['instruments'][0]['channels']['HHN'];r=rows[0]
    with pytest.raises(ValueError,match='缓冲区不完整'):
        _read_interval(rows,r['start_ns']-int(1e9),100,20,{},core_slice=(20,100))
    with pytest.raises(ValueError,match='分析窗内缺测'):
        _read_interval(rows,r['start_ns']-int(1e9),100,20,{},core_slice=(10,100))
