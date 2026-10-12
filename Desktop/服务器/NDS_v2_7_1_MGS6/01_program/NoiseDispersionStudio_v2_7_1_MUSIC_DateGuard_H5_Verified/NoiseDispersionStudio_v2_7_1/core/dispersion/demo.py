"""Known dispersive wavefield for end-to-end recovery and noise rejection."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.special import j0
from .common import csv,npz,write_json,DEFAULTS

def generate(directory,duration_s=14400,field="plane"):
    directory=Path(directory).resolve();directory.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(1901);fs=40.;n=round(duration_s*fs);f=np.fft.rfftfreq(n,1/fs)
    xy=np.array([[0,0],[1000,0],[2000,0],[3000,0]],float);distance=np.linalg.norm(xy[:,None]-xy[None,:],axis=-1)
    cphase=600.-20*f;use=(f>=.25)&(f<=8);fu=f[use];k=2*np.pi*fu/cphase[use]
    independent=(rng.normal(size=(len(fu),4))+1j*rng.normal(size=(len(fu),4)))/np.sqrt(2)
    if field=="plane":
        # Known propagation parallel to the line; its phase offset is exactly zero.
        common=(rng.normal(size=len(fu))+1j*rng.normal(size=len(fu)))/np.sqrt(2)
        fourier=common[:,None]*np.exp(-1j*k[:,None]*xy[None,:,0])+.05*independent
    elif field=="noise": fourier=independent
    elif field=="isotropic":
        cov=j0(k[:,None,None]*distance[None])+np.eye(4)[None]*1e-8
        fourier=np.einsum("fij,fj->fi",np.linalg.cholesky(cov),independent)
    else: raise ValueError("Unknown synthetic field")
    spectrum=np.zeros((4,len(f)),complex);spectrum[:,use]=fourier.T
    offsets=[0.,.007,.013,.021];stations=[]
    for i,label in enumerate(("D01","D02","D03","D04")):
        # Samples represent the true field at the station's asynchronous start.
        x=np.fft.irfft(spectrum[i]*np.exp(2j*np.pi*f*offsets[i]),n)*np.sqrt(n)
        valid=np.ones(n,bool)
        if i==1:
            x[round(1800*fs)]+=200*np.std(x)
            x[round(2400*fs):round(2402*fs)]=0
        if i==2: valid[round(3600*fs):round(3602*fs)]=False
        from obspy import UTCDateTime
        start=str(UTCDateTime("2026-01-01T00:00:00Z")+offsets[i])
        npz(directory/f"{label}.npz",data=x.astype(np.float32),valid=valid,fs_hz=fs,start_utc=start)
        stations.append(dict(station=label,file=label+".npz",seed_id=f"XX.{label}..HHZ",x_m=xy[i,0],y_m=xy[i,1],
                             response_xml="",precorrected=True,units="m/s^2",polarity=1,polarity_verified=True,clock_offset_s=0))
    csv(directory/"stations.csv",stations)
    config={"stations_csv":"stations.csv","data_root":".","output":"output","filter_guard_s":6.,
        "window_s":240.,"max_lag_s":12.,"bandpass_hz":[.5,7.],"min_clean_segment_s":40.,
        "response":{"pre_filt_hz":[.1,.3,7.5,9.]},
        "scan":{"min_hz":2.,"max_hz":5.,"step_hz":.05,"arrival_max_s":8.,"noise_lag_s":[9.,11.]},
        "cohorts":{"clock_groups":False,"states":0},
        "validation":{"block_s":1800.,"min_windows_per_block":2,"bootstrap_draws":99},
        "controls":{"trials":20},"phase":{"source_phase_rad":[0. if field=="plane" else np.pi/4],"source_phase_verified":field=="plane"}}
    write_json(directory/"config.json",config)
    write_json(directory/"truth.json",dict(type=field+"_scalar_random_field_not_a_layered_earth",
        phase_velocity_formula_m_s="600 - 20*f_hz",group_velocity_formula_m_s="(600-20*f_hz)^2 / 600",
        asynchronous_starts_s=offsets,inserted_single_sample_spike_station="D02",flatline_station="D02",gap_station="D03"))
    return directory/"config.json"

def check(directory):
    directory=Path(directory);curves=pd.read_csv(directory/"output/conditional_dispersion.csv")
    truth=json.loads((directory/"truth.json").read_text(encoding='utf-8-sig'))
    if truth["type"].startswith("noise"):
        assert curves.empty,"Independent noise was incorrectly retained as dispersion"
        write_json(directory/"output/synthetic_truth_check.json",dict(checks_passed=True,independent_noise_rejected=True))
        return pd.DataFrame([dict(independent_noise_rejected=True)])
    rows=[]
    for cid,q in curves.groupby("curve_id"):
        error=abs(q.phase_velocity_m_s/(600-20*q.frequency_hz)-1)
        rows.append(dict(curve_id=cid,pair_index=int(q.pair_index.iloc[0]),median_fraction=float(np.median(error)),maximum_fraction=float(error.max())))
    result=pd.DataFrame(rows)
    if result.empty: raise AssertionError("Synthetic dispersive field yielded no validated candidate")
    best=result.sort_values("median_fraction").groupby("pair_index").head(1)
    recovered=best[(best.median_fraction<.05)&(best.maximum_fraction<.1)]
    assert recovered.pair_index.nunique()>=2,"Known dispersion not recovered on >=2 pairs"
    audit=pd.read_csv(directory/"output/record_audit.csv")
    assert audit.repaired_native_samples.sum()>=1 and audit.masked_native_samples.sum()>0
    quality=pd.read_csv(directory/"output/window_quality.csv")
    assert ((quality.clean_fraction<.99)&quality.accepted).any(),"Local anomalies unnecessarily removed all affected windows"
    write_json(directory/"output/synthetic_truth_check.json",dict(checks_passed=True,matched_pairs=int(recovered.pair_index.nunique()),
        best_median_fraction=float(best.median_fraction.min()),truth_used_only_for_testing=True,
        no_layered_earth_or_unique_cycle_inference=True))
    csv(directory/"output/synthetic_closest_families.csv",best)
    return best
