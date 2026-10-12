"""Station-consistent amplitude-spectrum-preserving wrong-time controls."""
import numpy as np
import pandas as pd
from .common import ensure_ready,csv,write_json
from .prepare import normalized_fft,cc_from_fft,window_segments
from .scan import classify,envelopes,frequency_grid

def packet_statistic(cc,branch,reference,c,f):
    env=envelopes(cc,c,f)[branch];t=np.arange(env.shape[-1])/c["sample_rate_hz"]
    noise=(t>=c["scan"]["noise_lag_s"][0])&(t<=c["scan"]["noise_lag_s"][1]);values=[]
    for k,arrival in enumerate(reference):
        use=abs(t-arrival)<=max(.12,.1*arrival)
        nr=np.sqrt(np.mean(env[k,noise]**2))
        values.append(float(env[k,use].max()/max(nr,1e-20)) if use.any() else np.nan)
    return float(np.nanmedian(values))

def controls(c):
    out,cache,meta=ensure_ready(c,"scan_complete.json")
    catalog=pd.read_csv(out/"measured_segment_catalog.csv");phases=pd.read_csv(out/"measured_phase_segments.csv")
    if catalog.empty:
        csv(out/"timeshift_summary.csv",[],columns=["segment_id","observed_median_snr","null_q95","valid_trials","background_threshold_exceeded"])
        write_json(out/"controls_complete.json",dict(trials=c["controls"]["trials"],segments=0));return
    z=np.load(cache/"window_routes.npz");starts=z["starts"];times=z["timestamps"];training=z["training"];baseline=z["ram5"]
    assignments=pd.read_csv(out/"cohort_assignments.csv");names=catalog.cohort.unique()
    fixed={name:assignments[name].to_numpy(bool) for name in names}
    data=np.load(cache/"aligned_filtered.npy",mmap_mode="r");valid=np.load(cache/"aligned_valid.npy",mmap_mode="r")
    pairs=pd.read_csv(out/"pairs.csv");stacks=np.load(out/"cohort_stacks.npz");model=None
    if (out/"state_classifier.npz").exists() and c["cohorts"]["states"]:
        model=dict(np.load(out/"state_classifier.npz"))
    ref={sid:q.sort_values("frequency_hz") for sid,q in phases.groupby("segment_id")}
    observed={r.segment_id:packet_statistic(stacks[f"{r.cohort}__ram5__combined"][r.pair_index],r.branch,
              ref[r.segment_id].arrival_s.to_numpy(),c,ref[r.segment_id].frequency_hz.to_numpy()) for r in catalog.itertuples()}
    trials=c["controls"]["trials"];rng=np.random.default_rng(c["controls"]["seed"]);rows=[];nulls={sid:[] for sid in ref}
    baseline_error=0.;infeasible=0
    for trial in range(trials):
        null=np.zeros_like(baseline);validwindow=np.ones(len(starts),bool)
        for win,start in enumerate(starts):
            segments=window_segments(data,valid,int(start),c);total=sum(b-a for a,b in segments)
            rebuilt=np.zeros_like(baseline[win])
            for a,b in segments:
                ft,f,nf=normalized_fft(np.asarray(data[:,a:b],float),c,c["routes"][0])
                if trial==0: rebuilt+=cc_from_fft(ft,nf,b-a,c,pairs)*(b-a)/total
                # Delays are separated on the padded Fourier circle. This is a
                # diagnostic surrogate, not a calibrated surface-wave null model.
                period=nf/c["sample_rate_hz"];spacing=period/len(data)
                minimum=c["max_lag_s"]+c["scan"]["arrival_max_s"]+5
                if spacing<=minimum+2:
                    validwindow[win]=False;infeasible+=1;break
                jitter=min(10,(spacing-minimum)/4)
                delay=rng.permutation(len(data))*spacing+rng.uniform(-jitter,jitter,len(data))
                shifted=ft*np.exp(-2j*np.pi*delay[:,None]*f[None,:])
                null[win]+=cc_from_fft(shifted,nf,b-a,c,pairs)*(b-a)/total
            if trial==0 and validwindow[win]: baseline_error=max(baseline_error,float(np.max(abs(rebuilt-baseline[win]))))
        labels=classify(null,c,model) if model is not None else None
        trialmasks={}
        for name in names:
            selection=(labels==int(name[1:])) if name.startswith("G") and labels is not None else fixed[name]
            trialmasks[name]=selection&validwindow
        for r in catalog.itertuples():
            mask=trialmasks[r.cohort];count=int(mask.sum());stat=np.nan
            if count>=c["validation"]["min_windows_per_block"]:
                q=ref[r.segment_id];stat=packet_statistic(null[mask,r.pair_index].mean(axis=0),r.branch,q.arrival_s.to_numpy(),c,q.frequency_hz.to_numpy())
            nulls[r.segment_id].append(stat);rows.append(dict(trial=trial,segment_id=r.segment_id,cohort=r.cohort,windows=count,median_snr=stat))
        print("controls: wrong-time trial",trial+1,"/",trials,flush=True)
    if baseline_error>2e-6: raise ValueError("CCF reproduction failed before null controls")
    summary=[]
    for sid,values in nulls.items():
        x=np.asarray(values);x=x[np.isfinite(x)];obs=observed[sid];q=float(np.quantile(x,.95)) if len(x) else np.nan
        summary.append(dict(segment_id=sid,observed_median_snr=obs,null_q95=q,valid_trials=len(x),
            plus_one_exceedance_rank=(1+int(sum(x>=obs)))/(len(x)+1) if len(x) else np.nan,
            calibrated_p_value=False,minimum_valid_trials=20,background_threshold_exceeded=bool(len(x)>=20 and obs>q),
            classifier_reapplied_without_refitting=True,station_delay_closure_preserved=True))
    csv(out/"timeshift_trials.csv",rows,columns=["trial","segment_id","cohort","windows","median_snr"])
    csv(out/"timeshift_summary.csv",summary)
    write_json(out/"controls_complete.json",dict(trials=trials,segments=len(catalog),maximum_baseline_reproduction_error=baseline_error,
               infeasible_segment_trials=infeasible,calibrated_p_values=False,missing_controls_stored_as_nan=True))
