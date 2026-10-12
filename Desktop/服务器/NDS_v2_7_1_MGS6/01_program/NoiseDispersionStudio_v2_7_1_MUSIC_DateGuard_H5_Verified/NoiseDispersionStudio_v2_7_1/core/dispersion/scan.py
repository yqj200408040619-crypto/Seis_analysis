"""Development-only multi-peak FTAN tracking and fixed-cohort phase validation."""
import json
import numpy as np
import pandas as pd
from scipy import signal
from scipy.cluster.vq import kmeans2
from scipy.fft import next_fast_len
from .common import sha,csv,npz,write_json,ensure_ready,bands,circ,invalidate_models

def frequency_grid(c):
    s=c["scan"];return np.round(np.arange(s["min_hz"],s["max_hz"]+s["step_hz"]*.25,s["step_hz"]),9)

def envelopes(cc,c,f,width=None):
    n=len(cc);nf=max(8192,next_fast_len(2*n));ff=np.fft.fftfreq(nf,1/c["sample_rate_hz"])
    s=c["scan"];width=s["ftan_width"] if width is None else width
    sigma=np.maximum(s["minimum_sigma_hz"],width*f)
    kernel=np.exp(-.5*((abs(ff)[None,:]-f[:,None])/sigma[:,None])**2)
    y=np.fft.ifft(np.fft.fft(cc,nf)[None,:]*kernel,axis=-1).real[:,:n]
    env=abs(signal.hilbert(y,axis=-1));m=n//2
    return {"positive":env[:,m:],"negative":env[:,m::-1]}

def ftan_peaks(cc,distance,c,f,width=None):
    s=c["scan"];fs=c["sample_rate_hz"];envs=envelopes(cc,c,f,width);t=np.arange(len(cc)//2+1)/fs
    low,high=s["group_velocity_bounds_m_s"]
    search=(t>=max(s["arrival_min_s"],distance/high))&(t<=min(s["arrival_max_s"],distance/low))
    ix=np.flatnonzero(search);noise=(t>=s["noise_lag_s"][0])&(t<=s["noise_lag_s"][1]);result={}
    for branch,env in envs.items():
        candidates=[]
        for k,freq in enumerate(f):
            if len(ix)<3: candidates.append([]);continue
            ys=env[k,ix];nr=float(np.sqrt(np.mean(env[k,noise]**2)))
            pp,_=signal.find_peaks(ys,distance=max(2,int(max(.2,.4/freq)*fs)),prominence=max(nr*.5,float(ys.max())*.06))
            strongest=int(ys.argmax())
            if strongest not in pp: pp=np.r_[pp,strongest]
            pp=pp[np.argsort(ys[pp])[::-1]][:s["max_peaks"]];nodes=[]
            for rank,q in enumerate(pp):
                j=ix[q];a,mid,b=env[k,j-1:j+2];curv=a-2*mid+b
                delta=float(np.clip(.5*(a-b)/curv,-.5,.5)) if curv < -1e-20 else 0.
                height=float(mid-.25*(a-b)*delta);snr=height/max(nr,1e-20)
                if snr<3 or height<.08*float(ys.max()): continue
                nodes.append(dict(arrival=float(t[j]+delta/fs),snr=snr,rank=rank+1,
                                  height_fraction=height/max(float(ys.max()),1e-20),edge=bool(q in (0,len(ix)-1))))
            candidates.append(nodes)
        result[branch]=candidates
    return result

def track_peaks(candidates,f,c):
    ending={};allstates=[];s=c["scan"]
    for k,nodes in enumerate(candidates):
        for j,node in enumerate(nodes):
            if node["edge"]: continue
            reward=1+min(float(np.log(node["snr"])),4);options=[dict(path=[(k,j)],score=reward)]
            for gap in (1,2):
                previous=k-gap
                if previous<0: continue
                for jj,prev in enumerate(candidates[previous]):
                    change=np.log(node["arrival"]/prev["arrival"])
                    if abs(change)>.18*gap: continue
                    for state in ending.get((previous,jj),[]):
                        penalty=(change/(.07*gap))**2+1.5*(gap-1)
                        if len(state["path"])>=2:
                            kk,jjj=state["path"][-2];prior=candidates[kk][jjj]
                            old=np.log(prev["arrival"]/prior["arrival"])/(previous-kk)
                            penalty+=.2*((change/gap-old)/.04)**2
                        options.append(dict(path=state["path"]+[(k,j)],score=state["score"]+reward-penalty))
            options.sort(key=lambda x:x["score"],reverse=True);ending[k,j]=options[:3];allstates+=options[:3]
    allstates.sort(key=lambda x:x["score"],reverse=True);selected=[]
    for state in allstates:
        path=state["path"]
        if len(path)<s["min_points"] or f[path[-1][0]]-f[path[0][0]]<s["min_band_hz"]-1e-8: continue
        values={k:candidates[k][j]["arrival"] for k,j in path};duplicate=False
        for chosen in selected:
            common=sorted(set(values)&set(chosen["values"]))
            if len(common)<.6*min(len(values),len(chosen["values"])): continue
            if np.mean([abs(values[k]-chosen["values"][k])/values[k]<.05 for k in common])>=.65:
                duplicate=True;break
        if not duplicate:
            selected.append(dict(**state,values=values))
            if len(selected)>=s["max_tracks"]: break
    return selected

def peak_arrays(peaks,reference):
    n=len(reference);arrival=np.full(n,np.nan);snr=np.zeros(n);edge=np.ones(n,bool)
    for k in np.flatnonzero(np.isfinite(reference)):
        if not peaks[k]: continue
        p=min(peaks[k],key=lambda x:abs(x["arrival"]-reference[k]))
        if abs(p["arrival"]-reference[k])>max(.12,.15*reference[k]): continue
        arrival[k]=p["arrival"];snr[k]=p["snr"];edge[k]=p["edge"]
    return arrival,snr,edge

def fixed_gate(t,a,b):
    if b<=a: raise ValueError("Empty gate")
    w=np.zeros_like(t);w[(t>=a)&(t<=b)]=1;ramp=min(.5,(b-a)/4)
    lo=(t>=a)&(t<a+ramp);hi=(t>b-ramp)&(t<=b)
    w[lo]=np.sin(np.pi/2*(t[lo]-a)/ramp)**2;w[hi]=np.cos(np.pi/2*(t[hi]-(b-ramp))/ramp)**2
    return w

def spectrum(cc,branch,limits,c,f):
    m=len(cc)//2;y=cc[m:] if branch=="positive" else cc[m::-1];t=np.arange(len(y))/c["sample_rate_hz"]
    return np.exp(-2j*np.pi*f[:,None]*t[None,:])@(y*fixed_gate(t,*limits))/c["sample_rate_hz"]

def phase_delay(h,f,index):
    result=np.full(len(f),np.nan)
    result[index]=-np.gradient(np.unwrap(np.angle(h[index])),f[index])/(2*np.pi)
    return result

def classify(cc,c,model):
    fs=c["sample_rate_hz"];p=c["cohorts"]
    filt=signal.sosfiltfilt(signal.butter(3,p["feature_band_hz"],btype="bandpass",fs=fs,output="sos"),cc,axis=-1)
    t=(np.arange(cc.shape[-1])-cc.shape[-1]//2)/fs
    features=filt[:,:,abs(t)<=p["feature_max_lag_s"]][:,:,::2].copy()
    features/=np.maximum(np.linalg.norm(features,axis=-1,keepdims=True),1e-20)
    x=features.reshape(len(cc),-1)
    if model is None: return x
    scores=(x-model["feature_mean"])@model["components"].T
    return np.sum((scores[:,None]-model["centers"][None])**2,axis=-1).argmin(axis=1)

def cohort_masks(cc,times,training,c,out):
    masks={"all":np.ones(len(cc),bool)};p=c["cohorts"];model=None
    if p["clock_groups"]:
        hours=(times%86400)/3600;a,b=p["night_hours_utc"]
        night=(hours>=a)|(hours<b) if a>b else (hours>=a)&(hours<b)
        if night.any() and (~night).any(): masks.update(night=night,day=~night)
    if p["states"]:
        x=classify(cc,c,None);mean=x[training].mean(axis=0)
        _,singular,v=np.linalg.svd(x[training]-mean,full_matrices=False)
        comp=v[:min(p["pca_components"],len(v))];scores=(x-mean)@comp.T
        if training.sum()<=p["states"]: raise ValueError("Not enough training windows for state classifier")
        centers,_=kmeans2(scores[training],p["states"],iter=60,minit="++",seed=c["seed"])
        labels=np.sum((scores[:,None]-centers[None])**2,axis=-1).argmin(axis=1)
        order=sorted(range(p["states"]),key=lambda k:-np.sum(training&(labels==k)))
        model=dict(feature_mean=mean,components=comp,centers=centers[order])
        labels=classify(cc,c,model)
        for k in range(p["states"]): masks[f"G{k}"]=labels==k
        npz(out/"state_classifier.npz",**model,training=training)
    table=pd.DataFrame(dict(window=np.arange(len(times)),timestamp=times,training=training,**masks))
    csv(out/"cohort_assignments.csv",table)
    return masks,model

PHASE_COLUMNS=["segment_id","gate_id","cohort","pair_index","pair","branch","frequency_hz",
    "measured_unwrapped_phase_rad","group_velocity_m_s","arrival_s","relative_phase_q025_rad",
    "relative_phase_q975_rad","bootstrap_valid_draws","calibrated_confidence_interval",
    "instrument_response_corrected","absolute_cycle_resolved","rayleigh_mode_verified"]
CATALOG_COLUMNS=["segment_id","gate_id","cohort","pair_index","pair","branch","fmin_hz","fmax_hz","points",
    "group_velocity_median_m_s","eligible_blocks","gate_start_s","gate_end_s","absolute_cycle_resolved","rayleigh_mode_verified"]

def scan(c):
    out,cache,meta=ensure_ready(c,"prepare_complete.json")
    print(f"本次扫描参数：最小连续带宽 {c['scan']['min_band_hz']:g} Hz；最少连续频点 {c['scan']['min_points']}；扫描步长 {c['scan']['step_hz']:g} Hz",flush=True)
    invalidate_models(out)
    for name in ("scan_complete.json","controls_complete.json","families_complete.json","verification.json"):
        (out/name).unlink(missing_ok=True)
    with np.load(cache/"window_routes.npz") as z: arrays={k:z[k] for k in z.files}
    times=arrays.pop("timestamps");training=arrays.pop("training");starts=arrays.pop("starts")
    f=frequency_grid(c);s=c["scan"];v=c["validation"];fs=c["sample_rate_hz"]
    pairs=pd.read_csv(out/"pairs.csv");masks,_=cohort_masks(arrays["ram5"],times,training,c,out)
    block_id=np.floor(times/v["block_s"]).astype(np.int64);block_ids=np.unique(block_id)
    rng=np.random.default_rng(c["seed"]+1);draws=rng.integers(0,len(block_ids),size=(v["bootstrap_draws"],len(block_ids)))
    npz(out/"block_bootstrap_draws.npz",draws=draws,block_ids=block_ids)
    required=[r["name"] for r in c["routes"] if r.get("required",False)]
    points=[];routes=[];catalog=[];measured=[];gates=[];localpeaks=[];counts=[];saved={};gidnum=0
    for cname,mask in masks.items():
        defs={"train":mask&training,"validate":mask&~training,"combined":mask}
        eligible=[b for b in block_ids if np.sum(mask&(block_id==b))>=v["min_windows_per_block"]]
        counts.append(dict(cohort=cname,train_windows=int(sum(defs["train"])),validation_windows=int(sum(defs["validate"])),eligible_blocks=len(eligible)))
        if min(sum(defs["train"]),sum(defs["validate"]))<v["min_windows_each_split"] or len(eligible)<v["min_eligible_blocks"]: continue
        stacks={(r,k):arr[selected].mean(axis=0) for r,arr in arrays.items() for k,selected in defs.items()}
        blocks={b:arrays["ram5"][mask&(block_id==b)].mean(axis=0) for b in eligible}
        leaves={b:arrays["ram5"][mask&(block_id!=b)].mean(axis=0) for b in eligible}
        for (r,k),stack in stacks.items(): saved[f"{cname}__{r}__{k}"]=stack
        for p,distance in enumerate(pairs.distance_m):
            peaks={k:ftan_peaks(stack[p],distance,c,f) for k,stack in stacks.items()}
            narrow={k:ftan_peaks(stacks["ram5",k][p],distance,c,f,s["narrow_width"]) for k in defs}
            bp={b:ftan_peaks(stack[p],distance,c,f) for b,stack in blocks.items()}
            lp={b:ftan_peaks(stack[p],distance,c,f) for b,stack in leaves.items()}
            for branch in ("positive","negative"):
                dev=peaks["ram5","train"][branch];tracks=track_peaks(dev,f,c)
                print("scan:",cname,p,branch,"development tracks",len(tracks),flush=True)
                for k,nodes in enumerate(dev):
                    for node in nodes: localpeaks.append(dict(cohort=cname,pair_index=p,branch=branch,frequency_hz=f[k],**node))
                for track in tracks:
                    ref=np.full(len(f),np.nan);sn=np.zeros(len(f));edge=np.ones(len(f),bool)
                    for k,j in track["path"]: node=dev[k][j];ref[k]=node["arrival"];sn[k]=node["snr"];edge[k]=node["edge"]
                    alt=next((r for r in required if r!="ram5"),"ram5")
                    ar,_,ed=peak_arrays(peaks[alt,"train"][branch],ref);nar,_,ned=peak_arrays(narrow["train"][branch],ref)
                    with np.errstate(invalid="ignore",divide="ignore"):
                        development=np.isfinite(ref)&(sn>=v["snr_min"])&~edge&~ed&~ned
                        development&=(abs(ar-ref)/ref<=v["arrival_fraction_max"])&(abs(nar-ref)/ref<=v["arrival_fraction_max"])
                    for ix in bands(development,f,s["min_points"],s["min_band_hz"]):
                        gidnum+=1;gid=f"P{gidnum:05d}"
                        limits=(max(s["arrival_min_s"],float(min(ref[ix])-s["gate_margin_s"])),min(s["arrival_max_s"],float(max(ref[ix])+s["gate_margin_s"])))
                        h={k:spectrum(stack[p],branch,limits,c,f) for k,stack in stacks.items()};hbase=h["ram5","combined"]
                        arrivals={k:peak_arrays(pk[branch],ref) for k,pk in peaks.items()};basear,_,baseed=arrivals["ram5","combined"]
                        passing=np.ones(len(f),bool);phmax=np.zeros(len(f));splitmax=np.zeros(len(f));delaymax=np.zeros(len(f));arrmax=np.zeros(len(f))
                        for route in required:
                            ph=abs(circ(np.angle(h[route,"combined"]),np.angle(hbase)))
                            split=abs(circ(np.angle(h[route,"train"]),np.angle(h[route,"validate"])))
                            delay=phase_delay(h[route,"combined"],f,ix)
                            a1,s1,e1=arrivals[route,"train"];a2,s2,e2=arrivals[route,"validate"];ac,_,ec=arrivals[route,"combined"]
                            with np.errstate(invalid="ignore",divide="ignore"):
                                da=abs(a2-a1)/a1;dc=abs(delay-ac)/ac
                                ok=(s1>=v["snr_min"])&(s2>=v["snr_min"])&~e1&~e2&~ec
                                ok&=(da<=v["arrival_fraction_max"])&(abs(a1-ref)/ref<=v["arrival_fraction_max"])
                                ok&=(abs(ac-basear)/basear<=v["arrival_fraction_max"])&(dc<=v["phase_delay_fraction_max"])
                                ok&=(ph<=v["phase_route_max_rad"])&(split<=v["phase_split_max_rad"])
                            passing&=ok;phmax=np.maximum(phmax,ph);splitmax=np.maximum(splitmax,split)
                            delaymax=np.maximum(delaymax,np.nan_to_num(dc,nan=np.inf));arrmax=np.maximum(arrmax,np.nan_to_num(da,nan=np.inf))
                            for k in ix: routes.append(dict(gate_id=gid,cohort=cname,pair_index=p,branch=branch,route=route,
                                frequency_hz=f[k],phase_change_rad=ph[k],split_phase_difference_rad=split[k],snr_train=s1[k],snr_validate=s2[k],
                                split_arrival_fraction=da[k],phase_vs_packet_delay_fraction=dc[k],route_pass=bool(ok[k])))
                        expand,move=s["gate_change_s"];controls=[(max(.01,limits[0]-expand),min(c["max_lag_s"],limits[1]+expand)),
                            (limits[0]+move,limits[1]-move),(max(.01,limits[0]-move),limits[1]-move),
                            (limits[0]+move,min(c["max_lag_s"],limits[1]+move))]
                        changes=[]
                        for limit in controls:
                            changes.append(abs(circ(np.angle(spectrum(stacks["ram5","combined"][p],branch,limit,c,f)),np.angle(hbase))) if limit[1]>limit[0] else np.full(len(f),np.inf))
                        gatechange=np.max(changes,axis=0)
                        nc,_,nce=peak_arrays(narrow["combined"][branch],ref);n1,_,n1e=peak_arrays(narrow["train"][branch],ref)
                        with np.errstate(invalid="ignore",divide="ignore"): ndiff=abs(nc-basear)/basear;ndev=abs(n1-ref)/ref
                        phasevotes=[];packetvotes=[];leavephase=[];leavearrival=[];blockh=[];blockcounts=[]
                        for b in block_ids:
                            bm=mask&(block_id==b);num=int(bm.sum());blockcounts.append(num)
                            bh=spectrum(arrays["ram5"][bm,p].mean(axis=0),branch,limits,c,f) if num else np.zeros(len(f),complex)
                            blockh.append(bh)
                            if b not in eligible: continue
                            ba,bs,be=peak_arrays(bp[b][branch],ref)
                            with np.errstate(invalid="ignore",divide="ignore"): packet=(bs>=v["block_snr_min"])&~be&(abs(ba-ref)/ref<=v["block_arrival_fraction_max"])
                            packetvotes.append(packet);phasevotes.append(packet&(abs(circ(np.angle(bh),np.angle(h["ram5","train"])))<=v["block_phase_max_rad"]))
                            lh=spectrum(leaves[b][p],branch,limits,c,f);leavephase.append(abs(circ(np.angle(lh),np.angle(hbase))))
                            la,_,le=peak_arrays(lp[b][branch],ref)
                            with np.errstate(invalid="ignore",divide="ignore"): err=abs(la-basear)/basear
                            err[le]=np.inf;leavearrival.append(err)
                        pv=np.mean(phasevotes,axis=0);bv=np.mean(packetvotes,axis=0);lv=np.max(leavephase,axis=0);la=np.max(leavearrival,axis=0)
                        intrain=np.zeros(len(f),bool);intrain[ix]=True
                        good=passing&intrain&(gatechange<=v["gate_phase_max_rad"])&~nce&~n1e
                        good&=(ndiff<=v["arrival_fraction_max"])&(ndev<=v["arrival_fraction_max"])
                        good&=(pv>=v["block_fraction_min"])&(bv>=v["block_fraction_min"])&(lv<=v["leaveout_phase_max_rad"])&(la<=v["leaveout_arrival_fraction_max"])
                        # Reject undersampled phase wraps; denser scans still do not add physical resolution.
                        sampling=2*np.pi*s["step_hz"]*ref < .9*np.pi;good&=sampling
                        gates.append(dict(gate_id=gid,cohort=cname,pair_index=p,pair=pairs.pair.iloc[p],branch=branch,
                            gate_start_s=limits[0],gate_end_s=limits[1],training_min_hz=f[ix[0]],training_max_hz=f[ix[-1]],
                            eligible_blocks=len(eligible),train_windows=int(sum(defs["train"])),validation_windows=int(sum(defs["validate"])),
                            accepted_points=int(sum(good)),gate_learned_from_training_only=True))
                        for k in ix: points.append(dict(gate_id=gid,cohort=cname,pair_index=p,branch=branch,frequency_hz=f[k],
                            arrival_train_s=ref[k],arrival_combined_s=basear[k],group_velocity_m_s=distance/basear[k],
                            required_routes_pass=bool(passing[k]),maximum_route_phase_change_rad=phmax[k],maximum_split_phase_difference_rad=splitmax[k],
                            maximum_split_arrival_fraction=arrmax[k],maximum_phase_vs_packet_delay_fraction=delaymax[k],gate_change_rad=gatechange[k],
                            narrow_arrival_fraction=ndiff[k],narrow_development_fraction=ndev[k],block_phase_fraction=pv[k],block_packet_fraction=bv[k],
                            maximum_leaveout_phase_rad=lv[k],maximum_leaveout_arrival_fraction=la[k],phase_sampling_pass=bool(sampling[k]),robust_point=bool(good[k])))
                        for number,jx in enumerate(bands(good,f,s["min_points"],s["min_band_hz"]),1):
                            sid=f"{gid}_b{number}";ph=np.unwrap(np.angle(hbase[jx]));qlo=np.full(len(jx),np.nan);qhi=qlo.copy();validdraws=0
                            if len(draws):
                                weights=np.apply_along_axis(lambda q:np.bincount(q,minlength=len(block_ids)),1,draws)*np.asarray(blockcounts)
                                totals=weights.sum(axis=1);validdraws=int(sum(totals>0))
                                hb=np.divide(weights@np.asarray(blockh)[:,jx],totals[:,None],out=np.full((len(draws),len(jx)),np.nan+1j*np.nan),where=totals[:,None]>0)
                                btheta=np.unwrap(np.angle(hb),axis=-1)
                                use=totals>0
                                btheta[use]+=2*np.pi*np.rint(np.median((ph-btheta[use])/(2*np.pi),axis=-1))[:,None]
                                if use.any(): qlo,qhi=np.quantile(btheta[use]-ph,[.025,.975],axis=0)
                            catalog.append(dict(segment_id=sid,gate_id=gid,cohort=cname,pair_index=p,pair=pairs.pair.iloc[p],branch=branch,
                                fmin_hz=f[jx[0]],fmax_hz=f[jx[-1]],points=len(jx),group_velocity_median_m_s=float(np.median(distance/basear[jx])),
                                eligible_blocks=len(eligible),gate_start_s=limits[0],gate_end_s=limits[1],absolute_cycle_resolved=False,rayleigh_mode_verified=False))
                            for k,phase,lo,hi in zip(jx,ph,qlo,qhi): measured.append(dict(segment_id=sid,gate_id=gid,cohort=cname,pair_index=p,
                                pair=pairs.pair.iloc[p],branch=branch,frequency_hz=f[k],measured_unwrapped_phase_rad=phase,
                                group_velocity_m_s=distance/basear[k],arrival_s=basear[k],relative_phase_q025_rad=lo,relative_phase_q975_rad=hi,
                                bootstrap_valid_draws=validdraws,calibrated_confidence_interval=False,instrument_response_corrected=meta["response_corrected"],
                                absolute_cycle_resolved=False,rayleigh_mode_verified=False))
        # Atomic checkpoints retain all diagnostic frequencies, including rejected ones.
        csv(out/"candidate_gates.csv",gates);csv(out/"point_diagnostics.csv",points);csv(out/"route_diagnostics.csv",routes)
        csv(out/"measured_segment_catalog.csv",catalog,CATALOG_COLUMNS);csv(out/"measured_phase_segments.csv",measured,PHASE_COLUMNS)
        print("scan: cohort complete",cname,"accepted segments",len(catalog),flush=True)
    csv(out/"candidate_gates.csv",gates,columns=None if gates else ["gate_id","cohort","pair_index","branch","gate_start_s","gate_end_s"])
    csv(out/"point_diagnostics.csv",points,columns=None if points else ["gate_id","frequency_hz","robust_point"])
    csv(out/"route_diagnostics.csv",routes,columns=None if routes else ["gate_id","frequency_hz","route","route_pass"])
    csv(out/"measured_segment_catalog.csv",catalog,CATALOG_COLUMNS);csv(out/"measured_phase_segments.csv",measured,PHASE_COLUMNS)
    csv(out/"cohort_counts.csv",counts);csv(out/"local_packet_peaks.csv",localpeaks,columns=None if localpeaks else ["cohort","pair_index","branch","frequency_hz","arrival","snr"])
    npz(out/"cohort_stacks.npz",frequency_hz=f,lags_s=(np.arange(2*round(c["max_lag_s"]*fs)+1)-round(c["max_lag_s"]*fs))/fs,**saved)
    write_json(out/"scan_complete.json",dict(config=c,point_diagnostics_sha256=sha(out/"point_diagnostics.csv"),
               candidate_gates_sha256=sha(out/"candidate_gates.csv"),segments=len(catalog),measured_points=len(measured),gates=len(gates),
               window_sets_fixed_across_frequency=True,selection_uses_structural_model=False,validation_status=v["status"],
               frequency_step_hz=s["step_hz"],frequency_step_is_not_physical_resolution=True,bootstrap_is_calibrated_confidence_interval=False))
