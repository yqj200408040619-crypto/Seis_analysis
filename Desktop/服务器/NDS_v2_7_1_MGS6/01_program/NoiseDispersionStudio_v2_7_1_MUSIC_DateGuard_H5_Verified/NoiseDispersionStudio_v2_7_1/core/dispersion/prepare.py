"""UTC alignment, local quality masks, nominal response removal and linear CCFs."""
from pathlib import Path
from fractions import Fraction
import glob, itertools, json, math
import numpy as np
import pandas as pd
from scipy import signal, ndimage
from scipy.fft import next_fast_len
from .common import same_preparation_config,csv,npz,write_json,sha,intervals,paths,invalidate_models,parse_utc

def boolean(x):
    if pd.isna(x) or x=="": return False
    if str(x).lower() in ("true","1","yes"): return True
    if str(x).lower() in ("false","0","no"): return False
    raise ValueError("Boolean field must be true/false: "+str(x))

def station_table(c):
    from pyproj import Geod,CRS,Transformer
    table=pd.read_csv(c["stations_csv"],keep_default_na=False,dtype={"station":str,"file":str,"seed_id":str})
    if not {"station","file","seed_id"}<=set(table): raise ValueError("stations.csv requires station,file,seed_id")
    if len(table)<2 or table.station.duplicated().any(): raise ValueError("At least two uniquely named stations required")
    if table.seed_id.duplicated().any(): raise ValueError("Different stations must have different seed_id values")
    if any(not s or not all(a.isalnum() or a in "_-" for a in s) for s in table.station): raise ValueError("Station labels: letters, digits, _ and - only")
    for col,val in {"response_xml":"","precorrected":False,"units":"counts","polarity":1,
                    "polarity_verified":False,"clock_offset_s":0.0}.items():
        if col not in table: table[col]=val
    for col in ("precorrected","polarity_verified"): table[col]=table[col].map(boolean)
    table["polarity"]=pd.to_numeric(table.polarity);table["clock_offset_s"]=pd.to_numeric(table.clock_offset_s)
    if not table.polarity.isin([-1,1]).all(): raise ValueError("Polarity must be declared +1 or -1")
    if table["clock_offset_s"].isna().any(): raise ValueError("Missing clock offset")
    geo={"latitude","longitude"}<=set(table) and (table.latitude!="").all() and (table.longitude!="").all()
    if geo:
        table["latitude"]=pd.to_numeric(table.latitude);table["longitude"]=pd.to_numeric(table.longitude)
        crs=CRS.from_proj4(f'+proj=aeqd +lat_0={table.latitude.mean()} +lon_0={table.longitude.mean()} +datum=WGS84 +units=m')
        tr=Transformer.from_crs(4326,crs,always_xy=True)
        table["x_m"],table["y_m"]=tr.transform(table.longitude.to_numpy(),table.latitude.to_numpy())
    else:
        if not {"x_m","y_m"}<=set(table): raise ValueError("Supply latitude/longitude or local x_m/y_m")
        table["x_m"]=pd.to_numeric(table.x_m);table["y_m"]=pd.to_numeric(table.y_m)
    pairs=[];geod=Geod(ellps="WGS84")
    for i,j in itertools.combinations(range(len(table)),2):
        a,b=table.iloc[i],table.iloc[j]
        if geo: az,_,d=geod.inv(a.longitude,a.latitude,b.longitude,b.latitude)
        else:
            dx=b.x_m-a.x_m;dy=b.y_m-a.y_m;d=float(np.hypot(dx,dy));az=float(np.degrees(np.arctan2(dx,dy)))
        if not np.isfinite(d) or d<=0: raise ValueError("Zero/invalid pair distance")
        pairs.append(dict(pair_index=len(pairs),i=i,j=j,pair=a.station+"--"+b.station,distance_m=float(d),azimuth_deg=float(az)))
    centered=table[["x_m","y_m"]].to_numpy()-table[["x_m","y_m"]].to_numpy().mean(axis=0)
    sv=np.linalg.svd(centered,compute_uv=False)
    return table,pd.DataFrame(pairs),bool(len(sv)<2 or sv[1]<0.01*max(sv[0],1e-15))

def source_files(row,c):
    found=[]
    for pattern in row["file"].split(";"):
        found+=glob.glob(str(Path(c["data_root"])/pattern))
    found=sorted(set(str(Path(p).resolve()) for p in found))
    if not found: raise FileNotFoundError("No waveform file for "+row["station"])
    return found

def read_station(row,c):
    from obspy import read,Stream,Trace,UTCDateTime
    st=Stream();files=source_files(row,c)
    for path in files:
        if Path(path).suffix.lower()==".npz":
            with np.load(path,allow_pickle=False) as z:
                data=np.asarray(z["data"],float);fs=float(z["fs_hz"]);start=parse_utc(str(z["start_utc"].item()))
                if data.ndim!=1: raise ValueError("NPZ data must be one-dimensional")
                if "valid" in z: data=np.where(z["valid"],data,np.nan)
            ids=row["seed_id"].split(".")
            if len(ids)!=4: raise ValueError("seed_id requires NET.STA.LOC.CHA, including empty LOC")
            tr=Trace(data,dict(network=ids[0],station=ids[1],location=ids[2],channel=ids[3],sampling_rate=fs,starttime=start));st+=tr
        else: st+=read(path)
    if not row["seed_id"] or any(x in row["seed_id"] for x in "*?"):
        raise ValueError("Specify one exact seed_id per station")
    st=st.select(id=row["seed_id"]);st.sort(keys=["starttime"])
    if not st: raise ValueError("No matching seed_id: "+row["seed_id"])
    rates={round(tr.stats.sampling_rate,9) for tr in st}
    if len(rates)!=1: raise ValueError("Within-station sample-rate changes require separate input epochs")
    if any(b.stats.starttime<=a.stats.endtime for a,b in zip(st,st[1:])):
        raise ValueError("Overlapping waveform samples: resolve duplicates or mask the overlap explicitly")
    for tr in st:
        tr.stats.starttime-=float(row["clock_offset_s"])
        tr.data=np.asarray(np.ma.filled(tr.data,np.nan),float)
    # Join only exactly adjacent native grids. Rounding a fractional native-grid
    # restart across a gap changes UTC and hence phase; each run is instead
    # resampled from its own actual timestamp below. Concatenate each group once.
    groups=[];current=[]
    for tr in st:
        if current:
            prev=current[-1];residual=float(tr.stats.starttime-prev.stats.endtime)*tr.stats.sampling_rate-1
            if abs(residual)>1e-5: groups.append(current);current=[]
        current.append(tr)
    if current: groups.append(current)
    joined=Stream()
    for group in groups:
        tr=group[0].copy()
        if len(group)>1: tr.data=np.concatenate([q.data for q in group])
        joined+=tr
    st=joined
    return st,files

def repair_local(raw,fs,p):
    x=np.asarray(raw,float).copy();valid=np.isfinite(x);repaired=np.zeros(len(x),bool)
    events=[dict(sample_start=int(a),sample_end_exclusive=int(b),action="missing_mask") for a,b in intervals(~valid)]
    for a,b in intervals(valid):
        y=x[a:b]
        if len(y)<8: continue
        if p["enabled"]:
            baseline=ndimage.median_filter(y,size=3,mode="nearest")
            delta=y-baseline;dy=np.diff(y)
            scale=1.4826*np.median(abs(dy-np.median(dy)))/np.sqrt(2)
            floor=np.finfo(float).eps*max(float(np.max(abs(y))),1.0)*32
            med=np.median(y);amp=1.4826*np.median(abs(y-med))
            # Optional local robust scales follow changing ambient amplitudes.
            # The detector changes; waveform normalization/validation do not.
            local=max(0.,float(p.get("scale_window_s",0)))
            if local:
                block=max(32,round(local*fs));scales=[];amps=[];medians=[];centers=[]
                for lo in range(0,len(y),block):
                    hi=min(len(y),lo+block);q=y[lo:hi];dd=np.diff(q)
                    qm=float(np.median(q));qa=1.4826*float(np.median(abs(q-qm)))
                    qs=1.4826*float(np.median(abs(dd-np.median(dd))))/np.sqrt(2) if len(dd) else 0.
                    scales.append(max(qs,.05*scale,floor));amps.append(max(qa,.05*amp,floor));medians.append(qm);centers.append((lo+hi-1)/2)
                index=np.arange(len(y));scale=np.interp(index,centers,scales)
                amp=np.interp(index,centers,amps);med=np.interp(index,centers,medians)
            suspect=abs(delta)>p["threshold_mad"]*np.maximum(scale,floor)
            # A long elevated sample plateau may have small three-point residuals.
            # Detect extreme excursions as well, but only repair a bounded brief run.
            suspect|=abs(y-med)>p["threshold_mad"]*np.maximum(amp,floor)*2
            for i,j in intervals(suspect):
                left=a+i;right=a+j
                if (j-i)/fs<=p["max_repair_s"]+1e-12 and i>0 and j<len(y):
                    endpoints=(x[left-1],x[right])
                    endpoint_scale=float(np.max(scale[max(0,i-1):min(len(y),j+1)])) if np.ndim(scale) else scale
                    if abs(endpoints[1]-endpoints[0])<=8*max(endpoint_scale,floor):
                        x[left:right]=np.linspace(*endpoints,j-i+2)[1:-1]
                        repaired[left:right]=True;action="short_local_interpolation"
                    else: valid[left:right]=False;action="local_mask"
                else: valid[left:right]=False;action="local_mask"
                events.append(dict(sample_start=int(left),sample_end_exclusive=int(right),action=action))
        # Long exact flatlines/rails are unusable, not interpolated.
        for i,j in intervals(np.diff(y)==0):
            if (j-i+1)/fs>=p["flatline_min_s"]:
                valid[a+i:a+j+1]=False
                events.append(dict(sample_start=int(a+i),sample_end_exclusive=int(a+j+1),action="flatline_mask"))
    bad=~valid
    if p["mask_pad_s"] and bad.any():
        radius=int(np.ceil(p["mask_pad_s"]*fs))
        bad=ndimage.maximum_filter1d(bad,size=2*radius+1,mode="constant",cval=1)
    valid&=~bad
    return x,valid,repaired,events

def fractional_samples(y,fraction):
    """Constant fractional delay on a uniform grid, bandlimited windowed sinc."""
    if abs(fraction)<1e-9: return np.asarray(y)
    j=np.arange(-32,33);h=np.sinc(j+fraction)*np.kaiser(65,8.6);h/=h.sum()
    return signal.fftconvolve(y,h,mode="same")

def resample_align(x,fs_in,fs_out,start,origin,count):
    ratio=Fraction(fs_out/fs_in).limit_denominator(100000)
    if abs(float(ratio)-fs_out/fs_in)>1e-10: raise ValueError("Cannot represent resampling ratio accurately")
    y=signal.resample_poly(x,ratio.numerator,ratio.denominator,padtype="line")
    offset=float(start-origin)*fs_out
    first=max(0,int(np.ceil(offset-1e-8)))
    u=first-offset;integer=int(np.floor(u+1e-8));frac=u-integer
    if frac<0 and frac>-1e-7: frac=0.0
    z=fractional_samples(y,frac)
    end=min(count,first+len(z)-integer)
    return first,end,z[integer:integer+max(0,end-first)]

def cosine_band(freq,edges):
    a,b,c,d=edges;w=np.zeros_like(freq);w[(freq>=b)&(freq<=c)]=1
    ix=(freq>a)&(freq<b);w[ix]=.5-.5*np.cos(np.pi*(freq[ix]-a)/(b-a))
    ix=(freq>c)&(freq<d);w[ix]=.5+.5*np.cos(np.pi*(freq[ix]-c)/(d-c))
    return w

def response_table(row,st,c):
    from obspy import read_inventory
    if row["precorrected"]:
        expected={"ACC":"m/s^2","VEL":"m/s","DISP":"m"}[c["response"]["output"]]
        if row["units"]!=expected: raise ValueError("Precorrected units must match response.output: "+expected)
        return None,dict(response_status="physical_input_user_declared",nominal_dip_deg=None,nominal_azimuth_deg=None)
    if not row["response_xml"]:
        if c["response"]["required"]: raise ValueError("Missing response_xml for "+row["station"])
        return None,dict(response_status="uncorrected_provisional",nominal_dip_deg=None,nominal_azimuth_deg=None)
    inv=read_inventory(str(Path(c["stations_csv"]).parent/row["response_xml"]))
    net,sta,loc,cha=row["seed_id"].split(".")
    candidates=[ch for n in inv if n.code==net for s in n if s.code==sta for ch in s if ch.code==cha and ch.location_code==loc]
    start=st[0].stats.starttime;end=st[-1].stats.endtime
    candidates=[ch for ch in candidates if (ch.start_date is None or ch.start_date<=start) and (ch.end_date is None or ch.end_date>=end)]
    if len(candidates)!=1: raise ValueError("A unique response epoch must cover the full input; split inputs at response changes: "+row["station"])
    ch=candidates[0]
    if abs(ch.sample_rate-st[0].stats.sampling_rate)>1e-6: raise ValueError("Metadata/native sample-rate mismatch")
    step=c["response"]["table_step_hz"];ny=c["sample_rate_hz"]/2
    freq=np.arange(step,ny+step/2,step)
    h=ch.response.get_evalresp_response_for_frequencies(freq,output=c["response"]["output"])
    if not np.isfinite(h).all() or (abs(h)==0).any(): raise ValueError("Zero/nonfinite instrument response")
    test=np.linspace(c["bandpass_hz"][0],c["bandpass_hz"][1],151)+step*.173
    exact=ch.response.get_evalresp_response_for_frequencies(test,output=c["response"]["output"])
    approx=np.exp(np.interp(test,freq,np.log(abs(h)))+1j*np.interp(test,freq,np.unwrap(np.angle(h))))
    err=float(np.max(abs(np.angle(approx/exact))))
    if err>c["response"]["interpolation_phase_tolerance_rad"]: raise ValueError("Response interpolation too coarse; reduce table_step_hz")
    gains=[s.stage_gain for s in ch.response.response_stages]
    return (freq,h),dict(response_status="nominal_response_corrected",nominal_dip_deg=ch.dip,nominal_azimuth_deg=ch.azimuth,
                        sensitivity=float(ch.response.instrument_sensitivity.value),input_units=ch.response.instrument_sensitivity.input_units,
                        all_stage_gains_positive=bool(all(g is not None and g>0 for g in gains)),
                        response_interpolation_phase_error_rad=err,epoch_start=str(ch.start_date),epoch_end=str(ch.end_date))

def correct_response(x,table,c):
    if table is None: return x
    fs=c["sample_rate_hz"];pad=min(len(x)-1,int(60*fs))
    z=np.pad(signal.detrend(x,type="linear"),(pad,pad),mode="reflect");nf=next_fast_len(len(z))
    f=np.fft.rfftfreq(nf,1/fs);ff,h=table
    hh=np.exp(np.interp(f,ff,np.log(abs(h)))+1j*np.interp(f,ff,np.unwrap(np.angle(h))))
    result=np.fft.irfft(np.fft.rfft(z,nf)*cosine_band(f,c["response"]["pre_filt_hz"])/hh,nf)
    return result[pad:pad+len(x)]

def normalized_fft(x,c,route):
    fs=c["sample_rate_hz"];n=x.shape[-1];nf=next_fast_len(2*n-1)
    z=signal.detrend(x,axis=-1,type="linear")
    if "notch_hz" in route:
        b,a=signal.iirnotch(route["notch_hz"],route.get("notch_q",30),fs);z=signal.filtfilt(b,a,z,axis=-1)
    if route.get("normalization","ram")=="onebit": z=np.sign(z)
    else:
        scale=ndimage.uniform_filter1d(abs(z),max(1,int(route.get("ram_s",5)*fs)),axis=-1,mode="reflect")
        floor=np.maximum(np.median(scale,axis=-1,keepdims=True)*.05,np.finfo(float).tiny);z/=np.maximum(scale,floor)
    z*=signal.windows.tukey(n,.05)
    ft=np.fft.rfft(z,nf,axis=-1);f=np.fft.rfftfreq(nf,1/fs)
    smooth=ndimage.uniform_filter1d(abs(ft),max(3,int(route["whiten_hz"]/(fs/nf))),axis=-1,mode="nearest")
    lo,hi=c["bandpass_hz"];core=(f>=max(lo,.2))&(f<=min(hi,5))
    if not core.any(): core=(f>=lo)&(f<=hi)
    floor=np.maximum(np.median(smooth[:,core],axis=-1,keepdims=True)*.01,1e-12)
    edge_lo=min(.05,(hi-lo)/10);edge_hi=min(.3,(hi-lo)/10)
    taper=cosine_band(f,[lo,lo+edge_lo,hi-edge_hi,hi])
    return ft/np.maximum(smooth,floor)*taper,f,nf

def cc_from_fft(ft,nf,n,c,pairs):
    m=round(c["max_lag_s"]*c["sample_rate_hz"]);result=[]
    for r in pairs.itertuples():
        y=np.fft.irfft(np.conj(ft[r.i])*ft[r.j],nf)
        result.append(np.r_[y[-m:],y[:m+1]]/n)
    return np.asarray(result,np.float32)

def window_segments(data,valid,start,c):
    n=round(c["window_s"]*c["sample_rate_hz"]);minimum=round(c["min_clean_segment_s"]*c["sample_rate_hz"])
    joint=np.all(valid[:,start:start+n],axis=0)
    return [(start+i,start+j) for i,j in intervals(joint) if j-i>=minimum]

def prepare(c):
    from obspy import UTCDateTime
    out,cache=paths(c);table,pairs,line=station_table(c);records=[];sources=[]
    headers=[]
    for row in table.to_dict("records"):
        st,files=read_station(row,c);headers.append((st[0].stats.starttime,st[-1].stats.endtime+1/st[-1].stats.sampling_rate))
        for file in files: sources.append(dict(path=file,sha256=sha(file)))
        if row["response_xml"]:
            p=(Path(c["stations_csv"]).parent/row["response_xml"]).resolve();sources.append(dict(path=str(p),sha256=sha(p)))
        del st
    sources.append(dict(path=c["stations_csv"],sha256=sha(c["stations_csv"])))
    common_start=max(t[0] for t in headers);common_end=min(t[1] for t in headers)
    if c["start_utc"]: common_start=max(common_start,parse_utc(c["start_utc"]))
    if c["end_utc"]: common_end=min(common_end,parse_utc(c["end_utc"]))
    # One common trim after alignment, never 600 native seconds per station independently.
    origin=common_start+c["discard_first_s"];end=common_end
    fs=c["sample_rate_hz"];count=int(np.floor(float(end-origin)*fs+1e-7))
    if count<round(c["window_s"]*fs): raise ValueError("Insufficient common duration after initial discard")
    fingerprint=dict(config=c,sources=sources,code_sha256={p.name:sha(p) for p in sorted(Path(__file__).parent.glob("*.py"))})
    provenance=cache/"provenance.json"
    if provenance.exists():
        old=json.loads(provenance.read_text(encoding='utf-8-sig'))
        if same_preparation_config(old.get("config",{}),c) and all(old.get(k)==fingerprint[k] for k in ("sources","code_sha256")) and (out/"prepare_complete.json").exists():
            print("prepare: verified input/config/code cache reused",flush=True);return
    # Invalidate completion markers before producing a new generation.
    invalidate_models(out)
    for p in (out/"prepare_complete.json",out/"scan_complete.json",out/"controls_complete.json",out/"families_complete.json",out/"verification.json"): p.unlink(missing_ok=True)
    data=np.lib.format.open_memmap(cache/"aligned_filtered.npy",mode="w+",dtype="float32",shape=(len(table),count));data[:]=0
    valid=np.lib.format.open_memmap(cache/"aligned_valid.npy",mode="w+",dtype="bool",shape=data.shape);valid[:]=False
    eventrows=[];guard=max(c["filter_guard_s"],3/c["bandpass_hz"][0]);sos=signal.butter(4,c["bandpass_hz"],btype="bandpass",fs=fs,output="sos")
    for k,row in enumerate(table.to_dict("records")):
        st,_=read_station(row,c);nativefs=st[0].stats.sampling_rate
        if nativefs<fs: raise ValueError("Output sampling rate must not exceed any native rate")
        if nativefs/2<=c["response"]["pre_filt_hz"][-1]: raise ValueError("Native Nyquist below response pre-filter")
        resp,audit=response_table(row,st,c);repcount=badcount=0
        for tr in st:
            x,ok,repaired,events=repair_local(tr.data,nativefs,c["spikes"]);repcount+=int(repaired.sum());badcount+=int((~ok).sum())
            for event in events:
                eventrows.append(dict(station=row["station"],start_utc=str(tr.stats.starttime+event["sample_start"]/nativefs),
                                      end_exclusive_utc=str(tr.stats.starttime+event["sample_end_exclusive"]/nativefs),**event))
            for a,b in intervals(ok):
                if (b-a)/nativefs<2*guard+c["min_clean_segment_s"]: continue
                t0=tr.stats.starttime+a/nativefs
                lo,hi,y=resample_align(x[a:b],nativefs,fs,t0,origin,count)
                if hi<=lo: continue
                # Use full run as filter/response padding; crop only afterward.
                ratio=Fraction(fs/nativefs).limit_denominator(100000)
                full=signal.resample_poly(x[a:b],ratio.numerator,ratio.denominator,padtype="line")
                full=signal.sosfiltfilt(sos,correct_response(full,resp,c))*int(row["polarity"])
                off=float(t0-origin)*fs;u=lo-off;integer=int(np.floor(u+1e-8));fraction=u-integer
                if fraction<0 and fraction>-1e-7: fraction=0.
                aligned=fractional_samples(full,fraction)[integer:integer+hi-lo]
                goodlo=max(lo,int(np.ceil(float(t0+guard-origin)*fs)))
                goodhi=min(hi,int(np.floor(float(t0+(b-a)/nativefs-guard-origin)*fs)))
                if goodhi>goodlo:
                    data[k,goodlo:goodhi]=aligned[goodlo-lo:goodhi-lo]
                    valid[k,goodlo:goodhi]=True
        records.append(dict(station=row["station"],seed_id=row["seed_id"],native_sample_rate_hz=float(nativefs),
                            reported_start_utc=str(st[0].stats.starttime+float(row["clock_offset_s"])),
                            clock_corrected_start_utc=str(st[0].stats.starttime),clock_corrected_end_utc=str(st[-1].stats.endtime),
                            native_samples=int(sum(tr.stats.npts for tr in st)),
                            common_origin_native_offset_samples=float(origin-st[0].stats.starttime)*nativefs,
                            common_origin_utc=str(origin),initial_discard_s=c["discard_first_s"],
                            declared_clock_offset_s=float(row["clock_offset_s"]),declared_polarity=int(row["polarity"]),
                            polarity_verified=bool(row["polarity_verified"]),repaired_native_samples=repcount,
                            masked_native_samples=badcount,usable_aligned_samples=int(valid[k].sum()),**audit))
        print("prepare:",row["station"],"repaired",repcount,"masked",badcount,flush=True)
    data.flush();valid.flush()
    # These station-level diagnostics must survive a later window rejection.
    # Otherwise a real-data failure leaves no evidence identifying its cause.
    csv(out/"record_audit.csv",records)
    csv(out/"local_quality_events.csv",eventrows,columns=["station","start_utc","end_exclusive_utc","sample_start","sample_end_exclusive","action"])
    csv(out/"stations.csv",table);csv(out/"pairs.csv",pairs)
    # UTC-anchored non-overlapping windows; a calendar validation block is never split.
    w=c["window_s"];block=c["validation"]["block_s"];timestamps=[];starts=[]
    first=math.ceil(float(origin.timestamp)/w)*w
    for ts in np.arange(first,float(end.timestamp)-w+1e-5,w):
        if math.floor(ts/block)!=math.floor((ts+w-1e-6)/block): continue
        start=round(float(UTCDateTime(float(ts))-origin)*fs)
        if 0<=start and start+round(w*fs)<=count: timestamps.append(float(ts));starts.append(start)
    m=round(c["max_lag_s"]*fs);shape=(len(starts),len(pairs),2*m+1)
    routes={r["name"]:np.zeros(shape,np.float32) for r in c["routes"]};quality=[];accepted=[]
    for win,(start,ts) in enumerate(zip(starts,timestamps)):
        segments=window_segments(data,valid,start,c);usable=sum(b-a for a,b in segments);frac=usable/(w*fs)
        good=frac>=c["min_clean_fraction"] and usable>0
        if good:
            for a,b in segments:
                x=np.asarray(data[:,a:b],float)
                if np.min(np.std(x,axis=-1))<1e-15: good=False;break
                for r in c["routes"]:
                    ft,_,nf=normalized_fft(x,c,r);routes[r["name"]][win]+=cc_from_fft(ft,nf,b-a,c,pairs)*(b-a)/usable
        if not good:
            for arr in routes.values(): arr[win]=0
        quality.append(dict(window=win,start_utc=str(UTCDateTime(ts)),timestamp=ts,start_sample=start,
                            clean_fraction=float(frac),clean_segments=len(segments),effective_clean_seconds=usable/fs,accepted=bool(good)))
        accepted.append(good)
        if win%30==0: print("prepare: CCF windows",win+1,"/",len(starts),flush=True)
    quality=csv(out/"window_quality.csv",quality)
    selected=np.asarray(accepted,bool)
    if selected.sum()<c["minimum_windows"]: raise ValueError("Too few usable windows; inspect window_quality.csv before changing thresholds")
    times=np.asarray(timestamps)[selected];split=c["validation"]["split_utc"]
    if split is None:
        blocks=np.unique(np.floor(times/block).astype(np.int64));ntrain=max(1,min(len(blocks)-1,int(np.floor(len(blocks)*c["validation"]["train_fraction"]))))
        if len(blocks)<2: raise ValueError("At least two non-overlapping temporal blocks required")
        split=float(blocks[ntrain]*block)
    else: split=float(UTCDateTime(split).timestamp)
    train=times<split
    if min(train.sum(),(~train).sum())<c["validation"]["min_windows_each_split"]: raise ValueError("Too few development or validation windows")
    if np.intersect1d(np.floor(times[train]/block),np.floor(times[~train]/block)).size: raise ValueError("split_utc must fall on a validation block boundary")
    schedule=quality[selected].copy();schedule["training"]=train;schedule["block_id"]=np.floor(times/block).astype(np.int64)
    csv(out/"window_schedule.csv",schedule);csv(out/"stations.csv",table);csv(out/"pairs.csv",pairs);csv(out/"record_audit.csv",records)
    csv(out/"local_quality_events.csv",eventrows,columns=["station","start_utc","end_exclusive_utc","sample_start","sample_end_exclusive","action"])
    npz(cache/"window_routes.npz",timestamps=times,training=train,starts=np.asarray(starts)[selected],**{k:v[selected] for k,v in routes.items()})
    fingerprint.update(origin_utc=str(origin),end_exclusive_utc=str(end),split_utc=str(UTCDateTime(split)),
                       line_array=line,response_corrected=all(x["response_status"]!="uncorrected_provisional" for x in records),
                       polarity_verified=all(x["polarity_verified"] for x in records))
    write_json(provenance,fingerprint)
    write_json(out/"prepare_complete.json",dict(scheduled_windows=len(starts),accepted_windows=int(selected.sum()),
               discarded_common_first_s=c["discard_first_s"],filter_guard_s=guard,common_utc_aligned=True))
