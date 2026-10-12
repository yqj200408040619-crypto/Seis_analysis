from pathlib import Path
import hashlib, json, os, tempfile,re
import numpy as np
import pandas as pd

DEFAULTS = {
    "stations_csv": "stations.csv", "data_root": ".", "output": "output",
    "start_utc": None, "end_utc": None, "sample_rate_hz": 20.0,
    "discard_first_s": 600.0, "filter_guard_s": 30.0,
    "bandpass_hz": [0.1, 7.142857142857], "window_s": 600.0,
    "max_lag_s": 60.0, "min_clean_fraction": 0.7,
    "min_clean_segment_s": 150.0, "minimum_windows": 24,
    "response": {"required": True, "output": "ACC", "pre_filt_hz": [0.03, 0.075, 7.5, 9.0],
                 "table_step_hz": 0.005, "interpolation_phase_tolerance_rad": 0.001},
    "spikes": {"enabled": True, "threshold_mad": 18.0, "scale_window_s": 0.0, "max_repair_s": 0.04,
               "mask_pad_s": 0.1, "flatline_min_s": 1.0},
    "routes": [
        {"name": "ram5", "ram_s": 5.0, "whiten_hz": 0.05, "required": True},
        {"name": "ram20", "ram_s": 20.0, "whiten_hz": 0.05, "required": True},
        {"name": "white025", "ram_s": 5.0, "whiten_hz": 0.025, "required": True},
        {"name": "white100", "ram_s": 5.0, "whiten_hz": 0.1, "required": True},
        {"name": "onebit", "normalization": "onebit", "whiten_hz": 0.05, "required": False}],
    "scan": {"min_hz": 1.5, "max_hz": 5.0, "step_hz": 0.025,
             "ftan_width": 0.15, "narrow_width": 0.08, "minimum_sigma_hz": 0.04,
             "group_velocity_bounds_m_s": [200.0, 5000.0], "arrival_min_s": 0.2,
             "arrival_max_s": 25.0, "noise_lag_s": [35.0, 55.0],
             "min_band_hz": 0.25, "min_points": 8, "max_peaks": 8, "max_tracks": 12,
             "gate_margin_s": 0.9, "gate_change_s": [0.3, 0.2]},
    "cohorts": {"clock_groups": True, "night_hours_utc": [18, 6], "states": 2,
                "pca_components": 6, "feature_band_hz": [2.65, 4.8], "feature_max_lag_s": 10.0},
    "validation": {"block_s": 21600.0, "train_fraction": 0.5, "split_utc": None,
                   "status": "internal_chronological", "min_windows_each_split": 12,
                   "min_windows_per_block": 6, "min_eligible_blocks": 4,
                   "snr_min": 10.0, "phase_route_max_rad": 0.2,
                   "phase_split_max_rad": 0.35, "arrival_fraction_max": 0.1,
                   "phase_delay_fraction_max": 0.2, "gate_phase_max_rad": 0.35,
                   "block_snr_min": 5.0, "block_phase_max_rad": 0.5,
                   "block_arrival_fraction_max": 0.15, "block_fraction_min": 0.75,
                   "leaveout_phase_max_rad": 0.2, "leaveout_arrival_fraction_max": 0.1,
                   "bootstrap_draws": 499},
    "controls": {"trials": 40, "seed": 9074},
    "phase": {"source_phase_rad": [0.7853981633974483], "source_phase_verified": False,
              "cycle_min": 0, "cycle_max": 80, "velocity_bounds_m_s": [150.0, 3000.0],
              "link_median_fraction": 0.03, "link_max_fraction": 0.05,
              "link_group_median_fraction": 0.15, "max_direction_difference_deg":45.0},
    "seed": 9072
}

def merge(default, custom, context="config"):
    if not isinstance(custom, dict):
        raise ValueError(f"{context}: expected object")
    extra = set(custom) - set(default)
    if extra: raise ValueError(f"{context}: unknown keys {sorted(extra)}")
    out = json.loads(json.dumps(default))
    for k,v in custom.items():
        out[k] = merge(default[k], v, context+"."+k) if isinstance(default[k],dict) else v
    return out

def load_config(path):
    path = Path(path).resolve()
    c = merge(DEFAULTS, json.loads(path.read_text(encoding="utf-8-sig")))
    for k in ("stations_csv", "data_root", "output"):
        c[k] = str((path.parent / c[k]).resolve())
    validate_config(c)
    return c

def validate_config(c):
    fs=c["sample_rate_hz"]; s=c["scan"]; v=c["validation"]; r=c["response"]
    if not 0 < c["bandpass_hz"][0] < s["min_hz"] < s["max_hz"] < c["bandpass_hz"][1] < fs/2:
        raise ValueError("Require 0 < bandpass low < scan low < scan high < bandpass high < output Nyquist")
    pf=r["pre_filt_hz"]
    if len(pf)!=4 or not 0<pf[0]<pf[1]<=c["bandpass_hz"][0]<c["bandpass_hz"][1]<=pf[2]<pf[3]<fs/2:
        raise ValueError("Response pre-filter must enclose bandpass and stay below output Nyquist")
    if not 0<s["arrival_min_s"]<s["arrival_max_s"]<s["noise_lag_s"][0]<s["noise_lag_s"][1]<c["max_lag_s"]:
        raise ValueError("Arrival search, late-noise interval and max lag must be separated")
    if c["min_clean_segment_s"]<=c["max_lag_s"] or c["min_clean_segment_s"]>c["window_s"]:
        raise ValueError("Clean segment must be longer than max lag and no longer than a window")
    if c["window_s"]*fs % 1 or c["max_lag_s"]*fs % 1: raise ValueError("Window/lag must have integer samples")
    if c["discard_first_s"]<0 or c["filter_guard_s"]<0 or not 0<c["min_clean_fraction"]<=1:
        raise ValueError("Invalid trimming/quality settings")
    if not 0<v["train_fraction"]<1 or v["block_s"]<c["window_s"] or v["min_eligible_blocks"]<2:
        raise ValueError("Invalid validation block/split settings")
    if not 0<v["block_fraction_min"]<=1 or v["bootstrap_draws"]<0: raise ValueError("Invalid validation settings")
    if v["status"] not in ("internal_chronological", "internal_reuse", "external_user_declared"):
        raise ValueError("Unknown validation status")
    names=[x["name"] for x in c["routes"]]
    if not names or names[0]!="ram5" or len(names)!=len(set(names)): raise ValueError("Unique routes; first baseline must be ram5")
    if not c["routes"][0].get("required",False): raise ValueError("Baseline route is required")
    for route in c["routes"]:
        if set(route)-{"name","ram_s","whiten_hz","required","normalization","notch_hz","notch_q"}: raise ValueError("Unknown route key")
        if route.get("whiten_hz",0)<=0 or route.get("normalization","ram") not in ("ram","onebit"): raise ValueError("Invalid route")
        if route.get("ram_s",5)<=0: raise ValueError("RAM duration must be positive")
        if "notch_hz" in route and not 0<route["notch_hz"]<fs/2: raise ValueError("Invalid notch frequency")
    if not isinstance(s["min_points"],int) or s["min_points"]<3 or s["step_hz"]<=0 or s["min_band_hz"]<=0 or s["max_tracks"]<1:
        raise ValueError("Invalid frequency/track settings")
    if not 0<=c["phase"]["cycle_min"]<=c["phase"]["cycle_max"]: raise ValueError("Invalid cycle range")
    if any(not isinstance(c["phase"][k],int) for k in ("cycle_min","cycle_max")): raise ValueError("Cycle limits must be integers")
    if not c["phase"]["source_phase_rad"] or not np.isfinite(c["phase"]["source_phase_rad"]).all(): raise ValueError("Finite source phase hypotheses required")
    if not 0<c["phase"]["velocity_bounds_m_s"][0]<c["phase"]["velocity_bounds_m_s"][1]: raise ValueError("Invalid phase velocity bounds")
    if not 0<s["group_velocity_bounds_m_s"][0]<s["group_velocity_bounds_m_s"][1]: raise ValueError("Invalid group velocity bounds")
    if c["cohorts"]["states"]<0 or c["controls"]["trials"]<0: raise ValueError("Negative count")
    if c["cohorts"]["states"]:
        lo,hi=c["cohorts"]["feature_band_hz"]
        if not c["bandpass_hz"][0]<=lo<hi<=c["bandpass_hz"][1]: raise ValueError("State feature band outside preprocessing band")

def sha(path):
    with Path(path).open("rb") as h: return hashlib.file_digest(h,"sha256").hexdigest()

def parse_utc(value):
    """Preserve ISO nanoseconds instead of datetime's microsecond rounding."""
    from obspy import UTCDateTime
    if isinstance(value,UTCDateTime):return UTCDateTime(ns=value.ns,precision=9)
    if isinstance(value,str):
        match=re.fullmatch(r'(.*[T ]\d{2}:\d{2}:\d{2})[.,](\d+)(.*)',value)
        if match:
            if len(match[2])>9:raise ValueError('UTC supports at most 9 fractional digits')
            base=UTCDateTime(match[1]+match[3])
            return UTCDateTime(ns=base.ns+int(match[2].ljust(9,'0')),precision=9)
    return UTCDateTime(value,precision=9)

def write_json(path, data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+".part")
    tmp.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False),encoding="utf-8")
    os.replace(tmp,path)

def csv(path, records, columns=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    df=records if isinstance(records,pd.DataFrame) else pd.DataFrame(records,columns=columns)
    df=df.replace([np.inf,-np.inf],np.nan)
    tmp=path.with_name(path.name+".part");df.to_csv(tmp,index=False);os.replace(tmp,path)
    return df

def npz(path, **values):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.stem+".part.npz");np.savez_compressed(tmp,**values);os.replace(tmp,path)

def intervals(mask):
    d=np.diff(np.r_[False,np.asarray(mask,bool),False].astype(np.int8))
    return list(zip(np.flatnonzero(d==1),np.flatnonzero(d==-1)))

def circ(a,b): return np.angle(np.exp(1j*(a-b)))

def bands(mask,f,min_points,min_width):
    return [np.arange(i,j) for i,j in intervals(mask) if j-i>=min_points and f[j-1]-f[i]>=min_width-1e-8]

def paths(c):
    out=Path(c["output"]);out.mkdir(parents=True,exist_ok=True)
    cache=out/"cache";cache.mkdir(exist_ok=True)
    return out,cache

def invalidate_models(out):
    # Generated downstream files cannot survive a new measurement generation.
    for name in ("inversion_complete.json","conditional_models.json","selected_observed_dispersion.csv",
                 "model_predictions.csv","inversion_attempts.csv"):
        (Path(out)/name).unlink(missing_ok=True)
    (Path(out)/"figures/05_optional_conditional_models.png").unlink(missing_ok=True)

def same_preparation_config(a,b):
    """Only contiguous-band gates may change without recomputing waveforms."""
    left=json.loads(json.dumps(a));right=json.loads(json.dumps(b))
    for item in (left,right):
        for key in ('min_band_hz','min_points'):item.get('scan',{}).pop(key,None)
    return left==right


def ensure_ready(c, stage):
    out,cache=paths(c)
    meta=json.loads((cache/"provenance.json").read_text(encoding='utf-8-sig'))
    if not same_preparation_config(meta["config"],c):
        raise ValueError("预处理参数已改变；请先重新预处理，再扫描。Config changed; rerun prepare")
    here=Path(__file__).parent
    now={p.name:sha(p) for p in sorted(here.glob("*.py"))}
    if meta["code_sha256"]!=now: raise ValueError("程序版本已改变；请先重新预处理。Program changed; rerun prepare")
    for x in meta["sources"]:
        if sha(x["path"])!=x["sha256"]: raise ValueError("Source changed: "+x["path"])
    if not (out/stage).exists(): raise ValueError("Missing previous stage: "+stage)
    if stage!="prepare_complete.json":
        marker=out/"scan_complete.json"
        scanned=json.loads(marker.read_text(encoding='utf-8-sig')).get('config',meta['config']) if marker.exists() else meta['config']
        if scanned!=c:raise ValueError("扫描参数已改变；请先重新扫描，再运行后续阶段。Scan settings changed; rerun scan")
    return out,cache,meta
