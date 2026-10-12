#!/usr/bin/env python3
"""Follow-up diagnostics: overlapping MUSIC gates and synchronized low-band windows."""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal.windows import hann

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_validation import (BASE, FS, OUT, baseline_settings, estimate_one, load_actual,
                            local_peaks, session_for, simulate_independent_band_waves,
                            simulate_series, station_geometry, valid_windows)

from ndstudio.direction import REASONS, angular_distance, estimate, pair_gate
from ndstudio.phase_records import full_window_starts, phase_window_geometry

LOW_FREQUENCIES = np.round(np.arange(.08, .131, .01), 2)
PAIR = "BW_MGS05--BW_MGSBH"
SEED = 20261103
REPS = 999


def gate_flags(result, c):
    """Return independent gate flags; q!=1 residual is explicitly unmeasured."""
    q = np.asarray(result["music_source_count"], int)
    enough = np.asarray(result["station_count"], int) >= c["min_stations"]
    geometry = np.asarray(result["geometry_ratio"], float)
    geometry_ok = geometry >= c["min_geometry_ratio"]
    estimable = enough & geometry_ok
    q0 = estimable & (q == 0)
    q1 = estimable & (q == 1)
    q2 = estimable & (q >= 2)
    coherence = np.asarray(result["coherence"], float)
    rms = np.asarray(result["phase_rms_deg"], float)
    residual_evaluated = q1 & np.isfinite(rms)
    flags = {
        "few_stations": ~enough,
        "geometry": enough & ~geometry_ok,
        "no_signal": q0,
        "multiple_sources": q2,
        "low_coherence": np.isfinite(coherence) & (coherence < c["min_coherence"]),
        "phase_residual_fail": residual_evaluated & (rms > c["phase_rms_max_deg"]),
        "uncertain": np.isfinite(result["direction_uncertainty_deg"]) &
                     (np.asarray(result["direction_uncertainty_deg"]) > c["max_direction_uncertainty_deg"]),
        "ambiguous": np.isfinite(result["secondary_ratio"]) &
                     (np.asarray(result["secondary_ratio"]) >= c["max_secondary_peak_ratio"]),
        "boundary": np.isfinite(result["velocity_m_s"]) &
                    ((np.asarray(result["velocity_m_s"]) <= c["velocity_min_m_s"]*1.001) |
                     (np.asarray(result["velocity_m_s"]) >= c["velocity_max_m_s"]*.999)),
        "grid_resolution": np.isfinite(result["grid_half_cell_phase_rad"]) &
                           (np.asarray(result["grid_half_cell_phase_rad"]) > np.pi/2),
    }
    metric_names = ["low_coherence", "phase_residual_fail", "uncertain", "ambiguous", "boundary", "grid_resolution"]
    metric_failures = np.sum(np.column_stack([flags[k] for k in metric_names]), axis=1)
    return dict(flags=flags, q=q, estimable=estimable, q0=q0, q1=q1, q2=q2,
                residual_evaluated=residual_evaluated, metric_failure_count=metric_failures,
                q1_all_gates_pass=q1 & (metric_failures == 0))


def run_front_rejection_audit():
    arrays, _, names, xy, fs, _ = load_actual()
    c = baseline_settings()
    session = session_for(arrays, names, xy, fs, c)
    freq_file = BASE / "music_front/frequency_directions.csv"
    frequencies = np.sort(pd.read_csv(freq_file).frequency_hz.unique())
    rows, summary_rows, check_rows = [], [], []
    stored = pd.read_csv(BASE / "music_front/direction/frequency_diagnostics.csv")
    for fi, f in enumerate(frequencies, 1):
        geo = phase_window_geometry(float(f), fs, c["music_front_cycles"], c["music_front_cycles"])
        n, step = geo["window_samples"], geo["slide_samples"]
        starts = full_window_starts(len(next(iter(arrays.values()))), n, step, include_last=True)
        coeff = session.transform(starts, n, hann(n, sym=True), [float(f)], snapshots=True)
        result, _ = estimate(coeff, xy, np.array([f]), c)
        audit = gate_flags(result, c)
        q, flags = audit["q"], audit["flags"]
        reason = np.asarray(result["reason"], str)
        each = pd.DataFrame({
            "frequency_hz": float(f), "start_sample": starts,
            "station_count": result["station_count"], "geometry_ratio": result["geometry_ratio"],
            "source_count": np.where(audit["estimable"], q, -1),
            "source_count_evaluated": audit["estimable"],
            "reported_phase_rms_deg": result["phase_rms_deg"],
            "phase_rms_evaluated": audit["residual_evaluated"],
            "phase_rms_analysis_deg": np.where(audit["residual_evaluated"], result["phase_rms_deg"], np.nan),
            "coherence": result["coherence"], "azimuth_deg": result["azimuth_deg"],
            "uncertainty_deg": result["direction_uncertainty_deg"],
            "secondary_ratio": result["secondary_ratio"], "velocity_m_s": result["velocity_m_s"],
            "grid_half_cell_phase_rad": result["grid_half_cell_phase_rad"],
            "final_reason": reason, "accepted": result["quality"]})
        for key, val in flags.items():
            each[key + "_fail"] = val
        each["single_source_gate_failure_count"] = np.where(audit["q1"], audit["metric_failure_count"], np.nan)
        each["single_source_all_gates_pass"] = audit["q1_all_gates_pass"]
        rows.append(each)
        s = dict(frequency_hz=float(f), windows=len(starts),
                 not_estimable=int((~audit["estimable"]).sum()), q0=int(audit["q0"].sum()),
                 q1=int(audit["q1"].sum()), q2plus=int(audit["q2"].sum()),
                 accepted=int(result["quality"].sum()),
                 multi_gate_q1=int((audit["q1"] & (audit["metric_failure_count"] >= 2)).sum()),
                 q1_residual_evaluated=int(audit["residual_evaluated"].sum()))
        for key in flags:
            s[key + "_all_windows"] = int(flags[key].sum())
            s[key + "_q1"] = int((flags[key] & audit["q1"]).sum())
        summary_rows.append(s)
        actual_counts = Counter(reason)
        expected = stored[np.isclose(stored.frequency_hz, f)]
        expected_counts = {english: int(row.windows) for row in expected.itertuples()
                           for english, chinese in REASONS.items() if row.reason == chinese}
        for label in sorted(set(actual_counts) | set(expected_counts)):
            check_rows.append(dict(frequency_hz=float(f), reason=label,
                                   recomputed=actual_counts.get(label, 0),
                                   baseline=expected_counts.get(label, 0),
                                   delta=actual_counts.get(label, 0)-expected_counts.get(label, 0)))
        if fi % 5 == 0 or fi == len(frequencies):
            print(f"rejection audit {fi}/{len(frequencies)} frequencies; f={f:.3f} Hz; windows={len(starts)}", flush=True)
    pd.concat(rows, ignore_index=True).to_csv(OUT / "rejection_condition_windows.csv", index=False)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "rejection_condition_summary.csv", index=False)
    pd.DataFrame(check_rows).to_csv(OUT / "rejection_reason_reproduction.csv", index=False)
    print(f"rejection audit complete: {int(summary.windows.sum())} window-frequency rows; "
          f"max label delta={int(pd.DataFrame(check_rows).delta.abs().max())}", flush=True)


def longest_consistent_run(valid, azimuth, max_step):
    best = run = 0
    previous = None
    for ok, az in zip(valid, azimuth):
        if not ok or not np.isfinite(az):
            run, previous = 0, None
            continue
        if previous is None or angular_distance(float(az), previous) <= max_step:
            run += 1
        else:
            run = 1
        best = max(best, run)
        previous = float(az)
    return best


def run_real_coherence_counterfactual():
    """Recompute candidate counts/bands under lower coherence floors only."""
    q = pd.read_csv(OUT / "rejection_condition_windows.csv")
    q = q[q.source_count == 1].copy()
    other_cols = ["phase_residual_fail_fail", "uncertain_fail", "ambiguous_fail",
                  "boundary_fail", "grid_resolution_fail"]
    other_pass = ~q[other_cols].astype(bool).any(axis=1)
    per_frequency, summary = [], []
    freqs = np.sort(q.frequency_hz.unique())
    for threshold in [.40, .50, .60, .65, .70]:
        eligible = other_pass & (q.coherence >= threshold)
        valid_freqs = []
        for f, g in q.groupby("frequency_hz", sort=True):
            mask = (other_pass.loc[g.index] & (g.coherence >= threshold)).to_numpy(bool)
            count = int(mask.sum())
            valid_freqs.append(count > 0)
            per_frequency.append(dict(min_coherence=threshold, frequency_hz=float(f),
                                      single_source_windows=len(g),
                                      pass_other_gates=int(other_pass.loc[g.index].sum()),
                                      counterfactual_accepted_windows=count))
        # Existing scan is at 0.005 Hz; runs are exact contiguous grid indices.
        best_start = best_end = run_start = None
        best_len = run_len = 0
        for i, (f, ok) in enumerate(zip(freqs, valid_freqs)):
            if ok:
                if run_len == 0:
                    run_start = float(f)
                run_len += 1
                if run_len > best_len:
                    best_len, best_start, best_end = run_len, run_start, float(f)
            else:
                run_len = 0
        summary.append(dict(min_coherence=threshold, real_q1_candidates=len(q),
                            pass_all_other_gates=int(other_pass.sum()),
                            counterfactual_accepted_windows=int(eligible.sum()),
                            n_frequencies_with_any=int(np.sum(valid_freqs)),
                            longest_contiguous_points=best_len,
                            longest_band_min_hz=best_start,
                            longest_band_max_hz=best_end,
                            settings_promoted=False))
    pd.DataFrame(per_frequency).to_csv(OUT / "real_coherence_floor_counterfactual_by_frequency.csv", index=False)
    pd.DataFrame(summary).to_csv(OUT / "real_coherence_floor_counterfactual.csv", index=False)
    print("real-data coherence counterfactual complete", flush=True)


def run_common_low_frequency():
    arrays, masks, names, xy, fs, origin = load_actual()
    c = baseline_settings()
    session = session_for(arrays, names, xy, fs, c)
    origin_ts = pd.Timestamp(origin)
    if origin_ts.tzinfo is None:
        origin_ts = origin_ts.tz_localize("UTC")
    midnight = pd.Timestamp("2023-11-01T00:00:00Z")
    centers = [midnight + pd.Timedelta(minutes=30+60*i) for i in range(96)]
    center_samples = np.rint(np.array([(t-origin_ts).total_seconds()*fs for t in centers])).astype(np.int64)
    geometries = {float(f): phase_window_geometry(float(f), fs, c["music_front_cycles"],
                                                   c["music_front_cycles"]) for f in LOW_FREQUENCIES}
    starts_by_f = {f: center_samples - g["window_samples"]//2 for f,g in geometries.items()}
    common = np.ones(len(centers), bool)
    for f, starts in starts_by_f.items():
        common &= valid_windows(starts, geometries[f]["window_samples"], masks, names)
    centers = [t for t, ok in zip(centers, common) if ok]
    center_samples = center_samples[common]
    if not len(centers):
        raise RuntimeError("No UTC-centered windows are valid at all six frequencies")
    starts_by_f = {f: s[common] for f,s in starts_by_f.items()}

    station_ix = {name:i for i,name in enumerate(names)}
    a,b = PAIR.split("--")
    delta = xy[station_ix[b]] - xy[station_ix[a]]
    distance = float(np.linalg.norm(delta))
    tapers, results_rows, peak_rows = [], [], []
    az_by_f, usable_by_f, speed_by_f = {}, {}, {}
    for fi, f in enumerate(LOW_FREQUENCIES, 1):
        f=float(f); n=geometries[f]["window_samples"]; starts=starts_by_f[f]
        taper=hann(n,sym=True)
        dpss_coeff=session.transform(starts,n,taper,[f],snapshots=True)
        exact_coeff=session.transform(starts,n,taper,[f],snapshots=False)[0]
        metrics=[]; peaks_this=[]
        for wi, start in enumerate(starts):
            r, examples=estimate(dpss_coeff[:,:,wi:wi+1,:],xy,np.array([f]),c)
            q=int(r["music_source_count"][0]); reason=str(r["reason"][0])
            audit=gate_flags(r,c); az=float(r["azimuth_deg"][0]) if np.isfinite(r["azimuth_deg"][0]) else np.nan
            gate=pair_gate(r,delta[0],delta[1],c)
            proj=float(delta[0]*np.sin(np.radians(az))+delta[1]*np.cos(np.radians(az))) if np.isfinite(az) else np.nan
            pair_usable=bool(gate["direction_valid"][0] and np.isfinite(proj) and proj>distance*.1)
            cross=exact_coeff[station_ix[b],wi]*np.conj(exact_coeff[station_ix[a],wi])
            phase=float(np.angle(cross)) if np.isfinite(cross.real) and np.isfinite(cross.imag) and abs(cross)>0 else np.nan
            v0=(-2*np.pi*f*proj/phase) if pair_usable and np.isfinite(phase) and phase<0 else np.nan
            if np.isfinite(v0) and not (c["velocity_min_m_s"] < v0 < c["velocity_max_m_s"]):v0=np.nan
            if q>=1 and examples:
                peaks=local_peaks(examples[0]["power"],examples[0]["azimuth"],examples[0]["velocity"],max_peaks=4)
                for rank,p in enumerate(peaks,1):
                    peaks_this.append(dict(frequency_hz=f,center_utc=str(centers[wi]),start_sample=int(start),
                                           source_count=q,rank=rank,azimuth_deg=p["azimuth_deg"],
                                           velocity_m_s=p["velocity_m_s"],relative_power=p["relative_power"]))
            row=dict(frequency_hz=f,center_utc=str(centers[wi]),start_sample=int(start),window_samples=n,
                     source_count=q,final_reason=reason,accepted=bool(r["quality"][0]),
                     coherence=float(r["coherence"][0]),phase_rms_reported_deg=float(r["phase_rms_deg"][0]),
                     phase_rms_evaluated=bool(audit["residual_evaluated"][0]),
                     phase_rms_analysis_deg=float(r["phase_rms_deg"][0]) if audit["residual_evaluated"][0] else np.nan,
                     uncertainty_deg=float(r["direction_uncertainty_deg"][0]),
                     secondary_ratio=float(r["secondary_ratio"][0]),azimuth_deg=az,
                     pair_direction_usable=pair_usable,pair_angle_deg=float(gate["pair_angle_deg"][0]),
                     projected_distance_m=proj,raw_pair_phase_rad=phase,
                     k0_velocity_diagnostic_m_s=v0,velocity_cycle_assumption=0 if np.isfinite(v0) else np.nan,
                     q1_all_direction_gates_pass=bool(audit["q1_all_gates_pass"][0]))
            for key,val in audit["flags"].items():row[key+"_fail"]=bool(val[0])
            metrics.append(row)
        frame=pd.DataFrame(metrics)
        results_rows.append(frame)
        peak_rows.extend(peaks_this)
        az_by_f[f]=frame.azimuth_deg.to_numpy(float)
        usable_by_f[f]=(frame.source_count.to_numpy(int)==1)&frame.accepted.to_numpy(bool)&frame.pair_direction_usable.to_numpy(bool)
        speed_by_f[f]=frame.k0_velocity_diagnostic_m_s.to_numpy(float)
        print(f"common-window {fi}/{len(LOW_FREQUENCIES)}: {f:.2f} Hz, "
              f"centers={len(frame)}, q1 accepted={int(usable_by_f[f].sum())}, "
              f"k0 speeds={int(np.isfinite(speed_by_f[f]).sum())}",flush=True)

    windows=pd.concat(results_rows,ignore_index=True)
    peaks=pd.DataFrame(peak_rows)
    windows.to_csv(OUT/"low_frequency_common_window_results.csv",index=False)
    peaks.to_csv(OUT/"low_frequency_common_window_peaks.csv",index=False)

    summary=[]
    for f in LOW_FREQUENCIES:
        f=float(f); g=windows[np.isclose(windows.frequency_hz,f)]
        valid_speed=g.k0_velocity_diagnostic_m_s.dropna().to_numpy(float)
        primary=peaks[(np.isclose(peaks.frequency_hz,f))&(peaks["rank"]==1)] if len(peaks) else pd.DataFrame()
        modal_share=np.nan
        if len(primary):
            bins=np.floor((primary.azimuth_deg.to_numpy(float)+7.5)/15).astype(int)%24
            modal_share=float(Counter(bins).most_common(1)[0][1]/len(primary))
        summary.append(dict(frequency_hz=f,common_utc_centers=len(g),q0=int((g.source_count==0).sum()),
                            q1=int((g.source_count==1).sum()),q2plus=int((g.source_count>=2).sum()),
                            q1_direction_accepted=int(g.q1_all_direction_gates_pass.sum()),
                            pair_direction_usable=int(g.pair_direction_usable.sum()),
                            k0_velocity_count=len(valid_speed),k0_velocity_median_m_s=float(np.median(valid_speed)) if len(valid_speed) else np.nan,
                            k0_velocity_q025_m_s=float(np.quantile(valid_speed,.025)) if len(valid_speed) else np.nan,
                            k0_velocity_q975_m_s=float(np.quantile(valid_speed,.975)) if len(valid_speed) else np.nan,
                            primary_peak_count=len(primary),primary_direction_modal_share_15deg=modal_share))
    summary=pd.DataFrame(summary)
    summary.to_csv(OUT/"low_frequency_common_window_summary.csv",index=False)

    branch_rows=[]
    for i,center in enumerate(centers):
        az=np.array([az_by_f[float(f)][i] for f in LOW_FREQUENCIES],float)
        ok=np.array([usable_by_f[float(f)][i] for f in LOW_FREQUENCIES],bool)
        row=dict(center_utc=str(center),usable_single_source_frequencies=int(ok.sum()),
                 all_six_usable=bool(ok.all()),longest_branch_run_15deg=longest_consistent_run(ok,az,15.),
                 longest_branch_run_30deg=longest_consistent_run(ok,az,30.))
        row["full_six_branch_15deg"]=bool(ok.all() and np.all(angular_distance(az[1:],az[:-1])<=15.))
        row["full_six_branch_30deg"]=bool(ok.all() and np.all(angular_distance(az[1:],az[:-1])<=30.))
        for j,f in enumerate(LOW_FREQUENCIES):
            row[f"azimuth_{f:.2f}"]=az[j] if ok[j] else np.nan
            row[f"k0_velocity_{f:.2f}"]=speed_by_f[float(f)][i]
        speeds=np.array([speed_by_f[float(f)][i] for f in LOW_FREQUENCIES],float)
        row["all_six_k0_speed_valid"]=bool(np.isfinite(speeds).all())
        branch_rows.append(row)
    branch=pd.DataFrame(branch_rows)
    top_peak_rows=[]
    if len(peaks):
        prim=peaks[peaks["rank"]==1]
        for center in centers:
            g=prim[prim.center_utc==str(center)].set_index("frequency_hz").reindex(LOW_FREQUENCIES)
            az=g.azimuth_deg.to_numpy(float);ok=np.isfinite(az)
            row=dict(center_utc=str(center),primary_peak_count=int(ok.sum()),
                     primary_peak_longest_run_15deg=longest_consistent_run(ok,az,15.),
                     primary_peak_longest_run_30deg=longest_consistent_run(ok,az,30.),
                     primary_peak_full_six_15deg=bool(ok.all() and np.all(angular_distance(az[1:],az[:-1])<=15.)),
                     primary_peak_full_six_30deg=bool(ok.all() and np.all(angular_distance(az[1:],az[:-1])<=30.)))
            for j,f in enumerate(LOW_FREQUENCIES):row[f"primary_peak_azimuth_{f:.2f}"]=az[j]
            top_peak_rows.append(row)
    top_peaks=pd.DataFrame(top_peak_rows)
    if len(top_peaks):branch=branch.merge(top_peaks,on="center_utc",how="left")
    branch.to_csv(OUT/"low_frequency_common_window_branch_tracks.csv",index=False)

    peak_links=[]
    if len(peaks):
        for i in range(len(LOW_FREQUENCIES)-1):
            f0=float(LOW_FREQUENCIES[i]);f1=float(LOW_FREQUENCIES[i+1])
            p0=peaks[(np.isclose(peaks.frequency_hz,f0))&(peaks["rank"]==1)].set_index("center_utc")
            p1=peaks[(np.isclose(peaks.frequency_hz,f1))&(peaks["rank"]==1)].set_index("center_utc")
            matched=p0.join(p1,lsuffix="_low",rsuffix="_high",how="inner")
            dist=angular_distance(matched.azimuth_deg_low.to_numpy(float),matched.azimuth_deg_high.to_numpy(float)) if len(matched) else np.array([])
            peak_links.append(dict(frequency_low_hz=f0,frequency_high_hz=f1,centers_with_both_primary_peaks=len(matched),
                                   primary_peak_links_15deg=int(np.sum(dist<=15)),primary_peak_links_30deg=int(np.sum(dist<=30)),
                                   median_primary_peak_azimuth_shift_deg=float(np.median(dist)) if len(dist) else np.nan))
    pd.DataFrame(peak_links).to_csv(OUT/"low_frequency_common_peak_links.csv",index=False)

    # Paired adjacent-frequency differences at the same UTC centers, resampled
    # synchronously by 6-hour blocks. Marginal medians are not used as paired deltas.
    block_labels=np.array([pd.Timestamp(t).floor("6h").isoformat() for t in centers])
    unique_blocks=np.unique(block_labels)
    wide=branch.set_index("center_utc")[[f"k0_velocity_{f:.2f}" for f in LOW_FREQUENCIES]]
    wide=wide.reindex([str(t) for t in centers])
    rng=np.random.default_rng(SEED); boot_rows=[]
    vals=wide.to_numpy(float); paired_masks=[]; paired_deltas=[]
    for j in range(len(LOW_FREQUENCIES)-1):
        paired=np.isfinite(vals[:,j])&np.isfinite(vals[:,j+1])
        paired_masks.append(paired);paired_deltas.append(vals[paired,j+1]-vals[paired,j])
    for bi in range(REPS):
        chosen=rng.choice(unique_blocks,size=len(unique_blocks),replace=True)
        block_counts=Counter(chosen.tolist())
        weights=np.array([block_counts.get(b,0) for b in block_labels],int)
        row=dict(draw=bi,unique_blocks=int(np.sum(weights>0)))
        for j,(paired,delta) in enumerate(zip(paired_masks,paired_deltas)):
            take=paired&(weights>0)
            weighted=np.repeat(delta[np.isin(np.flatnonzero(paired),np.flatnonzero(take))],weights[take]) if take.any() else np.array([])
            name=f"delta_v_{LOW_FREQUENCIES[j]:.2f}_{LOW_FREQUENCIES[j+1]:.2f}"
            row[name]=float(np.median(weighted)) if len(weighted) else np.nan
        boot_rows.append(row)
    boot=pd.DataFrame(boot_rows);boot.to_csv(OUT/"low_frequency_common_window_block_bootstrap.csv",index=False)
    delta_rows=[]
    for j in range(len(LOW_FREQUENCIES)-1):
        name=f"delta_v_{LOW_FREQUENCIES[j]:.2f}_{LOW_FREQUENCIES[j+1]:.2f}"
        vals=boot[name].dropna().to_numpy(float)
        paired=paired_masks[j]; ncenters=int(paired.sum()); nblocks=int(np.unique(block_labels[paired]).size)
        enough=ncenters>=10 and nblocks>=4
        if not enough:vals=np.array([],float)
        raw_delta=paired_deltas[j]
        delta_rows.append(dict(frequency_low_hz=LOW_FREQUENCIES[j],frequency_high_hz=LOW_FREQUENCIES[j+1],
                               paired_centers=ncenters,paired_six_hour_blocks=nblocks,
                               bootstrap_estimable=bool(enough),observed_delta_median_m_s=float(np.median(raw_delta)) if len(raw_delta) else np.nan,
                               complete_bootstrap_draws=len(vals),delta_q025_m_s=float(np.quantile(vals,.025)) if len(vals) else np.nan,
                               delta_median_m_s=float(np.quantile(vals,.5)) if len(vals) else np.nan,
                               delta_q975_m_s=float(np.quantile(vals,.975)) if len(vals) else np.nan,
                               positive_fraction=float(np.mean(vals>0)) if len(vals) else np.nan))
    pd.DataFrame(delta_rows).to_csv(OUT/"low_frequency_common_window_block_deltas.csv",index=False)
    print(f"common-window analysis complete: {len(centers)} common UTC centers; "
          f"six-way q1+pair+phase speeds={int(branch.all_six_k0_speed_valid.sum())}; "
          f"all-six 15deg q1 branches={int(branch.full_six_branch_15deg.sum())}; "
          f"all-six 15deg primary peaks={int(branch.primary_peak_full_six_15deg.sum())}",flush=True)


def run_low_coherence_simulation():
    """Counterfactual lower coherence floors on matched known-answer trials."""
    names,xy=station_geometry();base=baseline_settings()
    thresholds=[.40,.50,.60,.65,.70]
    cases=[]
    for f,v in [(.10,3000.),(.45,1600.)]:
        for trial in range(100):
            arrays,n=simulate_series(names,xy,f,v,45.,400000+trial,snr_db=None,signal_enabled=False)
            cases.append(("pure_noise",f,np.nan,trial,arrays,n))
        for trial in range(24):
            arrays,n=simulate_independent_band_waves(names,xy,f,v,45.,620000+trial,snr_db=0)
            cases.append(("single_wave_snr0",f,0.,trial,arrays,n))
    rows=[]
    ctest=dict(base,min_coherence=min(thresholds))
    for i,(case,f,snr,trial,arrays,n) in enumerate(cases,1):
        z,r=estimate_one(arrays,names,xy,f,n,ctest)
        audit=gate_flags(r,ctest);q1=audit["q1"]
        other=["phase_residual_fail","uncertain","ambiguous","boundary","grid_resolution"]
        other_pass=q1.copy()
        for key in other:other_pass &= ~audit["flags"][key]
        coh=float(r["coherence"][0])
        for threshold in thresholds:
            rows.append(dict(case=case,frequency_hz=f,snr_db=snr,trial=trial,
                             min_coherence=threshold,source_count=int(r["music_source_count"][0]),
                             coherence=coh,other_gates_pass=bool(other_pass[0]),
                             accepted_counterfactual=bool(other_pass[0] and coh>=threshold)))
        if i%40==0 or i==len(cases):print(f"lower-coherence simulation {i}/{len(cases)}",flush=True)
    frame=pd.DataFrame(rows);frame.to_csv(OUT/"coherence_lower_floor_simulation_trials.csv",index=False)
    summary=frame.groupby(["case","frequency_hz","min_coherence"],as_index=False).agg(
        trials=("accepted_counterfactual","size"),accepted=("accepted_counterfactual","sum"),
        acceptance_rate=("accepted_counterfactual","mean"),q1_other_gates_pass=("other_gates_pass","sum"),
        median_coherence=("coherence","median"))
    summary.to_csv(OUT/"coherence_lower_floor_simulation_summary.csv",index=False)
    print(f"lower-coherence simulation complete: {len(cases)} paired cases",flush=True)


def main():
    if len(sys.argv)==1 or sys.argv[1]=="all":
        run_front_rejection_audit();run_common_low_frequency();run_low_coherence_simulation()
    elif sys.argv[1]=="rejection":run_front_rejection_audit()
    elif sys.argv[1]=="common":run_common_low_frequency()
    elif sys.argv[1]=="coherence":run_low_coherence_simulation()
    elif sys.argv[1]=="counterfactual":run_real_coherence_counterfactual()
    else:raise SystemExit("Usage: run_next_diagnostics.py [all|rejection|common|coherence|counterfactual]")


if __name__=="__main__":main()
