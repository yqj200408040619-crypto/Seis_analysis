#!/usr/bin/env python3
"""Isolated validation of the MGS6 MUSIC workflow; never edits baseline outputs."""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent / "project"
WORKSPACE = HERE.parents[1]
SOURCE = WORKSPACE / "source_2023" / "NoiseDispersionStudio_v2_7_1"
OUT = HERE
OUT.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("NUMBA_CACHE_DIR", str(OUT / "numba_cache"))
sys.path[:0] = [str(SOURCE), str(SOURCE / "core")]

import h5py
import numpy as np
import pandas as pd
from scipy.signal.windows import hann

from ndstudio.direction import DirectionSession, angular_distance, settings
from ndstudio.music import estimate
from ndstudio.phase_screening import screen_cross
from ndstudio.phase_stack import stack_stats
from ndstudio.phase_velocity import lift_band
from ndstudio.reconstruction import align_cross, phase_speed


BASE = PROJECT / "output"
FS = 10.0
BOOT_SEED = 20261102
BOOT_REPS = 499


def baseline_settings():
    doc = json.loads((BASE / "music_front/music_settings_used.json").read_text())
    return settings(doc["settings"]["direction"])


def station_geometry():
    table = pd.read_csv(BASE / "music_front/direction/station_geometry.csv")
    names = table.station.astype(str).tolist()
    xy = table.set_index("station").loc[names, ["x_m", "y_m"]].to_numpy(float)
    return names, xy


def session_for(arrays, names, xy, fs, c):
    s = DirectionSession.__new__(DirectionSession)
    s.arrays = arrays
    s.names = list(names)
    s.xy = np.asarray(xy, float)
    s.fs = float(fs)
    s.limits = {name: (0, len(arrays[name])) for name in names}
    s.stable = 0
    s.c = c
    return s


def simulate_series(names, xy, f, velocity, azimuth, seed, snr_db=None,
                    second=None, gains=None, signal_enabled=True):
    rng = np.random.default_rng(seed)
    n = int(round(200.0 / f * FS))
    t = np.arange(n, dtype=float) / FS
    theta = np.radians(azimuth)
    axis = np.array([np.sin(theta), np.cos(theta)])
    phase1 = rng.uniform(-np.pi, np.pi)
    delays1 = (xy @ axis) / velocity
    gains = np.ones(len(names)) if gains is None else np.asarray(gains, float)
    sig = np.stack([
        gains[i] * np.cos(2*np.pi*f*(t-delays1[i]) + phase1)
        for i in range(len(names))
    ])
    if second is not None:
        az2, v2, amp2 = second
        theta2 = np.radians(az2)
        delays2 = (xy @ np.array([np.sin(theta2), np.cos(theta2)])) / v2
        phase2 = rng.uniform(-np.pi, np.pi)
        sig += np.stack([
            gains[i] * amp2 * np.cos(2*np.pi*f*(t-delays2[i]) + phase2)
            for i in range(len(names))
        ])
    if not signal_enabled:
        sig.fill(0.0)
    if snr_db is None:
        noise_std = 1.0
    else:
        # Unit-amplitude sinusoid has RMS 1/sqrt(2).
        noise_std = 1.0 / (np.sqrt(2.0) * 10.0**(snr_db/20.0))
    noise = rng.normal(0.0, noise_std, size=sig.shape)
    return {name: sig[i] + noise[i] for i, name in enumerate(names)}, n


def narrowband_spectrum(f_center, bandwidth, n, rng):
    freq = np.fft.fftfreq(n, d=1.0/FS)
    positive = (freq > 0) & (np.abs(freq-f_center) <= bandwidth)
    spec = np.zeros(n, complex)
    pos_ix = np.flatnonzero(positive)
    spec[pos_ix] = rng.normal(size=len(pos_ix)) + 1j*rng.normal(size=len(pos_ix))
    spec[-pos_ix] = np.conj(spec[pos_ix])
    wave = np.fft.ifft(spec).real
    rms = np.sqrt(np.mean(wave*wave))
    return spec / max(rms, 1e-30), freq


def simulate_independent_band_waves(names, xy, f, velocity, azimuth, seed,
                                    snr_db=20, second=None):
    """Independent narrowband random wavefields delayed across the real array."""
    rng = np.random.default_rng(seed)
    n = int(round(200.0 / f * FS))
    duration = n/FS
    bw = 1.45 * 3.5 / duration
    waves = [(f, velocity, azimuth, 1.0)]
    if second is not None:
        az2, v2, amp2 = second
        waves = [(f-0.35*bw, velocity, azimuth, 1.0),
                 (f+0.35*bw, v2, az2, amp2)]
    freq = np.fft.fftfreq(n, d=1.0/FS)
    spectra = [narrowband_spectrum(fc, bw, n, rng)[0] for fc, *_ in waves]
    sigma = 1.0/(np.sqrt(2.0)*10.0**(snr_db/20.0))
    out = {}
    for i, name in enumerate(names):
        signal = np.zeros(n, float)
        for spec, (fc, vel, az, amp) in zip(spectra, waves):
            theta = np.radians(az)
            delay = np.dot(xy[i], [np.sin(theta), np.cos(theta)]) / vel
            signal += amp*np.fft.ifft(spec*np.exp(-2j*np.pi*freq*delay)).real
        out[name] = signal + rng.normal(0.0, sigma, n)
    return out, n


def estimate_one(arrays, names, xy, f, n, c):
    s = session_for(arrays, names, xy, FS, c)
    z = s.transform(np.array([0], dtype=np.int64), n, hann(n, sym=True), [f], snapshots=True)
    result, _ = estimate(z, xy, np.array([f]), c)
    return z, result


def covariance_summary(z):
    x = np.asarray(z[0, :, 0, :], complex)
    energy = np.mean(abs(x)**2, axis=1)
    x = x / np.sqrt(np.maximum(energy, 1e-300))[:, None]
    r = x @ x.conj().T / x.shape[1]
    ev = np.linalg.eigvalsh((r + r.conj().T) / 2)[::-1]
    return float(ev[0] / max(np.mean(ev[1:]), 1e-30)), ev.tolist()


def summarize_trials(frame, group_cols, out_path):
    rows = []
    for key, g in frame.groupby(group_cols, dropna=False, sort=True):
        key = key if isinstance(key, tuple) else (key,)
        row = dict(zip(group_cols, key))
        row.update(trials=len(g), accepted=int(g.accepted.sum()),
                   acceptance_rate=float(g.accepted.mean()),
                   q0=int((g.source_count == 0).sum()),
                   q1=int((g.source_count == 1).sum()),
                   q2plus=int((g.source_count >= 2).sum()))
        if "expected_sources" in g:
            row["source_count_correct"] = int((g.source_count == g.expected_sources).sum())
            row["source_count_accuracy"] = float((g.source_count == g.expected_sources).mean())
        if "azimuth_error_deg" in g:
            q = g[g.azimuth_error_deg.notna()]
            row["direction_scored"] = len(q)
            row["median_azimuth_error_deg"] = float(q.azimuth_error_deg.median()) if len(q) else None
            row["p90_azimuth_error_deg"] = float(q.azimuth_error_deg.quantile(.9)) if len(q) else None
            row["median_abs_velocity_error_fraction"] = float(q.abs_velocity_error_fraction.median()) if len(q) else None
        rows.append(row)
    pd.DataFrame(rows).to_csv(out_path, index=False)
    return rows


def run_simulation():
    names, xy = station_geometry()
    base = baseline_settings()
    rng = np.random.default_rng(7711)
    rows = []
    trial = 0

    for f, truth_v, truth_az in [(0.10, 3000.0, 45.0), (0.45, 1600.0, 45.0)]:
        n = int(round(200.0 / f * FS))
        for snr in [20, 10, 0, -5]:
            for k in range(8):
                trial += 1
                arrays, n = simulate_series(names, xy, f, truth_v, truth_az,
                                            100000 + trial, snr_db=snr)
                z, r = estimate_one(arrays, names, xy, f, n, base)
                q = int(r["music_source_count"][0])
                accepted = bool(r["quality"][0])
                azerr = float(angular_distance(r["azimuth_deg"][0], truth_az)) if q else np.nan
                verr = abs(float(r["velocity_m_s"][0]) / truth_v - 1) if q else np.nan
                eig_ratio, eig = covariance_summary(z)
                rows.append(dict(case="single_wave", frequency_hz=f, snr_db=snr,
                    gain_pattern="equal", trial=k, expected_sources=1, source_count=q,
                    reason=str(r["reason"][0]), accepted=accepted,
                    azimuth_error_deg=azerr, abs_velocity_error_fraction=verr,
                    coherence=float(r["coherence"][0]) if q else np.nan,
                    eigenvalue_ratio=eig_ratio, eigenvalues_json=json.dumps(eig)))

        gains = np.array([0.12, 0.22, 0.45, 0.9, 1.8, 3.4])
        for k in range(12):
            trial += 1
            arrays, n = simulate_series(names, xy, f, truth_v, truth_az,
                200000 + trial, snr_db=10, gains=gains)
            z, r = estimate_one(arrays, names, xy, f, n, base)
            q = int(r["music_source_count"][0])
            eig_ratio, eig = covariance_summary(z)
            rows.append(dict(case="single_wave_gain_imbalance", frequency_hz=f,
                snr_db=10, gain_pattern="0.12,0.22,0.45,0.9,1.8,3.4", trial=k,
                expected_sources=1, source_count=q, reason=str(r["reason"][0]),
                accepted=bool(r["quality"][0]),
                azimuth_error_deg=float(angular_distance(r["azimuth_deg"][0], truth_az)) if q else np.nan,
                abs_velocity_error_fraction=abs(float(r["velocity_m_s"][0])/truth_v-1) if q else np.nan,
                coherence=float(r["coherence"][0]) if q else np.nan,
                eigenvalue_ratio=eig_ratio, eigenvalues_json=json.dumps(eig)))

    for f, v in [(0.10, 3000.0), (0.45, 1600.0)]:
        for k in range(20):
            trial += 1
            arrays, n = simulate_series(names, xy, f, v, 45.0,
                300000 + trial, snr_db=20, second=(125.0, 2200.0 if f < .2 else 2100.0, .6))
            z, r = estimate_one(arrays, names, xy, f, n, base)
            q = int(r["music_source_count"][0])
            eig_ratio, eig = covariance_summary(z)
            rows.append(dict(case="two_coherent_tones", frequency_hz=f, snr_db=20,
                gain_pattern="equal", trial=k, expected_sources=2, source_count=q,
                reason=str(r["reason"][0]), accepted=bool(r["quality"][0]),
                azimuth_error_deg=np.nan, abs_velocity_error_fraction=np.nan,
                coherence=float(r["coherence"][0]) if q else np.nan,
                eigenvalue_ratio=eig_ratio, eigenvalues_json=json.dumps(eig)))

        for k in range(24):
            trial += 1
            arrays, n = simulate_independent_band_waves(names, xy, f, v, 45.0,
                500000 + trial, snr_db=20,
                second=(125.0, 2200.0 if f < .2 else 2100.0, .6))
            z, r = estimate_one(arrays, names, xy, f, n, base)
            q = int(r["music_source_count"][0])
            eig_ratio, eig = covariance_summary(z)
            rows.append(dict(case="two_independent_narrowband_waves", frequency_hz=f,
                snr_db=20, gain_pattern="equal", trial=k, expected_sources=2,
                source_count=q, reason=str(r["reason"][0]), accepted=bool(r["quality"][0]),
                azimuth_error_deg=np.nan, abs_velocity_error_fraction=np.nan,
                coherence=float(r["coherence"][0]) if q else np.nan,
                eigenvalue_ratio=eig_ratio, eigenvalues_json=json.dumps(eig)))

        for k in range(100):
            trial += 1
            arrays, n = simulate_series(names, xy, f, v, 45.0,
                                         400000 + trial, snr_db=None, signal_enabled=False)
            z, r = estimate_one(arrays, names, xy, f, n, base)
            q = int(r["music_source_count"][0])
            eig_ratio, eig = covariance_summary(z)
            rows.append(dict(case="pure_noise", frequency_hz=f, snr_db=np.nan,
                gain_pattern="equal", trial=k, expected_sources=0, source_count=q,
                reason=str(r["reason"][0]), accepted=bool(r["quality"][0]),
                azimuth_error_deg=np.nan, abs_velocity_error_fraction=np.nan,
                coherence=float(r["coherence"][0]) if q else np.nan,
                eigenvalue_ratio=eig_ratio, eigenvalues_json=json.dumps(eig)))

    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "simulation_trials.csv", index=False)
    summarize_trials(frame, ["case", "frequency_hz", "snr_db", "gain_pattern"],
                     OUT / "simulation_summary.csv")

    # The false-accept rate and signal recovery are compared across coherence cutoffs.
    sweep_rows = []
    for f, v in [(0.10, 3000.0), (0.45, 1600.0)]:
        for case, snr, count in [("pure_noise", None, 40), ("single_wave_snr0", 0, 24)]:
            for seed_ix in range(count):
                if case == "pure_noise":
                    arrays, n = simulate_series(names, xy, f, v, 45.0,
                        610000 + seed_ix, snr_db=None, signal_enabled=False)
                else:
                    arrays, n = simulate_independent_band_waves(names, xy, f, v, 45.0,
                        620000 + seed_ix, snr_db=snr)
                sess = session_for(arrays, names, xy, FS, base)
                z = sess.transform(np.array([0], dtype=np.int64), n, hann(n, sym=True),
                                   [f], snapshots=True)
                for threshold in [.70, .80, .90, .95, .98]:
                    c = settings(dict(base, min_coherence=threshold))
                    r, _ = estimate(z, xy, np.array([f]), c)
                    q = int(r["music_source_count"][0])
                    sweep_rows.append(dict(case=case, frequency_hz=f, trial=seed_ix,
                        min_coherence=threshold, source_count=q, reason=str(r["reason"][0]),
                        accepted=bool(r["quality"][0]),
                        coherence=float(r["coherence"][0]) if q else np.nan,
                        azimuth_error_deg=float(angular_distance(r["azimuth_deg"][0],45.0))
                            if q and case != "pure_noise" else np.nan,
                        abs_velocity_error_fraction=abs(float(r["velocity_m_s"][0])/v-1)
                            if q and case != "pure_noise" else np.nan))
    sweep = pd.DataFrame(sweep_rows)
    sweep.to_csv(OUT / "coherence_threshold_sweep_trials.csv", index=False)
    sweep_summary = sweep.groupby(["case","frequency_hz","min_coherence"],as_index=False).agg(
        trials=("accepted","size"), accepted=("accepted","sum"),
        acceptance_rate=("accepted","mean"), median_coherence=("coherence","median"),
        median_azimuth_error_deg=("azimuth_error_deg","median"),
        median_abs_velocity_error_fraction=("abs_velocity_error_fraction","median"))
    sweep_summary.to_csv(OUT / "coherence_threshold_sweep_summary.csv", index=False)
    # Evaluate a predeclared training/validation split and the attainable ROC of this gate.
    base_rows = sweep[sweep.min_coherence == .70].copy()
    trade_rows = []
    split_rows = []
    roc_thresholds = [.70,.90,.95,.98,.99,.995,.999,.9995,.9999]
    for f in sorted(base_rows.frequency_hz.unique()):
        noise = base_rows[(base_rows.frequency_hz == f) & (base_rows.case == "pure_noise")]
        wave = base_rows[(base_rows.frequency_hz == f) & (base_rows.case == "single_wave_snr0")]
        for threshold in roc_thresholds:
            nscore = noise.coherence.to_numpy(float)
            wscore = wave.coherence.to_numpy(float)
            fpr = float(np.mean(noise.accepted.to_numpy(bool) & (nscore >= threshold)))
            tpr = float(np.mean(wave.accepted.to_numpy(bool) & (wscore >= threshold)))
            trade_rows.append(dict(frequency_hz=f, min_coherence=threshold,
                end_to_end_pure_noise_false_accept_rate=fpr,
                end_to_end_single_wave_snr0_recall=tpr))
        ntrain = noise[noise.trial < 24]
        nvalid = noise[noise.trial >= 24]
        wtrain = wave[wave.trial < 14]
        wvalid = wave[wave.trial >= 14]
        eligible_thresholds = []
        for candidate in roc_thresholds:
            train_fpr = float(np.mean(ntrain.accepted.to_numpy(bool) &
                (ntrain.coherence.to_numpy(float) >= candidate)))
            train_tpr = float(np.mean(wtrain.accepted.to_numpy(bool) &
                (wtrain.coherence.to_numpy(float) >= candidate)))
            if train_fpr <= .05 and train_tpr >= .90:
                eligible_thresholds.append(candidate)
        threshold = max(eligible_thresholds) if eligible_thresholds else np.nan
        split_rows.append(dict(frequency_hz=f,
            calibration_rule="max cutoff meeting training false-accept<=5% and wave recall>=90%",
            training_noise_trials=len(ntrain), training_wave_trials=len(wtrain),
            heldout_noise_trials=len(nvalid), heldout_wave_trials=len(wvalid),
            selected_min_coherence=threshold,
            training_false_accept_rate=float(np.mean(ntrain.accepted.to_numpy(bool) &
                (ntrain.coherence.to_numpy(float) >= threshold))) if np.isfinite(threshold) else np.nan,
            training_single_wave_recall=float(np.mean(wtrain.accepted.to_numpy(bool) &
                (wtrain.coherence.to_numpy(float) >= threshold))) if np.isfinite(threshold) else np.nan,
            heldout_false_accept_rate=float(np.mean(nvalid.accepted.to_numpy(bool) &
                (nvalid.coherence.to_numpy(float) >= threshold))) if np.isfinite(threshold) else np.nan,
            heldout_single_wave_recall=float(np.mean(wvalid.accepted.to_numpy(bool) &
                (wvalid.coherence.to_numpy(float) >= threshold))) if np.isfinite(threshold) else np.nan))
    pd.DataFrame(trade_rows).to_csv(OUT / "coherence_threshold_tradeoff.csv", index=False)
    pd.DataFrame(split_rows).to_csv(OUT / "coherence_threshold_holdout.csv", index=False)

    # One factor changes at a time: K at fixed NW, then NW at fixed K.
    ablation = []
    configs = [(3.5, 4), (3.5, 6), (3.5, 7), (4.0, 6)]
    for f, v in [(0.10, 3000.0), (0.45, 1600.0)]:
        for snr in [10, 0]:
            for nw, kmax in configs:
                c = settings(dict(base, music_time_bandwidth=nw, music_snapshots=kmax))
                for seed_ix in range(10):
                    arrays, n = simulate_series(names, xy, f, v, 45.0,
                        900000 + seed_ix, snr_db=snr)
                    z, r = estimate_one(arrays, names, xy, f, n, c)
                    q = int(r["music_source_count"][0])
                    eratio, eig = covariance_summary(z)
                    ablation.append(dict(case="single_wave", frequency_hz=f, snr_db=snr, NW=nw,
                        snapshots=kmax, half_bandwidth_hz=nw/(200.0/f), trial=seed_ix,
                        expected_sources=1,
                        source_count=q, reason=str(r["reason"][0]), accepted=bool(r["quality"][0]),
                        azimuth_error_deg=float(angular_distance(r["azimuth_deg"][0],45.0)) if q else np.nan,
                        abs_velocity_error_fraction=abs(float(r["velocity_m_s"][0])/v-1) if q else np.nan,
                        coherence=float(r["coherence"][0]) if q else np.nan,
                        eigenvalue_ratio=eratio, eigenvalues_json=json.dumps(eig)))
                if snr == 10:
                    for seed_ix in range(12):
                        arrays, n = simulate_series(names, xy, f, v, 45.0,
                            930000 + seed_ix, snr_db=None, signal_enabled=False)
                        z, r = estimate_one(arrays, names, xy, f, n, c)
                        q = int(r["music_source_count"][0])
                        eratio, eig = covariance_summary(z)
                        ablation.append(dict(case="pure_noise", frequency_hz=f, snr_db=np.nan, NW=nw,
                            snapshots=kmax, half_bandwidth_hz=nw/(200.0/f), trial=seed_ix,
                            expected_sources=0, source_count=q, reason=str(r["reason"][0]),
                            accepted=bool(r["quality"][0]), azimuth_error_deg=np.nan,
                            abs_velocity_error_fraction=np.nan,
                            coherence=float(r["coherence"][0]) if q else np.nan,
                            eigenvalue_ratio=eratio, eigenvalues_json=json.dumps(eig)))
    abl = pd.DataFrame(ablation)
    abl.to_csv(OUT / "nw_snapshot_ablation_trials.csv", index=False)
    summarize_trials(abl, ["case", "frequency_hz", "snr_db", "NW", "snapshots"],
                     OUT / "nw_snapshot_ablation_summary.csv")
    print(f"synthetic trials={len(frame)}, ablation trials={len(abl)}")


def load_actual():
    path = BASE / "music_front/aligned_waveforms.h5"
    arrays = {}
    masks = {}
    with h5py.File(path, "r") as h:
        fs = float(h.attrs["fs_hz"])
        origin = str(h.attrs["start_utc"])
        for name in sorted(h["stations"].keys()):
            arr = h[f"stations/{name}/data"][...]
            valid = h[f"stations/{name}/valid"][...].astype(bool)
            arr[~valid] = np.nan
            arrays[name] = arr
            masks[name] = valid
    names, xy = station_geometry()
    return arrays, masks, names, xy, fs, origin


def get_direction_archive(target_f, folder="music_front/direction"):
    files = list((BASE / folder).glob("*.npz"))
    scored = []
    for path in files:
        with np.load(path, allow_pickle=False) as z:
            if "frequency_hz" not in z:
                continue
            f = float(z["frequency_hz"])
            scored.append((abs(f-target_f), f, path))
    if not scored:
        raise FileNotFoundError(f"No direction archive in {folder}")
    d, f, path = min(scored, key=lambda x: x[0])
    with np.load(path, allow_pickle=False) as z:
        data = {k: z[k] for k in z.files}
    return path, f, data


def valid_windows(starts, n, masks, names):
    good = np.ones(len(starts), bool)
    for name in names:
        m = masks[name]
        prefix = np.r_[0, np.cumsum(~m, dtype=np.int64)]
        good &= (prefix[starts+n] - prefix[starts]) == 0
    return good


def result_summary(label, result, names, drop=None):
    q = result["music_source_count"]
    goodq = q == 1
    acceptedq = goodq & result["quality"]
    quality = result["quality"]
    return dict(label=label, dropped_station=drop or "none", windows=len(q),
        accepted=int(quality.sum()), accepted_fraction=float(quality.mean()) if len(q) else np.nan,
        q0=int((q == 0).sum()), q1=int((q == 1).sum()), q2=int((q >= 2).sum()),
        low_coherence=int((result["reason"] == "low_coherence").sum()),
        uncertain=int((result["reason"] == "uncertain").sum()),
        ambiguous=int((result["reason"] == "ambiguous").sum()),
        median_coherence=float(np.nanmedian(result["coherence"])) if goodq.any() else np.nan,
        median_azimuth_deg=float(np.nanmedian(result["azimuth_deg"][acceptedq])) if acceptedq.any() else np.nan,
        median_velocity_m_s=float(np.nanmedian(result["velocity_m_s"][acceptedq])) if acceptedq.any() else np.nan,
        geometry_ratio=float(result["geometry_ratio"][0]) if len(q) else np.nan)


def local_peaks(power, az, vel, min_separation_deg=12.0, max_peaks=4):
    p = np.asarray(power, float)
    nv, na = p.shape
    cand = []
    for vi in range(1, nv-1):
        for ai in range(na):
            value = p[vi, ai]
            hood = p[vi-1:vi+2, :][:, [(ai-1) % na, ai, (ai+1) % na]]
            if value >= np.max(hood) and np.isfinite(value):
                cand.append((float(value), vi, ai))
    cand.sort(reverse=True)
    kept = []
    for value, vi, ai in cand:
        theta = float(az[ai])
        speed = float(vel[vi])
        if any(float(angular_distance(theta, q["azimuth_deg"])) < min_separation_deg and
               abs(np.log(speed/q["velocity_m_s"])) < .12 for q in kept):
            continue
        kept.append(dict(azimuth_deg=theta, velocity_m_s=speed,
                         relative_power=value/ max(cand[0][0],1e-30)))
        if len(kept) >= max_peaks:
            break
    return kept


def run_real():
    arrays, masks, names, xy, fs, origin = load_actual()
    base = baseline_settings()
    sess = session_for(arrays, names, xy, fs, base)
    rows = []
    coeff_by_f = {}
    for target, folder in [(.10, "phase_guided/direction"), (.45, "music_front/direction")]:
        archive, f, meta = get_direction_archive(target, folder)
        n = int(meta["window_samples"])
        all_starts = meta["start_sample"]
        ok = valid_windows(all_starts, n, masks, names)
        candidates = all_starts[ok]
        if len(candidates) < 20:
            raise RuntimeError(f"Only {len(candidates)} fully valid windows at {f} Hz")
        take = np.unique(np.linspace(0, len(candidates)-1, min(32, len(candidates)), dtype=int))
        starts = candidates[take]
        coeff = sess.transform(starts, n, hann(n, sym=True), [f], snapshots=True)
        coeff_by_f[f] = (coeff, starts, n)
        full, _ = estimate(coeff, xy, np.array([f]), base)
        rows.append(result_summary(f"{f:.3f}_all_six", full, names))
        common = full["music_source_count"] == 1
        for dropped in names:
            keep = np.array([j for j, name in enumerate(names) if name != dropped], int)
            sub, _ = estimate(coeff[:, keep, :, :], xy[keep], np.array([f]), base)
            row = result_summary(f"{f:.3f}_leave_one_out", sub, names, dropped)
            paired = common & (sub["music_source_count"] == 1)
            row["paired_q1_windows"] = int(paired.sum())
            row["paired_median_azimuth_shift_deg"] = float(np.nanmedian(angular_distance(
                sub["azimuth_deg"][paired], full["azimuth_deg"][paired]))) if paired.any() else np.nan
            row["paired_median_velocity_change_fraction"] = float(np.nanmedian(
                sub["velocity_m_s"][paired]/full["velocity_m_s"][paired]-1)) if paired.any() else np.nan
            rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "real_leave_one_station_out.csv", index=False)

    # Grid comparison uses identical windows and the same covariance snapshots.
    gridrows = []
    fine_c = settings(dict(base, direction_step_deg=1.0, velocity_points=256))
    for f, (coeff, starts, n) in coeff_by_f.items():
        grid_ix = np.unique(np.linspace(0, coeff.shape[2]-1, min(12, coeff.shape[2]), dtype=int))
        zgrid = coeff[:, :, grid_ix, :]
        coarse, _ = estimate(zgrid, xy, np.array([f]), base)
        fine, _ = estimate(zgrid, xy, np.array([f]), fine_c)
        paired = (coarse["music_source_count"] == 1) & (fine["music_source_count"] == 1)
        gridrows.append(dict(frequency_hz=f, windows=len(grid_ix), coarse_step_deg=base["direction_step_deg"],
            coarse_velocity_points=base["velocity_points"], fine_step_deg=fine_c["direction_step_deg"],
            fine_velocity_points=fine_c["velocity_points"], coarse_accepted=int(coarse["quality"].sum()),
            fine_accepted=int(fine["quality"].sum()), paired_q1_windows=int(paired.sum()),
            paired_median_azimuth_difference_deg=float(np.nanmedian(angular_distance(
                coarse["azimuth_deg"][paired], fine["azimuth_deg"][paired]))) if paired.any() else np.nan,
            paired_median_velocity_change_fraction=float(np.nanmedian(
                fine["velocity_m_s"][paired]/coarse["velocity_m_s"][paired]-1)) if paired.any() else np.nan))
    pd.DataFrame(gridrows).to_csv(OUT / "real_grid_refinement.csv", index=False)

    # Track local MUSIC peaks across the six reconstructed low frequencies.
    peak_rows = []
    low_freqs = np.arange(.08, .131, .01)
    for lf in low_freqs:
        path, actual_f, info = get_direction_archive(float(lf), "phase_guided/direction")
        ln = int(info["window_samples"])
        starts0 = info["start_sample"]
        valid = valid_windows(starts0, ln, masks, names)
        possible = starts0[valid]
        if not len(possible):
            continue
        # Spread windows over the four days, avoiding over-weighting one time block.
        pick = np.unique(np.linspace(0, len(possible)-1, min(5, len(possible)), dtype=int))
        use_starts = possible[pick]
        zall = sess.transform(use_starts, ln, hann(ln, sym=True), [actual_f], snapshots=True)
        cnow = settings(dict(base, frequency_half_width_hz=0., frequency_bins=1))
        for wi, start in enumerate(use_starts):
            res, ex = estimate(zall[:, :, wi:wi+1, :], xy, np.array([actual_f]), cnow)
            q = int(res["music_source_count"][0])
            if q == 0 or not ex:
                peak_rows.append(dict(frequency_hz=float(actual_f), start_sample=int(start),
                    source_count=q, accepted=bool(res["quality"][0]), peak_rank=0,
                    azimuth_deg=np.nan, velocity_m_s=np.nan, relative_power=np.nan,
                    reason=str(res["reason"][0])))
                continue
            peaks = local_peaks(ex[0]["power"], ex[0]["azimuth"], ex[0]["velocity"])
            for rank, peak in enumerate(peaks, 1):
                peak_rows.append(dict(frequency_hz=float(actual_f), start_sample=int(start),
                    source_count=q, accepted=bool(res["quality"][0]), peak_rank=rank,
                    azimuth_deg=peak["azimuth_deg"], velocity_m_s=peak["velocity_m_s"],
                    relative_power=peak["relative_power"], reason=str(res["reason"][0])))
    peaks_df = pd.DataFrame(peak_rows)
    peaks_df.to_csv(OUT / "real_multipeak_observations.csv", index=False)
    good_peaks = peaks_df[(peaks_df.peak_rank > 0) & peaks_df.azimuth_deg.notna()].copy()
    if len(good_peaks):
        good_peaks["direction_bin_15deg"] = (np.floor((good_peaks.azimuth_deg+7.5)/15).astype(int) % 24)
        summary = good_peaks.groupby(["frequency_hz", "direction_bin_15deg"], as_index=False).agg(
            windows=("start_sample", "nunique"), median_azimuth_deg=("azimuth_deg", "median"),
            median_velocity_m_s=("velocity_m_s", "median"), median_relative_power=("relative_power", "median"),
            rank1_windows=("peak_rank", lambda x: int((x == 1).sum())))
        summary.to_csv(OUT / "real_multipeak_tracks.csv", index=False)
        primary = good_peaks[good_peaks.peak_rank == 1].copy()
        if len(primary):
            primary["direction_bin_15deg"] = (np.floor((primary.azimuth_deg+7.5)/15).astype(int) % 24)
            modes = primary.groupby(["frequency_hz","direction_bin_15deg"],as_index=False).agg(
                primary_windows=("start_sample","nunique"), median_azimuth_deg=("azimuth_deg","median"))
            total = primary.groupby("frequency_hz").start_sample.nunique().rename("primary_windows_total")
            modes = modes.join(total,on="frequency_hz")
            modes["primary_share"] = modes.primary_windows/modes.primary_windows_total
            modes = modes.sort_values(["frequency_hz","primary_windows"],ascending=[True,False])
            modes.groupby("frequency_hz",as_index=False).head(1).to_csv(
                OUT / "real_multipeak_frequency_summary.csv", index=False)
    else:
        pd.DataFrame(columns=["frequency_hz","direction_bin_15deg","windows","median_azimuth_deg",
                              "median_velocity_m_s","median_relative_power","rank1_windows"]).to_csv(
                                  OUT / "real_multipeak_tracks.csv", index=False)
        pd.DataFrame().to_csv(OUT / "real_multipeak_frequency_summary.csv", index=False)
    print(f"real windows={len(starts)} at {f:g} Hz; low-frequency peak rows={len(peaks_df)}")


def circular_mean_deg(values):
    a = np.radians(np.asarray(values, float))
    return float(np.degrees(np.angle(np.mean(np.exp(1j*a)))) % 360.0)


def run_bootstrap():
    import re
    summary = pd.read_csv(BASE / "phase_guided/phase_stack_summary.csv")
    want = summary[(summary.curve_id == "P00038_b1_S0_N1") &
                   summary.frequency_hz.between(.079, .131)].sort_values("frequency_hz")
    if len(want) != 6:
        raise RuntimeError(f"Expected 6 low-frequency observations, got {len(want)}")
    geometry = pd.read_csv(BASE / "phase_guided/direction/station_geometry.csv").set_index("station")
    block_starts = np.arange(16, dtype=int)
    block_ids = [f"2023-11-{day:02d}T{hour:02d}" for day in range(1,5) for hour in [0,6,12,18]]
    window_data = []
    x0 = np.datetime64("2023-10-31T23:59:58.470", "ms")
    for row in want.to_dict("records"):
        p = BASE / "phase_guided" / row["windows_file"]
        with np.load(p, allow_pickle=False) as z:
            d = {k: z[k] for k in z.files}
        pair = str(d["pair"])
        a, b = pair.split("--")
        delta = (geometry.loc[b, ["x_m", "y_m"]].to_numpy(float) -
                 geometry.loc[a, ["x_m", "y_m"]].to_numpy(float))
        sec = d["start_sample"].astype(float) / float(d["fs_hz"])
        dates = x0 + (sec * 1000).astype("timedelta64[ms]")
        day_hour = []
        for dt in dates:
            stamp = pd.Timestamp(str(dt)).tz_localize("UTC")
            day_hour.append((stamp.date().isoformat(), int(stamp.hour//6)*6))
        index = np.full(len(sec), -1, int)
        block_map = {key:i for i,key in enumerate(block_ids)}
        for i, (date, hour) in enumerate(day_hour):
            key = (f"2023-11-{int(date[-2:]):02d}T{hour:02d}" if date.startswith("2023-11-") else "")
            if key in block_map:
                index[i] = block_map[key]
        d.update(delta=delta, block_index=index)
        window_data.append(d)

    rng = np.random.default_rng(BOOT_SEED)
    draws = rng.integers(0, len(block_ids), size=(BOOT_REPS, len(block_ids)), endpoint=False)
    counts = np.stack([np.bincount(r, minlength=len(block_ids)) for r in draws])
    boot_rows = []
    fvals = want.frequency_hz.to_numpy(float)
    original_projection = want.projected_distance_m.to_numpy(float)
    original_phases = np.radians(want.after_raw_phase_deg.to_numpy(float))
    original_anchors = want.predicted_unwrapped_phase_rad.to_numpy(float)
    original_velocity = phase_speed(fvals,
        lift_band(original_phases, original_anchors)[0], original_projection)
    original_az = want.confirmed_azimuth_deg.to_numpy(float)

    for di, mult in enumerate(counts):
        phases = np.full(6, np.nan)
        azs = np.full(6, np.nan)
        projections = np.full(6, np.nan)
        accepted_counts = np.zeros(6, int)
        unique_counts = np.zeros(6, int)
        anchors = np.full(6, np.nan)
        for fi, d in enumerate(window_data):
            block = d["block_index"]
            chosen = np.flatnonzero((block >= 0) & (mult[np.maximum(block,0)] > 0))
            if not len(chosen):
                continue
            weights = mult[block[chosen]]
            elig = d["direction_valid"][chosen].astype(bool)
            az = d["propagation_azimuth_deg"][chosen].astype(float)
            elig &= np.isfinite(az)
            if not elig.any():
                continue
            bins = (np.floor((az[elig] + 7.5)/15).astype(int) % 24)
            weighted_counts = np.bincount(bins, weights=weights[elig], minlength=24)
            vals = np.flatnonzero(weighted_counts)
            nvals = weighted_counts[vals]
            mode_bin = int(vals[np.argmax(nvals)])
            az_bins = np.full(len(az), -1, dtype=int)
            az_bins[elig] = np.floor((az[elig] + 7.5)/15).astype(int) % 24
            group = elig & (az_bins == mode_bin)
            theta = circular_mean_deg(np.repeat(az[group], weights[group]))
            dx, dy = d["delta"]
            axis = np.array([np.sin(np.radians(theta)), np.cos(np.radians(theta))])
            ref = float(np.dot(d["delta"], axis))
            if ref <= 0:
                continue
            projection_window = dx*np.sin(np.radians(az)) + dy*np.cos(np.radians(az))
            # Reapply the stored pair-angle decision, chosen directional group and branch sign.
            eligible = elig & group & (projection_window > np.linalg.norm(d["delta"])*.1)
            eligible &= np.isfinite(d["raw_cross"][chosen].real) & np.isfinite(d["raw_cross"][chosen].imag)
            if np.unique(chosen[eligible]).size < 3:
                continue
            f = float(d["frequency_hz"])
            prior_v = float(d["prior_true_velocity_m_s"])
            aligned, _ = align_cross(d["raw_cross"][chosen], f, projection_window, ref, prior_v)
            center = -2*np.pi*f*ref/prior_v
            low = float(d["phase_offset_low_rad"])
            high = float(d["phase_offset_high_rad"])
            scr = screen_cross(aligned, center, low, high, "curve", True, False, eligible)
            accepted = scr["accepted"]
            if not accepted.any():
                continue
            # Count original windows only; repeated blocks are weights, not new independent windows.
            unique_counts[fi] = int(np.unique(chosen[accepted]).size)
            accepted_counts[fi] = int(np.sum(weights[accepted]))
            if unique_counts[fi] < 3:
                continue
            phase = float(np.angle(np.sum(aligned[accepted] * weights[accepted])))
            phases[fi] = phase
            azs[fi] = theta
            projections[fi] = ref
            anchors[fi] = center
        if np.isfinite(phases).all() and np.isfinite(anchors).all() and np.isfinite(projections).all():
            try:
                unwrapped, shift = lift_band(phases, anchors)
                velocity = phase_speed(fvals, unwrapped, projections)
                if not np.isfinite(velocity).all() or (velocity <= 0).any():
                    velocity[:] = np.nan
            except Exception:
                velocity = np.full(6, np.nan)
        else:
            velocity = np.full(6, np.nan)
        monotone = bool(np.isfinite(velocity).all() and np.all(np.diff(velocity) <= 1e-8))
        row = dict(draw=di, all_six_valid=bool(np.isfinite(velocity).all()), monotone_nonincreasing=monotone,
                   accepted_counts_json=json.dumps(accepted_counts.tolist()),
                   unique_accepted_json=json.dumps(unique_counts.tolist()))
        for i, f in enumerate(fvals):
            row[f"azimuth_{f:.2f}"] = azs[i]
            row[f"projection_{f:.2f}"] = projections[i]
            row[f"phase_{f:.2f}"] = phases[i]
            row[f"velocity_{f:.2f}"] = velocity[i]
            row[f"n_unique_{f:.2f}"] = unique_counts[i]
        for i in range(5):
            row[f"delta_v_{fvals[i]:.2f}_{fvals[i+1]:.2f}"] = velocity[i+1]-velocity[i] if np.isfinite(velocity).all() else np.nan
        boot_rows.append(row)
    df = pd.DataFrame(boot_rows)
    np.savez_compressed(OUT / "bootstrap_common_blocks.npz", draws=draws,
        block_ids=np.asarray(block_ids), seed=BOOT_SEED, bootstrap_replicates=BOOT_REPS)
    df.to_csv(OUT / "low_frequency_joint_bootstrap_draws.csv", index=False)

    summary_rows = []
    full = df[df.all_six_valid]
    for i, row in enumerate(want.to_dict("records")):
        f = fvals[i]
        key = f"velocity_{f:.2f}"
        good = df[key].dropna().to_numpy(float)
        azkey = f"azimuth_{f:.2f}"
        pkey = f"projection_{f:.2f}"
        azdraw = df[azkey].dropna().to_numpy(float)
        pdraw = df[pkey].dropna().to_numpy(float)
        summary_rows.append(dict(frequency_hz=f, original_velocity_m_s=original_velocity[i],
            original_azimuth_deg=original_az[i], original_projected_distance_m=original_projection[i],
            valid_velocity_draws=len(good), velocity_q025_m_s=float(np.quantile(good,.025)) if len(good) else np.nan,
            velocity_q50_m_s=float(np.quantile(good,.5)) if len(good) else np.nan,
            velocity_q975_m_s=float(np.quantile(good,.975)) if len(good) else np.nan,
            azimuth_q025_deg=float(np.quantile(azdraw,.025)) if len(azdraw) else np.nan,
            azimuth_q975_deg=float(np.quantile(azdraw,.975)) if len(azdraw) else np.nan,
            projection_q025_m=float(np.quantile(pdraw,.025)) if len(pdraw) else np.nan,
            projection_q975_m=float(np.quantile(pdraw,.975)) if len(pdraw) else np.nan))
    pd.DataFrame(summary_rows).to_csv(OUT / "low_frequency_joint_bootstrap_summary.csv", index=False)
    diff_rows = []
    for i in range(5):
        key = f"delta_v_{fvals[i]:.2f}_{fvals[i+1]:.2f}"
        vals = df[key].dropna().to_numpy(float)
        diff_rows.append(dict(frequency_low_hz=fvals[i], frequency_high_hz=fvals[i+1],
            bootstrap_draws=len(vals), delta_velocity_q025_m_s=float(np.quantile(vals,.025)) if len(vals) else np.nan,
            delta_velocity_q50_m_s=float(np.quantile(vals,.5)) if len(vals) else np.nan,
            delta_velocity_q975_m_s=float(np.quantile(vals,.975)) if len(vals) else np.nan,
            observed_original_delta_m_s=float(original_velocity[i+1]-original_velocity[i]),
            bootstrap_positive_fraction=float(np.mean(vals>0)) if len(vals) else np.nan))
    pd.DataFrame(diff_rows).to_csv(OUT / "low_frequency_adjacent_differences.csv", index=False)
    print(f"joint block bootstrap draws={len(df)}, complete six-point draws={int(df.all_six_valid.sum())}")


def run_holdout():
    """Fixed train/holdout date comparison; no settings are selected on Nov 4."""
    arrays, masks, names, xy, fs, origin = load_actual()
    base = baseline_settings()
    origin_ts = pd.Timestamp(origin)
    variants = [(3.5, 4), (3.5, 6), (3.5, 7), (4.0, 6)]
    out = []
    for target, folder in [(.10, "phase_guided/direction"), (.45, "music_front/direction")]:
        _, f, meta = get_direction_archive(target, folder)
        n = int(meta["window_samples"])
        starts0 = meta["start_sample"]
        valid = valid_windows(starts0, n, masks, names)
        starts0 = starts0[valid]
        stamps = origin_ts + pd.to_timedelta(starts0/fs, unit="s")
        days = np.asarray(stamps.strftime("%Y-%m-%d"))
        splits = {"train_Nov01_03": np.isin(days, ["2023-11-01","2023-11-02","2023-11-03"]),
                  "holdout_Nov04": days == "2023-11-04"}
        for split, mask in splits.items():
            candidates = starts0[mask]
            if not len(candidates):
                continue
            pick = np.unique(np.linspace(0, len(candidates)-1, min(16,len(candidates)), dtype=int))
            chosen = candidates[pick]
            for nw, kmax in variants:
                c = settings(dict(base, music_time_bandwidth=nw, music_snapshots=kmax))
                sess = session_for(arrays, names, xy, fs, c)
                z = sess.transform(chosen, n, hann(n, sym=True), [f], snapshots=True)
                for threshold in [.70, .98]:
                    current = settings(dict(c, min_coherence=threshold))
                    result, _ = estimate(z, xy, np.array([f]), current)
                    q1 = result["music_source_count"] == 1
                    usable = q1 & result["quality"]
                    out.append(dict(frequency_hz=float(f), split=split,
                        NW=nw, snapshots=kmax, min_coherence=threshold,
                        half_bandwidth_hz=nw/(200.0/f), windows=len(chosen),
                        accepted=int(result["quality"].sum()),
                        accepted_fraction=float(result["quality"].mean()),
                        q0=int((result["music_source_count"]==0).sum()),
                        q1=int(q1.sum()), q2plus=int((result["music_source_count"]>=2).sum()),
                        median_coherence=float(np.nanmedian(result["coherence"][q1])) if q1.any() else np.nan,
                        median_accepted_azimuth_deg=float(np.nanmedian(result["azimuth_deg"][usable])) if usable.any() else np.nan,
                        median_accepted_velocity_m_s=float(np.nanmedian(result["velocity_m_s"][usable])) if usable.any() else np.nan,
                        settings_promoted=False))
    frame = pd.DataFrame(out)
    frame.to_csv(OUT / "holdout_date_validation.csv", index=False)
    print(f"holdout rows={len(frame)}; 2023-11-04 evaluated as untouched date")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=["simulation", "real", "bootstrap", "holdout", "all"], default="all")
    args = ap.parse_args()
    if args.part in ("simulation", "all"):
        run_simulation()
    if args.part in ("real", "all"):
        run_real()
    if args.part in ("bootstrap", "all"):
        run_bootstrap()
    if args.part in ("holdout", "all"):
        run_holdout()


if __name__ == "__main__":
    main()
