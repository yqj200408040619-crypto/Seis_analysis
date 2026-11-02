"""Step 4+5: predicted-arrival-gated short-window CCF stacking at 5-10 s.

Uses the fitted/extrapolated prior curve Vph(f) (step 3) as the PREDICTED arrival. For each
station pair and each target frequency (0.1-0.2 Hz), we:
  - cut short Goertzel windows (N cycles) across the whole record (program geometry),
  - form the cross-spectrum C = G_b * conj(G_a) per window,
  - predict the inter-station phase phi_pred = 2*pi*f*d/v_pred (both lag branches),
  - KEEP a window only if its observed cross-phase is within the tolerance that corresponds
    to a +/-10% velocity error at that f and distance (per-pair, per-frequency phase margin),
  - complex-stack the kept windows, then invert the stacked phase to a measured velocity,
    choosing the cycle integer nearest the predicted velocity.

This directly realises the user's design. The prior is a-priori (HF measured + base-model
anchor); the recovered low-frequency velocity is then compared to the prior/authoritative.

Reuses the program's response-corrected alignment (align_project) and exact Goertzel
(goertzel_windows); no new signal engine.
"""
import sys, json, os, argparse, itertools, traceback
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import ndstudio
from dispersion.common import load_config, parse_utc
from dispersion.prepare import station_table
from ndstudio.phase_stack import align_project, goertzel_windows
from ndstudio.phase_records import available_phase_extent, phase_window_geometry, full_window_starts
from scipy.signal.windows import hann
from pyproj import Geod

GEOD = Geod(ellps='WGS84')


def wrap(x):
    return np.angle(np.exp(1j*np.asarray(x)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--project', default=str(HERE/'runs'/'core_la_14d'/'project'))
    ap.add_argument('--prior', default=str(HERE/'runs'/'hf_core'/'fitted_prior_curve.csv'))
    ap.add_argument('--out', default=str(HERE/'runs'/'lf_gated_stack'))
    ap.add_argument('--fmin', type=float, default=0.10)
    ap.add_argument('--fmax', type=float, default=0.20)
    ap.add_argument('--fstep', type=float, default=0.01)
    ap.add_argument('--cycles', type=float, default=12.0)   # short windows (~N periods)
    ap.add_argument('--slide-cycles', type=float, default=4.0)
    ap.add_argument('--vel-tol', type=float, default=0.10)  # +/-10% velocity gate
    a = ap.parse_args()

    root = Path(a.project).resolve()
    out = Path(a.out).resolve(); out.mkdir(parents=True, exist_ok=True)
    prior = pd.read_csv(a.prior).sort_values('frequency_hz')
    pf = prior.frequency_hz.to_numpy(float); pv = prior.phase_velocity_m_s.to_numpy(float)
    def v_pred(f):
        return float(np.interp(f, pf, pv))

    c = load_config(root/'config.json'); fs = float(c['sample_rate_hz'])
    table, _, _ = station_table(c); xy = table.set_index('station')
    names = list(table.station)

    # Align all stations on the shared UTC grid (reuse existing aligned H5 if present).
    first, last, bounds = available_phase_extent(root, set(names))
    aligned_h5 = out/'aligned_waveforms.h5'
    import tempfile
    with tempfile.TemporaryDirectory(prefix='.lf_work_', dir=out) as work:
        resources = []
        reuse = root/'output'/'music_front'/'aligned_waveforms.h5'
        src = reuse if reuse.exists() else None
        if src is None:
            arrays, fsa, origin, meta = align_project(root, set(names), str(first), str(last), aligned_h5, work, resources)
        else:
            import h5py
            with h5py.File(src) as h:
                fsa = float(h.attrs['fs_hz']); origin = parse_utc(h.attrs['start_utc'])
                arrays = {s: np.where(h['stations/'+s+'/valid'][:], h['stations/'+s+'/data'][:], np.nan) for s in names}
        fs = fsa

        # station coords in local meters (use geometry x/y if present, else from lat/lon)
        latlon = {r.station: (r.latitude, r.longitude) for _, r in table.iterrows()}
        freqs = np.round(np.arange(a.fmin, a.fmax+1e-9, a.fstep), 6)
        pairs = list(itertools.combinations(names, 2))
        per_pair_rows = []; agg = {}
        for (sa, sb) in pairs:
            la, lo = latlon[sa]; lb, lob = latlon[sb]
            _, _, dist = GEOD.inv(lo, la, lob, lb)
            for f in freqs:
                vP = v_pred(float(f))
                geo = phase_window_geometry(float(f), fs, a.cycles, a.slide_cycles)
                n = geo['window_samples']; step = geo['slide_samples']
                L = len(arrays[sa])
                if n > L:
                    continue
                starts = full_window_starts(L, n, step)
                taper = hann(n, sym=True)
                ga = goertzel_windows(np.ascontiguousarray(arrays[sa]), starts, n, float(f), fs, taper)
                gb = goertzel_windows(np.ascontiguousarray(arrays[sb]), starts, n, float(f), fs, taper)
                cross = gb*np.conj(ga)
                good = np.isfinite(cross) & (np.abs(cross) > 0)
                if good.sum() < 10:
                    continue
                cross = cross[good]
                phase = np.angle(cross)
                # predicted inter-station phase magnitude at v_pred; both lag branches +/-
                phi = 2*np.pi*float(f)*dist/vP
                # phase margin for +/-vel_tol velocity: d(phi)/phi = -dv/v -> dphi = phi*vel_tol
                margin = phi*a.vel_tol
                margin = float(np.clip(margin, np.radians(5.), np.radians(60.)))  # keep sane bounds
                # keep windows whose wrapped phase is near +phi or -phi (either propagation sense)
                keep_pos = np.abs(wrap(phase - phi)) <= margin
                keep_neg = np.abs(wrap(phase + phi)) <= margin
                keep = keep_pos | keep_neg
                kept = int(keep.sum())
                ngood = int(good.sum())
                # --- UNBIASED coherence diagnostics (NOT gated by the prior) ---
                # (1) ungated complex stack: does ANY coherent inter-station phase exist?
                ung = cross.sum(); ung_res = abs(ung)/ngood; ung_phase = float(np.angle(ung))
                # velocity implied by the UNGATED stacked phase (cycle nearest prior)
                Nu = round((phi - ung_phase)/(2*np.pi)); phiu = ung_phase + 2*np.pi*Nu
                v_ungated = 2*np.pi*float(f)*dist/phiu if phiu > 0 else np.nan
                # unit-amplitude (phase-only) ungated resultant = phase concentration
                ung_unit = abs((cross/np.abs(cross)).sum())/ngood
                # (2) retention vs uniform-phase expectation (two arcs of half-width margin)
                exp_frac = 2*(2*margin)/(2*np.pi)
                ret_frac = kept/ngood
                if kept < 5:
                    per_pair_rows.append(dict(pair=f'{sa}--{sb}', frequency_hz=float(f), distance_m=dist,
                                              v_pred_m_s=vP, windows=ngood, kept=kept, v_meas_m_s=np.nan,
                                              resultant_length=np.nan, ungated_resultant=ung_res,
                                              ungated_unit_resultant=ung_unit, v_ungated_m_s=v_ungated,
                                              retention_fraction=ret_frac, uniform_expected_fraction=exp_frac,
                                              retention_ratio=ret_frac/exp_frac if exp_frac>0 else np.nan, note='too_few_kept'))
                    continue
                # gated stack (fold negative branch), also report its UNIT resultant (coherence)
                z = np.where(keep_neg & ~keep_pos, np.conj(cross), cross)
                stacked = z[keep].sum()
                stacked_phase = np.angle(stacked)
                resultant = abs(stacked)/kept
                unit_resultant = abs((z[keep]/np.abs(z[keep])).sum())/kept
                N = round((phi - stacked_phase)/(2*np.pi))
                phi_full = stacked_phase + 2*np.pi*N
                v_meas = 2*np.pi*float(f)*dist/phi_full if phi_full > 0 else np.nan
                per_pair_rows.append(dict(pair=f'{sa}--{sb}', frequency_hz=float(f), distance_m=dist,
                                          v_pred_m_s=vP, windows=ngood, kept=kept,
                                          v_meas_m_s=float(v_meas) if np.isfinite(v_meas) else np.nan,
                                          resultant_length=float(resultant), unit_resultant=float(unit_resultant),
                                          ungated_resultant=float(ung_res), ungated_unit_resultant=float(ung_unit),
                                          v_ungated_m_s=float(v_ungated) if np.isfinite(v_ungated) else np.nan,
                                          retention_fraction=float(ret_frac), uniform_expected_fraction=float(exp_frac),
                                          retention_ratio=float(ret_frac/exp_frac) if exp_frac>0 else np.nan,
                                          cycle_N=int(N), note=''))
        for item in reversed(resources):
            try: item._mmap.close()
            except Exception:
                try: item.close()
                except Exception: pass

    pp = pd.DataFrame(per_pair_rows)
    pp.to_csv(out/'gated_pair_frequency.csv', index=False)
    # UNBIASED recovery: use the UNGATED stack (not seeded by the prior). Keep pairs whose
    # raw per-window cross-phase is coherent (ungated unit resultant above a noise floor) and
    # resolvable (d/lambda in a sane band). Compare recovered v to BOTH prior and authoritative.
    vm = pp[pp.v_ungated_m_s.notna()].copy()
    vm['d_over_lambda'] = vm.distance_m*vm.frequency_hz/vm.v_ungated_m_s
    coh = vm[(vm.ungated_unit_resultant >= 0.15) & (vm.windows >= 50) &
             (vm.d_over_lambda >= 0.15) & (vm.d_over_lambda <= 6.0)].copy()
    recov = coh.groupby('frequency_hz').agg(
        v_ungated_median_m_s=('v_ungated_m_s', 'median'), v_ungated_std_m_s=('v_ungated_m_s', 'std'),
        ungated_unit_resultant_median=('ungated_unit_resultant', 'median'),
        v_pred_m_s=('v_pred_m_s', 'first'), n_pairs=('pair', 'nunique')).reset_index()
    recov['period_s'] = 1.0/recov.frequency_hz
    recov['rel_err_vs_prior'] = (recov.v_ungated_median_m_s - recov.v_pred_m_s)/recov.v_pred_m_s
    recov.to_csv(out/'recovered_lowfreq_velocity.csv', index=False)
    # gated (prior-seeded) result kept separately for contrast + to show the circularity guard
    gated = pp[(pp.v_meas_m_s.notna()) & (pp.kept >= 8)].groupby('frequency_hz').agg(
        v_gated_median_m_s=('v_meas_m_s', 'median'),
        gated_unit_resultant_median=('unit_resultant', 'median'),
        retention_ratio_median=('retention_ratio', 'median')).reset_index()
    gated.to_csv(out/'gated_prior_seeded_contrast.csv', index=False)
    summary = dict(
        target_band_hz=[a.fmin, a.fmax], cycles=a.cycles, vel_tol=a.vel_tol,
        pairs=len(pairs), frequencies=list(map(float, freqs)),
        recovered_UNGATED=recov.round(4).to_dict('records'),
        gated_prior_seeded_contrast=gated.round(4).to_dict('records'),
        interpretation=('UNGATED is the honest recovery (not seeded by the prior). '
                        'ungated_unit_resultant is the per-window cross-phase concentration '
                        '(coherence); near 0 => no coherent LF signal (prior-gated v is then circular). '
                        'retention_ratio_median ~1 => kept fraction only matches uniform-random => no signal.'),
        note='v_meas (gated) is prior-seeded and must NOT be read as independent recovery.')
    (out/'lf_gated_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=str), encoding='utf-8')
    print('done')


if __name__ == '__main__':
    main()
