#!/usr/bin/env python3
"""Isolated counterfactual rebuild of MUSIC-front references under coherence floors."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent / "project"
WORKSPACE = HERE.parents[1]
SOURCE = WORKSPACE / "source_2023" / "NoiseDispersionStudio_v2_7_1"
sys.path[:0] = [str(SOURCE), str(SOURCE / "core")]
from dispersion.common import load_config
from ndstudio.music_front import direction_settings, frequency_reference, correct_candidates
from ndstudio.phase_records import phase_window_geometry


def main():
    c = load_config(ROOT / "config.json")
    dc = direction_settings(ROOT)
    curves = pd.read_csv(ROOT / "output/conditional_dispersion.csv")
    pairs = pd.read_csv(ROOT / "output/pairs.csv").set_index("pair_index")
    stations = pd.read_csv(ROOT / "output/music_front/direction/station_geometry.csv").set_index("station")
    data = pd.read_csv(HERE / "rejection_condition_windows.csv")
    output, band_rows, ref_rows = [], [], []
    for threshold in [.40, .50, .60, .65, .70]:
        refs = {}
        for f, q in data.groupby("frequency_hz", sort=True):
            q = q.sort_values("start_sample")
            f = float(f)
            geometry = phase_window_geometry(f, c["sample_rate_hz"], dc["music_front_cycles"],
                                             dc["music_front_cycles"])
            n = int(geometry["window_samples"])
            other_pass = ~q[["phase_residual_fail_fail", "uncertain_fail", "ambiguous_fail",
                             "boundary_fail", "grid_resolution_fail"]].astype(bool).any(axis=1)
            quality = ((q.source_count.to_numpy(int) == 1) & other_pass.to_numpy(bool) &
                       (q.coherence.to_numpy(float) >= threshold))
            az = q.azimuth_deg.to_numpy(float)
            groups = np.floor((np.nan_to_num(az) + dc["direction_group_width_deg"] / 2) /
                              dc["direction_group_width_deg"]).astype(int) % int(round(360 / dc["direction_group_width_deg"]))
            result = dict(quality=quality, azimuth_deg=az,
                          direction_group=groups,
                          direction_uncertainty_deg=q.uncertainty_deg.to_numpy(float),
                          velocity_m_s=q.velocity_m_s.to_numpy(float))
            starts = q.start_sample.to_numpy(np.int64)
            for (pair, branch), _ in curves[np.isclose(curves.frequency_hz, f)].groupby(["pair", "branch"]):
                a, b = pair.split("--")
                delta = stations.loc[b, ["x_m", "y_m"]].to_numpy(float) - stations.loc[a, ["x_m", "y_m"]].to_numpy(float)
                ref = frequency_reference(result, starts, n, delta, dc, branch, angle_gate=False)
                refs[(pair, branch, round(f, 9))] = ref
                ref_rows.append(dict(min_coherence=threshold, frequency_hz=f, pair=pair, branch=branch,
                                     accepted=bool(ref), independent_windows=ref["independent_windows"] if ref else 0,
                                     selected_direction_group=ref["selected_direction_group"] if ref else np.nan,
                                     azimuth_deg=ref["azimuth_deg"] if ref else np.nan,
                                     projected_distance_m=ref["projected_distance_m"] if ref else np.nan,
                                     direction_uncertainty_deg=ref["direction_uncertainty_deg"] if ref else np.nan,
                                     music_velocity_m_s=ref["music_velocity_m_s"] if ref else np.nan))
        final, audit, catalog = correct_candidates(curves, pairs, refs, c["scan"])
        accepted_refs = sum(bool(x) for x in refs.values())
        band_rows.append(dict(min_coherence=threshold, pair_frequency_references=len(refs),
                              accepted_frequency_references=accepted_refs,
                              corrected_points=len(final), accepted_audit_points=int(audit.accepted.sum()),
                              candidate_bands=len(catalog), longest_band_points=int(catalog.points.max()) if len(catalog) else 0,
                              longest_band_min_hz=float(catalog.loc[catalog.points.idxmax(), "fmin_hz"]) if len(catalog) else np.nan,
                              longest_band_max_hz=float(catalog.loc[catalog.points.idxmax(), "fmax_hz"]) if len(catalog) else np.nan,
                              settings_promoted=False))
        final.to_csv(HERE / f"coherence_front_{threshold:.2f}_corrected.csv", index=False)
        catalog.to_csv(HERE / f"coherence_front_{threshold:.2f}_catalog.csv", index=False)
        print(f"coherence={threshold:.2f}: refs={accepted_refs}/{len(refs)}, corrected={len(final)}, bands={len(catalog)}", flush=True)
    pd.DataFrame(ref_rows).to_csv(HERE / "coherence_front_counterfactual_references.csv", index=False)
    pd.DataFrame(band_rows).to_csv(HERE / "coherence_front_counterfactual_summary.csv", index=False)


if __name__ == "__main__":
    main()
