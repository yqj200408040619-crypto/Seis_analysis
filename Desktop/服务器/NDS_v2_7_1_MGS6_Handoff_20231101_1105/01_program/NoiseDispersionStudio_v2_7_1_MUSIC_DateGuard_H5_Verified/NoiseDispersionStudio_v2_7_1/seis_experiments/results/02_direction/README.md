# 02 — Direction & velocity (extended aperture + dense core)

Per-frequency MUSIC direction results. The stable direction only emerged under the
USER-AUTHORIZED relaxed scan `validation` gate (see relaxed_validation_used.json); these are
non-default, non-calibrated estimates.

- `extended_geometry.json`           — extended 11-stn array, aperture 62.3 km.
- `relaxed_validation_used.json`     — the relaxed scan-gate thresholds that were applied.
- `extended_frequency_directions.csv`— extended array per-frequency MUSIC directions.
- `stable_direction_summary.json`    — extended stable band 0.10-0.115 Hz: az 32-46 deg, Vph 2.5-2.7 km/s.
- `core_frequency_directions.csv`    — dense 30 km core per-frequency MUSIC directions.
- `core_refine_summary.json`         — core refinement 0.115-0.125 Hz: az 27-31 deg, Vph 2.76-2.95 km/s.

Cross-validation: the 62 km array stabilizes the long periods (8.7-10 s); the 30 km core
refines toward 5 s. Both agree on NE propagation (~30-37 deg) and ~2.5-2.95 km/s Rayleigh Vph.
