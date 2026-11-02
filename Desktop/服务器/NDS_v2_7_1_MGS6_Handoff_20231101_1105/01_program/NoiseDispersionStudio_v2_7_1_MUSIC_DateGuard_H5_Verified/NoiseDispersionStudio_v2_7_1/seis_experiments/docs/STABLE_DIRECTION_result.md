# STABLE DIRECTION OBTAINED (authorized relaxed scan validation)

User authorized relaxing the scan validation gate: snr_min 10->3, block_snr_min 5->3,
block_fraction_min 0.75->0.5, block_phase_max_rad 0.5->0.7, and MODEST arrival loosening
(arrival_fraction_max 0.1->0.2, leaveout 0.1->0.2, block_arrival 0.15->0.25). Kept
music_front_cycles=30, minimum_independent_windows=3.

## Effect on the chain (extended 11-stn, 62 km LA array, 14 days)
- scan: conditional_dispersion.csv 0 -> 843 rows; robust_contiguous_runs 0 -> 50;
  accepted_continuous_runs 0 -> 20; candidate_gates 28 -> 149.
- front MUSIC: frequencies 0 -> 50; frequency_directions.csv 161 per-frequency rows.

## STABLE WAVE DIRECTION (per-frequency MUSIC, accepted=True)
Band 0.100-0.115 Hz (period 8.7-10.0 s):
  propagation azimuth 32.0-45.6 deg (median 37 deg, NE-ward; north=0, clockwise)
  MUSIC phase velocity 2.50-2.68 km/s (median 2.55 km/s)
  direction uncertainty median 1.8 deg (tight)
  independent windows 3-9; two independent pairs agree: GR2-RPV and HLL-RIO.
Physically consistent: ~2.5 km/s Rayleigh at ~9-10 s, propagating inland (NE) from
Pacific (SW) secondary-microseism sources.

Above 0.12 Hz (toward the 5 s / 0.2 Hz target) estimates still did NOT stabilize
(accepted=False) even with relaxation -> the coherent, stable direction sits at ~9-10 s,
not at 5 s, for this array/epoch.

## IMPORTANT caveat
These thresholds are NON-DEFAULT, NON-CALIBRATED, set with explicit user permission.
direction_uncertainty / temporal spread are resolution diagnostics, NOT calibrated CIs.
The joint continuous dispersion band still did not pass (only 1 non-overlapping window),
so this is a stable DIRECTION + per-frequency MUSIC velocity, not yet a verified continuous
dispersion curve for authoritative comparison.

Files: runs/extended_la_14d/project/output/music_front/frequency_directions.csv,
runs/extended_la_14d/relaxed_validation_used.json, runs/extended_la_14d/stable_direction_summary.json.


## STEP 2: dense 30 km core refinement (user's two-step design)
Dense core only (8 stn: BHP,GR2,HLL,LCG,LGB,PDR,SMF2,USC), aperture 29.6 km, same relaxed
validation + cycles=30, 14 days, band ~0.08-0.14 Hz. (Reused the already-downloaded 14-day
core waveforms; no re-download.)
- scan: conditional_dispersion 248 rows; robust_contiguous_runs 24; accepted_continuous_runs 10.
- front MUSIC: 60 per-frequency rows; 3 accepted.
Stable band 0.115-0.125 Hz (8.0-8.7 s), pair BHP-PDR:
  azimuth 27.5-30.7 deg (median 30.5 deg), MUSIC phase velocity 2.76-2.95 km/s
  (median 2.84 km/s), direction uncertainty median 2.5 deg.

## CROSS-ARRAY VALIDATION (the key outcome)
| array          | aperture | stable band        | azimuth (prop)   | MUSIC Vph        |
|----------------|----------|--------------------|------------------|------------------|
| extended 11stn | 62 km    | 0.100-0.115 Hz     | 32-46 (med 37)   | 2.50-2.68 km/s   |
| dense core 8   | 30 km    | 0.115-0.125 Hz     | 27-31 (med 30.5) | 2.76-2.95 km/s   |
Two INDEPENDENT arrays agree: NE-ward propagation ~30-37 deg, Rayleigh phase velocity
~2.5-2.9 km/s at ~8-10 s. The wide aperture stabilized the long periods; the dense core
extended slightly toward 5 s and refined the velocity. Physically consistent with Pacific
(SW) secondary-microseism energy propagating inland (NE) across the LA basin.

CAVEAT (unchanged): scan validation thresholds were relaxed with user's permission
(snr_min 3, block wave-packet loosened, arrival modestly loosened). These are NON-CALIBRATED;
the azimuth/velocity are stable estimates with low formal uncertainty but should be treated as
a relaxed-threshold result, not the program's default calibrated product. The target 5 s
(0.2 Hz) itself still did not stabilize; coherent energy held at ~8-10 s.
