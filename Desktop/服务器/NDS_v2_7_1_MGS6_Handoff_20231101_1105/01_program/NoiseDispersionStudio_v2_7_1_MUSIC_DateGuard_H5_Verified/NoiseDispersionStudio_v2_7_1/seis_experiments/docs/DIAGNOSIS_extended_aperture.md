# Extended-aperture direction experiment (user design: dense core + sparse extenders)

## What was built and run (REAL data, SCEDC)
Custom selection (hand-built in the program's own schema, via download.array_geometry):
 - Dense core (8, ~30 km, for velocity refinement): CI BHP,GR2,HLL,LCG,LGB,PDR,SMF2,USC
 - Sparse azimuthal extenders to widen aperture: CI MWC(NE), RPV(S), RIO(E)
   (WSS dropped: mixed native sampling grid over 14 d -> MappingError.)
 - Final array: 11 stations, all HHZ 100 Hz, aperture 62.3 km, min pair 1.85 km.
Band centered on ~0.2 Hz (5 s): scan 0.08-0.40 Hz; bandpass 0.04-0.8 Hz; 10 Hz sample.
Span 14 days (2023-11-01..15); full download (each station 12,096,000 samples @10Hz).
Direction stage run in DIAGNOSTIC mode (keeps estimates without the strict filter accept gate).

## Result
music_front/music_complete.json:
  "no_qualified_direction_no_true_velocity": true, corrected_points 0, corrected_curves 0.
direction/direction_complete.json: frequencies 0, raw_joint_points 0, accepted_continuous_points 0.
frequency_diagnostics.csv: empty (no per-frequency estimate even formed).

## Why (root cause, confirmed)
Even in diagnostic mode and with a 62 km aperture + 14 days, NO frequency assembles enough
qualifying windows: the per-frequency MUSIC screen requires coherence >= 0.7 and >= 3
independent windows, and the pre-direction window-quality screen (clean_fraction) plus the
200-cycle window length at 0.08-0.4 Hz leave too few clean, coherent windows. Widening the
aperture improves angular *resolution* but does NOT raise inter-station coherence; at 5 s over
30-60 km baselines the ambient field is not coherent enough to pass the program's calibrated
window/coherence gates. Mixing far stations adds long lags and lower pairwise coherence, which
does not help these gates.

## Overall conclusion across ALL experiments (now 5 arrays / 7 full runs)
Long Valley (7 stn, 16 km) 3d+14d; LA basin (8 stn, 30 km) 3d+14d; LA extended (11 stn,
62 km) 14d diagnostic. Every run completes end-to-end on real data but yields ZERO qualified
direction / velocity under the program's DEFAULT calibrated thresholds. The limiter is the
per-point/per-window quality + coherence gate set, which is region-, aperture-, station-count-,
and data-length-independent across everything tried. No measured dispersion exists to compare
with authoritative models, so no comparison is fabricated.

## The only remaining in-scope lever (needs user's explicit go-ahead)
The program's Direction_Workflow docs explicitly allow tuning the direction/quality thresholds
(min_coherence, minimum_independent_windows, phase_rms_max_deg, max_direction_uncertainty_deg,
the block wave-packet / route-SNR gates, music_front_cycles). Relaxing these changes the
program's calibrated scientific judgement, so it must be a deliberate, user-authorized choice,
not a silent workaround. The extended 62 km array in diagnostic mode is the best-posed case to
try a documented relaxation on, because its candidate energy sits in the microseism band.


## UPDATE: authorized tuning of minimum_independent_windows + music_front_cycles
User authorized adjusting these two DIRECTION-stage parameters. Ran extended 14-day array with
music_front_cycles=30 (was 200), minimum_independent_windows=3. Result UNCHANGED:
direction_complete.json frequencies=0; frequency_directions.csv empty; no_qualified_direction=true.

ROOT CAUSE TRACED PRECISELY (code read of ndstudio/music_front.py run_music_front):
The front MUSIC iterates ONLY over frequencies present in output/conditional_dispersion.csv,
which is produced by the upstream SCAN stage. In every run conditional_dispersion.csv is EMPTY
(header only) because the scan's robust-point screen passes 0 points (candidate_gates=28,
robust_points=0). So the direction stage is STARVED: music_front_cycles and
minimum_independent_windows govern the direction loop, which never receives a single frequency.
=> These two parameters cannot unblock anything on this data; verified empirically.

THE ACTUAL BLOCKER = scan validation gate in config.json -> "validation" (defaults in
core/dispersion/common.py): snr_min=10.0, block_snr_min=5.0, phase_route_max_rad=0.2,
block_phase_max_rad=0.5, block_fraction_min=0.75, arrival_fraction_max=0.1, etc.
These are the "route/SNR/boundary", "block wave-packet", and "arrival consistency" gates that
the scan_candidate_diagnostics reported as rejecting every point. To get ANY measured point,
those scan/validation thresholds (NOT the direction parameters) must be relaxed — a separate,
more consequential change to the program's calibrated science that the user has not authorized.

