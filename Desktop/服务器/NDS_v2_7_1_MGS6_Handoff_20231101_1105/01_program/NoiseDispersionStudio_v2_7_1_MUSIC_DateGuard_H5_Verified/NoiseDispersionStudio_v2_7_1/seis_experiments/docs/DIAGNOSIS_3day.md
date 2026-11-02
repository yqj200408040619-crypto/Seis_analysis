# Long Valley NC array — end-to-end run diagnosis (3-day window)

## Array (VERIFIED, live NCEDC metadata + waveforms)
- Region: Long Valley Caldera, Mammoth Lakes, CA. Network NC (served by NCEDC).
- 7 broadband HHZ stations, compact subset chosen by program's own select_array:
  MCB, MCY, MDPB, MDY, MINS, MLI, MMLB (screening later swapped MCB->MQ1P).
- Geometry: aperture (max pair) = 16.17 km; min pair = 1.20 km; mean NN = 2.97 km; 7 > 6.
- Sensors: STS-2, Trillium Compact, CMG-3ESPC, Cascadia; all 100 Hz; full coverage 2023-11-01..04.

## Aperture/frequency gate (matched to aperture)
- resolvable band f_min = c_lo/D, f_max = c_hi/(2 d_min), c in [1.5, 4.5] km/s.
- At selection: [0.093, 1.868] Hz. After screening (effective d_min changed): [0.161, 1.868] Hz.
- Spans the secondary microseism.

## Pipeline stages (3-day span 2023-11-01..04, NCEDC)
array-select OK -> fetch OK (real data, aligned_waveforms.h5 = 139 MB)
 -> audit-acquisition OK -> screen-array OK -> preprocess+scan OK -> front-MUSIC OK
 -> select-before: FAILED hard constraint.

## Result: NO stable direction/velocity passed strict default thresholds
- MUSIC front: 0 frequencies accepted; 0 joint raw points; 0 continuous-band points.
- Scan: 11 candidate dispersion gates formed (pairs incl. MCY-MMLB, MDPB-MMLB, MDY-MINS,
  MDPB-MINS, MLI-MMLB), all at 0.24-1.25 Hz, but robust_points = 0 for every gate.
- Dominant rejection bottlenecks (point_diagnostics via scan_candidate_diagnostics):
  route/SNR/boundary comprehensive check; block wave-packet fraction; block packet+phase
  joint; time-segment arrival. i.e. insufficient coherent SNR / too few robust windows.
- Scan needs >=8 contiguous points over >=0.05 Hz; longest run = 0 points. 0 robust runs.

## Interpretation
Signal fragments exist (11 gates) but none survive SNR/wave-packet/arrival-consistency gates,
and MUSIC direction accepts nothing. The candidate gates cluster ABOVE the microseism peak
(0.24-1.25 Hz) where 16 km inter-station coherence is weak and a 3-day span yields too few
independent 200-cycle windows. This is a scientific-threshold negative, not a code error:
with 3 days + default strict gates, this array does not yield a calibrated stable direction+velocity.

## Legitimate next lever (not threshold-weakening)
Longer time span -> more independent windows -> stabler ambient-noise stack. Standard remedy.
Re-run with a longer --start/--end (e.g. ~2 weeks) reusing the pipeline unchanged.
