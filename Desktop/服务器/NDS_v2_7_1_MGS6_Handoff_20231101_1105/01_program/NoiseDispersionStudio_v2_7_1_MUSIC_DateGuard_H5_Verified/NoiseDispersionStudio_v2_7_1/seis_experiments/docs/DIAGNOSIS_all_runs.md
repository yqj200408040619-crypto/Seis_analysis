# Dense-array direction/velocity/dispersion — consolidated diagnosis (4 full runs)

All runs used NDS v2.7.1 array_search_e2e (unified MUSIC per-frequency direction + phase
velocity), real FDSN waveforms, program's own gates. Evidence = files under runs/.

## Arrays actually run end-to-end on REAL data
| Region (network@provider)        | #stns | aperture | d_min | band gate (Hz)   | span   | h5 chunks |
|----------------------------------|-------|----------|-------|------------------|--------|-----------|
| Long Valley caldera (NC@NCEDC)   | 7     | 16.2 km  | 1.2km | 0.161–1.868      | 3 day  | real 139MB|
| Long Valley caldera (NC@NCEDC)   | 7     | 16.2 km  | 1.2km | 0.161–1.868      | 14 day | real      |
| LA basin (CI@SCEDC)              | 8     | 29.6 km  | 1.85km| 0.051–1.216      | 3 day  | 96        |
| LA basin (CI@SCEDC)              | 8     | 29.6 km  | 1.85km| 0.051–1.216      | 14 day | 448       |

All four: array-select -> fetch -> audit -> screen -> preprocess+scan -> front-MUSIC ->
review -> verify -> select-before ALL OK. Then select-before FAILS: "no directionally-
consistent candidate curves survived selection".

## Provider/federation fixes discovered (both verified)
- NC-network waveforms live at NCEDC (IRIS federates NC metadata only -> empty waveforms).
- CI-network waveforms live at SCEDC (IRIS serves only USC; 7/8 fail FDSNNoDataException).
Both candidates in array_search_e2e.py set to the archive that actually serves waveforms.

## Result (identical across all four runs)
- MUSIC front: frequencies=0, raw_joint_points=0, accepted_continuous_points=0.
- Scan candidate gates: LV 3d=11, LV 14d=25, LA 3d=9, LA 14d=29. ALL with robust_points=0.
- No measured_dispersion.csv / dispersion_reference.json / comparison/ in any run
  (pipeline never reached select-after). Verified by Test-Path.

## Where the LA-basin gates sit (physically correct)
LA-basin 14d: 29 gates, almost all in the secondary-microseism band 0.05–0.34 Hz, on
classic basin pairs (GR2-USC, BHP-USC, LCG-USC, PDR-SMF2, ...). The scan clearly detects
coherent microseism-band surface waves. The caldera gates (Long Valley) sat higher/more
scattered (0.16–1.43 Hz) as expected for a volcanic scattering medium.

## Binding constraint (triple-confirmed, region-independent)
The default per-point quality gates reject EVERY point:
  - 分块波包比例 / 分块波包与相位联合 (block wave-packet fraction / packet+phase joint)
  - 处理路线综合检验（含SNR/边界）(route/SNR/boundary comprehensive check)
  - 时段到时 / 相位时延·波包到时 (time-segment & packet arrival consistency)
plus MUSIC coherence>=0.7, direction uncertainty<=20deg, >=3 independent windows, and
scan needs >=8 contiguous robust points over >=0.05 Hz.

Increasing data (3->14 d) raised the gate count and moved LA gates into the microseism band,
but did NOT make any point pass the block wave-packet / SNR gates. Changing region
(volcanic caldera -> coherent sedimentary basin), raising station count (7->8), and enlarging
aperture (16->29.6 km) also did not. => The limiter is the strict DEFAULT threshold set,
not the array or the data volume.

## Honest verdict for the task
Under the program's default calibrated scientific thresholds, NO stable wave direction and
NO phase/group velocity survive on any tried array, so there is NO measured dispersion to
compare against an authoritative model. Fabricating a direction/velocity/dispersion from the
empty outputs would be wrong and is not done. Authoritative references were recorded in the
candidate definitions (western-US ANT: Lin/Moschetti/Ritzwoller 2008; LA basin: Lin et al.
2013 / Lee et al. 2014) for use IF a future run (with justified threshold relaxation and/or
much longer stacks) produces a surviving measured band.
