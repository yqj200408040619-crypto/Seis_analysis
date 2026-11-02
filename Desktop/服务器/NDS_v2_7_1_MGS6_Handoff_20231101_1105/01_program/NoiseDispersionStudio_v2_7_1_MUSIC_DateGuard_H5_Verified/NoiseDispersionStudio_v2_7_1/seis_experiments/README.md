# Seismic ambient-noise experiments on NoiseDispersionStudio (NDS) v2.7.1

This folder collects the analysis scripts, curated results, and documentation for a series of
ambient-noise experiments built **on top of** NoiseDispersionStudio v2.7.1 (the unified-MUSIC
per-frequency direction / phase-velocity program in the parent directory). Nothing here
replaces the program; every script reuses the program's own acquisition, alignment, scan
(traditional FTAN/whitening route), MUSIC direction, and Goertzel phase-stack machinery.

The work answers four successive requests:
1. Find a dense ~20 km array (>6 broadband stations), match the frequency band to the aperture,
   find directionally-consistent station pairs, recover a stable wave direction + velocity, and
   compare the measured dispersion to authoritative data.
2. Try other regions with slightly more stations and a slightly larger aperture (<30 km).
3. Extend the aperture with a few sparse (possibly cross-network) stations to sharpen the
   direction near 5 s; then refine with the dense 20-30 km sub-array.
4. Use the program's traditional route to get a HIGH-frequency dispersion curve from ~2 days,
   anchor it with an a-priori low-frequency base-model velocity point, fit + extrapolate to
   5-10 s, and use the predicted arrival to gate short CCF windows and stack toward a low-
   frequency velocity; compare to authoritative values.

## Environment (verified)
- Python 3.14 (`py -3.14`), numpy 2.4, scipy 1.16.3, pandas 2.3.3, h5py 3.16, obspy 1.5.1.
- Run from the PROGRAM ROOT with `PYTHONPATH` = program root and `PYTHONUTF8=1`.
  `ndstudio/__init__.py` auto-adds `../core`, so `import dispersion.*` works.
- Live FDSN required. IMPORTANT provider/federation facts discovered here:
  - NC-network (Long Valley) WAVEFORMS are at **NCEDC** (IRIS federates NC metadata only).
  - CI-network (LA basin / Southern California) WAVEFORMS are at **SCEDC** (IRIS serves only a
    few CI stations; most return FDSNNoDataException). Metadata is available on IRIS.
  Always fetch waveforms from the archive that actually holds them.

## Folder layout
```
seis_experiments/
  README.md                         this file
  docs/                             narrative diagnoses & results (read these)
    DIAGNOSIS_all_runs.md           arrays run, default-threshold negative, provider fixes
    DIAGNOSIS_extended_aperture.md  62 km extended-aperture experiment + trace to scan gate
    STABLE_DIRECTION_result.md      stable direction under authorized relaxed scan validation
    ANCHORED_EXTRAPOLATION_result.md  HF + anchor -> extrapolate -> gated stack -> LF recovery
    DIAGNOSIS_3day.md               first Long Valley 3-day diagnosis
  scripts/
    probing/        live-FDSN station/geometry probes (metadata only, cheap)
    launchers/      Windows .bat/.ps1 wrappers used to run the drivers
    direction/      run_extended_direction.py (extended + core MUSIC direction)
    anchored_lf/    run_hf_scan.py, run_fit_extrapolate.py, run_lf_gated_stack.py
  results/
    01_array_search/   array/aperture/geometry probe outputs + driver statuses
    02_direction/      stable per-frequency MUSIC directions (extended 62 km + core 30 km)
    03_anchored_extrapolation/  HF dispersion, fitted prior, LF recovery vs authoritative
```
Heavy artefacts (multi-GB `.h5` downloads, `.npz` window archives, caches) are deliberately
NOT committed; only compact JSON/CSV summaries and the narrative docs are kept. See
`run_*` scripts to regenerate the heavy intermediates from live FDSN.

## Headline results (all from program files; see docs for full evidence)
- Dense arrays verified against live metadata: Long Valley NC (7 stn, 16 km, NCEDC) and
  LA basin CI (8 stn, 30 km, SCEDC). Preset Parkfield candidates do NOT hold (only 2 BB stns).
- Under the program's DEFAULT calibrated thresholds, NO stable direction / continuous
  dispersion survives on any array tried (region-, aperture-, data-length-independent). The
  binding gate is the SCAN `validation` block (route-SNR / block wave-packet / arrival), not
  the direction parameters. No dispersion comparison was fabricated.
- With the USER-AUTHORIZED relaxation of the scan `validation` gate (snr_min 10->3, block
  wave-packet and arrival loosened), the extended 62 km array yields a stable direction at
  0.10-0.115 Hz (azimuth ~32-46 deg NE, Vph ~2.5-2.7 km/s, dir. uncertainty ~1.8 deg); the
  dense 30 km core refines it (0.115-0.125 Hz, ~27-31 deg, ~2.8-2.95 km/s). Two independent
  apertures agree. These are relaxed-threshold (non-default) estimates, clearly labelled.
- Anchored-extrapolation low-frequency recovery WORKS: a 2-day traditional-route HF dispersion
  (0.3-1.16 Hz, ~0.5 km/s basin sediments) + an a-priori 5-10 s base-model anchor, fitted and
  extrapolated, then an UNBIASED (non-circular) short-window cross-phase stack, independently
  recovers LA-basin Rayleigh phase velocity at 5-9 s to ~5% median of authoritative values.
  The prior is used ONLY to set the band and resolve the cycle integer; a symmetric
  predicted-phase gate is circular and is reported separately as such.

## Reproduce (from program root, PYTHONPATH=program root, PYTHONUTF8=1)
```
# 1. probe live FDSN for dense arrays (metadata only)
py -3.14 scripts/probing/run_newregion_probe.py
# 2. array search driver (downloads ~days of waveforms) - see ndstudio/array_search_e2e.py
py -3.14 -m ndstudio.array_search_e2e --workdir runs/array_search_labasin --only LABasin_CI --start 2023-11-01 --end 2023-11-04
# 3. extended-aperture + core MUSIC direction (authorized relaxed validation)
#    env: EXT_RELAX=1 EXT_SNR_MIN=3 EXT_FRONT_CYCLES=30 ; see docs/STABLE_DIRECTION_result.md
py -3.14 scripts/direction/run_extended_direction.py --workdir runs/extended_la_14d
# 4. anchored-extrapolation LF recovery
py -3.14 scripts/anchored_lf/run_hf_scan.py --workdir runs/hf_core
py -3.14 scripts/anchored_lf/run_fit_extrapolate.py
py -3.14 scripts/anchored_lf/run_lf_gated_stack.py --project runs/core_la_14d/project --prior runs/hf_core/fitted_prior_curve.csv --out runs/lf_gated_stack
```

## Honesty / caveats
- All "authoritative" comparisons use published Southern-California / LA-basin ambient-noise
  Rayleigh phase velocities (SCEC SoCal ANT; Lin/Moschetti/Ritzwoller 2008 for western US).
- Relaxed-threshold and anchored-prior results are explicitly non-default; uncertainties quoted
  are resolution diagnostics, not calibrated confidence intervals.
- The predicted-arrival symmetric gate is circular by construction; the reported low-frequency
  recovery uses the prior-independent ungated cross-phase, with the circular contrast kept for
  transparency (`results/03_anchored_extrapolation/gated_prior_seeded_contrast.csv`).
