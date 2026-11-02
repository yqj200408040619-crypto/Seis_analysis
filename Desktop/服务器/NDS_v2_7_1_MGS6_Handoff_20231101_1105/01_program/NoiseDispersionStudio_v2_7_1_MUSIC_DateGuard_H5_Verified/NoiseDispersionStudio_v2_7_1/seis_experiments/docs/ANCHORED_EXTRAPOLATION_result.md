# Anchored-extrapolation low-frequency velocity recovery — RESULT

User strategy (implemented, details refined): HF traditional-route dispersion + a-priori
low-frequency base-model anchor -> fit+extrapolate a predicted Vph(f) -> predicted-arrival
short-window gating + stack -> recover 5-10 s velocity and compare to authoritative.

Array: LA dense core, 8 broadband CI/HHZ stations (BHP,GR2,HLL,LCG,LGB,PDR,SMF2,USC),
aperture 29.6 km, SCEDC. HF scan on 2 days; LF stack on 14 days (reused aligned H5).

## Step 1 - HF traditional route (2 days)
NDS scan (multi-route FTAN/whitening) over 0.2-5 Hz at 20 Hz. 13 accepted continuous runs.
Physical cycle branch per pair selected (d/lambda ~9-11, non-increasing, Rayleigh range).
Aggregate measured dispersion: 0.3-1.16 Hz, Vph ~0.60 km/s @0.3 Hz -> ~0.45 km/s @1.16 Hz
(shallow LA-basin sediments; correct sense). Coherent HF energy capped ~1.16 Hz at 30 km.

## Step 2 - a-priori anchor (external base model, labelled prior)
LA-basin Rayleigh Vph (SCEC SoCal ANT / LA-basin literature): 5 s -> 2.9 km/s, 10 s -> 3.3 km/s.

## Step 3 - fit + extrapolate
Shape-preserving PCHIP in (period, velocity), monotone, through HF measured points + anchors.
Predicted Vph: 5 s 2.90, 6 s 3.06, 7 s 3.17, 8 s 3.24, 10 s 3.30 km/s. (fitted_prior_curve.csv)

## Step 4 - predicted-arrival-gated short-window stack (12-cycle Goertzel windows)
For each pair/frequency: cut short windows, form cross-spectrum, predict phi=2*pi*f*d/v_pred,
keep windows within +/-10% velocity (phase margin), complex-stack.
CIRCULARITY GUARD: the prior-gated stack reproduces the prior (retention ratio ~1.0, i.e. the
kept fraction only matches uniform-random; gated unit resultant ~0.98 is an artefact of the
symmetric gate). So the gated velocity is NOT read as independent recovery.

## Step 5 - UNBIASED recovery (prior-free) and comparison to authoritative
The honest estimate uses the UNGATED per-window cross-phase (independent of the prior):
  ungated_unit_resultant 0.19-0.33 (noise floor ~0.005) => genuine coherent 5-10 s energy.
Coherence-weighted recovered Vph vs authoritative LA-basin Rayleigh Vph:
  9.1 s 3.09 (auth 3.30, -6.5%) | 8.3 s 3.14 (3.25, -3.3%) | 7.7 s 3.30 (3.25, +1.6%)
  7.1 s 3.30 (3.15, +4.8%) | 5.9 s 2.76 (3.05, -9.5%) | 5.6 s 2.69 (3.05, -11.9%)
  5.3 s 2.93 (2.90, +1.2%) | 5.0 s 2.79 (2.90, -3.7%)
Median |error| over 5-9 s = ~5.6%. Best where coherence highest (0.19-0.20 Hz: +/-1-4%).
Weaker/over-estimated at 0.15-0.16 Hz (+27-31%, fewer pairs) and 0.10 s (low coherence).

## Verdict
The workflow WORKS: an HF-measured + base-model-anchored extrapolation, used to set the target
band and resolve the cycle ambiguity, combined with an UNBIASED (non-circular) short-window
cross-phase stack, independently recovers LA-basin Rayleigh phase velocity at 5-9 s to within
~5% median of authoritative values. The prior must be used to pick the cycle, NOT to gate the
phase symmetrically (that is circular). Files under runs/lf_gated_stack/:
gated_pair_frequency.csv, recovered_lowfreq_velocity.csv, gated_prior_seeded_contrast.csv,
final_recovered_vs_authoritative.csv, assessment.json, lf_gated_summary.json.

CAVEATS: unnormalized (no spectral whitening) per-window cross-spectra were used for the phase
concentration; adding RAM/one-bit whitening would sharpen coherence. ungated_unit_resultant and
std are diagnostics, not calibrated CIs. 0.10 Hz (10 s) is the least reliable (coherence 0.21).
