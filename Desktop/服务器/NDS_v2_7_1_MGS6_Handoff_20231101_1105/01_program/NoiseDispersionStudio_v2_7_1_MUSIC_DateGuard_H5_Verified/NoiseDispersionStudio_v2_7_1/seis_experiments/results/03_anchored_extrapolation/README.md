# 03 — Anchored-extrapolation low-frequency recovery

Pipeline: traditional-route HF dispersion (2 days) + a-priori base-model anchor ->
fit+extrapolate -> predicted-arrival short-window stack -> UNBIASED low-frequency recovery.

High-frequency measurement (NDS scan, 0.2-5 Hz @ 20 Hz, 2 days):
- `hf_summary.json`            — scan: 13 accepted continuous runs, conditional band 0.3-1.16 Hz.
- `hf_physical_dispersion.csv` — physical cycle-branch per pair (d/lambda ~9-11).
- `hf_dispersion_aggregate.csv`— median Vph(f): ~0.60 km/s @0.3 Hz -> ~0.45 km/s @1.16 Hz.
- `hf_physical_summary.json`   — per-pair selection summary.

Fit + extrapolation (PCHIP, monotone in period):
- `fit_summary.json`           — anchors 5 s=2.9, 10 s=3.3 km/s; predicted 5-10 s curve.
- `fitted_prior_curve.csv`     — dense predicted Vph(f) used as the prior.

Low-frequency gated stack (14 days, 0.1-0.2 Hz, 12-cycle windows):
- `recovered_lowfreq_velocity.csv`  — UNGATED (prior-free) recovered Vph per frequency.
- `gated_prior_seeded_contrast.csv` — prior-seeded gated result (CIRCULAR; retention ~1).
- `final_recovered_vs_authoritative.csv` — recovered vs authoritative LA-basin Vph.
- `assessment.json`            — median |error| ~5.6% over 5-9 s; best at 0.19-0.20 Hz (1-4%).
- `lf_gated_summary.json`      — full ungated + gated contrast with interpretation.

Verdict: anchored extrapolation sets the band and resolves the cycle; the unbiased short-window
cross-phase stack independently recovers 5-9 s Rayleigh Vph within ~5% of authoritative. The
symmetric predicted-phase gate is circular and is kept only as a transparency contrast.
