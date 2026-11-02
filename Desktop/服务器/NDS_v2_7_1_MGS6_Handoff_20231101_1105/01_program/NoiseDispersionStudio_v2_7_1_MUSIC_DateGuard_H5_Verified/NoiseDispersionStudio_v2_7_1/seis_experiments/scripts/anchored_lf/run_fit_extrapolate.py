"""Step 3: fit + bounded extrapolation of a Rayleigh phase-velocity curve.

Inputs:
  - measured high-frequency dispersion (aggregate median v(f), 0.3-1.16 Hz) from the
    NDS traditional-route scan (step 1), reflecting shallow LA-basin sediments (~0.45-0.6 km/s).
  - a-priori LOW-frequency anchor point(s) from an external base model (SCEC SoCal ANT /
    LA-basin literature), clearly labelled as priors, NOT measurements.

Output:
  - a smooth, shape-preserving (PCHIP) monotone curve Vph(period) through BOTH, evaluated on
    a dense frequency grid, with the 5-10 s extrapolated prediction used later to gate windows.
  - written as a prior curve CSV for the predicted-arrival-gated stack (step 4).

Monotonicity enforced in (period, velocity): velocity increases with period (decreases with
frequency), the physical Rayleigh sense for a basin over a faster crust.
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.interpolate import PchipInterpolator

HERE = Path(__file__).resolve().parent
RUN = HERE/'runs'/'hf_core'

# a-priori anchors (external base model; LABEL as prior). period_s -> Vph km/s
ANCHORS = [
    {'period_s': 5.0,  'velocity_km_s': 2.90, 'source': 'SCEC SoCal ANT / LA-basin lit (prior)'},
    {'period_s': 10.0, 'velocity_km_s': 3.30, 'source': 'SCEC SoCal ANT / LA-basin lit (prior)'},
]


def main():
    agg = pd.read_csv(RUN/'hf_dispersion_aggregate.csv')  # columns: frequency_hz, median, count, std
    agg = agg.rename(columns={'median': 'v_m_s'}).sort_values('frequency_hz')
    # Measured points in (period, velocity km/s)
    mf = agg.frequency_hz.to_numpy(float)
    mv = agg.v_m_s.to_numpy(float)/1000.0
    mp = 1.0/mf
    meas = pd.DataFrame({'period_s': mp, 'velocity_km_s': mv, 'kind': 'measured_hf'})

    anch = pd.DataFrame([{'period_s': a['period_s'], 'velocity_km_s': a['velocity_km_s'], 'kind': 'prior_anchor'} for a in ANCHORS])
    both = pd.concat([meas, anch], ignore_index=True).sort_values('period_s')

    # Collapse duplicate periods (median) and enforce strictly increasing period for PCHIP.
    both = both.groupby('period_s', as_index=False).agg(velocity_km_s=('velocity_km_s', 'median'),
                                                        kind=('kind', lambda s: 'both' if s.nunique() > 1 else s.iloc[0]))
    both = both.sort_values('period_s')
    p = both.period_s.to_numpy(float); v = both.velocity_km_s.to_numpy(float)

    # Enforce monotone-increasing velocity with period (cummax) to respect the physical sense
    # and give PCHIP a monotone input (shape-preserving -> monotone output).
    v_mono = np.maximum.accumulate(v)
    fit = PchipInterpolator(p, v_mono, extrapolate=True)

    # Dense prediction grid 0.08-1.2 Hz (period ~0.83-12.5 s)
    fgrid = np.round(np.arange(0.08, 1.2+1e-9, 0.005), 6)
    pgrid = 1.0/fgrid
    vgrid = fit(pgrid)
    # clip to measured range at the HF end, and to a sane Rayleigh ceiling
    vgrid = np.clip(vgrid, 0.15, 4.0)

    pred = pd.DataFrame({'frequency_hz': fgrid, 'period_s': pgrid, 'phase_velocity_m_s': vgrid*1000.0})
    pred['regime'] = np.where(pred.frequency_hz >= mf.min(), 'within_or_above_measured',
                              np.where(pred.frequency_hz >= 0.1, 'extrapolated_target_5_10s', 'extrapolated_beyond'))
    pred.to_csv(RUN/'fitted_prior_curve.csv', index=False)

    # Report the predicted velocity at the key target periods
    def at_period(T):
        return float(fit(T))
    summary = {
        'measured_hf_band_hz': [float(mf.min()), float(mf.max())],
        'measured_hf_v_km_s': [float(mv.min()), float(mv.max())],
        'anchors': ANCHORS,
        'fit': 'PCHIP shape-preserving in (period, velocity), cummax-monotone input',
        'predicted_vph_km_s': {
            '2s': round(at_period(2.0), 3), '3s': round(at_period(3.0), 3),
            '4s': round(at_period(4.0), 3), '5s': round(at_period(5.0), 3),
            '6s': round(at_period(6.0), 3), '7s': round(at_period(7.0), 3),
            '8s': round(at_period(8.0), 3), '10s': round(at_period(10.0), 3)},
        'output_csv': str(RUN/'fitted_prior_curve.csv'),
        'note': 'HF part measured; 5-10 s part is anchored extrapolation (prior), NOT a measurement.'}
    (RUN/'fit_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    print('done')


if __name__ == '__main__':
    main()
