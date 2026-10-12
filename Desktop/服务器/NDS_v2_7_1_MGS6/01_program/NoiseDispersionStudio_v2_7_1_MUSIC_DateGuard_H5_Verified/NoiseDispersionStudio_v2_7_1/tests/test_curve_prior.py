import numpy as np
import pandas as pd
from ndstudio.curve_prior import monotone_curves,screening_prior


def test_screening_only_prediction_and_increasing_curve_rejection(tmp_path):
    rows=[]
    for cid,vel in [('decreasing',[700.,650.,600.]),('increasing',[500.,520.,540.])]:
        for f,v in zip([2.,2.05,2.1],vel):
            rows.append(dict(curve_id=cid,frequency_hz=f,phase_velocity_m_s=v,branch='negative',
                source_phase_rad=np.pi/4,constant_cycle=2,distance_m=500.,is_model_prediction=False,
                background_discriminated=True,phase_offset_low_deg=-15.,phase_offset_high_deg=15.))
    q=screening_prior(pd.DataFrame(rows),tmp_path,dict(min_hz=1.9,max_hz=2.2,step_hz=.05),
        dict(enabled=True,extension_hz=.1),True)
    assert set(q.curve_id)=={'decreasing'}
    pred=q[q.is_prior_prediction];assert len(pred)==4
    assert pred.is_model_prediction.all() and not pred.background_discriminated.any()
    assert pred.measured_unwrapped_phase_rad.isna().all()
    measured=q[~q.is_prior_prediction];assert measured.frequency_hz.tolist()==[2.,2.05,2.1]
    assert (q.sort_values('frequency_hz').phase_velocity_m_s.diff().dropna()<=0).all()
    assert pd.read_csv(tmp_path/'input_curve_trend_screen.csv').accepted.tolist()==[True,False]


def test_cycle_aliases_keep_separate_extrapolated_phase_gates(tmp_path):
    rows=[]
    for cycle,velocities in [(2,[700.,650.,600.]),(3,[450.,430.,410.])]:
        for frequency,velocity in zip([2.,2.05,2.1],velocities):
            rows.append(dict(curve_id=f'alias_{cycle}',segment_id='observed_gate',
                frequency_hz=frequency,phase_velocity_m_s=velocity,branch='positive',
                source_phase_rad=np.pi/4,constant_cycle=cycle,distance_m=500.,
                is_model_prediction=False,background_discriminated=True))
    q=screening_prior(pd.DataFrame(rows),tmp_path,dict(min_hz=1.9,max_hz=2.2,step_hz=.05),
        dict(enabled=True,extension_hz=.1),True)
    observed=q[~q.is_prior_prediction]
    assert observed.segment_id.unique().tolist()==['observed_gate']
    extra=q[q.is_prior_prediction]
    assert extra.segment_id.nunique()==2
    assert extra.groupby('frequency_hz').segment_id.nunique().eq(2).all()
    assert extra.phase_equivalence_id.equals(extra.segment_id)
    theta=extra.source_phase_rad-2*np.pi*extra.frequency_hz*extra.distance_m/extra.phase_velocity_m_s+2*np.pi*extra.constant_cycle
    np.testing.assert_allclose(extra.predicted_unwrapped_phase_rad,theta)
