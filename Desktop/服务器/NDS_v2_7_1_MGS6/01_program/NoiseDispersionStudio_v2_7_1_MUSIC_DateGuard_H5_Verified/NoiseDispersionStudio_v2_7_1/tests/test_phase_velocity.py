"""Known-wave tests exercise cycle ambiguity, negative lags and wrap crossings."""
import unittest
import numpy as np
import pandas as pd
from ndstudio.phase_velocity import calculate


class PhaseVelocityTests(unittest.TestCase):
    def known_wave(self, branch):
        f=np.arange(1.0,4.001,0.05)
        distance=3000.; velocity=2000.; source=np.pi/4; cycle=3
        theta=source+2*np.pi*cycle-2*np.pi*f*distance/velocity
        sign=1 if branch=='positive' else -1
        raw=np.degrees(np.angle(np.exp(1j*sign*theta)))
        summary=pd.DataFrame(dict(curve_id='known_N3',frequency_hz=f,segment_id='known',pair='A--B',branch=branch,
            after_count=100,after_raw_phase_deg=raw,after_unit_phase_deg=raw,predicted_phase_deg=raw,
            phase_offset_low_deg=-10.,phase_offset_high_deg=10.))
        targets=pd.DataFrame(dict(curve_id='known_N3',frequency_hz=f,constant_cycle=cycle,source_phase_rad=source,
            distance_m=distance,predicted_unwrapped_phase_rad=sign*theta,phase_velocity_m_s=velocity))
        return summary,targets

    def test_multiple_wraps_positive(self):
        s,t=self.known_wave('positive');d,_=calculate(s,t)
        np.testing.assert_allclose(d.raw_velocity_m_s,2000.,rtol=1e-12)
        np.testing.assert_allclose(d.unit_velocity_m_s,2000.,rtol=1e-12)
        self.assertEqual(d.raw_global_unwrap_shift.nunique(),1)

    def test_multiple_wraps_negative(self):
        s,t=self.known_wave('negative');d,_=calculate(s,t)
        np.testing.assert_allclose(d.raw_velocity_m_s,2000.,rtol=1e-12)

    def test_cycle_alias_keeps_identical_observed_phase(self):
        s,t=self.known_wave('positive');d,_=calculate(s,t)
        t.constant_cycle=4
        alias,_=calculate(s,t)
        np.testing.assert_array_equal(d.final_raw_wrapped_phase_deg,alias.final_raw_wrapped_phase_deg)
        expected=360*s.frequency_hz*3000/(360*s.frequency_hz*3000/2000+360)
        np.testing.assert_allclose(alias.raw_velocity_m_s,expected,rtol=1e-12)
        self.assertTrue(np.all(alias.raw_velocity_m_s<d.raw_velocity_m_s))

    def test_empty_stack_must_not_invent_phase(self):
        s,t=self.known_wave('positive');s.loc[3,'after_count']=0
        with self.assertRaisesRegex(ValueError,'Empty stack'):calculate(s,t)


if __name__=='__main__':unittest.main()
