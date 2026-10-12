import unittest,sys,copy
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dispersion.common import DEFAULTS
from dispersion.prepare import repair_local

class LocalScale(unittest.TestCase):
    def test_smooth_amplitude_burst_is_kept_but_impulse_and_flatline_are_removed(self):
        fs=40.;t=np.arange(72000)/fs
        amplitude=1+99*.5*(np.tanh((t-300)/12)-np.tanh((t-480)/12))
        clean=amplitude*np.cos(2*np.pi*.4*t);x=clean.copy()
        impulse=int(360*fs)+7;x[impulse]+=20000
        x[int(1200*fs):int(1202*fs)]=500
        p=copy.deepcopy(DEFAULTS['spikes']);p.update(scale_window_s=60,max_repair_s=.075)
        fixed,valid,repaired,events=repair_local(x,fs,p)
        self.assertTrue(repaired[impulse]);self.assertTrue(valid[impulse])
        self.assertLess(abs(fixed[impulse]-clean[impulse]),.3)
        burst=(t>320)&(t<460);self.assertGreater(valid[burst].mean(),.999)
        self.assertFalse(valid[int(1201*fs)]);self.assertFalse(repaired[int(1201*fs)])
        unchanged=valid&~repaired;self.assertTrue(np.array_equal(fixed[unchanged],x[unchanged]))

if __name__=='__main__':unittest.main()
