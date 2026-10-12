import copy,sys,unittest
from pathlib import Path
import numpy as np
from scipy import signal
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dispersion.common import DEFAULTS,merge,validate_config,bands
from dispersion.prepare import repair_local,resample_align,cc_from_fft
from dispersion.scan import track_peaks,fixed_gate,spectrum

class NumericalTests(unittest.TestCase):
    def test_fractional_utc_alignment(self):
        from obspy import UTCDateTime
        t=np.arange(16000)/40.;base=UTCDateTime("2026-01-01");origin=base+2.
        # Different native grid phases must converge to the same physical UTC phase.
        for offset in (.007,.013,.021):
            x=np.cos(2*np.pi*4.2*(t+offset))+.2*np.sin(2*np.pi*2.3*(t+offset))
            lo,hi,y=resample_align(x,40.,20.,base+offset,origin,7800)
            expected=np.cos(2*np.pi*4.2*(2+np.arange(7800)/20))+.2*np.sin(2*np.pi*2.3*(2+np.arange(7800)/20))
            self.assertEqual(lo,0)
            self.assertLess(np.max(abs(y[40:-40]-expected[lo+40:hi-40])),.004)

    def test_positive_lag_sign_and_no_circular_wrap(self):
        import pandas as pd
        c=copy.deepcopy(DEFAULTS);c["sample_rate_hz"]=10.;c["max_lag_s"]=4.
        x=np.zeros((2,200));x[0,70]=1;x[1,90]=1;nf=400
        out=cc_from_fft(np.fft.rfft(x,nf),nf,200,c,pd.DataFrame([dict(i=0,j=1)]))
        self.assertEqual(out[0].argmax()-40,20)
        x[:]=0;x[0,2]=1;x[1,198]=1
        out=cc_from_fft(np.fft.rfft(x,nf),nf,200,c,pd.DataFrame([dict(i=0,j=1)]))
        self.assertLess(abs(out).max(),1e-8)

    def test_local_spike_and_gap_masks(self):
        fs=100.;t=np.arange(3000)/fs;x=np.sin(2*np.pi*2*t);x[501]+=100
        x[1200:1202]=np.nan;x[2000:2200]=50
        repaired,valid,mask,events=repair_local(x,fs,DEFAULTS["spikes"])
        self.assertTrue(mask[501]);self.assertTrue(valid[501])
        self.assertLess(abs(repaired[501]-np.sin(2*np.pi*2*t[501])),.02)
        self.assertFalse(valid[1200]);self.assertFalse(valid[2100]);self.assertTrue(valid[2500])
        self.assertFalse(mask[2100]);self.assertGreater(len(events),0)

    def test_multiple_peaks_are_preserved(self):
        c=copy.deepcopy(DEFAULTS);f=np.arange(1.5,3.01,.025);peaks=[]
        for k,_ in enumerate(f):
            peaks.append([dict(arrival=2.+.002*k,snr=30.,edge=False),dict(arrival=4.+.002*k,snr=15.,edge=False)])
        paths=track_peaks(peaks,f,c)
        arrivals=[np.median(list(p["values"].values())) for p in paths]
        self.assertTrue(any(x<3 for x in arrivals));self.assertTrue(any(x>3 for x in arrivals))

    def test_gaps_are_not_unwrapped_as_one_band(self):
        f=np.arange(1,3.01,.025);mask=np.ones(len(f),bool);mask[30:45]=False
        result=bands(mask,f,8,.25)
        self.assertEqual(len(result),2);self.assertTrue(all(not np.any((ix>=30)&(ix<45)) for ix in result))

    def test_branch_phase_and_cycle_reconstruction(self):
        c=copy.deepcopy(DEFAULTS);fs=20.;t=np.arange(1201)/fs;f=np.arange(2.,3.01,.025)
        y=np.cos(2*np.pi*2.5*(t-3))*np.exp(-.5*((t-3)/.5)**2)
        cc=np.r_[np.zeros(1200),y]
        h=spectrum(cc,"positive",(1,5),c,f)
        self.assertLess(np.median(abs(-np.gradient(np.unwrap(np.angle(h)),f)/(2*np.pi)-3)),.01)
        phi=np.pi/4;distance=1000.;velocity=500.+20*f;true=phi-2*np.pi*f*distance/velocity
        measured=np.unwrap(np.angle(np.exp(1j*true)));n=int(np.rint(np.median((measured-true)/(2*np.pi))))
        recovered=2*np.pi*f*distance/(phi-measured+2*np.pi*n)
        self.assertTrue(np.allclose(recovered,velocity))

    def test_invalid_config_and_unknown_keys_fail(self):
        with self.assertRaises(ValueError): merge(DEFAULTS,{"staions_csv":"typo"})
        c=copy.deepcopy(DEFAULTS);c["scan"]["noise_lag_s"]=[10,20]
        with self.assertRaises(ValueError): validate_config(c)
        validate_config(copy.deepcopy(DEFAULTS))

if __name__=="__main__": unittest.main()
