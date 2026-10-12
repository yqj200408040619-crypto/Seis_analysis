import tempfile,unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dispersion.prepare import read_station,resample_align

class RestartUTC(unittest.TestCase):
    def test_gap_restart_keeps_fractional_native_grid_and_physical_phase(self):
        from obspy import UTCDateTime
        base=UTCDateTime('2026-01-01');fs=40.;outfs=4.;f=.4
        with tempfile.TemporaryDirectory() as folder:
            for name,offset in [('a',.007),('b',320.012)]:
                t=offset+np.arange(12000)/fs
                np.savez(Path(folder)/(name+'.npz'),data=np.cos(2*np.pi*f*t),fs_hz=fs,start_utc=str(base+offset))
            st,_=read_station({'station':'X','file':'*.npz','seed_id':'XX.X..BHZ','clock_offset_s':0},{'data_root':folder})
            self.assertEqual(len(st),2)
            self.assertAlmostEqual(float(st[1].stats.starttime-base),320.012,places=8)
            for tr in st:
                lo,hi,y=resample_align(tr.data,fs,outfs,tr.stats.starttime,base,2600)
                expected=np.cos(2*np.pi*f*np.arange(lo,hi)/outfs)
                self.assertLess(np.max(abs(y[20:-20]-expected[20:-20])),.004)

if __name__=='__main__':unittest.main()
