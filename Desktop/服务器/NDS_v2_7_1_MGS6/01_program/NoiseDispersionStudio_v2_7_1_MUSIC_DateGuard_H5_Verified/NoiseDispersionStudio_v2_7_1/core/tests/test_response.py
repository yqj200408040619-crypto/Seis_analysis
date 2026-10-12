"""Response units/phase and epoch matching on a known analytic instrument."""
import copy,tempfile,sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dispersion.common import DEFAULTS
from dispersion.prepare import response_table,correct_response

class ResponseTests(unittest.TestCase):
    def test_nominal_response_removal_and_epoch_guard(self):
        from obspy import UTCDateTime,Trace,Stream
        from obspy.core.inventory import Inventory,Network,Station,Channel,Site
        from obspy.core.inventory.response import Response
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);date=UTCDateTime("2026-01-01");fs=20.
            response=Response.from_paz([],[-2*np.pi*3],stage_gain=1000.,input_units="M/S**2",output_units="COUNTS",
                normalization_factor=2*np.pi*np.sqrt(10.))
            ch=Channel(code="HGZ",location_code="",latitude=0.,longitude=0.,elevation=0.,depth=0.,azimuth=0.,dip=-90.,
                sample_rate=fs,start_date=date-1,end_date=date+86400,response=response)
            station=Station(code="TST",latitude=0.,longitude=0.,elevation=0.,site=Site(name="synthetic"),channels=[ch])
            Inventory([Network(code="XX",stations=[station])],source="synthetic response test").write(str(root/"response.xml"),format="STATIONXML")
            c=copy.deepcopy(DEFAULTS);c["stations_csv"]=str(root/"stations.csv")
            trace=Trace(np.zeros(6000),dict(network="XX",station="TST",location="",channel="HGZ",sampling_rate=fs,starttime=date))
            row=dict(station="TST",seed_id="XX.TST..HGZ",response_xml="response.xml",precorrected=False,units="counts")
            table,audit=response_table(row,Stream([trace]),c)
            self.assertEqual(audit["response_status"],"nominal_response_corrected")
            t=np.arange(6000)/fs;physical=np.sin(2*np.pi*2.3*t)+.3*np.cos(2*np.pi*4.1*t)
            freq=np.fft.rfftfreq(len(t),1/fs);h=response.get_evalresp_response_for_frequencies(np.maximum(freq,.00001),output="ACC")
            counts=np.fft.irfft(np.fft.rfft(physical)*h,len(t));recovered=correct_response(counts,table,c)
            self.assertLess(np.max(abs(recovered[1200:-1200]-physical[1200:-1200])),.002)
            trace.stats.starttime=date+172800
            with self.assertRaises(ValueError): response_table(row,Stream([trace]),c)

if __name__=="__main__": unittest.main()
