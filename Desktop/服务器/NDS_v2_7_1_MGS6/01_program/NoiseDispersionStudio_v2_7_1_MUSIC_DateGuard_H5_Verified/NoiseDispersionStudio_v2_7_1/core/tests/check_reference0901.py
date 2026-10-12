"""Numerical regression on the supplied stage7 fixed-cohort CCFs, not raw replay."""
from pathlib import Path
import sys,json,copy,argparse
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dispersion.common import DEFAULTS,write_json
from dispersion.scan import spectrum,ftan_peaks,peak_arrays

def check(stage7,destination):
    stage7=Path(stage7);c=copy.deepcopy(DEFAULTS)
    cat=pd.read_csv(stage7/"results/measured_segment_catalog.csv");ph=pd.read_csv(stage7/"results/measured_phase_segments.csv")
    pairs=pd.read_csv(stage7/"results/pairs.csv").set_index("pair_index");stacks=np.load(stage7/"results/cohort_stacks.npz")
    phase_error=0.;arrival_error=0.;points=0
    for r in cat.itertuples():
        q=ph[ph.segment_id==r.segment_id].sort_values("frequency_hz");f=q.frequency_hz.to_numpy();cc=stacks[f"{r.cohort}__ram5__combined"][r.pair_index]
        h=spectrum(cc,r.branch,(r.gate_start_s,r.gate_end_s),c,f);theta=np.unwrap(np.angle(h))
        theta+=2*np.pi*np.rint(np.median((q.measured_unwrapped_phase_rad-theta)/(2*np.pi)))
        phase_error=max(phase_error,float(np.max(abs(theta-q.measured_unwrapped_phase_rad))))
        peaks=ftan_peaks(cc,float(pairs.loc[r.pair_index,"distance_m"]),c,f)[r.branch]
        arrival,_,edge=peak_arrays(peaks,q.arrival_s.to_numpy());assert np.isfinite(arrival).all() and not edge.any()
        arrival_error=max(arrival_error,float(np.max(abs(arrival-q.arrival_s))));points+=len(q)
    assert phase_error<1e-6 and arrival_error<1e-6
    result=dict(checks_passed=True,reference_segments=len(cat),reference_phase_points=points,
        maximum_phase_error_rad=phase_error,maximum_arrival_error_s=arrival_error,
        raw_0901_end_to_end_replay=False,test_scope="generic FTAN and fixed-gate phase on previously validated CCFs")
    write_json(destination,result);print(json.dumps(result,indent=2))

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("stage7",type=Path);p.add_argument("--out",type=Path,required=True);a=p.parse_args();check(a.stage7,a.out)
