"""Independent all-window branch/MAD check plus scalar real-waveform recurrence."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
import h5py
import argparse
p=argparse.ArgumentParser(description='Independent check of uploaded MATLAB core on phase-stack MATLAB-mode output')
p.add_argument('--folder',required=True,type=Path);folder=p.parse_args().folder.resolve()
if json.loads((folder/'phase_stack_complete.json').read_text()).get('screening_mode')!='matlab':raise ValueError('Requires phase-stack screening_mode=matlab')
gates=pd.read_csv(folder/'phase_gate_summary.csv')
maximum=0.;checked=0;windows=0
with h5py.File(folder/'aligned_waveforms.h5') as h:
 for r in gates.itertuples():
  with np.load(folder/r.windows_file,allow_pickle=False) as z:
   c=z['cross'];amp=abs(c);vals=amp[np.isfinite(amp)];m=np.median(vals);s=1.4826*np.median(abs(vals-m));ath=m+5*s if s>0 else 5*m
   d=np.r_[0.,abs(np.diff(amp))];dv=d[1:][np.isfinite(d[1:])];m=np.median(dv) if len(dv) else 0.;s=1.4826*np.median(abs(dv-m)) if len(dv) else 0.;dth=m+5*s if s>0 else 5*m
   good=np.isfinite(c)&(amp>0)&(amp<=ath)&(d<=dth);phase=np.angle(c)
   for sign,side in [('positive',phase>=0),('negative',phase<0)]:
    preliminary=c[good&side].sum();center=np.angle(preliminary) if abs(preliminary)>0 else np.nan
    accepted=good&side&(abs(np.angle(c*np.exp(-1j*center)))<=np.pi/2)
    np.testing.assert_array_equal(accepted,z['accepted_'+sign]);actual=c[accepted].sum()
    expected=complex(getattr(r,'after_'+sign+'_raw_real'),getattr(r,'after_'+sign+'_raw_imag'))
    np.testing.assert_allclose(actual,expected,rtol=1e-12,atol=1e-30)
   a,b=r.pair.split('--');n=int(z['window_samples']);fs=float(z['fs_hz']);k=int(np.floor(r.frequency_hz*n/fs+.5));w=2*np.pi*k/n
   for pick in np.unique(np.linspace(0,len(z['start_sample'])-1,min(3,len(z['start_sample'])),dtype=int)):
    first=int(z['start_sample'][pick]);gs=[]
    for station in (a,b):
     x=h['stations/'+station+'/data'][first:first+n];valid=h['stations/'+station+'/valid'][first:first+n];x=np.where(valid,x,np.nan);x=x-x.mean();s1=s2=0.
     for value in x:s0=value+2*np.cos(w)*s1-s2;s2=s1;s1=s0
     gs.append(s1-np.exp(-1j*w)*s2)
    actual=gs[1]*np.conj(gs[0]);np.testing.assert_allclose(actual,c[pick],rtol=1e-10,atol=1e-30)
    maximum=max(maximum,float(abs(np.angle(actual*np.conj(c[pick])))));checked+=1
   windows+=len(c)
result=dict(passed=True,all_window_masks_checked=windows,gates_checked=len(gates),raw_scalar_cross_recomputed=checked,maximum_scalar_phase_error_rad=maximum,
 note='Independent translation of the uploaded MATLAB core; MATLAB runtime not executed.')
(folder/'matlab_reference_check.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
