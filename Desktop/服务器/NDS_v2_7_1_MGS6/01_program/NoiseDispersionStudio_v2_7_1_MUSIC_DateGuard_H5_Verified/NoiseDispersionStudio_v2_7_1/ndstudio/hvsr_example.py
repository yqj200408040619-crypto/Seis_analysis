"""Reproducible broadband spectral HVSR example, NOT a full wavefield simulation."""
from pathlib import Path
import json
import numpy as np
import h5py
from .hvsr import write_json,run_hvsr,DEFAULTS
from .project import new_project


def create_example(directory):
    directory=Path(directory).resolve();raw=directory/'raw';raw.mkdir(parents=True,exist_ok=True)
    fs=50.;window=60.;n=int(fs*window);rng=np.random.default_rng(20261010)
    f=np.fft.rfftfreq(n,1/fs);target=1+5*np.exp(-.5*(np.log(np.maximum(f,.001)/3.)/.2)**2)
    data=[[],[],[]]
    for w in range(30):
        for j in range(3):
            phases=rng.uniform(0,2*np.pi,len(f));amp=np.ones(len(f))*(target if j<2 else 1)
            coeff=amp*np.exp(1j*phases);coeff[0]=0;coeff[-1]=coeff[-1].real
            data[j].append(np.fft.irfft(coeff,n)*1e-6)
    data=np.asarray([np.concatenate(x) for x in data])
    # A missing interval and a transient are diagnosed; other windows remain intact.
    data[0,4*n+100:4*n+200]=np.nan;data[1,12*n+100]=.02
    with h5py.File(raw/'XX.DEMO..HH_3C.h5','w') as h:
        h.attrs['synthetic']=True
        for j,ch in enumerate(('HHN','HHE','HHZ')):
            ds=h.create_dataset('XX.DEMO..'+ch,data=data[j],compression='gzip')
            ds.attrs.update(role='waveform',seed_id='XX.DEMO..'+ch,channel=ch,network='XX',station='DEMO',location='',
                            sampling_rate_hz=fs,start_utc='2026-01-01T00:00:00Z',units='m/s',precorrected=True,
                            azimuth_deg=90. if j==1 else 0.,dip_deg=-90. if j==2 else 0.)
    root=new_project(directory/'project');spec=json.loads((root/'project.json').read_text(encoding='utf-8-sig'));spec['input_folder']='../raw';write_json(root/'project.json',spec)
    c=DEFAULTS|dict(fmin_hz=.3,fmax_hz=15.,vs_min_m_s=200.,vs_max_m_s=300.,bedrock_vs_m_s=1000.)
    write_json(root/'hvsr_settings.json',dict(folder='../raw',instrument_id='XX.DEMO..HH',settings=c))
    write_json(directory/'truth.json',dict(type='controlled broadband spectral example, not seismic full-wavefield simulation',designed_peak_hz=3.,
                                         conditional_example_vs_m_s=240.,conditional_example_depth_m=20.,injected_bad_windows=[4,12]))
    result=run_hvsr(raw,root/'output/hvsr','XX.DEMO..HH',c)
    selected=result['selected_peak']
    if selected is None or abs(selected['f0_hz']/3-1)>.06 or not selected['clear_peak']: raise AssertionError('Synthetic H/V peak recovery failed')
    return root
