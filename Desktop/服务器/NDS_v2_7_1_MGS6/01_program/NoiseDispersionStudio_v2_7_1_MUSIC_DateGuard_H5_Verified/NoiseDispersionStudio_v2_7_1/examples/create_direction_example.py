"""Create a reproducible 7-station plane wave and run the program's phase workflow.

Known propagation: north-clockwise 30 degrees, 500 m/s, 600 s at 20 Hz.
The illustrative CCF source/cycle aliases are conditional and differ from 500 m/s.
No download or external processing program is used.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
import pandas as pd
APP=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(APP),str(APP/'core')]
XY=np.array([[0,0],[67,12],[21,94],[133,48],[-72,173],[219,159],[87,-123]],float)

def create_example(root):
    import h5py
    from ndstudio.project import new_project
    from dispersion.common import write_json,sha
    root=new_project(root);out=root/'output';out.mkdir(exist_ok=True)
    fs=20.;duration=600.;t=np.arange(round(duration*fs))/fs;f=np.arange(1.3,1.81,.1)
    names=['S'+str(i) for i in range(len(XY))];az=30.;v=500.;slow=np.array([.5,np.sqrt(3)/2])/v
    # Broadband delayed stochastic wave; real-domain generation has independent spectral phases.
    rng=np.random.default_rng(510);freq=np.fft.rfftfreq(len(t),1/fs);spectrum=np.fft.rfft(rng.normal(size=len(t)))
    spectrum[(freq<.8)|(freq>2.5)]=0
    data=np.array([np.fft.irfft(spectrum*np.exp(-2j*np.pi*freq*(xy@slow)),n=len(t)) for xy in XY])
    stations=pd.DataFrame(dict(station=names,file=[s+'.npz' for s in names],seed_id=['XX.'+s+'..HHZ' for s in names],x_m=XY[:,0],y_m=XY[:,1],precorrected=True,units='m/s',polarity=1))
    stations.to_csv(root/'stations.csv',index=False)
    native=root/'canonical';native.mkdir(exist_ok=True)
    for name,x in zip(names,data):
        np.savez_compressed(native/(name+'.npz'),data=x,fs_hz=fs,start_utc='2020-01-01T00:00:00.000000Z')
    c=json.loads((root/'config.json').read_text());c['scan'].update(min_hz=1.,max_hz=2.,step_hz=.1,min_band_hz=.25,min_points=3);write_json(root/'config.json',c)
    curves=[];measured=[];pairs=[]
    for k,(a,b) in enumerate([(0,2),(0,1),(1,3)]):
        d=XY[b]-XY[a];distance=np.linalg.norm(d);pair=names[a]+'--'+names[b];pairs.append(dict(pair_index=k,pair=pair,distance_m=distance))
        for frequency in f:
            phase=-2*np.pi*frequency*(d@slow);cycle=1;source=np.pi/4
            conditional=2*np.pi*frequency*distance/(source-phase+2*np.pi*cycle)
            curves.append(dict(curve_id=pair,segment_id=pair,pair_index=k,pair=pair,branch='positive',frequency_hz=frequency,
                phase_velocity_m_s=conditional,group_velocity_m_s=500.,source_phase_rad=source,constant_cycle=cycle,
                cohort='synthetic',instrument_response_corrected=True,background_discriminated=True,source_phase_verified=False,is_model_prediction=False))
            measured.append(dict(segment_id=pair,frequency_hz=frequency,measured_unwrapped_phase_rad=phase,relative_phase_q025_rad=-.4,relative_phase_q975_rad=.4))
    pd.DataFrame(curves).to_csv(out/'conditional_dispersion.csv',index=False);pd.DataFrame(measured).to_csv(out/'measured_phase_segments.csv',index=False);pd.DataFrame(pairs).to_csv(out/'pairs.csv',index=False)
    path=root/'aligned.h5'
    with h5py.File(path,'w') as h:
        h.attrs.update(fs_hz=fs,start_utc='2020-01-01T00:00:00.000000Z',end_utc_exclusive='2020-01-01T00:10:00.000000Z',units='m/s',response_corrected=True,preprocessing_config_sha256=sha(root/'config.json'))
        for name,x in zip(names,data):
            g=h.create_group('stations/'+name);g.create_dataset('data',data=x);g.create_dataset('valid',data=np.ones(len(x),np.uint8))
    return root,path

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,default=APP/'examples/direction_plane/project')
    args=parser.parse_args()
    root,h5=create_example(args.directory)
    from ndstudio.phase_stack import run_phase_stack
    from ndstudio.phase_audit import audit
    from dispersion.common import write_json
    options=dict(stabilization_s=0.,repeat_cycles=40.,slide_cycles=2.,taper='hann',
        frequency_method='exact',time_mode='full',amplitude_rejection=False,tolerance_deg=20.,
        direction=dict(mode='filter',velocity_min_m_s=200.,velocity_max_m_s=1000.,velocity_points=41))
    write_json(root/'phase_settings.json',options)
    result=run_phase_stack(root,aligned_h5=h5,**options)
    # This example starts with prescribed target curves, rather than a CCF scan.
    write_json(root/'output/families_complete.json',dict(synthetic_prescribed_targets=True))
    checked=audit(root/'output/phase_guided')
    joint=pd.read_csv(root/'output/phase_guided/direction/joint_dispersion.csv')
    selected=joint[joint.accepted]
    summary=dict(synthetic=True,known_azimuth_deg=30.,known_velocity_m_s=500.,
        joint_velocity_median_m_s=float(selected.phase_velocity_m_s.median()),
        selected_curves=result['selected_curves'],accepted_window_evaluations=result['accepted_window_evaluations'],
        inversion_handoff_points=result['inversion_handoff_points'],audit=checked)
    write_json(root/'example_validation.json',summary)
    print(json.dumps(summary,ensure_ascii=False,indent=2))
    print('Report:',root/'output/phase_guided/phase_report.html')

if __name__=='__main__':main()
