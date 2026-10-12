#!/usr/bin/env python3
"""BW.MGS01–MGS05: invoke only Studio acquisition and analysis modules."""
import argparse
import json
import os
import sys
import time
from pathlib import Path

for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(key,'1')
SOURCE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(SOURCE))
import ndstudio
from ndstudio.acquisition import select_array,fetch_array,audit_acquisition
from ndstudio.project import new_project
from ndstudio.pipeline import execute
from dispersion.common import write_json

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workdir',type=Path,required=True)
    p.add_argument('--offline',action='store_true')
    p.add_argument('--start',default='2020-11-01')
    p.add_argument('--end',default='2020-11-05')
    p.add_argument('--minimum-wavelengths',type=float,default=3.)
    p.add_argument('--stages',default='import,prepare,scan,controls,families,review,verify,autoselect,invert,report')
    a=p.parse_args();root=a.workdir.resolve();root.mkdir(parents=True,exist_ok=True)
    selection=root/'array_selection.json';data=root/'download';project=root/'project'
    history=root/'e2e_steps.json'
    log=json.loads(history.read_text(encoding='utf-8-sig')).get('steps',[]) if history.exists() else []
    def step(name,action):
        start=time.time();print('BW MGS E2E:',name,flush=True)
        try:
            result=action();log.append({'stage':name,'status':'complete','elapsed_s':time.time()-start})
            write_json(root/'e2e_steps.json',{'complete':False,'steps':log})
            return result
        except Exception as exc:
            log.append({'stage':name,'status':'failed','error':str(exc),'elapsed_s':time.time()-start})
            write_json(root/'e2e_steps.json',{'complete':False,'steps':log});raise
    if not selection.exists() and not a.offline:
        step('array-select',lambda:select_array('LMU','BW',a.start,a.end,selection,size=5,
             min_stations=5,max_stations=100,station='MGS01,MGS02,MGS03,MGS04,MGS05',
             channel='HHZ',min_spacing_m=1000,max_spacing_m=15000,max_aperture_m=30000))
    selected=json.loads(selection.read_text(encoding='utf-8-sig'))
    expected={f'MGS{i:02d}' for i in range(1,6)}
    if {r['station'] for r in selected['stations']}!=expected:
        raise ValueError('Selection must contain exactly BW.MGS01–MGS05')
    if not a.offline:
        step('fetch',lambda:fetch_array(selection,data,chunk_seconds=21600,workers=4))
        step('acquisition-audit',lambda:audit_acquisition(data))
    if not (project/'project.json').exists():
        new_project(project)
        spec=json.loads((project/'project.json').read_text(encoding='utf-8-sig'))
        spec.update(input_folder=str(data),name='BW.MGS01–MGS05 real 96-hour test',automatic_geometry=False)
        write_json(project/'project.json',spec)
        c=json.loads((project/'config.json').read_text(encoding='utf-8-sig'))
        # Settings follow the actual kilometre-scale geometry; validation
        # thresholds remain unchanged. No measured points or models are supplied.
        c.update(sample_rate_hz=10.,bandpass_hz=[.05,2.5],max_lag_s=100.)
        c['response'].update(pre_filt_hz=[.01,.025,3.,4.],table_step_hz=.002)
        c['scan'].update(min_hz=.15,max_hz=1.5,step_hz=.005,minimum_sigma_hz=.015,
                         min_band_hz=.05,arrival_min_s=.5,arrival_max_s=45.,
                         noise_lag_s=[65.,95.],group_velocity_bounds_m_s=[500.,4500.],
                         gate_margin_s=3.,gate_change_s=[1.,.5])
        c['phase'].update(velocity_bounds_m_s=[500.,4000.])
        write_json(project/'config.json',c)
        write_json(project/'inversion_settings.json',{'thickness_bounds_m':[10.,10000.],
                   'vs_bounds_m_s':[200.,5000.],'require_independent_solver':False})
    spec=json.loads((project/'project.json').read_text(encoding='utf-8-sig'))
    spec['input_folder']=str(data)
    write_json(project/'project.json',spec)
    for name in a.stages.split(','):
        if name=='invert':
            chosen=json.loads((project/'selection.json').read_text(encoding='utf-8-sig'))
            if not chosen.get('cases'):
                reason='No measured dispersion group passed automatic inversion selection.'
                write_json(project/'output/inversion_skipped.json',{'reason':reason,'not_a_processing_error':True})
                log.append({'stage':'invert','status':'skipped','reason':reason})
                continue
        step(name,lambda name=name:execute(project,name,max_cases=12,modes=(0,),minimum_wavelengths=a.minimum_wavelengths))
    marker=project/'output/inversion_complete.json'
    models=json.loads(marker.read_text(encoding='utf-8-sig')).get('compatible_models',0) if marker.exists() else 0
    write_json(root/'e2e_steps.json',{'complete':a.stages.split(',')[-1]=='report','steps':log,'compatible_models':models,
               'acquisition_and_processing_only_through_studio':True})

if __name__=='__main__':main()
