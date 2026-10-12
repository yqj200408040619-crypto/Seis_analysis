#!/usr/bin/env python3
"""Reproduce the real FDSN → H5 → dispersion → structure test using Studio."""
import argparse
import json
import sys
import time
from pathlib import Path

SOURCE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(SOURCE))
import ndstudio  # registers the bundled numerical core
from ndstudio.acquisition import select_array,fetch_array,audit_acquisition,screen_array
from ndstudio.project import new_project
from ndstudio.pipeline import execute
from ndstudio.diagnostics import screen_prepared_array
from dispersion.common import write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workdir',type=Path,default=Path('DARE_real_e2e'))
    p.add_argument('--offline',action='store_true',help='Use the saved original H5 files in workdir/download')
    p.add_argument('--require-model',action='store_true',help='Return an error if processing finishes without an accepted conditional structure')
    p.add_argument('--max-cases',type=int,default=12)
    p.add_argument('--network',default='Y7');p.add_argument('--start',default='2020-10-12');p.add_argument('--end',default='2020-10-16')
    p.add_argument('--size',type=int,default=10);p.add_argument('--anchor')
    a=p.parse_args();root=a.workdir.resolve();root.mkdir(parents=True,exist_ok=True)
    log=[]
    def step(name,action):
        start=time.time();print('REAL E2E:',name,flush=True)
        try:
            result=action()
            log.append({'stage':name,'status':'complete','elapsed_s':time.time()-start})
            write_json(root/'e2e_steps.json',{'steps':log,'complete':False})
            return result
        except Exception as exc:
            log.append({'stage':name,'status':'failed','error':str(exc),'elapsed_s':time.time()-start})
            write_json(root/'e2e_steps.json',{'steps':log,'complete':False});raise
    selection=root/'array_selection.json';data=root/'download';project=root/'project'
    if not a.offline:
        step('array-select',lambda:select_array('GFZ',a.network,a.start,a.end,selection,
                   size=a.size,anchor=a.anchor or ('A00' if a.network=='Y7' else None),
                   min_stations=6 if a.network=='Y7' else 5,max_stations=49 if a.network=='Y7' else 100))
        step('fetch',lambda:fetch_array(selection,data,chunk_seconds=86400,workers=1,allow_missing_data=True))
    elif not selection.exists():
        selection=data/'array_selection.json'
    step('acquisition-audit',lambda:audit_acquisition(data))
    new_project(project)
    spec=json.loads((project/'project.json').read_text(encoding='utf-8-sig'))
    spec.update(input_folder=str(data),name=a.network+' 96h真实端到端测试')
    write_json(project/'project.json',spec)
    step('screen-array',lambda:screen_array(selection,data,project))
    # Numerical thresholds stay at the Studio defaults. Independent solver
    # verification is required for this reproducible real-data test.
    write_json(project/'inversion_settings.json',{'require_independent_solver':True})
    for name in ['import','prepare','scan','controls','families','review','verify','autoselect','invert','report']:
        try:
            step(name,lambda name=name:execute(project,name,max_cases=a.max_cases,modes=(0,)))
            if name=='scan':
                marker=json.loads((project/'output/scan_complete.json').read_text(encoding='utf-8-sig'))
                if not marker['segments']:
                    import pandas as pd
                    counts=pd.read_csv(project/'output/cohort_counts.csv')
                    config=json.loads((project/'config.json').read_text(encoding='utf-8-sig'))
                    if counts.eligible_blocks.max()<config['validation']['min_eligible_blocks']:
                        current=json.loads((data/'analysis_array_selection.json').read_text(encoding='utf-8-sig'))
                        minimum=current['constraints'].get('minimum_station_count',6)
                        if len(current['stations'])<=minimum:raise RuntimeError('Insufficient validation blocks even at the minimum station count')
                        step('screen-time-availability-subarray',lambda:screen_prepared_array(project,data,target_stations=minimum))
                        for again in ['import','prepare','scan']:
                            step(again+'-after-subarray',lambda again=again:execute(project,again))
        except ValueError as exc:
            if name!='prepare' or 'Too few usable windows' not in str(exc):raise
            step('screen-prepared',lambda:screen_prepared_array(project,data))
            step('import-after-availability-screen',lambda:execute(project,'import'))
            step('prepare-after-availability-screen',lambda:execute(project,'prepare'))
    marker=json.loads((project/'output/inversion_complete.json').read_text(encoding='utf-8-sig'))
    write_json(root/'e2e_steps.json',{'steps':log,'complete':True,'compatible_models':marker['compatible_models'],
                                  'structure_accepted':marker['compatible_models']>0,
                                  'acquisition_and_processing_only_through_studio':True})
    print('Report:',project/'output/studio_report.html',flush=True)
    if not marker['compatible_models']:
        print('Processing completed; no conditional structure passed. Rejection tables are preserved.',flush=True)
        if a.require_model:raise RuntimeError('No structure passed the unchanged acceptance conditions')


if __name__=='__main__':main()
