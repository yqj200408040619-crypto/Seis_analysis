#!/usr/bin/env python3
"""Desktop GUI or reproducible batch commands; run with --help for examples."""
import os
for var in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS','MKL_NUM_THREADS'):os.environ.setdefault(var,'1')
from pathlib import Path
import argparse,json,sys
from ndstudio import __version__

def main():
 p=argparse.ArgumentParser(description='Noise Dispersion Studio: HDF5 → ambient-noise dispersion → direction/phase correction → conditional inversion → geothermal/tomography handoff')
 p.add_argument('command',nargs='?',default='gui',choices=['gui','download','station-list','array-select','fetch','acquisition-audit','screen-array','screen-prepared','diagnose-prepare','init','inspect','import','all','prepare','scan','scan-diagnostics','controls','families','music','review','report','verify','autoselect','invert','example','invert-example','invert-csv','doctor','phase-export','phase-stack','phase-sensitivity','phase-audit','phase-velocity','phase-compare','hvsr','hvsr-inspect','hvsr-example','hvsr-fetch','tomography-export','tomography-invert','tomography-example','geothermal-validate','geothermal','geothermal-example'])
 p.add_argument('--project',type=Path);p.add_argument('--folder',type=Path);p.add_argument('--profile',type=Path);p.add_argument('--directory',type=Path,default=Path('examples/generated'))
 p.add_argument('--curves',type=Path);p.add_argument('--selection',type=Path);p.add_argument('--output',type=Path);p.add_argument('--auto-invert',action='store_true');p.add_argument('--max-cases',type=int,default=12);p.add_argument('--modes',default='0');p.add_argument('--field',choices=['plane','noise','isotropic'],default='plane');p.add_argument('--screenshot',type=Path)
 p.add_argument('--provider',default='GFZ');p.add_argument('--network');p.add_argument('--start');p.add_argument('--end');p.add_argument('--size',type=int,default=10);p.add_argument('--anchor');p.add_argument('--channel',help='Default HH? for hvsr-fetch, HHZ for other commands');p.add_argument('--min-spacing-m',type=float,default=200.);p.add_argument('--max-spacing-m',type=float,default=2000.);p.add_argument('--max-aperture-m',type=float,default=6000.);p.add_argument('--chunk-seconds',type=float,default=21600.);p.add_argument('--workers',type=int,default=4)
 p.add_argument('--station',default='*',help='Station codes/wildcards for station-list and array-select; comma-separated names supported')
 p.add_argument('--station-prefix',default='',help='Station initial letters, e.g. QS or QS,QT')
 p.add_argument('--min-latitude',type=float);p.add_argument('--max-latitude',type=float)
 p.add_argument('--min-longitude',type=float,help='West longitude; west > east crosses the dateline')
 p.add_argument('--max-longitude',type=float,help='East longitude')
 p.add_argument('--allow-missing-data',action='store_true')
 p.add_argument('--min-stations',type=int,default=6);p.add_argument('--max-stations',type=int,default=49)
 p.add_argument('--minimum-wavelengths',type=float,default=3.,help='Explicit automatic inversion-selection distance/wavelength floor (default 3)')
 p.add_argument('--diagnostic-top-k',type=int,default=8,help='Number of robust intervals and promising candidate gates printed after scan')
 p.add_argument('--curve-ids',help='Comma-separated curve IDs; phase-export always exports all curves')
 p.add_argument('--phase-tolerance-deg',type=float,default=15.)
 p.add_argument('--phase-range-source',choices=['fixed','measured'],default='fixed')
 p.add_argument('--repeat-cycles',type=float,default=20.)
 p.add_argument('--slide-cycles',type=float,default=1.)
 p.add_argument('--taper',choices=['none','hann'],default='none')
 p.add_argument('--phase-settings',type=Path)
 p.add_argument('--reconstruct-low',action='store_true')
 p.add_argument('--low-frequency-hz',type=float,default=.2)
 p.add_argument('--direction-mode',choices=['off','diagnostic','filter','joint'],help='Override saved direction mode')
 p.add_argument('--max-pair-angle-deg',type=float,help='Override saved maximum pair/wave-vector angle (default45)')
 p.add_argument('--stabilization-seconds',type=float,default=600.)
 p.add_argument('--goertzel-frequency',choices=['matlab_bin','exact'],default='matlab_bin')
 p.add_argument('--phase-screening',choices=['curve','matlab'],default='curve')
 p.add_argument('--keep-amplitude-anomalies',action='store_true')
 p.add_argument('--combine-phase-signs',action='store_true')
 p.add_argument('--phase-time-mode',choices=['full','custom'],help='Default: all imported records; explicit --start/--end selects a custom interval')
 p.add_argument('--aligned-h5',type=Path,help='Reuse response-corrected aligned H5 for offline reproduction')
 p.add_argument('--tolerances-deg',default='5,10,15,20')
 p.add_argument('--hvsr-settings',type=Path)
 p.add_argument('--instrument-id',help='NET.STA.LOC.channel-prefix; from hvsr-inspect')
 p.add_argument('--hvsr-tab',action='store_true')
 p.add_argument('--selected-only',action='store_true',help='For tomography-export, use only curve_ids in selection.json')
 a=p.parse_args()
 a.channel=a.channel or ('HH?' if a.command=='hvsr-fetch' else 'HHZ')
 station_filters=dict(station_prefix=a.station_prefix,minlatitude=a.min_latitude,maxlatitude=a.max_latitude,
                      minlongitude=a.min_longitude,maxlongitude=a.max_longitude)
 if a.command=='hvsr-fetch':
  if not all([a.network,a.start,a.end,a.output]):p.error('hvsr-fetch requires --network --station --start --end --output; use --channel HH? or BH?')
  from ndstudio.hvsr_acquisition import fetch_single_station
  result=fetch_single_station(a.provider,a.network,a.station,a.start,a.end,a.output,channel=a.channel,
    instrument_id=a.instrument_id,catalog_file=a.selection,chunk_seconds=a.chunk_seconds)
  print(json.dumps(result,ensure_ascii=False,indent=2));return 0
 if a.command=='hvsr-inspect':
  if not a.folder:p.error('hvsr-inspect requires --folder')
  from ndstudio.hvsr import discover
  print(json.dumps(discover(a.folder),ensure_ascii=False,indent=2));return 0
 if a.command=='hvsr-example':
  from ndstudio.hvsr_example import create_example
  print(create_example(a.directory));return 0
 if a.command=='hvsr':
  from ndstudio.hvsr import run_hvsr
  config_file=a.hvsr_settings or (a.project/'hvsr_settings.json' if a.project else None)
  data=json.loads(config_file.read_text(encoding='utf-8-sig')) if config_file else {}
  folder=a.folder or data.get('folder')
  if folder and not a.folder and config_file and not Path(folder).is_absolute():folder=config_file.parent/folder
  output=a.output or (a.project/'output/hvsr' if a.project else None)
  if not folder or not output:p.error('hvsr requires --folder --output or --project with hvsr_settings.json')
  config=dict(data.get('settings',data if 'folder' not in data else {}))
  if a.start:config['start']=a.start
  if a.end:config['end']=a.end
  result=run_hvsr(folder,output,a.instrument_id or data.get('instrument_id'),config)
  print(json.dumps({'instrument_id':result['instrument_id'],'peak':result['selected_peak'],'structure':result['structure'],'output':str(output)},ensure_ascii=False,indent=2));return 0
 if a.command=='geothermal-example':
  from ndstudio.geothermal import create_example,generate_geothermal
  root=create_example(a.directory);print(root);print(json.dumps(generate_geothermal(root),ensure_ascii=False,indent=2));return 0
 if a.command=='tomography-example':
  from ndstudio.tomography import create_tomography_example,run_phase_velocity_tomography
  root=create_tomography_example(a.directory);print(root);print(json.dumps(run_phase_velocity_tomography(root),ensure_ascii=False,indent=2));return 0
 if a.command in ('tomography-export','tomography-invert','geothermal-validate','geothermal'):
  if not a.project:p.error(a.command+' requires --project')
  if a.command=='tomography-export':
   from ndstudio.geothermal import export_tomography_paths
   print(json.dumps(export_tomography_paths(a.project,True if a.selected_only else None),ensure_ascii=False,indent=2));return 0
  if a.command=='tomography-invert':
   from ndstudio.tomography import run_phase_velocity_tomography
   print(json.dumps(run_phase_velocity_tomography(a.project),ensure_ascii=False,indent=2));return 0
  if a.command=='geothermal-validate':
   from ndstudio.geothermal import validate_spatial_vs
   print(json.dumps(validate_spatial_vs(a.project),ensure_ascii=False,indent=2));return 0
  from ndstudio.geothermal import generate_geothermal
  print(json.dumps(generate_geothermal(a.project),ensure_ascii=False,indent=2));return 0
 if a.command=='station-list':
  if not a.network or not a.output:p.error('station-list requires --network and --output (catalogue directory)')
  from ndstudio.acquisition import list_matching_stations
  r=list_matching_stations(a.provider,a.network,a.output,station=a.station,channel=a.channel,
                           start=a.start,end=a.end,**station_filters)
  print(json.dumps({'stations':r['matching_station_count'],'channel_epochs':r['matching_channel_epochs'],
                    'output':str(a.output)},ensure_ascii=False,indent=2));return 0
 if a.command=='phase-velocity':
  if not a.folder or not a.output:p.error('phase-velocity requires --folder and --output')
  from ndstudio.phase_velocity import run_phase_velocity
  focus=[x.strip() for x in a.curve_ids.split(',') if x.strip()] if a.curve_ids else None
  print(json.dumps(run_phase_velocity(a.folder,a.output,focus_curve_ids=focus),ensure_ascii=False,indent=2));return 0
 if a.command=='phase-audit':
  if not a.folder:p.error('phase-audit requires --folder pointing to phase output')
  from ndstudio.phase_audit import audit
  print(json.dumps(audit(a.folder),indent=2));return 0
 if a.command=='phase-sensitivity':
  if not a.folder:p.error('phase-sensitivity requires --folder pointing to phase output')
  from ndstudio.phase_stack import phase_sensitivity
  print(phase_sensitivity(a.folder,tuple(float(t) for t in a.tolerances_deg.split(','))).to_string(index=False));return 0
 if a.command=='phase-compare':
  if not a.project:p.error('--project is required')
  from ndstudio.reconstruction import paired_compare
  print(json.dumps(paired_compare(a.project),ensure_ascii=False,indent=2));return 0
 if a.command in ('phase-export','phase-stack'):
  if not a.project:p.error('--project is required')
  from ndstudio.phase_stack import export_targets,run_phase_stack
  settings_path=a.phase_settings or a.project/'phase_settings.json'
  settings=json.loads(settings_path.read_text(encoding='utf-8-sig')) if settings_path.exists() else {}
  if a.auto_invert and a.output and a.output.resolve()!=(a.project/'output/phase_guided').resolve():p.error('自动配对反演须使用标准相位目录，请省略 --output')
  if a.command=='phase-export':
   q,d=export_targets(a.project,a.output,settings.get('tolerance_deg',a.phase_tolerance_deg),settings.get('range_source',a.phase_range_source))
   print(json.dumps({'curves':int(q.curve_id.nunique()),'points':len(q),'output':str(d)},ensure_ascii=False));return 0
  r=run_phase_stack(a.project,a.output,settings.get('curve_ids',a.curve_ids),settings.get('start',a.start),settings.get('end',a.end),
    settings.get('repeat_cycles',a.repeat_cycles),settings.get('slide_cycles',a.slide_cycles),settings.get('tolerance_deg',a.phase_tolerance_deg),settings.get('range_source',a.phase_range_source),settings.get('taper',a.taper),a.aligned_h5,a.phase_time_mode or settings.get('time_mode'),
    settings.get('stabilization_s',a.stabilization_seconds),settings.get('frequency_method',a.goertzel_frequency),settings.get('amplitude_rejection',not a.keep_amplitude_anomalies),
    settings.get('split_branches',not a.combine_phase_signs),settings.get('screening_mode',a.phase_screening),
    (settings.get('direction',{'mode':'filter','method':'music'})|({'mode':a.direction_mode} if a.direction_mode else {})|({'max_pair_angle_deg':a.max_pair_angle_deg} if a.max_pair_angle_deg is not None else {})),settings.get('monotonic_only',False),settings.get('prior'),
    settings.get('reconstruction',{})|({'enabled':True,'lower_frequency_hz':a.low_frequency_hz} if a.reconstruct_low else {}))
  if a.auto_invert and r.get('reconstruction',{}).get('enabled'):
   from ndstudio.candidates import auto_selection
   from ndstudio.reconstruction import paired_compare
   auto_selection(a.project,max_cases=a.max_cases,modes=tuple(map(int,a.modes.split(','))),minimum_wavelengths=a.minimum_wavelengths)
   r['paired_comparison']=paired_compare(a.project)
  print(json.dumps({k:v for k,v in r.items() if k!='station_metadata'},ensure_ascii=False,indent=2));return 0
 if a.command=='array-select':
  if not all([a.network,a.start,a.end,a.output]):p.error('array-select requires --network --start --end --output')
  from ndstudio.acquisition import select_array
  r=select_array(a.provider,a.network,a.start,a.end,a.output,size=a.size,anchor=a.anchor,channel=a.channel,min_spacing_m=a.min_spacing_m,max_spacing_m=a.max_spacing_m,max_aperture_m=a.max_aperture_m,min_stations=a.min_stations,max_stations=a.max_stations,station=a.station,**station_filters)
  print(json.dumps({'stations':len(r['stations']),'geometry':r['geometry']},ensure_ascii=False,indent=2));return 0
 if a.command=='fetch':
  if not a.selection or not a.output:p.error('fetch requires --selection --output')
  from ndstudio.acquisition import fetch_array
  fetch_array(a.selection,a.output,chunk_seconds=a.chunk_seconds,workers=a.workers,allow_missing_data=a.allow_missing_data);return 0
 if a.command=='acquisition-audit':
  if not a.folder:p.error('acquisition-audit requires --folder')
  from ndstudio.acquisition import audit_acquisition
  audit_acquisition(a.folder);return 0
 if a.command=='screen-array':
  if not a.selection or not a.folder:p.error('screen-array requires --selection --folder; optional --project applies exclusions')
  from ndstudio.acquisition import screen_array
  r=screen_array(a.selection,a.folder,a.project)
  print(json.dumps({'stations':r['screened_stations'],'exclusions':r['exclusions'],'geometry':r['geometry']},ensure_ascii=False,indent=2));return 0
 if a.command=='download':
  from h5data_gui import main as download_gui
  download_gui();return 0
 if a.command=='gui':
  from ndstudio.gui import launch
  return launch(a.project,a.screenshot,a.hvsr_tab)
 if a.command=='doctor':
  import importlib
  for name in ['numpy','scipy','pandas','h5py','obspy','pyproj','matplotlib','disba','PySide6','pysurf96','geopandas','pygmt','xarray','netCDF4']:
   try:m=importlib.import_module(name);print(name,getattr(m,'__version__','available'))
   except ImportError:print(name,'NOT INSTALLED'+(' (optional independent solver)' if name=='pysurf96' else ' (optional geothermal)' if name in ('geopandas','pygmt','xarray','netCDF4') else ''))
  return 0
 if a.command=='example':
  from ndstudio.examples import waveform_example
  waveform_example(a.directory,a.field);return 0
 if a.command=='invert-example':
  from ndstudio.examples import inversion_example
  inversion_example(a.directory);return 0
 if a.command=='invert-csv':
  if not all([a.curves,a.selection,a.output]):p.error('invert-csv requires --curves --selection --output')
  import pandas as pd
  from ndstudio.inversion import invert_tables
  invert_tables(pd.read_csv(a.curves),json.loads(a.selection.read_text(encoding='utf-8-sig')),a.output);return 0
 if a.command=='inspect':
  if not a.folder:p.error('inspect requires --folder')
  from ndstudio.hdf import inspect_folder
  d=inspect_folder(a.folder,json.loads(a.profile.read_text(encoding='utf-8-sig')) if a.profile else {})
  print(json.dumps(d,ensure_ascii=False,indent=2));return 0 if d['ready'] else 2
 if not a.project:p.error('--project is required')
 if a.command=='diagnose-prepare':
  from ndstudio.diagnostics import diagnose_preparation
  diagnose_preparation(a.project);return 0
 if a.command=='screen-prepared':
  if not a.folder:p.error('screen-prepared requires --project --folder')
  from ndstudio.diagnostics import screen_prepared_array
  screen_prepared_array(a.project,a.folder,target_stations=a.size);return 0
 from ndstudio.project import new_project
 if a.command=='init':
  root=new_project(a.project)
  if a.folder:
   from dispersion.common import write_json
   spec=json.loads((root/'project.json').read_text(encoding='utf-8-sig'));spec['input_folder']=str(a.folder.resolve());write_json(root/'project.json',spec)
  print(root);return 0
 if a.folder:
  from dispersion.common import write_json
  path=a.project/'project.json';spec=json.loads(path.read_text(encoding='utf-8-sig'));spec['input_folder']=str(a.folder.resolve());write_json(path,spec)
 from ndstudio.pipeline import execute
 execute(a.project,a.command,auto_invert=a.auto_invert,max_cases=a.max_cases,modes=tuple(int(s) for s in a.modes.split(',')),minimum_wavelengths=a.minimum_wavelengths,scan_top_k=a.diagnostic_top_k)
 return 0

if __name__=='__main__':
 try:sys.exit(main())
 except Exception as exc:
  import traceback
  traceback.print_exc();sys.exit(1)
