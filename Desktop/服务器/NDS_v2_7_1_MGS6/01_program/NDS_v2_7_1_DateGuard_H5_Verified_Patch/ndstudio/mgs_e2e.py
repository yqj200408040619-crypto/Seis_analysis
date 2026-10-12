"""Native, resumable LMU BW MGS01–05 + MGSBH workflow."""
from pathlib import Path
import argparse,json,os
for name in ('OPENBLAS_NUM_THREADS','NUMBA_NUM_THREADS','OMP_NUM_THREADS'):os.environ.setdefault(name,'1')
import pandas as pd
from dispersion.common import write_json,sha
from .acquisition import select_array,fetch_array
from .project import new_project
from .real_compare import preprocess,completed_phase
from .pipeline import execute
from .candidates import auto_selection
from .phase_stack import run_phase_stack
from .phase_audit import audit
from .reconstruction import paired_compare

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workdir',type=Path,required=True);p.add_argument('--reuse-download',action='store_true');p.add_argument('--start',default='2022-11-01',help='UTC start date/time');p.add_argument('--end',default='2022-11-05',help='UTC end date/time');a=p.parse_args()
 root=a.workdir.resolve();root.mkdir(parents=True,exist_ok=True);data=root/'download';project=root/'project';selected=root/'array_selection.json';state=root/'mgs_e2e_status.json';stages=[]
 def norm_utc(value):
  t=pd.Timestamp(value)
  t=t.tz_localize('UTC') if t.tzinfo is None else t.tz_convert('UTC')
  return t.strftime('%Y-%m-%dT%H:%M:%S.%fZ')
 start_utc,end_utc=norm_utc(a.start),norm_utc(a.end)
 if pd.Timestamp(end_utc)<=pd.Timestamp(start_utc):raise ValueError('--end必须晚于--start')
 from .date_guard import bind_workdir,validate_download
 bind_workdir(root,start_utc,end_utc)
 run_meta=dict(program_version='2.7.1',requested_utc=[start_utc,end_utc])
 print(f'【MGS六站】目标时段 {start_utc} -> {end_utc}',flush=True)
 def run(name,fn):
  print('【MGS六站】'+name,flush=True)
  try:r=fn()
  except Exception as e:write_json(state,dict(run_meta,complete=False,stages=stages,failed_stage=name,error=str(e)));raise
  stages.append(name);write_json(state,dict(run_meta,complete=False,stages=stages));return r
 expected={'MGS01','MGS02','MGS03','MGS04','MGS05','MGSBH'}
 if selected.exists():
  old=json.loads(selected.read_text(encoding='utf-8-sig'))
  old_window=[norm_utc(old.get('start_utc')),norm_utc(old.get('end_utc'))]
  if old_window!=[start_utc,end_utc]:
   if (data.exists() and any(data.rglob('*'))) or (project/'project.json').exists():
    raise ValueError(f'工作目录已有其他时段数据 {old_window}；为避免混入旧H5，请改用新的--workdir')
   selected.unlink()
 if not selected.exists():run('array-select',lambda:select_array('https://erde.geophysik.uni-muenchen.de','BW',a.start,a.end,selected,station='MGS*',channel='HHZ',size=6,min_spacing_m=1000,max_spacing_m=15000,max_aperture_m=30000))
 s=json.loads(selected.read_text(encoding='utf-8-sig'))
 if [norm_utc(s.get('start_utc')),norm_utc(s.get('end_utc'))]!=[start_utc,end_utc]:raise ValueError('选台文件时段与本次请求不一致')
 if {r['station'] for r in s['stations']}!=expected:raise ValueError('六台物理仪器名单不一致')
 if not a.reuse_download:run('fetch',lambda:fetch_array(selected,data,workers=4,chunk_seconds=21600,allow_missing_data=True))
 d=json.loads((data/'acquisition_status.json').read_text(encoding='utf-8-sig'))
 if not d.get('complete') or {r['seed_id'].split('.')[1] for r in d['results'] if r.get('existing') or r.get('h5_files')}!=expected:raise ValueError('必须下载完成且六台均有波形')
 run('date-guard',lambda:write_json(root/'date_guard_audit.json',validate_download(data,start_utc,end_utc)))
 write_json(root/'download_gaps.json',{r['seed_id']:r['skipped'] for r in d['results']})
 if not (project/'project.json').exists():
  new_project(project);spec=json.loads((project/'project.json').read_text(encoding='utf-8-sig'));spec.update(input_folder=str(data),automatic_geometry=True,exclude_increasing_velocity=True);write_json(project/'project.json',spec)
  c=json.loads((project/'config.json').read_text(encoding='utf-8-sig'));c.update(sample_rate_hz=10.,bandpass_hz=[.05,2.5],max_lag_s=100.);c['response'].update(pre_filt_hz=[.01,.025,3.,4.],table_step_hz=.002)
  c['scan'].update(min_hz=.15,max_hz=1.5,step_hz=.005,minimum_sigma_hz=.015,min_band_hz=.05,arrival_min_s=.5,arrival_max_s=45.,noise_lag_s=[65.,95.],group_velocity_bounds_m_s=[500.,4500.],gate_margin_s=3.,gate_change_s=[1.,.5]);c['phase'].update(velocity_bounds_m_s=[500.,4000.]);write_json(project/'config.json',c)
  write_json(project/'inversion_settings.json',dict(thickness_bounds_m=[10.,10000.],vs_bounds_m_s=[200.,5000.],require_independent_solver=True))
 spec=json.loads((project/'project.json').read_text(encoding='utf-8-sig'))
 if Path(spec['input_folder']).resolve()!=data or spec.get('excluded_stations'):raise ValueError('项目输入或排除名单与六站测试不一致')
 spec['requested_utc']=[start_utc,end_utc];write_json(project/'project.json',spec)
 phase_settings_path=project/'phase_settings.json'
 phase_options=json.loads(phase_settings_path.read_text(encoding='utf-8-sig')) if phase_settings_path.exists() else {}
 phase_options['direction']=dict(mode='filter',max_pair_angle_deg=45.,direction_step_deg=2.,velocity_points=128,velocity_min_m_s=200.,velocity_max_m_s=5000.,frequency_half_width_hz=0.,frequency_bins=1,music_front_cycles=200.)
 write_json(phase_settings_path,phase_options)
 # An old phase marker must not block review of the frozen BEFORE input.
 old_phase=project/'output/phase_guided/phase_stack_complete.json'
 if old_phase.exists():
  recorded=json.loads(old_phase.read_text(encoding='utf-8-sig'))
  fingerprints={'module_sha256':'phase_stack.py','direction_sha256':'direction.py','music_module_sha256':'music.py','reconstruction_module_sha256':'reconstruction.py'}
  if any(recorded.get(key)!=sha(Path(__file__).with_name(name)) for key,name in fingerprints.items()):
   old_phase.replace(old_phase.with_name('phase_stack_complete.previous.json'))
   print('【相位续行】程序更新，旧完成标记已保留，后段自动重新处理。',flush=True)
 run('preprocess and scan',lambda:preprocess(root,expected_stations=expected))
 run('front MUSIC direction and velocity',lambda:execute(project,'music'))
 front=json.loads((project/'output/music_front/music_complete.json').read_text(encoding='utf-8-sig'))
 if not front['corrected_points']:print('【筛窗先验】前段无合格连续MUSIC频带，后段只用原CCF候选作先验，重新实测低频。',flush=True)
 for name in ('review','verify'):run(name,lambda name=name:execute(project,name))
 frozen=root/'before_selection.json';fingerprint=root/'before_selection_sources.json';source={str(p.relative_to(project)):sha(p) for p in (project/'config.json',project/'output/music_front/dispersion.csv',project/'output/conditional_dispersion.csv',project/'inversion_settings.json')}
 if frozen.exists() and fingerprint.exists() and json.loads(fingerprint.read_text(encoding='utf-8-sig'))==source:before=json.loads(frozen.read_text(encoding='utf-8-sig'))
 else:
  if (project/'output/phase_guided/phase_stack_complete.json').exists():raise ValueError('修正前输入变化，请使用新目录')
  before=run('select before',lambda:auto_selection(project,max_cases=3,modes=(0,),minimum_wavelengths=3.))
  if not before['cases']:
   before=run('select screening priors only',lambda:auto_selection(project,max_cases=3,modes=(0,),minimum_wavelengths=0.,screening_priors_only=not bool(front['corrected_points'])));before['selection_notes'].update(screening_prior_only=True,ccf_direct_inversion_qualified=False,before_structure_is_conditional_comparison=True);write_json(root/'prior_only_selection.json',before)
  write_json(frozen,before);write_json(fingerprint,source)
 if not before['cases']:raise ValueError('无合格非增速先验曲线')
 ids=sorted({b['curve_id'] for case in before['cases'] for b in case['bands']});from .candidates import read_curves
 q=pd.read_csv(project/'output/conditional_dispersion.csv') if before['selection_notes'].get('screening_priors_only') else read_curves(project)[0];q[q.curve_id.isin(ids)].to_csv(root/'before_selected_curves.csv',index=False)
 options=dict(curve_ids=ids,repeat_cycles=200.,slide_cycles=200.,taper='hann',frequency_method='exact',screening_mode='curve',tolerance_deg=15.,monotonic_only=True,direction=dict(mode='filter',max_pair_angle_deg=45.,direction_step_deg=2.,velocity_points=128,velocity_min_m_s=200.,velocity_max_m_s=5000.,frequency_half_width_hz=0.,frequency_bins=1,music_front_cycles=200.),reconstruction=dict(enabled=True,lower_frequency_hz=.08,step_hz=.01))
 write_json(project/'phase_settings.json',options);phase=project/'output/phase_guided'
 if not completed_phase(phase,options):run('Goertzel low-frequency reconstruction',lambda:run_phase_stack(project,**options))
 run('phase-audit',lambda:audit(phase))
 after=run('select after',lambda:auto_selection(project,max_cases=3,modes=(0,),minimum_wavelengths=3.))
 if not after['cases']:
  run('report',lambda:execute(project,'report'))
  write_json(state,dict(run_meta,complete=True,stages=stages,scientific_outcome='no_qualified_corrected_bands',inversion_skipped=True,physical_dispersion_and_structure_obtained=False,report=str(phase/'phase_report.html')))
  print('【科学结果】后段没有合格连续实测频带；处理与审计完成，自动反演跳过。',flush=True);return
 result=run('paired inversion',lambda:paired_compare(project,after,root/'comparison'));write_json(state,dict(run_meta,complete=True,stages=stages,comparison=result,extrapolated_predictions_in_inversion=False))
 print('【完成】'+str(root/'comparison/report.html'),flush=True)
if __name__=='__main__':main()
