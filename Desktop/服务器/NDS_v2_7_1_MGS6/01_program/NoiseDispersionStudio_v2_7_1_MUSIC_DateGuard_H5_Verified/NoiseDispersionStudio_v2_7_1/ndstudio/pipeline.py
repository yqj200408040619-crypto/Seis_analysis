"""CLI/GUI worker stages. Every completed stage has a validated marker."""
from pathlib import Path
import json,time,traceback
from dispersion.common import load_config,write_json,invalidate_models
from .project import import_hdf,suggest_geometry
from .candidates import review,auto_selection
from .inversion import invert_project

STAGES=['import','prepare','scan','controls','families','music','phase-stack','review','report','verify']
MARKERS={'prepare':'prepare_complete.json','scan':'scan_complete.json','controls':'controls_complete.json','families':'families_complete.json','verify':'verification.json'}

def event(stage,status,**more):
 print('@@EVENT '+json.dumps({'stage':stage,'status':status,**more},ensure_ascii=False),flush=True)

def execute(project,stage='all',auto_invert=False,max_cases=12,modes=(0,),minimum_wavelengths=3,scan_top_k=8):
 from .trace import Trace
 root=Path(project).resolve()
 phase=json.loads((root/'phase_settings.json').read_text(encoding='utf-8-sig')) if (root/'phase_settings.json').exists() else {}
 used={'analysis':load_config(root/'config.json'),'phase':phase,'requested_stage':stage,'auto_invert':auto_invert,'minimum_wavelengths':minimum_wavelengths,'max_cases':max_cases,'modes':list(modes)}
 with Trace(root/'output',used,'workflow') as trace:
  return _execute(project,stage,auto_invert,max_cases,modes,minimum_wavelengths,scan_top_k,trace)

def _execute(project,stage,auto_invert,max_cases,modes,minimum_wavelengths,scan_top_k,trace):
 root=Path(project).resolve();out=root/'output';out.mkdir(exist_ok=True)
 initial=load_config(root/'config.json')
 for key,relative in [('output','output'),('data_root','canonical'),('stations_csv','stations.csv')]:
  if Path(initial[key])!=root/relative:raise ValueError('Studio项目的'+key+'须固定为'+relative+'；请在项目内使用标准路径')
 spec=json.loads((root/'project.json').read_text(encoding='utf-8-sig'))
 if stage=='all' and spec.get('analysis_mode')=='curves_only':
  stages=['review','autoselect','invert','report'] if auto_invert else ['review','report']
 else:
  stages=STAGES if stage=='all' else [stage]
  if auto_invert:stages=stages+['autoselect','invert','report']
 log=[]
 for step in stages:
  t=time.time()
  from .trace import stage_parameters
  with trace.stage(step,stage_parameters(root,step,initial)) as stage_result:
   try:
    if step=='import':
     spec=json.loads((root/'project.json').read_text(encoding='utf-8-sig'))
     if spec.get('analysis_mode')=='curves_only':raise ValueError('当前是曲线反演项目；波形分析请新建项目')
     invalidate_models(out);(out/'structure_layers.csv').unlink(missing_ok=True)
     import_hdf(root)
    elif step=='autoselect':auto_selection(root,max_cases=max_cases,modes=modes,minimum_wavelengths=minimum_wavelengths)
    elif step=='invert':
     chosen=json.loads((root/'selection.json').read_text(encoding='utf-8-sig'))
     if auto_invert and not chosen.get('cases'):
      invalidate_models(out);reason='没有通过背景与距离/波长要求的自动反演曲线组'
      write_json(out/'inversion_skipped.json',{'reason':reason,'not_a_processing_error':True})
      log.append({'stage':step,'status':'skipped','reason':reason,'elapsed_s':time.time()-t});event(step,'skipped',message=reason)
      stage_result.update(status='skipped',reason=reason)
      write_json(root/'last_run.json',{'stages':log,'all_requested_stages_complete':len(log)==len(stages)});continue
     (out/'inversion_skipped.json').unlink(missing_ok=True);invert_project(root)
     marker=out/'phase_guided/phase_stack_complete.json'
     if marker.exists() and json.loads(marker.read_text(encoding='utf-8-sig')).get('reconstruction',{}).get('enabled'):
      from .reconstruction import paired_compare
      stage_result['paired_comparison']=paired_compare(root,chosen)
    elif step=='phase-stack':
     from .phase_stack import run_phase_stack
     import pandas as pd
     if not len(pd.read_csv(out/'conditional_dispersion.csv')):
      stage_result.update(status='skipped',reason='没有候选频散曲线；请查看扫描诊断并调参')
     else:
      from .music_front import run_music_front
      run_music_front(root)
      phase=json.loads((root/'phase_settings.json').read_text(encoding='utf-8-sig')) if (root/'phase_settings.json').exists() else {}
      allowed={'curve_ids','start','end','repeat_cycles','slide_cycles','tolerance_deg','range_source','taper','time_mode','stabilization_s','frequency_method','amplitude_rejection','split_branches','screening_mode','direction','monotonic_only','prior','reconstruction'}
      phase.setdefault('direction',{'mode':'filter','method':'music'})
      if set(phase)-allowed:raise ValueError('未知phase_settings参数: '+str(sorted(set(phase)-allowed)))
      r=run_phase_stack(root,**phase);stage_result.update(accepted_window_evaluations=r['accepted_window_evaluations'],handoff_points=r.get('inversion_handoff_points',0),trace=str(out/'phase_guided/phase_trace.json'))
    elif step=='music':
     from .music_front import run_music_front
     result=run_music_front(root)
     stage_result.update({key:result[key] for key in ('complete','method','corrected_points','corrected_curves','front_pair_angle_limit_applied')})
    elif step=='review':review(root)
    else:
     c=load_config(root/'config.json')
     if step in ('prepare','scan','controls','families'):
      invalidate_models(out)
      (out/'structure_layers.csv').unlink(missing_ok=True)
      for name in ('10_structure_models.png','11_selected_dispersion.png','12_vs_structure.png'):(out/'figures'/name).unlink(missing_ok=True)
      (out/'inversion_skipped.json').unlink(missing_ok=True)
      # Explicitly invalidate rerun markers before a worker can be cancelled.
      downstream=['prepare','scan','controls','families']
      for name in downstream[downstream.index(step):]:
       (out/MARKERS[name]).unlink(missing_ok=True)
      (out/'music_front/music_complete.json').unlink(missing_ok=True)
      for name in ['verification.json','inversion_complete.json']:(out/name).unlink(missing_ok=True)
      for name in ['phase_stack_complete.json','phase_report.html','inversion_curves.csv','inversion_catalog.csv']:(out/'phase_guided'/name).unlink(missing_ok=True)
     if step=='prepare':
      from dispersion.prepare import prepare
      prepare(c);write_json(root/'geometry_advice.json',suggest_geometry(root))
     elif step=='scan':
      from dispersion.scan import scan
      scan(c)
      from .scan_diagnostics import write_scan_diagnostics
      write_scan_diagnostics(c,scan_top_k)
     elif step=='scan-diagnostics':
      from .scan_diagnostics import write_scan_diagnostics
      write_scan_diagnostics(c,scan_top_k)
     elif step=='controls':
      from dispersion.controls import controls
      controls(c)
     elif step=='families':
      from dispersion.families import families
      families(c)
      from .music_front import run_music_front
      run_music_front(root);review(root)
     elif step=='report':
      spec=json.loads((root/'project.json').read_text(encoding='utf-8-sig'))
      if spec.get('analysis_mode')!='curves_only':
       from dispersion.products import report
       report(c)
      from .reports import studio_report
      studio_report(root)
     elif step=='verify':
      from dispersion.products import verify
      verify(c)
     else:raise ValueError('未知阶段: '+step)
    from .trace import file_summary
    stage_result.update(file_summary(root,step))
    log.append({'stage':step,'status':stage_result.get('status','complete'),'elapsed_s':time.time()-t})
    write_json(root/'last_run.json',{'stages':log,'all_requested_stages_complete':len(log)==len(stages)})
   except Exception as exc:
    write_json(root/'last_run.json',{'stages':log,'failed_stage':step,'error':str(exc),'all_requested_stages_complete':False});raise
 return log
