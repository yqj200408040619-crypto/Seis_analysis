"""Multigroup, multimode conditional Rayleigh inversion with explicit rejected attempts."""
from pathlib import Path
import json,warnings,hashlib
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from .candidates import INVERSION_DEFAULTS
from dispersion.inversion import model_from_log,forward,independent_phase
from dispersion.common import write_json,csv

class SelectionError(ValueError):pass

def selection_digest(selection):
 return hashlib.sha256(json.dumps(selection,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()

def code_digest():
 digest=hashlib.sha256()
 base=Path(__file__).resolve().parent.parent
 files=list((base/'ndstudio').glob('*.py'))+list((base/'core/dispersion').glob('*.py'))+[base/'run.py']
 for p in sorted(files):digest.update(p.relative_to(base).as_posix().encode());digest.update(p.read_bytes())
 return digest.hexdigest()

def inversion_current(project):
 root=Path(project);marker=root/'output/inversion_complete.json';selection=root/'selection.json'
 if not marker.exists() or not selection.exists():return False
 try:
  m=json.loads(marker.read_text(encoding='utf-8-sig'));s=json.loads(selection.read_text(encoding='utf-8-sig'))
  if m.get('selection_sha256')!=selection_digest(s) or m.get('studio_code_sha256')!=code_digest():return False
  if 'input_curve_file' in m:
   from dispersion.common import sha
   if sha(root/m['input_curve_file'])!=m['input_curve_sha256']:return False
   if sha(root/'config.json')!=m.get('project_config_sha256'):return False
   if m.get('phase_settings_sha256') and sha(root/'phase_settings.json')!=m['phase_settings_sha256']:return False
   if m.get('phase_completion_sha256') and sha(root/'output/phase_guided/phase_stack_complete.json')!=m['phase_completion_sha256']:return False
  return True
 except (OSError,ValueError):return False

def settings_from(custom):
 extra=set(custom)-set(INVERSION_DEFAULTS)
 if extra:raise SelectionError('未知反演设置: '+str(sorted(extra)))
 settings=dict(INVERSION_DEFAULTS);settings.update(custom)
 if not settings['layers'] or any(not isinstance(n,int) or not 1<=n<=4 for n in settings['layers']):raise SelectionError('layers需为1–4个有限层的列表')
 if settings['starts']<1 or settings['max_evaluations']<1 or settings['max_fit_points_per_band']<8:raise SelectionError('多起点数、迭代数和拟合采样数无效')
 if settings['vp_vs_ratio']<=np.sqrt(4/3) or settings['density_g_cm3']<=0:raise SelectionError('弹性Vp/Vs或密度无效')
 for k in ['thickness_bounds_m','vs_bounds_m_s']:
  b=settings[k]
  if len(b)!=2 or not 0<b[0]<b[1]:raise SelectionError(k+'需为两个递增的正数')
 for k in ['phase_rms_max_fraction','phase_point_max_fraction','group_median_max_fraction','independent_max_difference_m_s','finite_difference_log_step']:
  if settings[k]<=0:raise SelectionError(k+'必须为正数')
 if settings['group_fit_weight']<0:raise SelectionError('group_fit_weight不能为负')
 return settings

def validate_cases(curves,selection,minimum_points=8,minimum_band_hz=.25):
 from .hdf import truth
 required={'curve_id','cohort','pair_index','frequency_hz','phase_velocity_m_s','group_velocity_m_s','constant_cycle','source_phase_rad','instrument_response_corrected','background_discriminated','source_phase_verified','is_model_prediction'}
 if not required<=set(curves):raise SelectionError('曲线表缺少字段: '+str(sorted(required-set(curves))))
 curves=curves.copy()
 for column in ['instrument_response_corrected','background_discriminated','source_phase_verified','is_model_prediction']:
  curves[column]=curves[column].map(truth)
 settings=settings_from(selection.get('settings',{}));cases=selection.get('cases',[])
 if not cases:raise SelectionError('没有反演假设组，请先自动组合或人工选择')
 result=[];names=set()
 for case in cases:
  name=str(case.get('name',''))
  if not name or name in names:raise SelectionError('假设组名称不能为空或重复')
  names.add(name);bands=[];old=[];phis=[];cohorts=set()
  if not case.get('bands'):raise SelectionError(name+'没有曲线')
  for band in case['bands']:
   if 'mode' not in band or not isinstance(band['mode'],int) or isinstance(band['mode'],bool) or band['mode']<0:raise SelectionError('每条曲线必须显式指定Rayleigh模态，0为基阶')
   q=curves[curves.curve_id==band['curve_id']].copy().sort_values('frequency_hz')
   if 'fmin_hz' in band:q=q[q.frequency_hz>=band['fmin_hz']]
   if 'fmax_hz' in band:q=q[q.frequency_hz<=band['fmax_hz']]
   if len(q)<minimum_points or q.frequency_hz.max()-q.frequency_hz.min()<minimum_band_hz-1e-8:raise SelectionError('曲线不存在或所选频段太短: '+band['curve_id'])
   if q.frequency_hz.duplicated().any():raise SelectionError('每条曲线不能有重复频点')
   if q.pair_index.nunique()!=1 or q.cohort.nunique()!=1:raise SelectionError('一条曲线须有固定台对和时间组')
   if not np.isfinite(q[['frequency_hz','phase_velocity_m_s','group_velocity_m_s']]).all().all() or (q[['frequency_hz','phase_velocity_m_s','group_velocity_m_s']]<=0).any().any():raise SelectionError('反演观测中存在缺测、非正或非有限值')
   if not q.instrument_response_corrected.all():raise SelectionError('结构反演要求响应校正或明确声明物理量输入')
   if not settings['allow_background_unverified'] and not q.background_discriminated.all():raise SelectionError('所选曲线未通过错时背景分离')
   if q.is_model_prediction.any():raise SelectionError('模型预测不能混入观测反演')
   if q.constant_cycle.nunique()!=1 or q.source_phase_rad.nunique()!=1:raise SelectionError('每个连续段必须固定整数周期和源相位')
   cycle=float(q.constant_cycle.iloc[0])
   if not np.isfinite(cycle) or cycle<0 or cycle!=int(cycle):raise SelectionError('周期N必须为非负整数')
   for prior in old:
    if prior.pair_index.iloc[0]==q.pair_index.iloc[0] and set(np.round(prior.frequency_hz,9))&set(np.round(q.frequency_hz,9)):
     raise SelectionError('同一台对的重叠候选不能在同一假设组重复计权；请分别建组')
   old.append(q);phis.append(float(q.source_phase_rad.iloc[0]));cohorts.add(q.cohort.iloc[0])
   sample=np.unique(np.linspace(0,len(q)-1,min(len(q),settings['max_fit_points_per_band'])).astype(int))
   bands.append({'table':q,'fit':q.iloc[sample],'mode':band['mode'],'curve_id':band['curve_id']})
  if len(cohorts)!=1:raise SelectionError('联合组须来自同一固定时间组；全天和夜间可能共用观测')
  if not np.allclose(phis,phis[0]):raise SelectionError('联合组的源相位假设不同，应分别建组')
  result.append({'name':name,'bands':bands,'assumptions':case.get('assumptions',{'common_1d_structure':True}),'origin':case.get('origin','explicit')})
 return settings,result

def invert_tables(curves,selection,output,minimum_points=8,minimum_band_hz=.25):
 out=Path(output);out.mkdir(parents=True,exist_ok=True)
 (out/'inversion_complete.json').unlink(missing_ok=True)
 for name in ['05_optional_conditional_models.png','10_structure_models.png','11_selected_dispersion.png','12_vs_structure.png']:(out/'figures'/name).unlink(missing_ok=True)
 settings,cases=validate_cases(curves,selection,minimum_points,minimum_band_hz)
 independent_errors=(ValueError,RuntimeError,ZeroDivisionError)
 try:
  import pysurf96
  from pysurf96.wrapper import Surf96Error
  independent_errors=independent_errors+(Surf96Error,);independent_available=True
 except ImportError:independent_available=False
 if settings['require_independent_solver'] and not independent_available:raise ImportError('要求独立求解器，但pysurf96未安装；请安装可选依赖或明确采用未独立复验的探索性设置')
 models=[];attempts=[];observed=[];predictions=[];rng=np.random.default_rng(settings['seed'])
 # A killed attempt must not leave an old successful completion marker.
 (out/'inversion_complete.json').unlink(missing_ok=True)
 def checkpoint():
  write_json(out/'conditional_models.json',{'settings':settings,'models':models,'independent_solver_available':independent_available,
      'unique_structure_identified':False,'posterior_probabilities_calibrated':False})
  csv(out/'inversion_attempts.csv',attempts,columns=None if attempts else ['case','finite_layers','attempt','accepted','reason'])
 for case in cases:
  bands=case['bands'];base=np.concatenate([b['fit'].phase_velocity_m_s.to_numpy() for b in bands]);nb=len(bands)
  for nl in settings['layers']:
   hlo,hhi=settings['thickness_bounds_m'];vlo,vhi=settings['vs_bounds_m_s'];lower=np.log([hlo]*nl+[vlo]*(nl+1));upper=np.log([hhi]*nl+[vhi]*(nl+1))
   yphase=np.concatenate([b['fit'].phase_velocity_m_s.to_numpy() for b in bands]);ygroup=np.concatenate([b['fit'].group_velocity_m_s.to_numpy() for b in bands])
   weights=np.concatenate([np.full(len(b['fit']),1/np.sqrt(nb*len(b['fit']))) for b in bands]);gw=settings['group_fit_weight']
   size=len(yphase)*(2 if gw else 1)
   def residual(x):
    try:
     model=model_from_log(x,nl,settings['vp_vs_ratio'],settings['density_g_cm3'])
     pp=np.concatenate([forward(model,b['fit'].frequency_hz.to_numpy(),b['mode']) for b in bands])
     res=(pp/yphase-1)*weights
     if gw:
      gg=np.concatenate([forward(model,b['fit'].frequency_hz.to_numpy(),b['mode'],True) for b in bands]);res=np.r_[res,(gg/ygroup-1)*weights*gw]
     return res
    except (ValueError,ZeroDivisionError,RuntimeError):return np.full(size,10./np.sqrt(size))
   for attempt in range(settings['starts']):
    thick=np.exp(rng.uniform(np.log(hlo),np.log(hhi),nl));vs=np.median(base)*rng.uniform(.65,1.6,nl+1)
    vs[-1]=max(float(base.max())*rng.uniform(1.05,1.5),vs[-1]);x0=np.clip(np.log(np.r_[thick,vs]),lower+1e-6,upper-1e-6)
    fit=least_squares(residual,x0,bounds=(lower,upper),max_nfev=settings['max_evaluations'],diff_step=settings['finite_difference_log_step'],ftol=2e-7,xtol=2e-7,gtol=2e-7)
    diag={'case':case['name'],'finite_layers':nl,'attempt':attempt,'optimizer_status':int(fit.status),'evaluations':int(fit.nfev),'accepted':False,'reason':'not_evaluated'};attempts.append(diag)
    if not fit.success:
     diag['reason']='optimizer_not_converged';checkpoint();continue
    model=model_from_log(fit.x,nl,settings['vp_vs_ratio'],settings['density_g_cm3']);pps=[];ggs=[];errs=[];gerrs=[]
    try:
     for band in bands:
      q=band['table'];freq=q.frequency_hz.to_numpy();pp=forward(model,freq,band['mode']);gg=forward(model,freq,band['mode'],True)
      pps.append(pp);ggs.append(gg);errs.append(pp/q.phase_velocity_m_s.to_numpy()-1);gerrs.append(float(np.median(abs(gg/q.group_velocity_m_s.to_numpy()-1))))
    except (ValueError,ZeroDivisionError,RuntimeError) as exc:
     diag.update(reason='missing_modal_root',failure=str(exc));checkpoint();continue
    rms=float(np.sqrt(np.mean([np.mean(e*e) for e in errs])));maximum=float(max(np.max(abs(e)) for e in errs));gmaximum=float(max(gerrs))
    diag.update(phase_weighted_rms_fraction=rms,phase_max_fraction=maximum,maximum_group_median_fraction=gmaximum)
    if rms>settings['phase_rms_max_fraction'] or maximum>settings['phase_point_max_fraction']:
     diag['reason']='phase_misfit';checkpoint();continue
    if gmaximum>settings['group_median_max_fraction']:
     diag['reason']='group_misfit';checkpoint();continue
    independent_error=None;indok=False
    if independent_available:
     try:
      independent_error=float(max(np.max(abs(independent_phase(model,b['table'].frequency_hz.to_numpy(),b['mode'])-pp)) for b,pp in zip(bands,pps)))
      indok=independent_error<=settings['independent_max_difference_m_s']
     except independent_errors as exc:diag['independent_failure']=str(exc)
     if not indok:
      diag.update(reason='independent_solver_mismatch_or_missing_root',independent_max_difference_m_s=independent_error);checkpoint();continue
    elif settings['require_independent_solver']:
     diag['reason']='independent_solver_unavailable';checkpoint();continue
    if any(m['case']==case['name'] and m['finite_layers']==nl and np.max(abs(np.log(np.asarray(m['parameters'])/np.exp(fit.x))))<.025 for m in models):
     diag['reason']='duplicate_structure';checkpoint();continue
    mid=f'M{len(models)+1:05d}';diag.update(accepted=True,reason='conditional_compatible',model_id=mid,independent_max_difference_m_s=independent_error)
    bound_distance=np.minimum(fit.x-lower,upper-fit.x);active=bool((bound_distance<.02).any())
    models.append({'model_id':mid,'case':case['name'],'finite_layers':nl,'parameters':np.exp(fit.x).tolist(),
      'thickness_m':(model[:-1,0]*1000).tolist(),'vs_m_s':(model[:,2]*1000).tolist(),
      'vp_vs_ratio':settings['vp_vs_ratio'],'density_g_cm3':settings['density_g_cm3'],
      'phase_weighted_rms_fraction':rms,'phase_max_fraction':maximum,'maximum_group_median_fraction':gmaximum,
      'independent_solver_available':independent_available,'independent_solver_pass':indok,'independent_max_difference_m_s':independent_error,
      'parameters_near_search_bound':active,'optimizer_converged':bool(fit.success),
      'curve_ids':[b['curve_id'] for b in bands],'assumed_modes':[b['mode'] for b in bands],'assumptions':case['assumptions'],
      'source_phase_verified':bool(all(b['table'].source_phase_verified.all() for b in bands)),
      'absolute_cycles_or_modes_verified':False,'unique_structure_identified':False,'posterior_probability_calibrated':False})
    for b,pp,gg in zip(bands,pps,ggs):
     for record,pvalue,gvalue in zip(b['table'].to_dict('records'),pp,gg):
      observed.append({'model_id':mid,'case':case['name'],'assumed_mode':b['mode'],**record,
        'predicted_phase_velocity_m_s':float(pvalue),'predicted_group_velocity_m_s':float(gvalue),'is_observation':True})
     # Predictions extend only to the selected band's endpoints; never fill a measured gap as data.
     q=b['table'];grid=np.linspace(q.frequency_hz.min(),q.frequency_hz.max(),81)
     for freq in grid:
      value=None
      try:value=float(forward(model,np.asarray([freq]),b['mode'])[0])
      except (ValueError,RuntimeError,ZeroDivisionError):pass
      predictions.append({'model_id':mid,'case':case['name'],'curve_id':b['curve_id'],'assumed_mode':b['mode'],'frequency_hz':float(freq),
        'phase_velocity_m_s':value,'root_available':value is not None,'is_observation':False})
    checkpoint()
   print('invert:',case['name'],'finite layers',nl,'conditional models',len(models),flush=True)
 checkpoint()
 csv(out/'selected_observed_dispersion.csv',observed,columns=None if observed else ['model_id','case','curve_id','frequency_hz','is_observation'])
 csv(out/'model_predictions.csv',predictions,columns=['model_id','case','curve_id','assumed_mode','frequency_hz','phase_velocity_m_s','root_available','is_observation'])
 layers=[]
 for m in models:
  top=0.
  for i,vs in enumerate(m['vs_m_s']):
   thickness=m['thickness_m'][i] if i<len(m['thickness_m']) else None
   bottom=top+thickness if thickness is not None else None
   layers.append({'model_id':m['model_id'],'case':m['case'],'layer_index':i+1,
                  'layer_type':'finite' if thickness is not None else 'halfspace',
                  'top_depth_m':top,'bottom_depth_m':bottom,'thickness_m':thickness,
                  'vs_m_s':vs,'vp_m_s':vs*m['vp_vs_ratio'],'density_g_cm3':m['density_g_cm3'],
                  'independent_solver_pass':m['independent_solver_pass'],
                  'unique_structure_identified':False,'depth_resolution_established':False})
   if bottom is not None:top=bottom
 csv(out/'structure_layers.csv',layers,columns=['model_id','case','layer_index','layer_type','top_depth_m','bottom_depth_m',
        'thickness_m','vs_m_s','vp_m_s','density_g_cm3','independent_solver_pass','unique_structure_identified','depth_resolution_established'])
 write_json(out/'inversion_complete.json',{'cases':len(cases),'compatible_models':len(models),'independent_solver_available':independent_available,
   'unique_structure_identified':False,'observed_and_predicted_points_separate':True,'missing_roots_not_filled':True,'frequency_samples_not_treated_as_independent':True,
   'selection_sha256':selection_digest(selection),'studio_code_sha256':code_digest()})
 return models

def invert_project(project):
 from .candidates import read_curves
 root=Path(project);curves,_=read_curves(root);selection=json.loads((root/'selection.json').read_text(encoding='utf-8-sig'))
 if selection.get('selection_notes',{}).get('screening_priors_only'):raise SelectionError('仅供筛窗的先验不能直接反演；先完成MUSIC实测叠加并重新选择')
 from dispersion.common import load_config
 c=load_config(root/'config.json')
 result=invert_tables(curves,selection,root/'output',c['scan']['min_points'],c['scan']['min_band_hz'])
 from dispersion.common import sha
 spec=json.loads((root/'project.json').read_text(encoding='utf-8-sig'));curvefile=spec['external_curves_file'] if spec.get('analysis_mode')=='curves_only' else 'output/conditional_dispersion.csv'
 curvefile=curves.attrs.get('source_file',curvefile)
 path=root/'output/inversion_complete.json';m=json.loads(path.read_text(encoding='utf-8-sig'));m.update(input_curve_file=curvefile,input_curve_sha256=sha(root/curvefile),project_config_sha256=sha(root/'config.json'));write_json(path,m)
 if curvefile=='output/phase_guided/inversion_curves.csv':
  m['phase_completion_sha256']=sha(root/'output/phase_guided/phase_stack_complete.json')
  m['phase_settings_sha256']=sha(root/'phase_settings.json') if (root/'phase_settings.json').exists() else None
  write_json(path,m)
 return result
