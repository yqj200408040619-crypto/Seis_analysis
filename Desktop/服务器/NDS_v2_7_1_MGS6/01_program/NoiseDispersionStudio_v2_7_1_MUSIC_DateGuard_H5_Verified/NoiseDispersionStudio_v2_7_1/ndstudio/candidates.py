"""Review measured candidate families and build explicit conditional inversion cases."""
from pathlib import Path
import json,itertools
import numpy as np
import pandas as pd
from dispersion.common import csv,write_json

INVERSION_DEFAULTS={
 'layers':[1,2],'starts':8,'max_evaluations':180,'vp_vs_ratio':2.,'density_g_cm3':2.,
 'thickness_bounds_m':[2.,1500.],'vs_bounds_m_s':[120.,3500.],
 'phase_rms_max_fraction':.03,'phase_point_max_fraction':.05,'group_median_max_fraction':.15,
 'group_fit_weight':.25,'max_fit_points_per_band':32,'finite_difference_log_step':.003,
 'require_independent_solver':False,'independent_max_difference_m_s':.5,
 'allow_background_unverified':False,'seed':9075}

def project_settings(root):
 path=Path(root)/'inversion_settings.json'
 return dict(INVERSION_DEFAULTS)| (json.loads(path.read_text(encoding='utf-8-sig')) if path.exists() else {})

def read_curves(project):
 root=Path(project);out=root/'output'
 spec=json.loads((root/'project.json').read_text(encoding='utf-8-sig'))
 if spec.get('analysis_mode')=='curves_only':
  from dispersion.common import sha
  path=root/spec['external_curves_file']
  if sha(path)!=spec['external_curve_sha256']:raise ValueError('导入的曲线文件已改变，请重新导入')
  return pd.read_csv(path),pd.read_csv(root/'external_catalog.csv')
 if not (out/'families_complete.json').exists():raise ValueError('候选频散阶段未完成，请先完成分析')
 from dispersion.common import load_config,ensure_ready
 phase=out/'phase_guided'; marker=phase/'phase_stack_complete.json'
 if marker.exists():
  from dispersion.common import sha
  m=json.loads(marker.read_text(encoding='utf-8-sig'));export=json.loads((phase/'phase_export.json').read_text(encoding='utf-8-sig'))
  current=m.get('settings_sha256')==sha(root/'config.json') and export.get('source_files',{}).get('conditional_dispersion.csv')==sha(out/'conditional_dispersion.csv')
  if m.get('reconstruction',{}).get('enabled'):current &= m.get('reconstruction_module_sha256')==sha(Path(__file__).with_name('reconstruction.py')) and m.get('module_sha256')==sha(Path(__file__).with_name('phase_stack.py')) and m.get('stations_csv_sha256')==sha(root/'stations.csv')
  if m.get('music_module_sha256'):current &= m['music_module_sha256']==sha(Path(__file__).with_name('music.py'))
  if m.get('phase_settings_sha256') and (root/'phase_settings.json').exists():current &= m['phase_settings_sha256']==sha(root/'phase_settings.json')
  if not current and m.get('direction',{}).get('mode') in ('filter','joint'):raise ValueError('方向/相位结果已过期；参数或输入改变后，请重跑方向与相位阶段再进行反演。')
  if current and m.get('direction',{}).get('mode') in ('filter','joint') and (phase/'inversion_curves.csv').exists():
   print('【候选交接】使用方向约束后的连续叠加频散：'+str(phase/'inversion_curves.csv'),flush=True)
   curves=pd.read_csv(phase/'inversion_curves.csv');curves.attrs['source_file']='output/phase_guided/inversion_curves.csv'
   return curves,pd.read_csv(phase/'inversion_catalog.csv')
 ensure_ready(load_config(root/'config.json'),'families_complete.json')
 if (out/'music_front/music_complete.json').exists():
  from .music_front import current
  if not current(root):raise ValueError('前段MUSIC结果已过期，请重新运行music阶段')
  curves=pd.read_csv(out/'music_front/dispersion.csv');curves.attrs['source_file']='output/music_front/dispersion.csv'
  return curves,pd.read_csv(out/'music_front/catalog.csv')
 return pd.read_csv(out/'conditional_dispersion.csv'),pd.read_csv(out/'curve_catalog.csv')

def import_curve_csv(project,path):
 from .hdf import truth
 from dispersion.common import sha
 root=Path(project);frame=pd.read_csv(path)
 required=['curve_id','segment_id','cohort','pair_index','pair','branch','frequency_hz','phase_velocity_m_s','group_velocity_m_s','constant_cycle','source_phase_rad','distance_over_wavelength','instrument_response_corrected','background_discriminated','source_phase_verified','is_model_prediction']
 if not set(required)<=set(frame):raise ValueError('曲线CSV缺少字段: '+str(sorted(set(required)-set(frame))))
 for col in ['instrument_response_corrected','background_discriminated','source_phase_verified','is_model_prediction']:frame[col]=frame[col].map(truth)
 rows=[]
 for cid,q in frame.groupby('curve_id',sort=False):
  for field in ['segment_id','cohort','pair_index','pair','branch','constant_cycle','source_phase_rad']:
   if q[field].nunique()!=1:raise ValueError(cid+'存在不一致的'+field)
  if q.frequency_hz.duplicated().any():raise ValueError(cid+'有重复频点')
  rows.append(dict(curve_id=cid,segment_id=q.segment_id.iloc[0],cohort=q.cohort.iloc[0],pair_index=int(q.pair_index.iloc[0]),pair=q.pair.iloc[0],branch=q.branch.iloc[0],
    constant_cycle=int(q.constant_cycle.iloc[0]),source_phase_rad=float(q.source_phase_rad.iloc[0]),fmin_hz=float(q.frequency_hz.min()),fmax_hz=float(q.frequency_hz.max()),points=len(q),
    background_discriminated=bool(q.background_discriminated.all()),compatible_other_pairs=0))
 csv(root/'external_curves.csv',frame);csv(root/'external_catalog.csv',rows)
 spec=json.loads((root/'project.json').read_text(encoding='utf-8-sig'));spec.update(analysis_mode='curves_only',external_curves_file='external_curves.csv',external_curve_sha256=sha(root/'external_curves.csv'))
 write_json(root/'project.json',spec)
 return review(root)

def review(project):
 root=Path(project);out=root/'output';curves,cat=read_curves(root)
 cat=rank_candidates(curves,cat)
 csv(out/'candidate_review.csv',cat)
 source=root/'hdf_import_audit.json'
 if source.exists():
  info=json.loads(source.read_text(encoding='utf-8-sig'));write_json(out/'input_source_summary.json',{'raw_hdf_records':len(info['records']),'sources_unchanged_at_import':True,
   'source_files_opened_read_only':True,'utc_source_not_file_mtime':True})
 return cat

def rank_candidates(curves,cat):
 if len(curves):
  stats=curves.groupby('curve_id').agg(minimum_phase_velocity_m_s=('phase_velocity_m_s','min'),maximum_phase_velocity_m_s=('phase_velocity_m_s','max'),minimum_distance_over_wavelength=('distance_over_wavelength','min'))
  cat=cat.merge(stats,left_on='curve_id',right_index=True,validate='one_to_one')
  if 'pair' not in cat:cat=cat.merge(curves.groupby('curve_id').pair.first(),left_on='curve_id',right_index=True,validate='one_to_one')
 else:
  for col in ('minimum_phase_velocity_m_s','maximum_phase_velocity_m_s','minimum_distance_over_wavelength'):cat[col]=pd.Series(dtype=float)
 cat['bandwidth_hz']=cat.fmax_hz-cat.fmin_hz
 cat['review_priority']=cat.background_discriminated & (cat.minimum_distance_over_wavelength>=3) & (cat.compatible_other_pairs>=2)
 cat['unique_physical_dispersion_verified']=False
 cat=cat.sort_values(['review_priority','bandwidth_hz','compatible_other_pairs','curve_id'],ascending=[False,False,False,True])
 return cat


def auto_selection(project,max_cases=12,max_bands=3,modes=(0,),minimum_wavelengths=3,screening_priors_only=False):
 if not np.isfinite(minimum_wavelengths) or minimum_wavelengths<0:raise ValueError('距离/波长最低要求须为有限非负数')
 root=Path(project);out=root/'output'
 if screening_priors_only:
  curves=pd.read_csv(out/'conditional_dispersion.csv');catalog=pd.read_csv(out/'curve_catalog.csv');ranked=rank_candidates(curves,catalog)
 else:
  curves,catalog=read_curves(root);ranked=review(root)
 projected=ranked.get('measurement_method',pd.Series('',index=ranked.index)).eq('direction_projected_goertzel') & ranked.get('direction_quality_verified',pd.Series(False,index=ranked.index)).eq(True)
 ranked=ranked[ranked.background_discriminated & ((ranked.minimum_distance_over_wavelength>=minimum_wavelengths)|projected)]
 spec=json.loads((root/'project.json').read_text(encoding='utf-8-sig'))
 if spec.get('exclude_increasing_velocity',False):
  from .curve_prior import monotone_curves
  _,trend=monotone_curves(curves)
  csv(out/'automatic_curve_trend_screen.csv',trend)
  ranked=ranked[ranked.curve_id.isin(trend.loc[trend.accepted,'curve_id'])]
 if max_cases<1 or max_bands<1 or not modes:raise ValueError('组合数量、每组曲线数和模态列表不能为空')
 link_path=out/'music_front/cross_pair_compatible_curves.csv' if not screening_priors_only and (out/'music_front/music_complete.json').exists() else out/'cross_pair_compatible_curves.csv'
 links=pd.read_csv(link_path) if link_path.exists() else pd.DataFrame(columns=['curve_a','curve_b'])
 if 'parent_curve_id' in ranked and len(links):
  mapping={p:g.curve_id.tolist() for p,g in ranked.groupby('parent_curve_id')}
  front=out/'music_front/catalog.csv'
  if front.exists():
   fc=pd.read_csv(front)
   for parent,g in fc.groupby('parent_curve_id'):
    mapping[parent]=[cid for middle in g.curve_id for cid in mapping.get(middle,[middle]) if cid in set(ranked.curve_id)]
  expanded=[]
  for r in links.itertuples():
   for a,b in itertools.product(mapping.get(r.curve_a,[]),mapping.get(r.curve_b,[])):expanded.append(dict(curve_a=a,curve_b=b))
  links=pd.DataFrame(expanded,columns=['curve_a','curve_b'])
 valid=set(ranked.curve_id);adj={k:set() for k in valid}
 for r in links.itertuples():
  if r.curve_a in valid and r.curve_b in valid:adj[r.curve_a].add(r.curve_b);adj[r.curve_b].add(r.curve_a)
 meta=ranked.set_index('curve_id');rank={cid:i for i,cid in enumerate(ranked.curve_id)};cases=[];seen=set()
 for cid in ranked.curve_id:
  group=[cid];pairs={int(meta.loc[cid,'pair_index'])}
  for other in sorted(adj[cid],key=lambda k:rank[k]):
   if len(group)>=max_bands:break
   row=meta.loc[other]
   if int(row.pair_index) in pairs:continue
   if not all(other in adj[g] for g in group):continue
   if row.cohort!=meta.loc[cid,'cohort'] or not np.isclose(row.source_phase_rad,meta.loc[cid,'source_phase_rad']):continue
   group.append(other);pairs.add(int(row.pair_index))
   if len(group)>=max_bands:break
  key=tuple(sorted(group))
  if key in seen:continue
  seen.add(key)
  for mode in modes:
   if not isinstance(mode,int) or mode<0:raise ValueError('Rayleigh模态必须为非负整数，0代表基阶')
   cases.append({'name':f'AUTO_{len(cases)+1:03d}_M{mode}','bands':[{'curve_id':g,'mode':mode} for g in group],
     'assumptions':{'common_1d_structure':True,'rayleigh_modes_assumed_not_verified':True,'source_phase_and_cycles_conditional':True},
     'origin':'automatic_pairwise_compatible_group'})
   if len(cases)>=max_cases:break
  if len(cases)>=max_cases:break
 result={'settings':project_settings(root),'cases':cases,'selection_notes':{
   'automatic_selection_is_exploratory':True,'screening_priors_only':bool(screening_priors_only),'max_cases':max_cases,'max_bands_per_case':max_bands,
   'minimum_distance_over_wavelength':minimum_wavelengths,'all_other_candidates_retained':True,
   'projected_goertzel_uses_array_direction_independent_windows_and_shifted_null_instead_of_CCF_wavelength_floor':True,
   'groups_are_pairwise_cliques_not_transitive_components':True,'rank_is_not_a_posterior_probability':True}}
 write_json(root/'selection.json',result)
 return result

def manual_selection(project,curve_ids,mode=0,joint=False):
 root=Path(project);_,catalog=read_curves(root);valid=set(catalog.curve_id)
 ids=list(dict.fromkeys(curve_ids))
 if not ids or any(k not in valid for k in ids):raise ValueError('请选择存在的候选曲线')
 path=root/'selection.json';selection=json.loads(path.read_text(encoding='utf-8-sig')) if path.exists() else {'settings':project_settings(root),'cases':[]}
 groups=[ids] if joint else [[k] for k in ids]
 for group in groups:
  selection['cases'].append({'name':f'MANUAL_{len(selection["cases"])+1:03d}_M{mode}',
    'bands':[{'curve_id':k,'mode':int(mode)} for k in group],
    'assumptions':{'common_1d_structure':True,'rayleigh_modes_assumed_not_verified':True,'source_phase_and_cycles_conditional':True},'origin':'manual'})
 write_json(path,selection);return selection
