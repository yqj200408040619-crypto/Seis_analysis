"""Portable project setup, HDF5 import audit, and reproducible configuration."""
from pathlib import Path
import json,os,hashlib,re
import numpy as np
import pandas as pd
import h5py
from .hdf import inspect_folder,resolve,override_for,read_waveform,MappingError,truth
from dispersion.common import DEFAULTS,merge,write_json,csv,sha,validate_config

PROFILE={'dataset_pattern':'','defaults':{},'rules':[],'overrides':{}}

def new_project(directory):
 root=Path(directory).resolve();root.mkdir(parents=True,exist_ok=True)
 if (root/'project.json').exists():return root
 config=merge(DEFAULTS,{'data_root':'canonical','output':'output','sample_rate_hz':20.,'window_s':1200.,'min_clean_segment_s':600.,
   'response':{'output':'VEL'},'spikes':{'scale_window_s':60.},'scan':{'step_hz':.01},'cohorts':{'states':0},'phase':{'cycle_max':200}})
 write_json(root/'config.json',config);write_json(root/'mapping.json',PROFILE)
 write_json(root/'phase_settings.json',{'direction':{'mode':'filter','method':'music'}})
 write_json(root/'project.json',{'version':2,'input_folder':'','excluded_stations':[],'name':root.name,'automatic_geometry':True})
 return root

def resource_estimate(scan):
 good=[x['resolved']|{'shape':x['shape']} for x in scan['items'] if x.get('enabled') and 'resolved' in x]
 native=sum(x['shape'][x['sample_axis']]*8 for x in good)
 return {'native_float_bytes':int(native),'conversion_temporary_disk_bytes':int(native),
         'comment':'Estimate before compression; aligned caches and correlation arrays require additional space.'}

def import_hdf(project,folder=None):
 root=Path(project).resolve();spec=json.loads((root/'project.json').read_text(encoding='utf-8-sig'));profile=json.loads((root/'mapping.json').read_text(encoding='utf-8-sig'))
 folder=Path(folder or spec['input_folder'])
 if not folder.is_absolute():folder=root/folder
 folder=folder.resolve();scan=inspect_folder(folder,profile)
 write_json(root/'hdf_inventory.json',scan)
 if not scan['ready']:
  problems=[x['file']+'::'+x['dataset']+' — '+';'.join(x['issues']) for x in scan['items'] if x.get('enabled',True) and x['issues']]
  problems += [x['file']+' — '+x['error'] for x in scan['file_errors']]
  raise MappingError('请先补齐HDF5映射：\n'+'\n'.join(problems[:12]))
 selected=[x for x in scan['items'] if x.get('enabled',True) and x['resolved']['station'] not in spec.get('excluded_stations',[])]
 if len({x['resolved']['station'] for x in selected})<2:raise MappingError('至少需要两个未排除的台站')
 common_geo=all(all(k in x['resolved'] for k in ('latitude','longitude')) for x in selected)
 common_local=all(all(k in x['resolved'] for k in ('x_m','y_m')) for x in selected)
 if not common_geo and not common_local:raise MappingError('所有台站须统一使用经纬度或本地米制坐标')
 source_hash={x['file']:sha(x['file']) for x in selected};rows=[];audit=[];qc=[];groups={};response_cache={};channel_metadata={};dest=root/'canonical';dest.mkdir(exist_ok=True)
 for index,item in enumerate(selected):
  m=item['resolved'];station=m['station'];row={k:m.get(k,'') for k in ('station','seed_id','response_xml','precorrected','units','polarity','polarity_verified','clock_offset_s')}
  for k in ('latitude','longitude') if common_geo else ('x_m','y_m'):row[k]=m[k]
  if not row['precorrected']:
   embedded=m.get('response_embedded_dataset')
   if embedded:
    # Self-contained downloaded H5 remains portable even after moving just one
    # file. Extract identical embedded StationXML into the project, never write
    # to the source archive or depend on a sidecar relative to its old folder.
    with h5py.File(item['file'],'r') as h:
     if embedded not in h:raise MappingError('内嵌StationXML不存在: '+str(embedded))
     payload=h[embedded][()]
     expected_xml_sha=h.attrs.get('stationxml_sha256')
    if isinstance(payload,str):payload=payload.encode('utf-8')
    elif not isinstance(payload,bytes):raise MappingError('内嵌StationXML必须为UTF-8标量字符串')
    response_dir=root/'metadata';response_dir.mkdir(exist_ok=True)
    digest=hashlib.sha256(payload).hexdigest()
    if expected_xml_sha and digest!=expected_xml_sha:raise MappingError('内嵌StationXML校验失败: '+item['file'])
    response=response_dir/(digest+'.xml')
    if not response.exists():response.write_bytes(payload)
    # Validate embedded XML identity and complete response before accepting it.
    from obspy import read_inventory,Trace
    from dispersion.common import parse_utc
    from h5data_handoff import matching_channel
    ids=m['seed_id'].split('.');n=item['shape'][m['sample_axis']]
    tr=Trace(np.empty(0),dict(network=ids[0],station=ids[1],location=ids[2],channel=ids[3],
      starttime=parse_utc(m['start_utc']),sampling_rate=m['sampling_rate_hz'],npts=n))
    # Stats derived endtime requires npts after constructing an empty Trace.
    tr.stats.npts=n
    if digest not in response_cache:response_cache[digest]=read_inventory(str(response))
    channel_metadata[str(response.relative_to(root)),m['seed_id']]=matching_channel(response_cache[digest],tr)
    row['response_xml']=str(response.relative_to(root))
    item['stationxml_sha256']=digest
    item['original_project_response_xml']=row['response_xml']
    groups.setdefault(station,[]).append((item,m,row));continue
   response=Path(str(row['response_xml']))
   if not row['response_xml']:raise MappingError(station+'为counts数据，必须提供完整StationXML响应')
   if not response.is_absolute():response=Path(folder)/response
   if not response.is_file():raise MappingError('响应文件不存在: '+str(response))
   # Keep a portable copy; the StationXML epoch is verified later against every actual record.
   response_dir=root/'metadata';response_dir.mkdir(exist_ok=True);copy=response_dir/(sha(response)[:12]+'_'+response.name)
   if not copy.exists():copy.write_bytes(response.read_bytes())
   row['response_xml']=str(copy.relative_to(root))
  groups.setdefault(station,[]).append((item,m,row))
 # Verify station identity/metadata across all files, then resolve only identical native-grid overlap.
 for station,records in sorted(groups.items()):
  from dispersion.common import parse_utc
  records.sort(key=lambda x:parse_utc(x[1]['start_utc']).ns);template=records[0][2];kept=[];paths=[]
  # MiniSEED record headers have microsecond quantisation. Independently
  # decoded FDSN chunks can differ by 1–2 us despite identical overlapping
  # integer counts. Admit only this bounded residual for verified native
  # archives, require exact sample equality, and retain each original UTC.
  native_archive=all(item['attributes'].get('h5data_schema')=='h5data.ndstudio/1.0' for item,_,_ in records)
  time_tolerance_s=2e-6 if native_archive else 0.
  for idx,(item,m,row) in enumerate(records):
   # Canonical archives retain native integer counts. Floating conversion
   # belongs to numerical processing; promoting an entire multiday archive
   # doubles temporary memory and unnecessarily slows lossless compression.
   x=read_waveform(item,m,preserve_dtype=True)
   if row!=template:
    # FDSN XML creation/source headers may change between resumed requests.
    # Preserve every original XML, but use one canonical response only when
    # the entire selected Channel object (epoch, orientation, fs and every
    # response stage) is exactly equal; waveform metadata must also agree.
    a=channel_metadata.get((template['response_xml'],template['seed_id']));b=channel_metadata.get((row['response_xml'],row['seed_id']))
    normalized=dict(row,response_xml=template['response_xml'])
    if normalized==template and a is not None and b is not None and a==b:
     row=normalized
    else:raise MappingError(station+'跨文件的台站信息、响应或单位不同；请按明确历元拆分项目')
   start=parse_utc(m['start_utc']);fs=m['sampling_rate_hz'];duplicates=0;maximum_overlap_time_residual_ns=0
   if kept:
    pt,pfs,px=kept[-1]
    if fs!=pfs:raise MappingError(station+'的原始采样率跨文件改变；请按采样率分组')
    end_ns=pt.ns+round((len(px)-1)/fs*1e9)
    if start.ns<=end_ns:
     offset=(start.ns-pt.ns)/1e9*fs;grid_offset=round(offset)
     tolerance_samples=max(1e-5,time_tolerance_s*fs)+fs*1e-9
     if abs(offset-grid_offset)>tolerance_samples:
      raise MappingError(station+'存在不同原生网格的重叠；不能自动合并')
     overlap=len(px)-grid_offset
     overlap=min(overlap,len(x))
     covered=0
     for kt,kfs,kx in kept:
      offset=(kt.ns-start.ns)/1e9*fs;off=round(offset)
      lo=max(0,off);hi=min(overlap,off+len(kx))
      if hi<=lo:continue
      if kfs!=fs or abs(offset-off)>tolerance_samples or lo>covered or not np.array_equal(kx[lo-off:hi-off],x[lo:hi],equal_nan=True):
       raise MappingError(station+'存在冲突或不同原生网格的重叠；请人工确认来源后局部掩码，不自动拼接')
      maximum_overlap_time_residual_ns=max(maximum_overlap_time_residual_ns,round(abs(offset-off)/fs*1e9))
      covered=max(covered,hi)
     if covered!=overlap:raise MappingError(station+'重叠覆盖不完整，不能自动删除样点')
     x=x[overlap:]
     from obspy import UTCDateTime
     start=UTCDateTime(ns=start.ns+round(overlap/fs*1e9),precision=9);duplicates=overlap
   if not len(x):
    audit.append({'station':station,'file':item['file'],'dataset':item['dataset'],'identical_overlap_removed_samples':duplicates,'retained_samples':0});continue
   filename=f'{station}_{idx:05d}.npz';path=dest/filename;tmp=dest/(filename+'.part.npz')
   reused=False
   if path.exists():
    try:
     with np.load(path,allow_pickle=False) as z:
      reused=(float(z['fs_hz'])==fs and z['start_utc'].item()==str(start) and np.array_equal(z['data'],x,equal_nan=True))
    except (OSError,ValueError,KeyError):pass
   if not reused:
    np.savez_compressed(tmp,data=x,fs_hz=fs,start_utc=str(start));os.replace(tmp,path)
   # Check UTC and values once, including any NaN masks. Originals are never modified.
   with np.load(path,allow_pickle=False) as z:
    assert z['start_utc'].item()==str(start) and np.array_equal(z['data'],x,equal_nan=True)
   paths.append(filename);kept.append((start,fs,x))
   finite=x[np.isfinite(x)]
   med=float(np.median(finite)) if len(finite) else None
   mad=float(1.4826*np.median(abs(finite-med))) if len(finite) else None
   warning='no_finite_samples' if not len(finite) else 'constant_or_nearly_constant_record' if mad==0 else 'large_DC_relative_to_fluctuation_check_ADC_rail_or_offset' if m['units']=='counts' and abs(med)/mad>1e4 else ''
   qc.append({'station':station,'file':Path(item['file']).name,'dataset':item['dataset'],'median':med,'robust_scale':mad,'warning':warning,'automatic_station_exclusion':False})
   if warning:print('INPUT_QC:',station,warning,'check before interpreting; GUI can exclude station',flush=True)
   audit.append({'station':station,'file':item['file'],'dataset':item['dataset'],
     'embedded_stationxml_sha256':item.get('stationxml_sha256'),
     'original_project_response_xml':item.get('original_project_response_xml'),
     'canonical_response_xml':row['response_xml'],
     'source_sha256':source_hash[item['file']],'canonical_file':filename,'canonical_sha256':sha(path),
     'native_start_utc':str(start),'native_fs_hz':fs,'retained_samples':len(x),'identical_overlap_removed_samples':duplicates,
     'identical_overlap_time_tolerance_s':time_tolerance_s,'maximum_identical_overlap_time_residual_ns':maximum_overlap_time_residual_ns,
     'original_native_utc_shifted':False,'canonical_cache_reused':reused,
     'nonfinite_samples':int((~np.isfinite(x)).sum()),'finite_min':float(finite.min()) if len(finite) else None,
     'finite_max':float(finite.max()) if len(finite) else None,'retained_utc_and_values_preserved':True})
   print('import:',station,item['dataset'],'samples',len(x),'duplicate overlap',duplicates,flush=True)
  if not paths:raise MappingError(station+'没有保留样点')
  rows.append(dict(template,file=';'.join(paths)))
 for path,digest in source_hash.items():
  if sha(path)!=digest:raise MappingError('导入期间原始HDF5改变，请在停止写入后重试: '+path)
 csv(root/'stations.csv',rows);csv(root/'input_station_qc.csv',qc);write_json(root/'hdf_import_audit.json',{'records':audit,'sources':[{'path':k,'sha256':v} for k,v in source_hash.items()],
   'resource_estimate':resource_estimate(scan),'source_files_opened_read_only':True,'unknown_clock_offsets_not_estimated':True})
 spec['input_folder']=str(Path(folder).resolve());write_json(root/'project.json',spec)
 if spec.get('automatic_geometry',True):adapt_geometry_settings(root)
 return root/'config.json'

def adapt_geometry_settings(project):
 root=Path(project);c=json.loads((root/'config.json').read_text(encoding='utf-8-sig'));advice=suggest_geometry(root)
 arrival=float(np.ceil(advice['suggested_arrival_max_s']))
 c['scan']['arrival_max_s']=arrival
 noise=[max(c['scan']['noise_lag_s'][0],arrival+10),max(c['scan']['noise_lag_s'][1],arrival+30)]
 lag=float(np.ceil(max(c['max_lag_s'],noise[1]+5)))
 c['scan']['noise_lag_s']=noise;c['max_lag_s']=lag
 n=len(pd.read_csv(root/'stations.csv'));minimum=float(np.ceil((n*(lag+arrival+7)/2+2)/60)*60)
 c['min_clean_segment_s']=max(c['min_clean_segment_s'],minimum)
 c['window_s']=float(np.ceil(max(c['window_s'],c['min_clean_segment_s']/c['min_clean_fraction'])/60)*60)
 if c['validation']['block_s']/c['window_s']<c['validation']['min_windows_per_block']:
  raise MappingError('当前台数所需连续窗太长，验证块容纳不了最低窗数；请在高级配置中调整窗/块设置')
 write_json(root/'config.json',c);write_json(root/'geometry_advice.json',dict(advice,applied_before_window_selection=True))
 print('geometry: arrival',arrival,'lag',lag,'min clean',c['min_clean_segment_s'],flush=True)

def suggest_geometry(project):
 from dispersion.common import load_config
 from dispersion.prepare import station_table
 root=Path(project);c=load_config(root/'config.json');_,pairs,_=station_table(c);s=c['scan'];low,high=s['group_velocity_bounds_m_s'];rmax=float(pairs.distance_m.max());rmin=float(pairs.distance_m.min())
 suggested_arrival=max(s['arrival_max_s'],1.1*rmax/2000.)
 suggested_lag=max(c['max_lag_s'],suggested_arrival+35.)
 n=pd.read_csv(root/'stations.csv').shape[0]
 # Fourier wrong-time controls need delay spacing beyond lag+arrival+7 seconds.
 required_clean=np.ceil((n*(suggested_lag+suggested_arrival+7.)/2+2.)/60)*60
 return {'pair_distance_m':[rmin,rmax],'expected_arrival_s':[rmin/high,rmax/low],
   'arrival_reference_velocity_m_s':2000.,'arrival_search_covers_slowest_prior':rmax/low<=s['arrival_max_s'],
   'suggested_arrival_max_s':suggested_arrival,'suggested_max_lag_s':suggested_lag,
   'suggested_min_clean_segment_s':required_clean,'unknown_geometry_does_not_change_thresholds':True}
