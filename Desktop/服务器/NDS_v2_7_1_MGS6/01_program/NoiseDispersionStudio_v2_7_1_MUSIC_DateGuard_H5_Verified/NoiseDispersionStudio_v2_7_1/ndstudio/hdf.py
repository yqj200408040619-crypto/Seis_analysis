"""HDF5 discovery and explicit waveform mapping; originals are opened read-only."""
from pathlib import Path
import fnmatch,json,hashlib,os,re
from datetime import datetime,timezone
import numpy as np
import h5py

ALIASES={
 'sampling_rate_hz':['sampling_rate_hz','sample_rate_hz','sampling_rate','sample_rate','fs_hz','fs'],
 'start_utc':['start_utc','starttime','start_time','start_datetime','utc_start','start_timestamp'],
 'seed_id':['seed_id','nslc','trace_id'], 'station':['station','station_code','station_id','sta'],
 'network':['network','network_code','net'], 'station_label':['station_label'], 'location':['location','location_code','loc'],
 'channel':['channel','channel_code','component'], 'units':['units','unit','physical_units'],
 'latitude':['latitude','lat'], 'longitude':['longitude','lon','long'],
 'x_m':['x_m','east_m'],'y_m':['y_m','north_m'],
 'components':['components','channels','component_order'],
 'response_xml':['response_xml','stationxml'], 'response_embedded_dataset':['response_embedded_dataset'],
 'precorrected':['precorrected','response_corrected'],
 'polarity':['polarity'], 'polarity_verified':['polarity_verified'],'clock_offset_s':['clock_offset_s'],
 'valid_dataset':['valid_dataset'],'sample_axis':['sample_axis'],'component_index':['component_index'],
 'time_format':['time_format'],'vertical_confirmed':['vertical_confirmed']}

class MappingError(ValueError):pass

def integer(value,name):
 try:
  result=int(value)
  if float(value)!=result:raise ValueError()
  return result
 except (ValueError,TypeError,OverflowError):raise MappingError(name+'必须为整数，不能截断小数')

def plain(value):
 if isinstance(value,bytes):return value.decode('utf-8')
 if isinstance(value,np.ndarray):return [plain(x) for x in value.tolist()] if value.ndim else plain(value.item())
 if isinstance(value,np.generic):return plain(value.item())
 if isinstance(value,float) and not np.isfinite(value):return None
 if isinstance(value,(list,tuple)):return [plain(x) for x in value]
 if isinstance(value,complex):return str(value)
 return value

def truth(value):
 if isinstance(value,bool):return value
 if isinstance(value,(int,float)) and value in (0,1):return bool(value)
 if str(value).lower() in ('true','yes','1'):return True
 if str(value).lower() in ('false','no','0',''):return False
 raise MappingError('布尔字段必须为 true/false: '+str(value))

def utc(value,fmt='iso_utc'):
 if value is None or value=='':raise MappingError('缺少起始UTC；不能使用文件创建时间代替')
 if fmt in ('unix_s','unix_ms','unix_us','unix_ns'):
  scale={'unix_s':1,'unix_ms':1000,'unix_us':1000000,'unix_ns':1000000000}[fmt]
  from obspy import UTCDateTime
  if fmt=='unix_ns':return str(UTCDateTime(ns=int(value),precision=9))
  return str(UTCDateTime(float(value)/scale,precision=9))
 if fmt!='iso_utc':raise MappingError('未知时间格式: '+fmt)
 if not isinstance(value,str) or re.fullmatch(r'[-+\d.eE]+',value):raise MappingError('数值时间戳必须显式选择 unix_s/ms/us/ns 格式')
 from obspy import UTCDateTime
 # ISO timestamps without a suffix are interpreted as UTC only under this explicit mapping.
 try:
  match=re.fullmatch(r'(.*[T ]\d{2}:\d{2}:\d{2})[.,](\d+)(.*)',value)
  if match:
   if len(match[2])>9:raise ValueError('UTC小数最多支持9位')
   base=UTCDateTime(match[1]+match[3]);return str(UTCDateTime(ns=base.ns+int(match[2].ljust(9,'0')),precision=9))
  return str(UTCDateTime(value,precision=9))
 except Exception as exc:raise MappingError('UTC时间无法解析: '+str(value)) from exc

def attrs_for(ds):
 nodes=[];node=ds
 while True:
  nodes.append(node)
  if node.name=='/':break
  node=node.parent
 merged={}
 for node in reversed(nodes):
  for key,value in node.attrs.items():merged[str(key)]=plain(value)
 # Scalar siblings are useful for simple files storing fs/start/coordinates as datasets.
 for key,node in ds.parent.items():
  if isinstance(node,h5py.Dataset) and node.shape==():merged.setdefault(key,plain(node[()]))
 return merged

def extract_metadata(ds):
 attrs=attrs_for(ds);lower={k.lower():v for k,v in attrs.items()};meta={}
 for field,aliases in ALIASES.items():
  for name in aliases:
   if name in lower:meta[field]=lower[name];break
 if 'sampling_rate_hz' not in meta:
  for name in ('dt','sample_interval_s','delta'):
   if name in lower:
    if float(lower[name])<=0:raise MappingError('采样间隔必须大于0')
    meta['sampling_rate_hz']=1/float(lower[name]);break
 return meta,attrs

def describe(path,ds):
 meta,attrs=extract_metadata(ds)
 return {'file':str(Path(path).resolve()),'dataset':ds.name,'shape':list(ds.shape),'dtype':str(ds.dtype),
   'metadata':meta,'attributes':attrs,'nbytes':int(ds.size*ds.dtype.itemsize),'issues':[]}

def resolve(item,override=None):
 meta={k:v for k,v in item['metadata'].items() if v is not None};override=dict(override or {})
 # Attribute mappings are evaluated against the raw merged attributes, before explicit overrides.
 for field,key in override.get('attribute_map',{}).items():
  if key not in item['attributes']:raise MappingError('HDF5属性不存在: '+key)
  meta[field]=item['attributes'][key]
 for key,value in override.items():
  if key=='attribute_map':continue
  if value is None or value=='':meta.pop(key,None)
  else:meta[key]=value
 if 'components' in meta and isinstance(meta['components'],str):meta['components']=re.split(r'[,;\s]+',meta['components'].strip())
 shape=item['shape'];nd=len(shape)
 if nd==1:axis=0;component=None
 elif nd==2:
  axis=meta.get('sample_axis')
  if axis is None:
   short=[i for i,n in enumerate(shape) if n<=3];long=[i for i,n in enumerate(shape) if n>3]
   if len(short)!=1 or len(long)!=1:raise MappingError('二维数组需要指定sample_axis和component_index')
   axis=long[0]
  axis=integer(axis,'sample_axis')
  if axis not in (0,1):raise MappingError('sample_axis必须为0或1')
  component=meta.get('component_index')
  if component is None and 'components' in meta:
   comps=meta['components'];vertical=[i for i,c in enumerate(comps) if str(c).upper().endswith('Z')]
   if len(vertical)==1:component=vertical[0]
  if component is None and shape[1-axis]==1:component=0
  if component is None:raise MappingError('三分量顺序未声明；请选择Z分量索引')
  component=integer(component,'component_index')
  if not 0<=component<shape[1-axis]:raise MappingError('分量索引超出数组范围')
  if 'components' in meta:
   if len(meta['components'])!=shape[1-axis]:raise MappingError('components数量与数组分量维度不符')
   meta['channel']=meta['components'][component]
 else:raise MappingError('仅支持一维波形或二维分量数组；更高维需先拆分为明确台站记录')
 fs=float(meta.get('sampling_rate_hz',0))
 if not np.isfinite(fs) or fs<=0:raise MappingError('缺少有效原始采样率')
 start=utc(meta.get('start_utc'),meta.get('time_format','iso_utc'))
 sid=meta.get('seed_id')
 if sid:
  parts=str(sid).split('.')
  if len(parts)!=4 or not parts[1] or not parts[3]:raise MappingError('seed_id必须为NET.STA.LOC.CHA（空LOC保留两个点）')
  if 'channel' in meta and parts[3]!=str(meta['channel']):
   if nd==2:parts[3]=str(meta['channel']);sid='.'.join(parts)
   else:raise MappingError('seed_id通道与channel冲突')
  meta['station']=parts[1];meta['channel']=parts[3]
 else:
  if not meta.get('station') or not meta.get('channel'):raise MappingError('缺少台站或垂直分量通道；请填写映射')
  sid='.'.join([str(meta.get('network','XX')),str(meta['station']),str(meta.get('location','')),str(meta['channel'])])
 if not str(meta['channel']).upper().endswith('Z') and not truth(meta.get('vertical_confirmed',False)):
  raise MappingError('所选通道不以Z结尾；确认它是垂直分量后填写vertical_confirmed=true')
 station=str(meta.get('station_label',meta['station']))
 if not re.fullmatch(r'[A-Za-z0-9_-]+',station):raise MappingError('台站代码只能含字母、数字、下划线和连字符')
 units=str(meta.get('units','')).lower().replace(' ','')
 unitmap={'m/s':'m/s','m/s^2':'m/s^2','m/s2':'m/s^2','m':'m','counts':'counts','count':'counts','adc':'counts'}
 if units not in unitmap:raise MappingError('必须明确units=counts、m/s、m/s^2或m；其他单位请先换算')
 units=unitmap[units];pre=truth(meta.get('precorrected',units!='counts'))
 if pre and units=='counts':raise MappingError('counts不能声明为已响应校正的物理量')
 if not pre and units!='counts':raise MappingError('物理单位数据须明确precorrected=true')
 for field in ('latitude','longitude','x_m','y_m'):
  if field in meta:meta[field]=float(meta[field])
 geo=all(k in meta for k in ('latitude','longitude'));local=all(k in meta for k in ('x_m','y_m'))
 if not geo and not local:raise MappingError('缺少台站经纬度或本地米制坐标')
 if geo and (not -90<=meta['latitude']<=90 or not -180<=meta['longitude']<=180):raise MappingError('经纬度超出范围')
 if not all(np.isfinite(meta[k]) for k in ('latitude','longitude','x_m','y_m') if k in meta):raise MappingError('坐标不是有限数值')
 meta.update(seed_id=str(sid),station=station,start_utc=start,sampling_rate_hz=fs,sample_axis=axis,component_index=component,
   units=units,precorrected=pre,polarity=integer(meta.get('polarity',1),'polarity'),polarity_verified=truth(meta.get('polarity_verified',False)),
   clock_offset_s=float(meta.get('clock_offset_s',0)))
 if meta['polarity'] not in (-1,1):raise MappingError('polarity必须为+1或-1')
 if not np.isfinite(meta['clock_offset_s']):raise MappingError('clock_offset_s必须为有限数值')
 return meta

def inspect_folder(folder,profile=None):
 folder=Path(folder).resolve();profile=profile or {};items=[];errors=[]
 files=sorted(p for p in folder.rglob('*') if p.is_file() and p.suffix.lower() in ('.h5','.hdf5','.hdf'))
 if not files:raise MappingError('文件夹中没有.h5、.hdf5或.hdf文件')
 for path in files:
  try:
   with h5py.File(path,'r') as h:
    # Validate the entire self-contained archive, including files whose only
    # waveform is a short boundary fragment (<8 samples) excluded from analysis.
    expected_xml_sha=h.attrs.get('stationxml_sha256')
    if h.attrs.get('h5data_schema')=='h5data.ndstudio/1.0' and expected_xml_sha:
     if '/metadata/stationxml' not in h:raise MappingError('内嵌StationXML不存在: '+str(path))
     payload=h['/metadata/stationxml'][()]
     if isinstance(payload,str):payload=payload.encode('utf-8')
     if not isinstance(payload,bytes) or hashlib.sha256(payload).hexdigest()!=plain(expected_xml_sha):
      raise MappingError('内嵌StationXML校验失败: '+str(path))
    def visit(name,ds):
     if not isinstance(ds,h5py.Dataset) or ds.shape is None or ds.ndim not in (1,2) or ds.size<8:return
     if ds.dtype.kind not in ('i','u','f'):return
     if h.attrs.get('h5data_schema')=='h5data.ndstudio/1.0' and ds.attrs.get('role')!='waveform':return
     pattern=profile.get('dataset_pattern') or '*'
     if not fnmatch.fnmatch(ds.name,pattern):return
     if not profile.get('dataset_pattern') and ds.name.rsplit('/',1)[-1].lower() in ('time','times','timestamp','timestamps','valid','mask','latitude','longitude','coordinates'):return
     item=describe(path,ds);override=override_for(item,profile)
     # Download archives retain horizontal channels, but the Rayleigh workflow
     # selects Z by default. An explicit enabled=true mapping can opt one in.
     if h.attrs.get('h5data_schema')=='h5data.ndstudio/1.0' and not str(item['metadata'].get('channel','')).endswith('Z') and 'enabled' not in override:
      item['enabled']=False;item['exclusion_reason']='下载档案保留的非Z通道；默认不用于垂直分量频散';items.append(item);return
     if override.get('enabled',True) is False:item['enabled']=False;items.append(item);return
     item['enabled']=True
     try:item['resolved']=resolve(item,override)
     except (ValueError,TypeError,IndexError) as exc:item['issues']=[str(exc)]
     items.append(item)
    h.visititems(visit)
  except (OSError,ValueError) as exc:errors.append({'file':str(path),'error':str(exc)})
 return {'folder':str(folder),'files':len(files),'items':items,'file_errors':errors,
   'ready':not errors and bool(items) and all(not x['issues'] for x in items if x.get('enabled',True)) and any(x.get('enabled',True) for x in items)}

def override_for(item,profile):
 result=dict(profile.get('defaults',{}))
 rel=Path(item['file']).name
 for rule in profile.get('rules',[]):
  if fnmatch.fnmatch(rel,rule.get('file_pattern','*')) and fnmatch.fnmatch(item['dataset'],rule.get('dataset_pattern','*')):
   result.update(rule.get('values',{}))
 key=item['file']+'::'+item['dataset'];result.update(profile.get('overrides',{}).get(key,{}))
 return result

def read_waveform(item,meta,preserve_dtype=False):
 with h5py.File(item['file'],'r') as h:
  ds=h[item['dataset']]
  sel=slice(None) if ds.ndim==1 else ((slice(None),meta['component_index']) if meta['sample_axis']==0 else (meta['component_index'],slice(None)))
  values=np.asarray(ds[sel],dtype=None if preserve_dtype else float)
  valid_name=meta.get('valid_dataset')
  if valid_name:
   valid=h[valid_name];vv=np.asarray(valid[sel] if valid.ndim==ds.ndim else valid[:],bool)
   if vv.shape!=values.shape:raise MappingError('有效样点掩码与选定波形长度不一致')
   values=np.where(vv,values,np.nan)
 return values
