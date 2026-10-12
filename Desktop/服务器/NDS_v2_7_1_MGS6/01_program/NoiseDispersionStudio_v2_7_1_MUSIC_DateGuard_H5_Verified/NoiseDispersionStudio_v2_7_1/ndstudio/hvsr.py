"""Single-instrument HVSR on native UTC grids; no interpolation or whitening.

Inputs are downloaded H5 archives or explicitly calibrated one-dimensional
waveform datasets. Raw source files are always opened read-only. Each window
retains a decision and reason; no peak-driven window rejection is performed.
"""
from pathlib import Path
from io import BytesIO
import hashlib
import json
import math
import os
import numpy as np
import pandas as pd
import h5py
from scipy import signal
from obspy import Trace, UTCDateTime, read_inventory
from .hdf import attrs_for, extract_metadata, plain, truth, utc

DEFAULTS = dict(fmin_hz=0.2, fmax_hz=20., window_s=60., overlap=0.,
                ko_bandwidth=40., frequency_points=240, horizontal="rms",
                start="", end="", max_crest=12., max_block_rms_ratio=4.,
                reject_transients=True, water_level=None, response_padding_cycles=3.,
                vs_min_m_s=150., vs_max_m_s=500., bedrock_vs_m_s=1000.,
                peak_frequency_hz=None, rayleigh_assumption=False,
                ellipticity_trials=160, ellipticity_seed=42,
                thickness_min_m=1., thickness_max_m=300.,
                vp_vs=1.9, bedrock_vp_vs=1.8, density_g_cm3=1.8,
                bedrock_density_g_cm3=2.2)


def write_json(path, data):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    os.replace(tmp, path)


def settings(value=None):
    value = value or {}
    unknown = set(value)-set(DEFAULTS)
    if unknown: raise ValueError("未知 H/V 参数: " + ", ".join(sorted(unknown)))
    c = DEFAULTS | value
    for k in ("fmin_hz", "fmax_hz", "window_s", "ko_bandwidth", "max_crest",
              "max_block_rms_ratio", "response_padding_cycles", "vs_min_m_s", "vs_max_m_s",
              "bedrock_vs_m_s", "thickness_min_m", "thickness_max_m", "vp_vs", "bedrock_vp_vs",
              "density_g_cm3", "bedrock_density_g_cm3"):
        if not np.isfinite(c[k]) or c[k] <= 0: raise ValueError(k+" 必须为有限正数")
    if not 0 <= c['overlap'] < 1: raise ValueError("overlap 必须在 [0,1)")
    if not c['fmin_hz'] < c['fmax_hz']: raise ValueError("频率上限必须大于下限")
    if not c['vs_min_m_s'] < c['vs_max_m_s']: raise ValueError("Vs 上限必须大于下限")
    if not c['thickness_min_m'] < c['thickness_max_m']: raise ValueError("厚度上限必须大于下限")
    if min(c['vp_vs'], c['bedrock_vp_vs']) <= np.sqrt(4/3): raise ValueError("Vp/Vs 不满足正体积模量条件")
    for k, lo, hi in (("frequency_points",40,2000),("ellipticity_trials",20,10000)):
        if int(c[k]) != c[k] or not lo <= c[k] <= hi: raise ValueError(k+" 超出可用范围")
    if c['horizontal'] not in ('rms','geometric','total'): raise ValueError("horizontal 必须为 rms/geometric/total")
    if c['peak_frequency_hz'] is not None and not c['fmin_hz'] < c['peak_frequency_hz'] < c['fmax_hz']:
        raise ValueError("指定峰频率必须在分析频段内")
    if c['start'] and c['end'] and UTCDateTime(c['start']) >= UTCDateTime(c['end']): raise ValueError("结束UTC必须晚于开始UTC")
    return c


def discover(folder):
    """Metadata-only inventory, grouped by NET.STA.LOC.channel-family.

    Never mix instrument locations/families. 1/2 components require orientation
    metadata. N/E/Z suffixes are an explicit conventional fallback and recorded.
    """
    folder = Path(folder).resolve()
    paths = [folder] if folder.is_file() else sorted(p for p in folder.rglob('*') if p.suffix.lower() in ('.h5','.hdf5','.hdf'))
    groups = {}; errors = []
    for path in paths:
        try:
            with h5py.File(path, 'r') as h:
                def visit(name, ds):
                    if not isinstance(ds,h5py.Dataset) or ds.ndim != 1 or not np.issubdtype(ds.dtype,np.number): return
                    a = attrs_for(ds)
                    if a.get('role') in ('metadata','binary_archive','mask'): return
                    m, a = extract_metadata(ds)
                    if not m.get('seed_id') and not (m.get('station') and m.get('channel')): return
                    try:
                        sid = str(m.get('seed_id') or '.'.join([str(m.get('network','XX')),str(m['station']),str(m.get('location','')),str(m['channel'])]))
                        ids=sid.split('.')
                        if len(ids)!=4 or len(ids[3])<2: raise ValueError("NSLC格式错误")
                        fs=float(m['sampling_rate_hz'])
                        if not np.isfinite(fs) or fs<=0 or len(ds)<2: raise ValueError("采样率/样点数无效")
                        t=UTCDateTime(utc(m['start_utc'],m.get('time_format','iso_utc')))
                        orientation='metadata'; az=a.get('azimuth_deg',a.get('azimuth')); dip=a.get('dip_deg',a.get('dip'))
                        if az is None or dip is None:
                            suffix=ids[3][-1].upper()
                            if suffix not in 'NEZ': raise ValueError("1/2或自定义分量必须有azimuth_deg和dip_deg")
                            az=90. if suffix=='E' else 0.; dip=-90. if suffix=='Z' else 0.
                            orientation='inferred_from_channel_suffix'
                        az=float(az);dip=float(dip)
                        if not np.isfinite([az,dip]).all() or not -90<=dip<=90: raise ValueError("分量方向无效")
                        units=str(m.get('units','')).lower()
                        if units not in ('counts','m/s','m/s^2','m'): raise ValueError("明确声明单位 counts/m/s/m/s^2/m")
                        calibrated=truth(m.get('precorrected',False))
                        if units=='counts' and calibrated: raise ValueError("counts不能声明为已校正")
                        if units!='counts' and not calibrated: raise ValueError("物理单位须声明precorrected=true")
                        embedded=m.get('response_embedded_dataset'); response=str(m.get('response_xml') or '')
                        if units=='counts' and not embedded and not response: raise ValueError("counts缺少StationXML响应")
                        if embedded and embedded not in h: raise ValueError("内嵌StationXML路径不存在")
                        if response and not Path(response).is_absolute(): response=str((path.parent/response).resolve())
                        row=dict(file=str(path),dataset=name,seed_id=sid,channel=ids[3],fs_hz=fs,
                                 start_ns=t.ns,end_ns=t.ns+round((len(ds)-1)/fs*1e9),npts=len(ds),units=units,
                                 azimuth_deg=az,dip_deg=dip,orientation_source=orientation,
                                 response_embedded_dataset=embedded,response_xml=response,
                                 valid_dataset=m.get('valid_dataset'),clock_offset_s=float(m.get('clock_offset_s') or 0.),
                                 stationxml_sha256=a.get('stationxml_sha256'))
                        key='.'.join(ids[:3]+[ids[3][:-1]])
                        groups.setdefault(key,{}).setdefault(ids[3],[]).append(row)
                    except Exception as exc: errors.append(dict(file=str(path),dataset=name,error=str(exc)))
                h.visititems(visit)
        except Exception as exc: errors.append(dict(file=str(path),error=str(exc)))
    output=[]
    for key,channels in sorted(groups.items()):
        issues=[]
        if len(channels)!=3: issues.append(f"需要同一仪器3个分量，当前{len(channels)}个")
        rates={r['fs_hz'] for rows in channels.values() for r in rows}
        if len(rates)!=1: issues.append("采样率不同；请明确预先重采样，本入口不自动插值")
        output.append(dict(instrument_id=key,channels=channels,ready=not issues,issues=issues,
                           start_utc=str(UTCDateTime(ns=min(r['start_ns'] for rows in channels.values() for r in rows))),
                           end_utc=str(UTCDateTime(ns=max(r['end_ns'] for rows in channels.values() for r in rows)))))
    return dict(folder=str(folder),instruments=output,errors=errors)


def _read_interval(records, t_ns, n, fs, inventories, core_slice=None):
    """Read only one window plus response padding, preserve gaps and overlaps."""
    out=np.full(n,np.nan); assigned=np.zeros(n,bool); template=None; calibration=None
    for row in records:
        effective=row['start_ns']-round(row['clock_offset_s']*1e9)
        offset=(effective-t_ns)/1e9*fs; off=round(offset)
        lo=max(0,off);hi=min(n,off+row['npts'])
        if lo>=hi: continue
        tolerance=max(1e-5,2e-6*fs)+fs*1e-9
        if abs(offset-off)>tolerance: raise ValueError("三分量时间网格不齐；不能静默对齐或插值")
        with h5py.File(row['file'],'r') as h:
            x=np.asarray(h[row['dataset']][lo-off:hi-off],dtype=float)
            if row['valid_dataset']:
                x[~np.asarray(h[row['valid_dataset']][lo-off:hi-off],dtype=bool)]=np.nan
            inv=None
            if row['units']=='counts':
                cachekey=(row['file'],row['response_embedded_dataset'] or row['response_xml'])
                if cachekey not in inventories:
                    if row['response_embedded_dataset']:
                        xml=h[row['response_embedded_dataset']][()]
                        if isinstance(xml,str): xml=xml.encode('utf-8')
                        expected=h.attrs.get('stationxml_sha256')
                        if expected and hashlib.sha256(xml).hexdigest()!=expected: raise ValueError("StationXML校验失败")
                        inventories[cachekey]=read_inventory(BytesIO(xml))
                    else: inventories[cachekey]=read_inventory(row['response_xml'])
                inv=inventories[cachekey]
        old=assigned[lo:hi]
        if old.any() and not np.array_equal(out[lo:hi][old],x[old],equal_nan=True): raise ValueError("重叠波形冲突，不能自动拼接")
        current=None
        if inv is not None:
            from h5data_handoff import matching_channel
            ids=row['seed_id'].split('.')
            tr=Trace(np.zeros(0),dict(network=ids[0],station=ids[1],location=ids[2],channel=ids[3],sampling_rate=fs,starttime=UTCDateTime(ns=row['start_ns'])))
            tr.stats.npts=row['npts']; current=matching_channel(inv,tr)
        identity=(row['units'],row['azimuth_deg'],row['dip_deg'])
        if template and (identity!=calibration or current!=template.get('_channel')): raise ValueError("窗内响应或分量方向历元改变，须拆分分析")
        if template is None: template=dict(row,_channel=current); calibration=identity
        out[lo:hi]=x;assigned[lo:hi]=True
    valid=assigned & np.isfinite(out)
    if not valid.all():
        core=slice(*core_slice) if core_slice is not None else slice(None)
        if not valid[core].all(): raise ValueError("分析窗内缺测或非有限样点")
        raise ValueError("响应校正缓冲区不完整（分析窗内数据完整）")
    return out,template


def _correct(x,row,fs,c):
    x=signal.detrend(x,type='linear')
    if row['units']=='counts':
        ids=row['seed_id'].split('.')
        tr=Trace(x.copy(),dict(network=ids[0],station=ids[1],location=ids[2],channel=ids[3],sampling_rate=fs))
        tr.stats.response=row['_channel'].response
        # Frequency band lies strictly inside this prefilter. No sensitivity-only shortcut.
        tr.remove_response(output='VEL',water_level=c['water_level'],
                           pre_filt=(c['fmin_hz']*.5,c['fmin_hz']*.8,c['fmax_hz']*1.05,min(c['fmax_hz']*1.2,fs*.499)),
                           zero_mean=False,taper=True,taper_fraction=.05)
        return tr.data
    # All calibrated components must have the same physical unit. No unit conversion needed for a ratio.
    return x


def orientation_matrix(rows):
    az=np.radians([r['azimuth_deg'] for r in rows]);dip=np.radians([r['dip_deg'] for r in rows])
    m=np.column_stack([np.cos(dip)*np.cos(az),np.cos(dip)*np.sin(az),-np.sin(dip)])
    if np.linalg.cond(m)>10: raise ValueError("三分量方向矩阵退化或过度倾斜，不能作为单台三分量")
    return m


def ko_matrix(fft_f,frequencies,b):
    positive=fft_f>0
    z=b*np.log10(fft_f[positive][None,:]/frequencies[:,None])
    w=np.sinc(z/np.pi)**4
    w/=w.sum(axis=1,keepdims=True)
    return positive,w


def spectral_ratio(x,fs,c,matrix=None):
    """Smooth component Fourier amplitudes before forming the ratio."""
    if matrix is not None: x=np.linalg.solve(matrix,x)
    n=x.shape[1]; taper=signal.windows.tukey(n,.1)
    fa=np.abs(np.fft.rfft(signal.detrend(x,axis=1)*taper,axis=1))*np.sqrt(2/(fs*np.sum(taper*taper)))
    fft_f=np.fft.rfftfreq(n,1/fs);f=np.geomspace(c['fmin_hz'],c['fmax_hz'],int(c['frequency_points']))
    positive,w=matrix_cache(n,fs,c)
    smooth=fa[:,positive]@w.T; hn,he,v=smooth
    if np.any(v<=max(float(v.max())*1e-10,np.finfo(float).tiny)): raise ValueError("垂直分量谱接近零，H/V无定义")
    # Combine unsmoothed horizontal amplitudes, then smooth H. RMS/total retain
    # rotation invariance; smoothing each |FFT| before RMS would break it.
    hraw=np.sqrt(fa[0]*fa[1]) if c['horizontal']=='geometric' else np.sqrt(fa[0]**2+fa[1]**2)/(np.sqrt(2) if c['horizontal']=='rms' else 1.)
    h=hraw[positive]@w.T
    return f,h/v,smooth,hn/v,he/v


_KO_CACHE={}
def matrix_cache(n,fs,c):
    key=(n,fs,c['fmin_hz'],c['fmax_hz'],c['frequency_points'],c['ko_bandwidth'])
    if key not in _KO_CACHE:
        _KO_CACHE.clear()  # One station/job: bounded memory, no dense FFT-square smoothing matrix.
        _KO_CACHE[key]=ko_matrix(np.fft.rfftfreq(n,1/fs),np.geomspace(c['fmin_hz'],c['fmax_hz'],int(c['frequency_points'])),c['ko_bandwidth'])
    return _KO_CACHE[key]


def quality(x,fs,c):
    centered=x-np.median(x,axis=1,keepdims=True)
    rms=np.sqrt(np.mean(centered**2,axis=1))
    if not np.isfinite(rms).all() or np.any(rms<=np.finfo(float).tiny): return 'constant_or_nonfinite_component',{}
    crest=float(np.max(np.max(np.abs(centered),axis=1)/rms))
    k=max(2,int(round(fs)))
    blocks=[np.sqrt(np.mean(b*b,axis=1)) for b in np.array_split(centered,max(1,centered.shape[1]//k),axis=1)]
    block=np.asarray(blocks); ratio=float(np.max(block/np.maximum(np.median(block,axis=0),np.finfo(float).tiny)))
    # Repeated extrema on both rails indicate clipping; one transient extreme is not sufficient.
    rail=float(max(min(np.mean(row==row.min()),np.mean(row==row.max())) for row in x))
    metrics=dict(crest=crest,block_rms_ratio=ratio,rail_fraction=rail)
    if rail>.005: return 'clipping_or_quantized_rails',metrics
    if c['reject_transients'] and (crest>c['max_crest'] or ratio>c['max_block_rms_ratio']): return 'transient',metrics
    return '',metrics


def _independent_indices(starts,window_s):
    chosen=[];last=-math.inf
    for i,s in enumerate(starts):
        if s>=last+window_s-1e-6: chosen.append(i);last=s
    return np.asarray(chosen,dtype=int)


def peak_table(f,mean,std,curves,window_s,independent_count,selected=None,minimum_cycles=200):
    # Ignore numerical ripples or near-flat maxima: at least 10% log prominence.
    indices=signal.find_peaks(np.log(mean),prominence=np.log(1.1))[0]
    indices=[int(i) for i in indices if mean[i]>1.2]
    if selected is not None:
        near=[i for i in indices if abs(np.log(f[i]/selected))<.1]
        if not near: raise ValueError("指定频率附近±10%没有局部峰；请先查看曲线")
        selected=min(near,key=lambda i:abs(np.log(f[i]/selected)))
    rows=[]
    for i in indices:
        f0=float(f[i]);a0=float(mean[i]);local=(f>=f0*.8)&(f<=f0*1.2)
        # Frequency scatter follows the same local branch; peak support separately checks that a local maximum exists.
        frequencies=f[local]; individual=frequencies[np.argmax(curves[:,local],axis=1)]
        support=float(np.mean([any(.8*f0<=f[j]<=1.2*f0 for j in signal.find_peaks(y)[0]) for y in curves]))
        sigmaf=float(np.std(individual,ddof=1)) if len(individual)>1 else float('inf')
        eps,theta=(.25,3.) if f0<.2 else (.20,2.5) if f0<.5 else (.15,2.) if f0<1 else (.10,1.78) if f0<2 else (.05,1.58)
        left=(f>=f0/4)&(f<f0);right=(f>f0)&(f<=4*f0); reliability_band=(f>=.5*f0)&(f<=2*f0)
        coverage=bool(f[0]<=.5*f0 and f[-1]>=2*f0)
        peak_pm=[]
        for a in (mean/std,mean*std):
            ii=signal.find_peaks(a)[0];ii=[j for j in ii if .8*f0<=f[j]<=1.2*f0]
            peak_pm.append(bool(ii) and abs(f[max(ii,key=lambda j:a[j])]/f0-1)<=.05)
        checks=dict(left_drop=bool(left.any() and (a0/mean[left]).max()>2),right_drop=bool(right.any() and (a0/mean[right]).max()>2),
                    amplitude_gt_2=a0>2,plus_minus_peak_stable=all(peak_pm),frequency_scatter=bool(sigmaf<eps*f0),amplitude_scatter=bool(std[i]<theta))
        rchecks=dict(ten_cycles=bool(f0*window_s>10),significant_cycle_count_pass=bool(f0*window_s*independent_count>minimum_cycles),
                     band_scatter=bool(coverage and np.all(std[reliability_band]<(3. if f0<.5 else 2.))))
        half=signal.peak_widths(mean,[i],rel_height=.5)
        width=float(np.interp(half[3][0],np.arange(len(f)),f)-np.interp(half[2][0],np.arange(len(f)),f))
        count=sum(checks.values()); reliable=all(rchecks.values())
        rows.append(dict(f0_hz=f0,amplitude=a0,period_s=1/f0,sigma_f_hz=sigmaf if np.isfinite(sigmaf) else None,
                         amplitude_std_factor=float(std[i]),cycles_per_window=f0*window_s,independent_cycles=f0*window_s*independent_count,
                         reliability_pass=reliable,clarity_pass_count=count,clear_peak=bool(reliable and count>=5 and support>=.8),
                         individual_peak_support_fraction=support,minimum_significant_cycles=minimum_cycles,
                         full_scatter_band_covered=coverage,relative_width=width/f0,narrow_peak_check=bool(width/f0<.08),
                         selected=bool(i==selected),**rchecks,**checks))
    rows.sort(key=lambda r:(not r['clear_peak'],not r['reliability_pass'],-r['clarity_pass_count'],-r['amplitude']))
    if selected is None and rows: rows[0]['selected']=True
    return rows


def sha_file(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def run_hvsr(folder,output,instrument_id=None,config=None):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    # A failed/restarted job cannot leave a previous completion marker looking current.
    (output/'hvsr_complete.json').unlink(missing_ok=True)
    c=settings(config)
    # Remove only known generated structure products from a previous successful
    # run; failed or flat-curve runs must not expose them as current estimates.
    for name in ('hvsr_summary.json','shallow_structure.json','depth_vs_tradeoff.csv','ellipticity_trials.csv',
                 'ellipticity_fit.csv','ellipticity_models.json','hvsr_report.html'):
        (output/name).unlink(missing_ok=True)
    for name in ('shallow_structure.png','rayleigh_ellipticity.png'):
        (output/'figures'/name).unlink(missing_ok=True)
    inv=discover(folder);write_json(output/'inventory.json',inv)
    instruments=[i for i in inv['instruments'] if i['ready'] and (not instrument_id or i['instrument_id']==instrument_id)]
    if len(instruments)!=1: raise ValueError("请选择一个完整三分量仪器（NET.STA.LOC.通道前缀）；可用仪器数="+str(len(instruments)))
    item=instruments[0]; channels=[sorted(rows,key=lambda r:r['start_ns']) for _,rows in sorted(item['channels'].items())]
    fs=channels[0][0]['fs_hz'];n=round(c['window_s']*fs)
    if n<32 or abs(n/fs-c['window_s'])>1e-6: raise ValueError("窗长必须是整数个原生样点，且至少32点")
    if c['fmin_hz']<2/(n/fs): raise ValueError("频率下限不足2个FFT频率间隔；增大窗长")
    if c['fmax_hz']>=fs*.45: raise ValueError(f"频率上限须小于0.45×采样率={fs*.45:g} Hz，留响应预滤波余量")
    step=max(1,round(n*(1-c['overlap'])))
    allrows=[r for rows in channels for r in rows]
    unitset={r['units'] for r in allrows}
    if len(unitset)!=1: raise ValueError("三分量的单位/校正状态不同，不能混合counts与物理单位")
    start=max(min(r['start_ns']-round(r['clock_offset_s']*1e9) for r in rows) for rows in channels)
    end=min(max(r['end_ns']-round(r['clock_offset_s']*1e9) for r in rows) for rows in channels)
    origin=start
    if c['start']: start=max(start,UTCDateTime(c['start']).ns)
    if c['end']: end=min(end,UTCDateTime(c['end']).ns-1) # end is exclusive
    start=origin+round(math.ceil((start-origin)*fs/1e9-1e-5)/fs*1e9)
    count=math.floor((end-start)/1e9*fs)+1
    if count<n: raise ValueError("三分量共同记录不足一个完整时窗")
    padding=math.ceil(c['response_padding_cycles']/c['fmin_hz']*fs) if unitset=={'counts'} else 0
    inventories={};audit=[];ratios=[];spectra=[];nv=[];ev=[];starts=[]
    for offset in range(0,count-n+1,step):
        t=start+round(offset/fs*1e9);record=dict(window_id=len(audit),start_utc=str(UTCDateTime(ns=t)),accepted=False,reason='')
        try:
            raw=[];meta=[]
            for rows in channels:
                x,m=_read_interval(rows,t-round(padding/fs*1e9),n+2*padding,fs,inventories,
                                   core_slice=(padding,padding+n));raw.append(x);meta.append(m)
            raw=np.asarray(raw);reason,metrics=quality(raw[:,padding:padding+n],fs,c);record.update(metrics)
            if reason: raise ValueError(reason)
            matrix=orientation_matrix(meta)
            corrected=np.asarray([_correct(x,m,fs,c) for x,m in zip(raw,meta)])[:,padding:padding+n]
            f,ratio,amp,nratio,eratio=spectral_ratio(corrected,fs,c,matrix)
            if not np.isfinite(ratio).all() or np.any(ratio<=0): raise ValueError("H/V非有限或非正数")
            ratios.append(ratio);spectra.append(amp);nv.append(nratio);ev.append(eratio);starts.append(t/1e9)
            record['accepted']=True
        except ValueError as exc: record['reason']=str(exc)
        audit.append(record)
        if len(audit)%25==0: print(f"H/V: {len(audit)} windows, {len(ratios)} accepted",flush=True)
    pd.DataFrame(audit).to_csv(output/'window_qc.csv',index=False,encoding='utf-8-sig')
    if len(ratios)<3: raise ValueError(f"H/V至少需要3个有效时窗；当前{len(ratios)}个，查看window_qc.csv")
    curves=np.asarray(ratios);ind=_independent_indices(starts,n/fs)
    # Overlapping windows are displayed but statistics/SESAME counts use non-overlapping windows only.
    statistical=curves[ind];log=np.log(statistical);mean=np.exp(log.mean(axis=0))
    std=np.exp(log.std(axis=0,ddof=1)) if len(ind)>1 else np.full(len(f),np.inf)
    if not np.isfinite(std).all(): raise ValueError("至少需要2个互不重叠的有效时窗来计算离散度")
    peaks=peak_table(f,mean,std,statistical,n/fs,len(ind),c['peak_frequency_hz'],200 if c['reject_transients'] else 1000)
    peakframe=pd.DataFrame(peaks) if peaks else pd.DataFrame(columns=['f0_hz','amplitude','reliability_pass','clarity_pass_count','clear_peak','selected'])
    peakframe.to_csv(output/'peaks.csv',index=False,encoding='utf-8-sig')
    asp=np.exp(np.log(np.asarray(spectra)[ind]).mean(axis=0))
    table=pd.DataFrame(dict(frequency_hz=f,hv=mean,hv_std_factor=std,hv_minus_sigma=mean/std,hv_plus_sigma=mean*std,
                            n_over_v=np.exp(np.log(np.asarray(nv)[ind]).mean(axis=0)),e_over_v=np.exp(np.log(np.asarray(ev)[ind]).mean(axis=0)),
                            n_amplitude=asp[0],e_amplitude=asp[1],z_amplitude=asp[2]))
    table.to_csv(output/'hvsr_curve.csv',index=False,encoding='utf-8-sig')
    with h5py.File(output/'hvsr_windows.h5','w') as h:
        h.attrs.update(schema='ndstudio.hvsr/1',instrument_id=item['instrument_id'],horizontal=c['horizontal'])
        for k,v in dict(frequency_hz=f,hv=curves,component_amplitude=np.asarray(spectra),start_unix_s=starts,statistical_window_indices=ind).items(): h.create_dataset(k,data=v,compression='gzip')
    warnings=[]
    if any(r['orientation_source']!='metadata' for r in allrows): warnings.append("部分方向由N/E/Z命名推定，请核实仪器方位与倾斜")
    if c['window_s']*c['fmin_hz']<=10: warnings.append("频段低端每窗不超过10周期，增大窗长后复核")
    if len([p for p in peaks if p['clear_peak']])>1: warnings.append("多个清晰峰，不能自动认定最高峰为基频")
    selected=next((p for p in peaks if p['selected']),None)
    if not selected or not selected['clear_peak']: warnings.append("所选峰未通过可靠性与清晰度门限，暂停自动结构估计")
    if selected and selected['narrow_peak_check']: warnings.append("窄峰须检查机械/工业噪声来源")
    warnings.append("H/V峰幅值不等于地震放大倍数；响应校正不保证低于仪器噪声底的频段可信")
    summary=dict(schema='ndstudio.hvsr/1',instrument_id=item['instrument_id'],settings=c,windows=len(audit),
                 accepted_windows=len(curves),independent_windows=len(ind),statistics='lognormal; non-overlapping windows',
                 selected_peak=selected,peaks=peaks,warnings=warnings,source_files_opened_read_only=True,
                 source_sha256={p:sha_file(p) for p in sorted({r['file'] for r in allrows})})
    from .hvsr_structure import analyze_structure
    ambiguous=len([p for p in peaks if p['clear_peak']])>1 and c['peak_frequency_hz'] is None
    summary['structure']=analyze_structure(f,mean,std,None if ambiguous else selected,c,output)
    if summary['structure']['status']=='incomplete_peak_shape':
        summary['warnings'].append(summary['structure']['explanation'])
    if ambiguous:
        summary['structure'].update(status='ambiguous_multiple_peaks',explanation='多个清晰峰，须结合谱图/地质信息明确基频后填写指定峰，再计算结构。')
        write_json(output/'shallow_structure.json',summary['structure'])
    write_json(output/'hvsr_summary.json',summary)
    from .hvsr_reports import report
    report(output,table,curves,np.asarray(starts),summary)
    write_json(output/'hvsr_complete.json',dict(instrument_id=item['instrument_id'],summary_sha256=hashlib.sha256((output/'hvsr_summary.json').read_bytes()).hexdigest()))
    print(f"H/V complete: {item['instrument_id']}, {len(curves)} valid / {len(ind)} independent windows",flush=True)
    if selected: print(f"H/V peak: {selected['f0_hz']:.4g} Hz; reliability={selected['reliability_pass']}; clarity={selected['clarity_pass_count']}/6",flush=True)
    return summary
