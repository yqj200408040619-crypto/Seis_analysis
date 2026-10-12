"""Curve-guided exact-frequency Goertzel screening. No amplitude rejection by default.
This is conditional selection, not independent validation of an input dispersion curve.
"""
from pathlib import Path
import json, hashlib, math, os, zipfile, tempfile
import numpy as np
import pandas as pd
from scipy import signal
from numba import njit
from dispersion.common import csv, write_json, load_config, parse_utc, intervals, sha
from dispersion.prepare import station_table, read_station, response_table, correct_response, resample_align
from .phase_screening import screen_cross, amplitude_flags, geometric_fraction, accumulation_series
from .phase_records import available_phase_extent, NativePhaseReader, record_headers, phase_window_geometry, full_window_starts


def wrap(x):
    return np.angle(np.exp(1j*np.asarray(x)))


def calibration_signature(root,c):
    relevant={k:c[k] for k in ('sample_rate_hz','filter_guard_s','bandpass_hz','response')}
    table=Path(root)/'stations.csv';relevant['stations_sha256']=sha(table) if table.exists() else None
    relevant['response_code_sha256']=sha(Path(__file__).parent.parent/'core/dispersion/prepare.py')
    return hashlib.sha256(json.dumps(relevant,sort_keys=True).encode()).hexdigest()


def phase_interval_mask(phase, center, low, high):
    """Closed circular arc [center+low, center+high], offsets are unwrapped radians."""
    if not np.isfinite([center, low, high]).all() or high < low:
        raise ValueError('Invalid phase interval')
    if high-low >= 2*np.pi-1e-12:
        return np.isfinite(phase)
    midpoint=center+(low+high)/2
    return np.isfinite(phase) & (abs(wrap(np.asarray(phase)-midpoint)) <= (high-low)/2+1e-12)


@njit(cache=True)
def goertzel_windows(x, starts, n, f, fs, taper):
    """Exact requested frequency, standard DFT sign, demean, preserve NaN invalidity."""
    w=2*np.pi*f/fs;coeff=2*np.cos(w)
    out=np.empty(len(starts),np.complex128)
    rotation=np.exp(-1j*w*(n-1))
    for k in range(len(starts)):
        a=starts[k];mean=0.;good=True
        for j in range(n):
            v=x[a+j]
            if not np.isfinite(v):good=False
            mean+=v
        if not good:
            out[k]=complex(np.nan,np.nan);continue
        mean/=n;s1=0.;s2=0.
        for j in range(n):
            s0=(x[a+j]-mean)*taper[j]+coeff*s1-s2;s2=s1;s1=s0
        out[k]=(s1-np.exp(-1j*w)*s2)*rotation
    return out


def export_targets(project, output=None, tolerance_deg=15., range_source='fixed',screening_priors_only=False):
    root=Path(project).resolve();src=root/'output';dest=Path(output or src/'phase_guided');dest.mkdir(parents=True,exist_ok=True)
    if not np.isfinite(tolerance_deg) or not 0<=tolerance_deg<180:
        raise ValueError('Phase margin must be finite, >=0 and <180 degrees')
    if range_source not in ('fixed','measured'):raise ValueError('range_source must be fixed or measured')
    source_hashes={n:sha(src/n) for n in ('conditional_dispersion.csv','pairs.csv','measured_phase_segments.csv')}
    curves=pd.read_csv(src/'conditional_dispersion.csv');pairs=pd.read_csv(src/'pairs.csv')
    front=src/'music_front/music_complete.json'
    if front.exists():
        from .music_front import current
        if not current(root):raise ValueError('前段MUSIC结果已过期，请重新运行music阶段')
        corrected=pd.read_csv(src/'music_front/dispersion.csv')
        if corrected.empty and screening_priors_only:
            curves['target_role']='apparent_screening_prior_only';curves['origin_ccf_is_model_prediction']=curves.is_model_prediction;curves['is_model_prediction']=True
        else:curves=corrected
        source_hashes['music_front/dispersion.csv']=sha(src/'music_front/dispersion.csv')
    old_export=dest/'phase_export.json'
    if old_export.exists():
        old=json.loads(old_export.read_text(encoding='utf-8'))
        if old.get('source_files')!=source_hashes or old.get('tolerance_deg')!=tolerance_deg or old.get('range_source')!=range_source or old.get('screening_priors_only',False)!=screening_priors_only:
            (dest/'phase_stack_complete.json').unlink(missing_ok=True)
            (dest/'phase_report.html').unlink(missing_ok=True)
            from .phase_velocity import clear_comparison
            clear_comparison(dest/'velocity_comparison')
    measured=pd.read_csv(src/'measured_phase_segments.csv')
    if curves.empty:raise ValueError('No conditional curves to export')
    keys=['segment_id','frequency_hz']
    if measured.duplicated(keys).any():raise ValueError('Duplicate measured segment/frequency')
    q=curves.copy() if 'distance_m' in curves else curves.merge(pairs[['pair_index','distance_m']],on='pair_index',validate='many_to_one')
    q=q.merge(measured[keys+['measured_unwrapped_phase_rad','relative_phase_q025_rad','relative_phase_q975_rad']],on=keys,validate='many_to_one')
    if len(q)!=len(curves):raise ValueError('Every curve point must have corresponding measurement and pair')
    f=q.frequency_hz.to_numpy();v=q.phase_velocity_m_s.to_numpy();d=(abs(q.projected_distance_m).to_numpy() if 'projection_applied' in q and q.projection_applied.all() else q.distance_m.to_numpy())
    if not np.isfinite(np.c_[f,v,d]).all() or (f<=0).any() or (v<=0).any() or (d<=0).any():raise ValueError('Invalid frequency/velocity/distance')
    # Existing families: k*d = source_phase - theta + 2*pi*N.
    # Negative-lag CCF was reversed before spectrum(); return to raw C=Gj*conj(Gi).
    branch_sign=q.branch.map({'positive':1.,'negative':-1.}).to_numpy()
    if not np.isfinite(branch_sign).all():raise ValueError('Unknown CCF lag branch')
    theta=q.source_phase_rad.to_numpy()-2*np.pi*f*d/v+2*np.pi*q.constant_cycle.to_numpy()
    err=abs(wrap(theta-q.measured_unwrapped_phase_rad.to_numpy()))
    if err.max()>1e-6:raise ValueError('Curve and measured phase have incompatible source/cycle convention')
    q['predicted_unwrapped_phase_rad']=branch_sign*theta
    q['predicted_phase_deg']=np.degrees(wrap(branch_sign*theta))
    margin=np.radians(tolerance_deg)
    if range_source=='measured':
        lo=q.relative_phase_q025_rad.to_numpy();hi=q.relative_phase_q975_rad.to_numpy()
        if not np.isfinite(np.c_[lo,hi]).all() or (hi<lo).any():raise ValueError('Invalid measured phase bounds')
        rawlo=branch_sign*lo;rawhi=branch_sign*hi
        q['phase_offset_low_deg']=np.degrees(np.minimum(rawlo,rawhi)-margin)
        q['phase_offset_high_deg']=np.degrees(np.maximum(rawlo,rawhi)+margin)
    else:
        q['phase_offset_low_deg']=-tolerance_deg;q['phase_offset_high_deg']=tolerance_deg
    q['range_source']=range_source;q['phase_margin_deg']=tolerance_deg
    q['velocity_type']='wave_vector_phase_velocity' if 'projection_applied' in q and q.projection_applied.all() else 'conditional_pair_velocity';q['cross_convention']='Gtarget*conj(Greference)'
    # These bootstrap offsets describe a stacked gated CCF, not confidence bounds for a single short window.
    q['calibrated_single_window_interval']=False
    q['phase_equivalence_id']=q.segment_id
    q['prediction_reconstruction_error_rad']=err
    csv(dest/'curve_velocity_phase_table.csv',q)
    # Velocity bounds implied by the original segment's relative phase offsets; conditional only.
    kr=2*np.pi*f*d/v;lowden=kr-q.relative_phase_q975_rad.to_numpy();highden=kr-q.relative_phase_q025_rad.to_numpy()
    q['velocity_min_m_s']=np.where(highden>0,2*np.pi*f*d/highden,np.nan)
    q['velocity_max_m_s']=np.where(lowden>0,2*np.pi*f*d/lowden,np.nan)
    csv(dest/'curve_velocity_table.csv',q[['curve_id','segment_id','pair','branch','frequency_hz','phase_velocity_m_s','velocity_min_m_s','velocity_max_m_s','source_phase_rad','constant_cycle','velocity_type','range_source','calibrated_single_window_interval']])
    write_json(dest/'phase_export.json',dict(curves=int(q.curve_id.nunique()),points=len(q),phase_equivalence_groups=int(q.segment_id.nunique()),
        screening_priors_only=bool(screening_priors_only),tolerance_deg=tolerance_deg,range_source=range_source,max_reconstruction_error_rad=float(err.max()),
        source_files=source_hashes,
        interpretation='Conditional source phase retained; negative lag branch conjugated to raw cross-spectrum. Cycle aliases cannot be resolved on the originating pair.'))
    return q,dest


def align_project(root, station_names, start, end, aligned_path, workdir=None, resources=None):
    """Correct native data in bounded blocks; retain gaps on the full UTC grid."""
    import h5py
    c=load_config(Path(root)/'config.json');table,_,_=station_table(c)
    table=table[table.station.isin(station_names)]
    if set(table.station)!=set(station_names):raise ValueError('Missing requested station')
    origin=parse_utc(start);finish=parse_utc(end);fs=float(c['sample_rate_hz'])
    count=int(math.floor(float(finish-origin)*fs+1e-7))
    if count<2:raise ValueError('End must be after start')
    arrays={};metadata=[];guard=max(float(c['filter_guard_s']),60.)
    sos=signal.butter(4,c['bandpass_hz'],fs=fs,btype='bandpass',output='sos')
    block_samples=max(2,int(21600*fs))
    for _,row in table.iterrows():
        print('phase align',row.station,flush=True)
        if workdir is None:
            # Keep the established public in-memory API for small callers.
            stream,files=read_station(row,c);records=None;y=np.full(count,np.nan)
        else:
            records=record_headers(row,c);files=sorted({r['path'] for r in records})
            y=np.memmap(Path(workdir)/(row.station+'.float64'),dtype='float64',mode='w+',shape=(count,));y[:]=np.nan
            reader=NativePhaseReader(row,c,records,workdir)
            if resources is not None:resources.extend([y,reader])
        info=None
        for lo in range(0,count,block_samples):
            hi=min(count,lo+block_samples);block_start=origin+lo/fs;block_end=origin+hi/fs
            if records is None:
                st=stream.copy()
            else:
                st=reader.read(block_start-guard-2,block_end+guard+2)
            st.trim(block_start-guard-2,block_end+guard+2)
            if not st:continue
            response,current=response_table(row,st,c);info=current
            if response is None and not row['precorrected']:raise ValueError('Nominal response correction required')
            for tr in st:
                x=np.asarray(tr.data,float);valid=np.isfinite(x)
                for a,b in intervals(np.diff(x)==0):
                    if (b-a+1)/tr.stats.sampling_rate>=1.:valid[a:b+1]=False
                for a,b in intervals(valid):
                    if (b-a)/tr.stats.sampling_rate<=2*guard+2:continue
                    native_c=dict(c,sample_rate_hz=tr.stats.sampling_rate)
                    run=signal.detrend(correct_response(x[a:b]*float(row.polarity),response,native_c),type='linear')
                    from fractions import Fraction
                    ratio=Fraction(fs/tr.stats.sampling_rate).limit_denominator(100000)
                    filtered=signal.resample_poly(run,ratio.numerator,ratio.denominator,padtype='line')
                    filtered=signal.sosfiltfilt(sos,filtered);trim=int(math.ceil(guard*fs))+33
                    if len(filtered)<=2*trim:continue
                    filtered=filtered[trim:-trim]
                    first,last,z=resample_align(filtered,fs,fs,tr.stats.starttime+a/tr.stats.sampling_rate+trim/fs,block_start,hi-lo)
                    if last<=first:continue
                    if np.isfinite(y[lo+first:lo+last]).any():raise ValueError('Unresolved overlapping UTC runs')
                    y[lo+first:lo+last]=z
            print(f'phase align {row.station}: {hi}/{count} samples',flush=True)
        if workdir:reader.close()
        if info is None:raise ValueError('No requested-time samples: '+row.station)
        arrays[row.station]=y
        valid_samples=sum(int(np.isfinite(y[i:i+block_samples]).sum()) for i in range(0,count,block_samples))
        metadata.append(dict(station=row.station,seed_id=row.seed_id,valid_samples=valid_samples,total_samples=count,
                             response_status=info['response_status'],source_files=[dict(path=p,sha256=sha(p)) for p in files]))
    aligned_path=Path(aligned_path);aligned_path.parent.mkdir(parents=True,exist_ok=True)
    with h5py.File(aligned_path,'w') as h:
        h.attrs.update(fs_hz=fs,start_utc=str(origin),end_utc_exclusive=str(origin+count/fs),units={'VEL':'m/s','ACC':'m/s^2','DISP':'m'}[c['response']['output']],response_corrected=True)
        h.attrs['station_metadata_json']=json.dumps(metadata,ensure_ascii=False)
        h.attrs['preprocessing_config_sha256']=sha(Path(root)/'config.json')
        h.attrs['calibration_signature']=calibration_signature(root,c)
        for name,y in arrays.items():
            g=h.create_group('stations/'+name)
            records=record_headers(table[table.station==name].iloc[0],c) if workdir else None
            if records:
                g.attrs['record_start_sample']=int(np.ceil(float(records[0]['start']-origin)*fs-1e-7))
                g.attrs['record_end_sample_exclusive']=int(np.ceil(float(max(r['end'] for r in records)-origin)*fs-1e-7))
            data=g.create_dataset('data',shape=(count,),dtype='float64',compression='gzip')
            mask=g.create_dataset('valid',shape=(count,),dtype='uint8',compression='gzip')
            for i in range(0,count,block_samples):
                chunk=np.asarray(y[i:i+block_samples]);data[i:i+len(chunk)]=np.nan_to_num(chunk);mask[i:i+len(chunk)]=np.isfinite(chunk).astype(np.uint8)
    return arrays,fs,origin,metadata


def robust_flags(amplitude):
    return amplitude_flags(amplitude)[:2]


def stack_stats(cross, mask):
    z=cross[mask];n=len(z)
    if not n:return dict(count=0,raw_real=0.,raw_imag=0.,raw_phase_deg=None,raw_mean_amplitude=None,unit_phase_deg=None,unit_resultant_length=None)
    raw=z.sum();unit=(z/abs(z)).sum()
    return dict(count=n,raw_real=float(raw.real),raw_imag=float(raw.imag),raw_phase_deg=float(np.degrees(np.angle(raw))) if abs(raw)>0 else None,
                raw_mean_amplitude=float(abs(raw)/n),unit_phase_deg=float(np.degrees(np.angle(unit))) if abs(unit)>0 else None,unit_resultant_length=float(abs(unit)/n))


def save_windows(path, **arrays):
    """Publish only a closed, complete window archive; preserve any old file on failure."""
    path=Path(path);temporary=path.with_suffix('.npz.part')
    try:
        with temporary.open('wb') as handle:
            np.savez_compressed(handle,**arrays)
            handle.flush();os.fsync(handle.fileno())
        with zipfile.ZipFile(temporary) as archive:
            if archive.testzip() is not None:raise ValueError('Window archive CRC check failed')
        os.replace(temporary,path)
    finally:
        temporary.unlink(missing_ok=True)


def run_phase_stack(project,output=None,curve_ids=None,start=None,end=None,repeat_cycles=20.,slide_cycles=1.,tolerance_deg=15.,range_source='fixed',taper='none',aligned_h5=None,time_mode=None,stabilization_s=600.,frequency_method='matlab_bin',amplitude_rejection=True,split_branches=True,screening_mode='curve',direction=None,monotonic_only=False,prior=None,reconstruction=None):
    from .acquisition import acquisition_lock
    from .direction import settings
    from .trace import Trace
    options={k:v for k,v in locals().copy().items() if k not in ('settings','Trace','acquisition_lock')}
    direction=settings(direction or {'mode':'off'});options['direction']=direction
    root_path=Path(project)
    if (root_path/'output/music_front/music_complete.json').exists():
        from .music_front import current,run_music_front
        if not current(root_path):run_music_front(root_path)
    from .reconstruction import settings as reconstruction_settings
    reconstruction=reconstruction_settings(reconstruction);options['reconstruction']=reconstruction
    if reconstruction['enabled']:
        if direction['mode'] not in ('filter','joint'):raise ValueError('低频投影重建必须启用方向筛窗（filter/joint）')
        if frequency_method!='exact' or screening_mode!='curve':raise ValueError('低频重建须选择精确Goertzel频率及曲线引导门限')
        if not 0<tolerance_deg<90:raise ValueError('低频重建相位半宽须在0–90度之间')
    options={k:str(v) if isinstance(v,Path) else v for k,v in options.items()}
    options['project_config_sha256']=sha(Path(project)/'config.json')
    options['module_sha256']=sha(__file__)
    dest=Path(output or Path(project)/'output/phase_guided')
    dest.mkdir(parents=True,exist_ok=True)
    with acquisition_lock(dest), tempfile.TemporaryDirectory(prefix='.phase_work_',dir=dest) as workdir, Trace(dest,options,'phase') as trace:
        resources=[]
        try:
            return _run_phase_stack(project,output,curve_ids,start,end,repeat_cycles,slide_cycles,tolerance_deg,range_source,taper,aligned_h5,time_mode,workdir,stabilization_s,frequency_method,amplitude_rejection,split_branches,screening_mode,resources,direction,trace,monotonic_only,prior,reconstruction)
        except Exception:
            for name in ('phase_stack_complete.json','phase_report.html','inversion_curves.csv','inversion_catalog.csv'):(dest/name).unlink(missing_ok=True)
            raise
        finally:
            # Close every mapping before TemporaryDirectory removal, also on failure/Windows.
            for item in reversed(resources):
                if isinstance(item,np.memmap):
                    if not item._mmap.closed:item._mmap.close()
                else:item.close()


def _run_phase_stack(project,output=None,curve_ids=None,start=None,end=None,repeat_cycles=20.,slide_cycles=1.,tolerance_deg=15.,range_source='fixed',taper='none',aligned_h5=None,time_mode=None,workdir=None,stabilization_s=600.,frequency_method='matlab_bin',amplitude_rejection=True,split_branches=True,screening_mode='curve',resources=None,direction=None,trace=None,monotonic_only=False,prior=None,reconstruction=None):
    root=Path(project).resolve();c=load_config(root/'config.json')
    with trace.stage('候选曲线与相位门限',dict(tolerance_deg=tolerance_deg,range_source=range_source)) as log:
        targets,dest=export_targets(root,output,tolerance_deg,range_source,screening_priors_only=bool(reconstruction and reconstruction['enabled']))
        log.update(curves=int(targets.curve_id.nunique()),frequency_points=len(targets),output=str(dest/'curve_velocity_phase_table.csv'))
    from .phase_velocity import clear_comparison,run_phase_velocity
    clear_comparison(dest/'velocity_comparison')
    marker=dest/'phase_stack_complete.json';marker.unlink(missing_ok=True);(dest/'phase_report.html').unlink(missing_ok=True)
    for name in ('qualified_phase_summary.csv','inversion_curves.csv','inversion_catalog.csv','phase_continuous_bands.csv'):(dest/name).unlink(missing_ok=True)
    import shutil
    if (dest/'direction').exists():shutil.rmtree(dest/'direction')
    time_mode=time_mode or ('custom' if start or end else 'full')
    if time_mode not in ('full','custom'):raise ValueError('time_mode must be full or custom')
    if time_mode=='custom' and (not start or not end):raise ValueError('Custom interval requires both UTC start and end')
    if time_mode=='full':start=end=None
    record_bounds=[];pair_ranges=[]
    if not np.isfinite(stabilization_s) or stabilization_s<0:raise ValueError('Stabilization seconds must be finite and >=0')
    if frequency_method not in ('exact','matlab_bin'):raise ValueError('frequency_method must be exact or matlab_bin')
    if screening_mode not in ('curve','matlab'):raise ValueError('screening_mode must be curve or matlab')
    if not np.isfinite([repeat_cycles,slide_cycles]).all() or min(repeat_cycles,slide_cycles)<=0:raise ValueError('Cycles must be finite and positive')
    if taper not in ('none','hann'):raise ValueError('taper must be none or hann')
    if curve_ids:
        ids=[s.strip() for s in curve_ids.split(',') if s.strip()] if isinstance(curve_ids,str) else list(curve_ids)
        missing=set(ids)-set(targets.curve_id)
        if missing:raise ValueError('Unknown curves: '+','.join(sorted(missing)))
        targets=targets[targets.curve_id.isin(ids)].copy()
    if targets.empty:raise ValueError('Select at least one curve')
    from .curve_prior import screening_prior
    targets=screening_prior(targets,dest,c['scan'],None if reconstruction and reconstruction['enabled'] else prior,monotonic_only or bool(reconstruction and reconstruction['enabled']),trace)
    csv(dest/'selected_phase_targets.csv',targets)
    names=set(s for pair in targets.pair for s in pair.split('--'))
    xy_table=None
    if direction['mode']!='off':
        table,_,linear=station_table(c);xy_table=table.set_index('station');names=set(table.station)
        trace.note('阵列几何',f'{len(names)}台；近直线阵列={linear}；传播方向：北0°、顺时针；最大台对夹角={direction["max_pair_angle_deg"]:g}°',stations=sorted(names))
    trace.note('波形校正与UTC对齐',f'输入：{"校正H5" if aligned_h5 else "全部原始记录"}；采样率={c["sample_rate_hz"]:g} Hz；滤波={c["bandpass_hz"]}；响应输出={c["response"]["output"]}')
    front_aligned=root/'output/music_front/aligned_waveforms.h5'
    if not aligned_h5 and time_mode=='full' and (root/'output/music_front/music_complete.json').exists() and front_aligned.exists():
        from .music_front import current
        if current(root):aligned_h5=front_aligned;trace.note('复用前段MUSIC校正波形','沿用同一UTC时间轴、缺口及响应校正；逐频重算MUSIC和相位筛窗。')
    if not aligned_h5 and (dest/'aligned_waveforms.h5').exists():
        import h5py
        available_start,available_end,record_bounds=available_phase_extent(root,names)
        requested_start,requested_end=(str(available_start),str(available_end)) if time_mode=='full' else (start,end)
        target_origin=parse_utc(requested_start);target_finish=parse_utc(requested_end)
        target_count=int(math.floor(float(target_finish-target_origin)*c['sample_rate_hz']+1e-7))
        requested_end=str(target_origin+target_count/c['sample_rate_hz'])
        candidate=dest/'aligned_waveforms.h5'
        try:
            with h5py.File(candidate) as h:
                match=(h.attrs.get('calibration_signature')==calibration_signature(root,c) and parse_utc(h.attrs['start_utc']).ns==parse_utc(requested_start).ns and
                    parse_utc(h.attrs['end_utc_exclusive']).ns==parse_utc(requested_end).ns and names<=set(h['stations']))
                cached_meta=json.loads(h.attrs.get('station_metadata_json','[]'))
                if match:
                    for info in cached_meta:
                        if info['station'] in names:
                            for source in info['source_files']:
                                if not Path(source['path']).exists() or sha(source['path'])!=source['sha256']:match=False;break
            if match:
                aligned_h5=candidate;trace.note('复用校正波形','UTC、校正参数、台站表与原始文件校验均一致；重算方向、筛选与频散。')
        except (OSError,KeyError,ValueError):pass
    if aligned_h5:
        import h5py,shutil
        path=Path(aligned_h5).resolve();c=load_config(root/'config.json')
        with h5py.File(path,'r') as h:
            if not h.attrs.get('response_corrected',False):raise ValueError('Response-corrected aligned H5 required')
            expected_unit={'VEL':'m/s','ACC':'m/s^2','DISP':'m'}[c['response']['output']]
            if h.attrs.get('units')!=expected_unit:raise ValueError('Aligned H5 units do not match project')
            if 'calibration_signature' in h.attrs:
                if h.attrs['calibration_signature']!=calibration_signature(root,c):raise ValueError('Aligned H5 calibration differs from project')
            elif 'preprocessing_config_sha256' in h.attrs and h.attrs['preprocessing_config_sha256']!=sha(root/'config.json'):raise ValueError('Aligned H5 preprocessing config differs from project')
            fs=float(h.attrs['fs_hz']);origin=parse_utc(h.attrs['start_utc']);finish=parse_utc(h.attrs['end_utc_exclusive'])
            if time_mode=='full':start,end=str(origin),str(finish)
            if abs(fs-c['sample_rate_hz'])>1e-9 or abs(origin-parse_utc(start))>1e-9 or abs(finish-parse_utc(end))>1e-9:raise ValueError('Aligned H5 must match requested UTC interval and sample rate exactly')
            arrays={};station_limits={}
            for name in names:
                g=h['stations/'+name];data=g['data'];valid=g['valid']
                if data.ndim!=1 or valid.shape!=data.shape:raise ValueError('Invalid aligned data/mask')
                x=np.memmap(Path(workdir)/(name+'.float64'),dtype='float64',mode='w+',shape=data.shape) if workdir else np.empty(data.shape)
                for i in range(0,len(x),1000000):x[i:i+1000000]=np.where(valid[i:i+1000000],data[i:i+1000000],np.nan)
                arrays[name]=x
                if isinstance(x,np.memmap) and resources is not None:resources.append(x)
                station_limits[name]=(int(g.attrs.get('record_start_sample',0)),int(g.attrs.get('record_end_sample_exclusive',len(x))))
            if len({len(y) for y in arrays.values()})!=1 or len(next(iter(arrays.values())))!=round(float(finish-origin)*fs):raise ValueError('Aligned H5 sample count mismatch')
            meta=json.loads(h.attrs.get('station_metadata_json','[]'))
        if path!=(dest/'aligned_waveforms.h5').resolve():shutil.copyfile(path,dest/'aligned_waveforms.h5')
    else:
        available_start,available_end,record_bounds=available_phase_extent(root,names)
        if time_mode=='full':start,end=str(available_start),str(available_end)
        arrays,fs,origin,meta=align_project(root,names,start,end,dest/'aligned_waveforms.h5',workdir,resources)
        station_limits={r['station']:(int(np.ceil(float(parse_utc(r['start_utc'])-origin)*fs-1e-7)),int(np.ceil(float(parse_utc(r['end_utc_exclusive'])-origin)*fs-1e-7))) for r in record_bounds}
    print(f'相位处理范围：{time_mode} · {start} 至 {end}（结束不含）；逐频窗长 T=N/f',flush=True)
    trace.note('波形校正结果',f'{len(arrays)}台；{len(next(iter(arrays.values()))):,}样本/台；UTC {start} → {end}',output=str(dest/'aligned_waveforms.h5'),station_metadata=meta)
    session=None
    if direction['mode']!='off':
        from .direction import DirectionSession
        session=DirectionSession(arrays,xy_table,fs,origin,station_limits,dest,workdir,direction,stabilization_s,trace,amplitude_rejection)
    if reconstruction and reconstruction['enabled']:
        from .reconstruction import run
        return run(root,dest,targets,session,xy_table,c,reconstruction,repeat_cycles,slide_cycles,taper,amplitude_rejection,trace,meta)
    windows=dest/'windows';windows.mkdir(exist_ok=True);summaries=[];gate_rows=[];illustrations=[]
    for old in windows.glob('*.npz'):old.unlink()
    sample_count=len(next(iter(arrays.values())))
    # Pair-local starts are important: an unrelated station must not shift another pair's windows.
    for pair,pf in targets.groupby('pair',sort=True):
        a,b=pair.split('--');raw_first=max(station_limits[a][0],station_limits[b][0]);raw_last=min(station_limits[a][1],station_limits[b][1])
        first=max(0,raw_first+int(np.floor(stabilization_s*fs+.5)));last=min(sample_count,raw_last)
        pair_ranges.append(dict(pair=pair,common_start_utc=str(origin+raw_first/fs),common_end_utc_exclusive=str(origin+raw_last/fs),
                                processing_start_utc=str(origin+first/fs),processing_end_utc_exclusive=str(origin+last/fs),usable_duration_s=max(0.,(last-first)/fs),unavailable_reason='共同记录不足起点稳定时间' if first>=last else ''))
        for frequency,tf in pf.groupby('frequency_hz',sort=True):
            geometry=phase_window_geometry(float(frequency),fs,repeat_cycles,slide_cycles)
            n=geometry['window_samples'];step=geometry['slide_samples']
            # An unavailable pair/frequency remains explicit, rather than aborting other pairs.
            starts=full_window_starts(last-first,n,step)+first if last-first>=n else np.empty(0,np.int64)
            actual_frequency=float(frequency) if frequency_method=='exact' else int(np.floor(float(frequency)*n/fs+.5))*fs/n
            if len(starts) and not 0<actual_frequency<fs/2:raise ValueError('Selected MATLAB Goertzel bin is zero or outside Nyquist')
            print(f'phase {pair} {frequency:g} Hz: Goertzel={actual_frequency:g} Hz, window={n/fs:g} s, hop={step/fs:g} s, windows={len(starts)}',flush=True)
            window_taper=np.ones(n) if taper=='none' else signal.windows.hann(n,sym=True);gs={};direction_data={};eligible=None;selected_group=-1
            if session:
                from .direction import pair_gate,independent_indices
                all_direction=session.get(float(frequency),actual_frequency,n,step,window_taper)
                chosen=(all_direction['start_sample']>=first)&(all_direction['start_sample']+n<=last)
                starts=all_direction['start_sample'][chosen]
                local={k:v[chosen] for k,v in all_direction.items() if k not in ('coefficients','start_sample')}
                gs={station:all_direction['coefficients'][session.names.index(station),chosen] for station in (a,b)}
                delta=xy_table.loc[b,['x_m','y_m']].to_numpy(float)-xy_table.loc[a,['x_m','y_m']].to_numpy(float)
                direction_data=pair_gate(local,*delta,direction)
                direction_data.update(propagation_azimuth_deg=local['azimuth_deg'],direction_uncertainty_deg=local['direction_uncertainty_deg'],direction_group=local['direction_group'],music_velocity_m_s=local['velocity_m_s'])
                candidate=direction_data['direction_valid'].copy();cross0=gs[b]*np.conj(gs[a]);candidate&=np.isfinite(cross0)&(abs(cross0)>0)
                if amplitude_rejection:
                    large,jump,*_=amplitude_flags(abs(cross0));candidate&=~(large|jump)
                if candidate.any():
                    groups,counts=np.unique(local['direction_group'][candidate],return_counts=True);selected_group=int(groups[np.argmax(counts)])
                eligible=direction_data['direction_valid']&(local['direction_group']==selected_group)
                direction_data['direction_eligible']=eligible
                direction_data['direction_reason'][(direction_data['direction_reason']=='accepted')&~eligible]='other_direction'
                from .direction import REASONS
                trace.note('台对角度筛选',f'{pair} · {frequency:g} Hz · {len(starts)}窗 → 几何合格{int(direction_data["direction_valid"].sum())} → 主方向组{selected_group}保留{int(eligible.sum())}',
                    reasons={REASONS[k]:int((direction_data['direction_reason']==k).sum()) for k in np.unique(direction_data['direction_reason'])})
                if direction['mode']=='diagnostic':eligible=None
            for station in (() if session else (a,b)):
                out=np.memmap(Path(workdir)/(station+'.complex128'),dtype='complex128',mode='w+',shape=(len(starts),)) if workdir and len(starts) else np.empty(len(starts),complex)
                for i in range(0,len(starts),65536):
                    ix=starts[i:i+65536];low=int(ix[0]);high=int(ix[-1])+n
                    out[i:i+len(ix)]=goertzel_windows(arrays[station][low:high],ix-low,n,actual_frequency,fs,window_taper)
                gs[station]=out
                if isinstance(out,np.memmap) and resources is not None:resources.append(out)
            cross=gs[b]*np.conj(gs[a]);phase=np.angle(cross)
            for segment,q in tf.groupby('curve_id' if session and direction['mode'] in ('filter','joint') else 'segment_id',sort=False):
                r=q.iloc[0];center=np.radians(r.predicted_phase_deg);lo=np.radians(r.phase_offset_low_deg);hi=np.radians(r.phase_offset_high_deg)
                if np.max(abs(wrap(np.radians(q.predicted_phase_deg.to_numpy())-center)))>1e-8 or not np.allclose(q.phase_offset_low_deg,r.phase_offset_low_deg) or not np.allclose(q.phase_offset_high_deg,r.phase_offset_high_deg):raise ValueError('Inconsistent phase-equivalent curve gates')
                projection_fields={}
                cross=gs[b]*np.conj(gs[a])
                if session and direction['mode'] in ('filter','joint'):
                    from .music_front import frequency_reference
                    from .reconstruction import align_cross
                    ref=frequency_reference(local,starts,n,delta,direction,r.branch,angle_gate=True)
                    if ref:
                        reference=ref['projected_distance_m']
                        v=float(r.phase_velocity_m_s) if r.get('projection_applied',False) else float(r.phase_velocity_m_s)*abs(reference)/float(r.distance_m)
                        eligible=ref['eligible'];selected_group=ref['selected_direction_group']
                        # Keep the CCF source-phase/cycle gate, while rotating each measured window's projection.
                        aligned,rotation=align_cross(cross,actual_frequency,direction_data['projected_distance_m'],reference,v)
                        projection_fields=dict(raw_cross=cross.copy(),alignment_rotation_rad=rotation,window_projected_distance_m=direction_data['projected_distance_m'],
                            reference_projected_distance_m=reference,prior_true_velocity_m_s=v,confirmed_azimuth_deg=ref['azimuth_deg'])
                        cross=aligned
                        q=q.copy();q['phase_velocity_m_s']=v;q['projected_distance_m']=reference;q['projection_applied']=True
                        q['velocity_type']='wave_vector_phase_velocity';q['direction_method']='music'
                        targets.loc[q.index,['phase_velocity_m_s','projected_distance_m','projection_applied','velocity_type','direction_method']]=q[['phase_velocity_m_s','projected_distance_m','projection_applied','velocity_type','direction_method']]
                    else:eligible=np.zeros(len(starts),bool)
                    direction_data['direction_eligible']=eligible
                    split_branches_used=False
                else:split_branches_used=split_branches
                phase=np.angle(cross)
                masks=screen_cross(cross,center,lo,hi,screening_mode,amplitude_rejection,split_branches_used,eligible)
                branch=masks['phase_sign_branch'];base=masks['valid'];accepted=masks['accepted'];big=masks['amplitude_large'];jump=masks['amplitude_jump']
                if screening_mode=='matlab':
                    center=masks['preliminary_positive_phase_rad'] if branch=='nonnegative' else masks['preliminary_negative_phase_rad'];lo=-np.pi/2;hi=np.pi/2
                control_base=masks['amplitude_valid']
                if eligible is not None:control_base=control_base&eligible
                rotated90=control_base&phase_interval_mask(phase,center+np.pi/2,lo,hi) if np.isfinite(center) else np.zeros(len(phase),bool)
                rotated180=control_base&phase_interval_mask(phase,center+np.pi,lo,hi) if np.isfinite(center) else np.zeros(len(phase),bool)
                if split_branches_used and np.isfinite(center):
                    rotated90 &= (phase>=0) if float(wrap(center+np.pi/2))>=0 else (phase<0)
                    rotated180 &= (phase>=0) if float(wrap(center+np.pi))>=0 else (phase<0)
                gate_id=hashlib.sha256(f'{segment}|{pair}|{frequency:.9f}'.encode()).hexdigest()[:16];filename=gate_id+'.npz'
                save_windows(windows/filename,start_sample=starts,cross=cross,**projection_fields,**masks,**direction_data,direction_applied=eligible is not None,frequency_hz=frequency,actual_frequency_hz=actual_frequency,
                    window_samples=n,slide_samples=step,fs_hz=fs,start_utc=str(origin),pair=pair,segment_id=segment,
                    predicted_phase_rad=center,phase_offset_low_rad=lo,phase_offset_high_rad=hi,curve_ids=q.curve_id.to_numpy(dtype=str),
                    screening_mode=screening_mode,split_branches=bool(split_branches_used),amplitude_rejection=bool(amplitude_rejection),
                    processing_start_sample=first,processing_end_sample_exclusive=last,phase_residual_rad=wrap(phase-center))
                raw=stack_stats(cross,base);selected=stack_stats(cross,accepted);before=wrap(phase[base]-center);after=wrap(phase[accepted]-center)
                row=dict(gate_id=gate_id,segment_id=r.segment_id,pair=pair,branch=r.branch,phase_sign_branch=branch,frequency_hz=frequency,actual_frequency_hz=actual_frequency,
                    predicted_phase_deg=r.predicted_phase_deg,screening_center_deg=float(np.degrees(center)),phase_offset_low_deg=r.phase_offset_low_deg,phase_offset_high_deg=r.phase_offset_high_deg,
                    window_samples=n,slide_samples=step,window_duration_s=n/fs,slide_duration_s=step/fs,
                    processing_start_utc=str(origin+first/fs),processing_end_utc_exclusive=str(origin+last/fs),unwindowed_tail_samples=int(last-(starts[-1]+n)) if len(starts) else max(0,last-first),
                    unavailable_reason='处理时长不足一个完整逐频窗口' if not len(starts) else '',total_windows=len(starts),valid_windows=int(base.sum()),accepted_windows=int(accepted.sum()),
                    accepted_positive_windows=int(masks['accepted_positive'].sum()),accepted_negative_windows=int(masks['accepted_negative'].sum()),
                    amplitude_valid_windows=int(masks['amplitude_valid'].sum()),amplitude_direction_valid_windows=int(control_base.sum()),amplitude_rejected_windows=int((base&~masks['amplitude_valid']).sum()),
                    accepted_fraction=float(accepted.sum()/base.sum()) if base.any() else np.nan,
                    rotated90_accepted_windows=int(rotated90.sum()),rotated180_accepted_windows=int(rotated180.sum()),
                    uniform_phase_expected_fraction=geometric_fraction(center,lo,hi,branch) if np.isfinite(center) else np.nan,
                    amplitude_large_marked=int((big&base).sum()),amplitude_jump_marked=int((jump&base).sum()),
                    amplitude_marked_but_accepted=int(((big|jump)&accepted).sum()),windows_file='windows/'+filename,
                    amplitude_threshold=masks['amplitude_threshold'],delta_threshold=masks['delta_threshold'],
                    maximum_accepted_arc_residual_deg=float(np.degrees(abs(wrap(phase[accepted]-(center+(lo+hi)/2))).max())) if accepted.any() else np.nan,
                    median_abs_residual_before_deg=float(np.degrees(np.median(abs(before)))) if len(before) else np.nan,
                    median_abs_residual_after_deg=float(np.degrees(np.median(abs(after)))) if len(after) else np.nan)
                if session:
                    row.update(direction_method='music',confirmed_azimuth_deg=ref['azimuth_deg'] if projection_fields else np.nan,projected_distance_m=float(q.projected_distance_m.iloc[0]) if 'projected_distance_m' in q else np.nan,selected_direction_group=selected_group,direction_valid_windows=int(direction_data['direction_valid'].sum()),
                        direction_eligible_windows=int(direction_data['direction_eligible'].sum()),
                        direction_independent_windows=len(independent_indices(starts,n,accepted)),
                        direction_rejection_reasons=json.dumps({k:int((direction_data['direction_reason']==k).sum()) for k in np.unique(direction_data['direction_reason'])},ensure_ascii=False),
                        median_pair_angle_deg=float(np.nanmedian(direction_data['pair_angle_deg'])) if np.isfinite(direction_data['pair_angle_deg']).any() else np.nan)
                trace.note('相位筛选与叠加',f'{segment} · {frequency:g} Hz · 有效{int(base.sum())} → 幅值/方向合格{int(control_base.sum())} → 接受{int(accepted.sum())}；正{int(masks["accepted_positive"].sum())}/负{int(masks["accepted_negative"].sum())}；最终相位={selected["raw_phase_deg"]}',windows_file=str(windows/filename),amplitude_threshold=masks['amplitude_threshold'],delta_threshold=masks['delta_threshold'])
                row.update({'before_'+k:v for k,v in raw.items()});row.update({'after_'+k:v for k,v in selected.items()})
                for sign in ('positive','negative'):row.update({'after_'+sign+'_'+k:v for k,v in stack_stats(cross,masks['accepted_'+sign]).items()})
                gate_rows.append(row)
                for cid in q.curve_id:summaries.append(dict(curve_id=cid,**row))
                if len(starts) and len(illustrations)<6 and all(x['segment']!=segment for x in illustrations):
                    pick=np.unique(np.linspace(0,len(starts)-1,min(len(starts),20000),dtype=int))
                    series={sign:accumulation_series(cross,masks['accepted_'+sign]) for sign in ('positive','negative')} if not illustrations else {}
                    illustrations.append(dict(accumulation=series,segment=segment,pair=pair,f=frequency,phase=phase[pick],base=base[pick],accepted=accepted[pick],cross=cross[pick],center=center,lo=lo,hi=hi,sampled=len(pick)<len(starts)))
            for z in gs.values():
                if isinstance(z,np.memmap):z.flush();z._mmap.close()
            gs.clear()
            csv(dest/'phase_stack_summary.csv',summaries);csv(dest/'phase_gate_summary.csv',gate_rows)
            write_json(dest/'phase_progress.json',dict(complete=False,processed_gates=len(gate_rows),
                last_pair=pair,last_frequency_hz=float(frequency),direction_frequencies=len(session.summaries) if session else 0))
    frame=csv(dest/'phase_stack_summary.csv',summaries);gates=csv(dest/'phase_gate_summary.csv',gate_rows)
    branch_rows=[]
    for r in gate_rows:
        for sign in ('positive','negative'):
            stats={k.removeprefix('after_'+sign+'_'):v for k,v in r.items() if k.startswith('after_'+sign+'_')}
            branch_rows.append(dict(gate_id=r['gate_id'],pair=r['pair'],segment_id=r['segment_id'],frequency_hz=r['frequency_hz'],actual_frequency_hz=r['actual_frequency_hz'],phase_sign_branch=sign,**stats))
    csv(dest/'phase_branch_summary.csv',branch_rows)
    csv(dest/'pair_processing_intervals.csv',pair_ranges)
    csv(dest/'selected_phase_targets.csv',targets)
    observed_targets=targets[~targets.is_prior_prediction].copy()
    observed_summary=frame.merge(observed_targets[['curve_id','frequency_hz']],on=['curve_id','frequency_hz'])
    if session:session.finish(np.sort(observed_targets.frequency_hz.unique()),c['scan'])
    make_plots(dest,gates,illustrations)
    result=dict(complete=True,direction_method='music' if session else None,music_module_sha256=sha(Path(__file__).with_name('music.py')),
        station_geometry_sha256=sha(dest/'direction/station_geometry.csv') if session else None,selected_curves=int(targets.curve_id.nunique()),segments=int(targets.segment_id.nunique()),unique_pair_frequency_gates=len(gates),
        curve_frequency_rows=len(frame),valid_window_evaluations=int(gates.valid_windows.sum()),accepted_window_evaluations=int(gates.accepted_windows.sum()),
        time_mode=time_mode,time_extent_source='aligned_h5' if aligned_h5 else 'imported_station_records',station_record_bounds=record_bounds,
        start_utc=str(origin),end_utc_exclusive=str(origin+sample_count/fs),record_duration_s=sample_count/fs,sample_count=sample_count,sample_rate_hz=fs,
        window_duration_min_s=float(gates.window_duration_s.min()),window_duration_max_s=float(gates.window_duration_s.max()),end_anchored_windows=bool(session),stabilization_s=stabilization_s,frequency_method=frequency_method,screening_mode=screening_mode,split_phase_branches=bool(split_branches),pair_processing_intervals=pair_ranges,
        phase_screening_sha256=sha(Path(__file__).with_name('phase_screening.py')),
        plots_subsampled=any(e['sampled'] for e in illustrations),phase_records_sha256=sha(Path(__file__).with_name('phase_records.py')),repeat_cycles=repeat_cycles,slide_cycles=slide_cycles,taper=taper,tolerance_deg=tolerance_deg,range_source=range_source,
        direction=direction,direction_sha256=sha(Path(__file__).with_name('direction.py')),shared_utc_windows=bool(session),monotonic_only=bool(monotonic_only),prior=prior,
        phase_settings_sha256=sha(root/'phase_settings.json') if (root/'phase_settings.json').exists() else None,
        amplitude_rejection=bool(amplitude_rejection),response_corrected=True,station_metadata=meta,settings_sha256=sha(root/'config.json'),module_sha256=sha(__file__),aligned_h5_sha256=sha(dest/'aligned_waveforms.h5'),
        project_name=json.loads((root/'project.json').read_text(encoding='utf-8-sig')).get('name',root.name),
        warnings=['Input curves guide selection; this test is not independent validation.',
                  'Short overlapping windows are dependent; evaluation counts are not independent sample counts.',
                  'Gated CCF phase is transferred conditionally to raw window cross-spectra; selection cannot establish the source model.',
                  'Cycle aliases on the originating pair have identical gates and cannot be resolved by this test.'])
    if screening_mode=='curve':
        from .phase_bands import handoff
        result['inversion_handoff_points']=handoff(dest,observed_targets,observed_summary,c['scan'],direction,trace,monotonic_only)
    write_json(marker,result)
    write_json(dest/'phase_progress.json',dict(complete=True,processed_gates=len(gate_rows)))
    if screening_mode=='curve' and (gates.after_count>0).any():
        with trace.stage('叠加后频散与原曲线对比') as log:
            run_phase_velocity(dest,dest/'velocity_comparison',qualified=bool(session) or monotonic_only or bool(prior and prior.get('enabled')))
            log.update(output=str(dest/'velocity_comparison/report.html'))
    else:
        print('phase velocity comparison unavailable: MATLAB preliminary-phase mode or empty accepted stack; no velocities invented.',flush=True)
    make_report(dest,result,gates)
    trace.note('结果索引',str(dest/'phase_report.html'),direction_report=str(dest/'direction/report.html') if session else None,
        handoff=str(dest/'inversion_curves.csv'),trace=str(dest/'phase_trace.json'))
    return result


def phase_sensitivity(folder,tolerances=(5.,10.,15.,20.)):
    dest=Path(folder);gates=pd.read_csv(dest/'phase_gate_summary.csv');rows=[]
    for r in gates.itertuples():
        with np.load(dest/r.windows_file,allow_pickle=False) as w:
            if 'screening_mode' in w and str(w['screening_mode'])!='curve':raise ValueError('Tolerance sensitivity requires curve-guided mode')
            cross=w['cross'];base=w['valid'];valid=base.copy();phase=np.angle(cross);center=float(w['predicted_phase_rad'])
            if 'amplitude_valid' in w:valid=valid&w['amplitude_valid']
            if 'direction_applied' in w and bool(w['direction_applied']):valid &= w['direction_eligible']
            branch=str(w['phase_sign_branch']) if 'phase_sign_branch' in w else 'combined'
            if branch=='nonnegative':valid=valid&(phase>=0)
            elif branch=='negative':valid=valid&(phase<0)
            for tol in tolerances:
                if not 0<float(tol)<180:raise ValueError('Sensitivity tolerances must be in (0,180)')
                mask=valid&phase_interval_mask(phase,center,-np.radians(tol),np.radians(tol));n=int(mask.sum());base_count=int(base.sum())
                rows.append(dict(segment_id=r.segment_id,pair=r.pair,frequency_hz=r.frequency_hz,tolerance_deg=tol,valid_windows=base_count,accepted_windows=n,
                    accepted_fraction=n/max(1,base_count),uniform_phase_expected_fraction=geometric_fraction(center,-np.radians(tol),np.radians(tol),branch),unit_resultant_length=stack_stats(cross,mask)['unit_resultant_length']))
    frame=csv(dest/'phase_tolerance_sensitivity.csv',rows)
    summary=frame.groupby('tolerance_deg',as_index=False).agg(valid_windows=('valid_windows','sum'),accepted_windows=('accepted_windows','sum'))
    summary['accepted_fraction']=summary.accepted_windows/summary.valid_windows;summary['uniform_phase_expected_fraction']=[float(np.average(frame[frame.tolerance_deg==tol].uniform_phase_expected_fraction,weights=frame[frame.tolerance_deg==tol].valid_windows)) if frame[frame.tolerance_deg==tol].valid_windows.sum() else np.nan for tol in summary.tolerance_deg]
    csv(dest/'phase_tolerance_sensitivity_summary.csv',summary)
    marker=dest/'phase_stack_complete.json'
    if marker.exists():make_report(dest,json.loads(marker.read_text(encoding='utf-8')),gates)
    return summary


def make_plots(dest,gates,examples):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(13,7),constrained_layout=True)
    for ax,e in zip(axes.flat,examples):
        residual=np.degrees(wrap(e['phase']-e['center']));ax.hist(residual[e['base']],bins=np.arange(-180,181,6),color='#b8c4ce',label='before')
        ax.hist(residual[e['accepted']],bins=np.arange(-180,181,6),color='#168c93',label='accepted')
        ax.axvline(np.degrees(e['lo']),color='#c74a48');ax.axvline(np.degrees(e['hi']),color='#c74a48')
        ax.set(title=f"{e['segment']} | {e['f']:.2f} Hz",xlabel='Circular phase residual (deg)',ylabel='Windows');ax.legend(fontsize=8)
    for ax in list(axes.flat)[len(examples):]:ax.set_visible(False)
    fig.savefig(dest/'phase_distributions.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(11,5),constrained_layout=True)
    for sid,q in gates.groupby('segment_id'):ax.plot(q.frequency_hz,q.accepted_fraction,label=sid)
    geometric=float((gates.uniform_phase_expected_fraction*gates.valid_windows).sum()/max(1,gates.valid_windows.sum()))
    ax.axhline(geometric,ls='--',color='black',label='uniform-phase mean geometric fraction')
    ax.set(xlabel='Frequency (Hz)',ylabel='Accepted / valid windows',ylim=(0,1));ax.legend(fontsize=7,ncol=3);fig.savefig(dest/'retention_by_frequency.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(11,7),constrained_layout=True)
    if examples:
        e=examples[0]
        for sign,series in e.get('accumulation',{}).items():
            n=series['count'];z=series['raw']
            axes[0].plot(n,abs(z),label=sign+' phase branch')
            axes[1].plot(n,np.degrees(np.angle(z)),label=sign+' phase branch')
        axes[0].set(ylabel='Absolute cumulative complex cross-spectrum',title=f"{e['segment']} | {e['f']:.2f} Hz");axes[0].legend()
        axes[1].set(xlabel='All accepted overlapping windows',ylabel='Cumulative raw-stack phase (deg)');axes[1].legend()
    fig.savefig(dest/'accumulation.png',dpi=150);plt.close(fig)


def make_report(dest,result,gates):
    import html,base64
    rows=[]
    for sid,q in gates.groupby('segment_id'):
        rows.append(dict(segment=sid,pair=q.pair.iloc[0],branch=q.branch.iloc[0],frequencies=len(q),valid=int(q.valid_windows.sum()),accepted=int(q.accepted_windows.sum()),
            fraction=q.accepted_windows.sum()/max(1,q.valid_windows.sum()),median_before=q.before_unit_resultant_length.median(),median_after=q.after_unit_resultant_length.median(),
            rotated90_fraction=q.rotated90_accepted_windows.sum()/max(1,q.valid_windows.sum()),rotated180_fraction=q.rotated180_accepted_windows.sum()/max(1,q.valid_windows.sum())))
    body='<html><meta charset="utf-8"><title>曲线相位筛选测试</title><style>body{font:16px sans-serif;max-width:1100px;margin:35px auto;line-height:1.6}table{border-collapse:collapse;font-size:14px}td,th{border:1px solid #ccd;padding:7px}img{width:100%}</style><h1>'+html.escape(result.get('project_name','曲线引导'))+'：Goertzel 相位筛选</h1>'
    scope='全部已导入记录' if result.get('time_mode')=='full' else '自定义UTC区间'
    if result.get('time_extent_source')=='aligned_h5':scope+='（限所提供校正H5的全部时长）'
    body+=f'<p>数据范围：{scope}。UTC {html.escape(result["start_utc"])} 至 {html.escape(result["end_utc_exclusive"])}。{result["selected_curves"]} 条候选曲线，{result["segments"]} 个相位等价组，{result["unique_pair_frequency_gates"]} 个逐频筛选任务。窗口 {result["repeat_cycles"]} 周期，滑动 {result["slide_cycles"]} 周期；范围来源 {result["range_source"]}，容差/额外余量 {result["tolerance_deg"]}°。</p>'
    body+='<p><a href="phase_trace.json">逐阶段追溯JSON</a> · <a href="phase_processing.log">完整处理日志</a> · <a href="phase_settings_used.json">实际参数快照</a> · <a href="phase_continuous_bands.csv">重新检查的连续频带</a></p>'
    if result.get('shared_utc_windows'):
        body+='<p>已启用共同UTC窗口。<a href="direction/report.html">查看方向、夹角剔除与平面波联合频散</a>。曲线叠加选取曲线筛选之前确定的主传播方向组，未用平均角度直接修正CCF速度。</p>'
    body+=f'<p>起点稳定时间 {result.get("stabilization_s",0):g} 秒；频率方式 {html.escape(result.get("frequency_method","exact"))}；筛选模式 {html.escape(result.get("screening_mode","curve"))}；MAD幅值/突变剔除：{result.get("amplitude_rejection",False)}；正负相位分开叠加：{result.get("split_phase_branches",False)}。</p>'
    body+='<p>按各台站对的完整共同记录自动遍历，起点仅顺延一次稳定时间；按曲线频点 f 计算 T=N/f、Δt=M/f，并取整到采样点。去均值后执行 Goertzel，互谱 C=G目标·conj(G参考)。最近DFT频点模式与附件一致，actual_frequency_hz记录实际变换频率；精确频率模式直接计算曲线频点。先筛幅值/突变，再按 ≥0° / &lt;0° 分支筛相位、复数累加，不共轭、不反号。主结果使用预测相位同号分支，两分支完整结果见 phase_branch_summary.csv。缺测、滤波边缘和无定义相位不参与叠加。</p>'
    body+='<p>曲线引导模式的门限来自逐频速度、台距、源相位及CCF时延分支；MATLAB原门限模式则使用各相位分支初叠加相位 ±90°。相位正负与CCF正负时延是不同分类。最终相位 arg(ΣC)，另存等权相位 arg(ΣC/|C|)，不直接平均角度。</p>'
    if result.get('pair_processing_intervals'):body+='<h2>每对台站实际处理的全部共同记录</h2>'+pd.DataFrame(result['pair_processing_intervals']).to_html(index=False,float_format=lambda x:f'{x:.3f}')
    if result.get('plots_subsampled'):body+='<p>示例分布图最多取20000个等距窗口展示；累积图对全部保留窗口求和后抽取显示点。CSV、窗口档案及最终叠加使用全部窗口。</p>'
    body+='<p>本次是条件筛选。筛选后相位集中是门限作用，不能作为独立验证。窗口高度重叠，计数不可当作独立样本数。周期别名在原台站对上无法区分。</p>'
    if (dest/'velocity_comparison/report.html').exists():body+='<p><a href="velocity_comparison/report.html">查看叠加后波速与原候选曲线对比及连续频带检查结果</a>。</p>'
    else:body+='<p>本模式或当前保留窗不足，未生成曲线引导的速度对比。</p>'
    geometry=gates[['frequency_hz','actual_frequency_hz','window_duration_s','slide_duration_s']].drop_duplicates().sort_values('frequency_hz') if 'window_duration_s' in gates else None
    if geometry is not None:body+='<h2>逐频实际窗长与滑动步长</h2>'+geometry.to_html(index=False,float_format=lambda x:f'{x:.6f}')
    body+='<h2>两相位分支最终复数叠加</h2>'+pd.read_csv(dest/'phase_branch_summary.csv').to_html(index=False,float_format=lambda x:f'{x:.6g}') if (dest/'phase_branch_summary.csv').exists() else ''
    body+=pd.DataFrame(rows).to_html(index=False,float_format=lambda x:f'{x:.4f}')
    actual=float(gates.accepted_windows.sum()/max(1,gates.valid_windows.sum()))
    uniform=float((gates.uniform_phase_expected_fraction*gates.valid_windows).sum()/max(1,gates.valid_windows.sum()))
    body+=f'<p>按有效窗口数加权，总保留比例 {actual:.4%}；均匀相位的几何比例 {uniform:.4%}。旋转中心90°/180°的保留比例列只是分布对照，不是显著性检验。幅值标记但仍通过相位筛选的贡献数：{int(gates.amplitude_marked_but_accepted.sum())}。</p>'
    sensitivity=dest/'phase_tolerance_sensitivity_summary.csv'
    if sensitivity.exists():
        body+='<h2>固定相位容差敏感性检查</h2>'+pd.read_csv(sensitivity).to_html(index=False,float_format=lambda x:f'{x:.6f}')
    for name in ('phase_distributions.png','retention_by_frequency.png','accumulation.png'):
        body+='<img src="data:image/png;base64,'+base64.b64encode((dest/name).read_bytes()).decode()+'">'
    body+='<p>phase_stack_summary.csv：逐曲线逐频结果；phase_gate_summary.csv：去除周期重复的结果；windows/*.npz：全部窗口互谱、UTC起点索引、筛选与幅值标记；aligned_waveforms.h5：响应校正且保留缺测掩码的共同UTC波形。MATLAB可直接读取该H5及selected_phase_targets.csv。</p></html>'
    (dest/'phase_report.html').write_text(body,encoding='utf-8')
