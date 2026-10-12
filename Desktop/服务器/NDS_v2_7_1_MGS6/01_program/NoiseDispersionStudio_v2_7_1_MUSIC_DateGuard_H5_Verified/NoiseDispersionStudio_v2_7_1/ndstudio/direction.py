"""Shared, signed MUSIC direction and phase projection on shared UTC windows.

Angles are propagation azimuths: north=0, clockwise. The high-score envelope
is a resolution diagnostic, NOT a calibrated statistical confidence interval.
No input dispersion curve, CCF source phase or cycle enters the MUSIC estimate.
"""
from pathlib import Path
import hashlib, json, itertools
import numpy as np
import pandas as pd
from numba import njit
from dispersion.common import csv, write_json, sha

DEFAULTS=dict(mode='filter', max_pair_angle_deg=45., direction_step_deg=5., velocity_min_m_s=150.,
    velocity_max_m_s=3500., velocity_points=32, min_stations=3, min_geometry_ratio=.05,
    min_coherence=.7, envelope_fraction=.95, max_direction_uncertainty_deg=20.,
    max_secondary_peak_ratio=.95, direction_group_width_deg=15., conservative_angle=True,
    max_angle_velocity_error_fraction=.15, frequency_half_width_hz=0., frequency_bins=3,
    phase_rms_max_deg=25., minimum_independent_windows=3, minimum_accepted_windows=3, method='music', music_snapshots=6, music_time_bandwidth=3.5,
    music_sources=0, music_max_sources=2, music_diagonal_loading=1e-5, music_min_source_fraction=.1, music_front_cycles=20.)

REASONS={'accepted':'通过','missing':'缺测或无定义频谱','few_stations':'有效台站不足',
    'geometry':'子阵列近直线，二维方向不可辨','low_coherence':'阵列相干度不足',
    'boundary':'最优速度在搜索边界','ambiguous':'存在近等强的分离峰',
    'uncertain':'方向高分候选范围过宽','angle':'台对夹角超过阈值',
    'angle_boundary':'角度候选范围越过阈值','angle_error':'角度造成的速度误差预算超限',
    'other_direction':'属于另一传播方向组','grid_resolution':'阵列孔径与频率要求更细扫描网格',
    'no_signal':'MUSIC未识别出信号子空间','multiple_sources':'多波源窗不能唯一修正单一方向和速度'}


def settings(custom=None):
    supplied=custom or {}; extra=set(supplied)-set(DEFAULTS)
    if extra: raise ValueError('未知方向参数: '+str(sorted(extra)))
    c=dict(DEFAULTS)|supplied
    if c['method']!='music':raise ValueError('仅支持MUSIC定向；旧定向方法已停用')
    for key,lo,hi in [('music_snapshots',3,16),('music_sources',0,8),('music_max_sources',1,8)]:
        if isinstance(c[key],bool) or int(c[key])!=c[key] or not lo<=c[key]<=hi:raise ValueError(key+'整数范围无效')
        c[key]=int(c[key])
    if not 2<=c['music_front_cycles']<=10000:raise ValueError('前段MUSIC周期数须在2–10000之间')
    if not 0<c['music_min_source_fraction']<.5:raise ValueError('MUSIC次波源最小能量比例须在0–0.5之间')
    if c['music_sources']>c['music_max_sources']:raise ValueError('MUSIC信号数不能大于上限')
    if not .5<c['music_time_bandwidth']<10 or not 0<c['music_diagonal_loading']<.1:raise ValueError('MUSIC快照带宽/加载量无效')
    if c['music_snapshots']>2*c['music_time_bandwidth']:raise ValueError('DPSS快照数须不超过2NW')
    if c['mode'] not in ('off','diagnostic','filter','joint'): raise ValueError('方向模式无效')
    numeric=[v for k,v in c.items() if k not in ('mode','conservative_angle','method')]
    if not np.isfinite(numeric).all(): raise ValueError('方向参数必须为有限数值')
    if not 0<c['max_pair_angle_deg']<90: raise ValueError('最大台对夹角须在0–90°之间')
    if not 0<c['velocity_min_m_s']<c['velocity_max_m_s']: raise ValueError('方向速度范围无效')
    for k,lo,hi in [('velocity_points',3,256),('min_stations',3,100),('frequency_bins',1,15),
                    ('minimum_independent_windows',1,100000),('minimum_accepted_windows',1,100000)]:
        if isinstance(c[k],bool) or int(c[k])!=c[k] or not lo<=c[k]<=hi: raise ValueError(k+'整数范围无效')
        c[k]=int(c[k])
    if c['frequency_bins']%2!=1: raise ValueError('邻近频点数量须为奇数')
    if not 1<=c['direction_step_deg']<=30 or not np.isclose(360/c['direction_step_deg'],round(360/c['direction_step_deg'])):
        raise ValueError('方向步长须在1–30°且整除360°')
    if not 1<=c['direction_group_width_deg']<=90 or not np.isclose(360/c['direction_group_width_deg'],round(360/c['direction_group_width_deg'])):
        raise ValueError('方向分组宽度须在1–90°且整除360°')
    for k in ('min_geometry_ratio','min_coherence','envelope_fraction','max_secondary_peak_ratio'):
        if not 0<c[k]<1: raise ValueError(k+'须在0–1之间')
    if not 0<c['max_direction_uncertainty_deg']<90 or not 0<c['phase_rms_max_deg']<180: raise ValueError('方向或相位误差范围无效')
    if c['frequency_half_width_hz']<0 or c['max_angle_velocity_error_fraction']<=0: raise ValueError('频带半宽/角度误差预算无效')
    if not isinstance(c['conservative_angle'],bool): raise ValueError('conservative_angle须为布尔值')
    return c


def angular_distance(a,b): return abs((np.asarray(a)-np.asarray(b)+180)%360-180)


def geometry_ratio(xy):
    if len(xy)<3: return 0.
    sv=np.linalg.svd(xy-xy.mean(axis=0),compute_uv=False)
    return float(sv[1]/sv[0]) if len(sv)>1 and sv[0]>0 else 0.


def grid(c):
    # Uniform inverse velocity is a slowness grid, avoiding uneven resolution.
    slow=np.linspace(1/c['velocity_max_m_s'],1/c['velocity_min_m_s'],c['velocity_points'])
    az=np.arange(0.,360.,c['direction_step_deg']); a=np.radians(az)
    vectors=np.stack([slow[:,None]*np.sin(a),slow[:,None]*np.cos(a)],axis=-1).reshape(-1,2)
    return vectors,az,slow


@njit(cache=True)
def peak_diagnostics(power,az,slow,theta,velocity,best,envelope_fraction,separation_angle):
    """Fuse grid-envelope and eight-neighbor peak checks without grid-sized temporaries."""
    nw,ng=power.shape;na=len(az);nv=len(slow)
    uncertainty=np.zeros(nw);vmin=np.full(nw,np.inf);vmax=np.full(nw,-np.inf);secondary=np.zeros(nw)
    for w in range(nw):
        floor=best[w]*envelope_fraction
        for vi in range(nv):
            cv=1/slow[vi];separated_velocity=abs(np.log(cv/velocity[w]))>.25
            for ai in range(na):
                k=vi*na+ai;value=power[w,k];angle=abs((az[ai]-theta[w]+180)%360-180)
                if value>=floor:
                    uncertainty[w]=max(uncertainty[w],angle);vmin[w]=min(vmin[w],cv);vmax[w]=max(vmax[w],cv)
                if value<=secondary[w] or not (separated_velocity or angle>separation_angle):continue
                peak=True
                for dv in range(-1,2):
                    vj=vi+dv
                    if vj<0 or vj>=nv:continue
                    for da in range(-1,2):
                        if not dv and not da:continue
                        aj=(ai+da)%na
                        if value<power[w,vj*na+aj]:peak=False;break
                    if not peak:break
                if peak:secondary[w]=value
        secondary[w]/=max(best[w],1e-15)
    return uncertainty,vmin,vmax,secondary


def estimate(coefficients, xy, frequencies, c):
    from .music import estimate as music_estimate
    return music_estimate(coefficients,xy,frequencies,c)


def pair_gate(result, dx, dy, c):
    az=float(np.degrees(np.arctan2(dx,dy))%360)
    angle=angular_distance(result['azimuth_deg'],az); angle=np.minimum(angle,180-angle)
    upper=np.minimum(90.,angle+result['direction_uncertainty_deg'])
    lower=np.maximum(0.,angle-result['direction_uncertainty_deg'])
    denominator=np.cos(np.radians(angle))
    with np.errstate(divide='ignore',invalid='ignore'):
        error=np.maximum(abs(np.cos(np.radians(lower))/denominator-1),abs(np.cos(np.radians(upper))/denominator-1))
    reason=result['reason'].copy(); good=result['quality'].copy()
    reason[good&(angle>c['max_pair_angle_deg']+1e-9)]='angle'
    if c['conservative_angle']: reason[good&(angle<=c['max_pair_angle_deg']+1e-9)&(upper>c['max_pair_angle_deg']+1e-9)]='angle_boundary'
    reason[(reason=='accepted')&(error>c['max_angle_velocity_error_fraction'])]='angle_error'
    accepted=good&(reason=='accepted')
    return dict(direction_valid=accepted,direction_reason=reason,pair_angle_deg=angle,pair_angle_upper_deg=upper,
        angle_velocity_error_fraction=error,projected_distance_m=(dx*result['slow_e_s_m']+dy*result['slow_n_s_m'])*result['velocity_m_s'])


def independent_indices(starts,n,mask):
    selected=[]; end=-1
    for i in np.flatnonzero(mask):
        if starts[i]>=end: selected.append(i); end=int(starts[i])+n
    return np.asarray(selected,dtype=int)


def joint_estimate(coeff, xy, starts, n, f, result, c, amplitude_rejection=True):
    """Resolve pair phase wraps around independent MUSIC candidates; no curve gate."""
    from .phase_screening import amplitude_flags
    pairs=list(itertools.combinations(range(len(xy)),2)); d=np.array([xy[j]-xy[i] for i,j in pairs]); nw=len(starts)
    if not pairs or not nw:return [],[]
    cross=np.array([coeff[j]*np.conj(coeff[i]) for i,j in pairs]); eligible=[]
    for k,(dx,dy) in enumerate(d):
        m=pair_gate(result,dx,dy,c)['direction_valid']
        valid=np.isfinite(cross[k])&(abs(cross[k])>0); m=m&valid
        if amplitude_rejection:
            large,jump,*_=amplitude_flags(abs(cross[k])); m&=~(large|jump)
        eligible.append(m)
    w=np.asarray(eligible,float); measured=np.angle(cross)
    prediction=-2*np.pi*f*(d@np.array([result['slow_e_s_m'],result['slow_n_s_m']]))
    cycles=np.rint((prediction-measured)/(2*np.pi)); lifted=measured+2*np.pi*cycles
    rhs=np.nan_to_num(-lifted/(2*np.pi*f));dx=d[:,0,None];dy=d[:,1,None]
    # Direction is ALWAYS MUSIC. Fit only scalar slowness along that axis;
    # no legacy two-dimensional phase-line direction fit is retained.
    axis=np.array([np.sin(np.radians(result['azimuth_deg'])),np.cos(np.radians(result['azimuth_deg']))])
    projection=d@axis
    denominator=(w*np.nan_to_num(projection)**2).sum(axis=0)
    with np.errstate(divide='ignore',invalid='ignore'):
        scalar=(w*np.nan_to_num(projection)*rhs).sum(axis=0)/denominator
        se=axis[0]*scalar;sn=axis[1]*scalar;v=1/scalar
        pred=-2*np.pi*f*projection*scalar;residual=np.angle(np.exp(1j*(pred-measured)))
    rms=np.degrees(np.sqrt((w*np.nan_to_num(residual)**2).sum(axis=0)/np.maximum(w.sum(axis=0),1)))
    valid=(w.sum(axis=0)>=3)&(denominator>0)&np.isfinite(v)&(v>c['velocity_min_m_s'])&(v<c['velocity_max_m_s'])&(rms<=c['phase_rms_max_deg'])&result['quality']
    rows=[]; pair_rows=[]
    for group in np.unique(result['direction_group'][valid]):
        mask=valid&(result['direction_group']==group); ix=independent_indices(starts,n,mask)
        if not len(ix):continue
        mean_axis=np.mean(axis[:,ix],axis=1);mean_axis/=np.linalg.norm(mean_axis)
        reference=mean_axis*np.median(scalar[ix]);accum=[];design=[]
        for k,(i,j) in enumerate(pairs):
            selected=ix[eligible[k][ix]]
            if not len(selected):continue
            z=cross[k,selected]; window_pred=-2*np.pi*f*(d[k,0]*se[selected]+d[k,1]*sn[selected])
            ref_pred=-2*np.pi*f*(d[k]@reference)
            aligned=z/abs(z)*np.exp(1j*(ref_pred-window_pred)); total=aligned.sum(); phase=float(np.angle(total))
            cycle=int(np.rint((ref_pred-phase)/(2*np.pi)))
            accum.append(-(phase+2*np.pi*cycle)/(2*np.pi*f));design.append(d[k])
            pair_rows.append(dict(direction_group=int(group),i=i,j=j,independent_windows=len(selected),
                reference_projected_distance_m=float(d[k]@reference/np.linalg.norm(reference)),
                final_phase_deg=np.degrees(phase),phase_cycle=cycle,unit_resultant_length=abs(total)/len(selected)))
        if len(design)<3 or geometry_ratio(np.vstack([np.zeros(2),design]))<c['min_geometry_ratio']:continue
        projected=np.asarray(design)@mean_axis
        final_scalar=float(projected@np.asarray(accum)/(projected@projected))
        if not np.isfinite(final_scalar) or final_scalar<=0:continue
        final=mean_axis*final_scalar;velocity=1/final_scalar
        accepted=len(ix)>=c['minimum_independent_windows'] and c['velocity_min_m_s']<velocity<c['velocity_max_m_s']
        rows.append(dict(direction_method='music',direction_group=int(group),propagation_azimuth_deg=float(np.degrees(np.arctan2(final[0],final[1]))%360),
            phase_velocity_m_s=float(velocity),slow_e_s_m=float(final[0]),slow_n_s_m=float(final[1]),
            velocity_temporal_q025_m_s=float(np.quantile(v[ix],.025)),velocity_temporal_q975_m_s=float(np.quantile(v[ix],.975)),
            overlapping_windows=int(mask.sum()),independent_windows=len(ix),joint_pairs=len(design),
            median_phase_rms_deg=float(np.median(rms[ix])),accepted=bool(accepted),
            reason='accepted' if accepted else 'independent_windows_or_velocity',
            uncertainty_kind='nonoverlapping-window temporal spread, not a confidence interval'))
    return rows,pair_rows


class DirectionSession:
    def __init__(self, arrays, xy_table, fs, origin, station_limits, dest, workdir, c, stabilization_s, trace, amplitude_rejection):
        self.arrays=arrays; self.names=sorted(arrays); self.xy=np.array([xy_table.loc[s,['x_m','y_m']].to_numpy(float) for s in self.names])
        self.fs=fs;self.origin=origin;self.limits=station_limits;self.dest=Path(dest)/'direction';self.dest.mkdir(exist_ok=True)
        self.work=Path(workdir);self.c=c;self.stable=int(np.floor(stabilization_s*fs+.5));self.trace=trace
        self.cache={};self.summaries=[];self.joint=[];self.joint_pairs=[];self.examples=[];self.amp=amplitude_rejection
        csv(self.dest/'station_geometry.csv',pd.DataFrame(dict(station=self.names,x_m=self.xy[:,0],y_m=self.xy[:,1])))
        excluded={'mode','max_pair_angle_deg','direction_group_width_deg','conservative_angle','max_angle_velocity_error_fraction','minimum_independent_windows','minimum_accepted_windows'}
        compute={k:v for k,v in c.items() if k not in excluded}
        fingerprint=dict(settings=compute,waveforms=sha(Path(dest)/'aligned_waveforms.h5'),xy=self.xy.tolist(),stations=self.names,fs=fs,origin=str(origin),limits=station_limits,stable=self.stable,code=sha(__file__),music_code=sha(Path(__file__).with_name('music.py')),goertzel_code=sha(Path(__file__).with_name('phase_stack.py')))
        digest=hashlib.sha256(json.dumps(fingerprint,sort_keys=True).encode()).hexdigest()
        self.persistent=Path(dest)/'direction_cache'/digest;self.persistent.mkdir(parents=True,exist_ok=True)

    def get(self, nominal_f, actual_f, n, step, taper):
        from .phase_stack import goertzel_windows,save_windows
        from .phase_records import full_window_starts
        taper_sha=hashlib.sha256(np.ascontiguousarray(taper).tobytes()).hexdigest()
        key=(float(nominal_f),float(actual_f),n,step,taper_sha)
        if key in self.cache:
            with np.load(self.cache[key],allow_pickle=False) as z:
                data={k:z[k] for k in z.files}
            data['direction_group']=self.groups(data['azimuth_deg'])
            return data
        token=hashlib.sha256(repr(key).encode()).hexdigest()[:16];path=self.persistent/(token+'.npz')
        duration=len(next(iter(self.arrays.values())))-self.stable
        starts=full_window_starts(duration,n,step,include_last=True)+self.stable if duration>=n else np.empty(0,np.int64)
        half=self.c['frequency_half_width_hz'];bins=self.c['frequency_bins'] if half>0 else 1
        freqs=actual_f+np.linspace(-half,half,bins) if bins>1 else np.array([actual_f])
        if np.any(freqs<=0) or np.any(freqs>=self.fs/2):raise ValueError('MUSIC邻近频点超出0–Nyquist范围')
        reused=path.exists()
        with self.trace.stage(f'方向估计 {nominal_f:g} Hz',dict(actual_frequency_hz=actual_f,windows=len(starts),stations=self.names,
            window_s=n/self.fs,hop_s=step/self.fs,frequency_bins_hz=freqs.tolist(),cache_reused=reused)) as log:
            if reused:
                with np.load(path,allow_pickle=False) as z:data={k:z[k] for k in z.files}
                if not np.array_equal(starts,data.pop('start_sample')):raise ValueError('方向缓存的UTC窗口不一致')
                central=data.pop('coefficients');result=data;result['direction_group']=self.groups(result['azimuth_deg'])
                examples=[]
                # Bounded diagnostic recomputation; the full scan is reused.
                if len(self.examples)<3 and len(starts):
                    ix=starts[:1];preview=self.transform(ix,n,taper,freqs,snapshots=True)
                    _,examples=estimate(preview,self.xy,freqs,self.c)
            else:
                coefficients=self.transform(starts,n,taper,freqs,snapshots=True)
                result,examples=estimate(coefficients,self.xy,freqs,self.c)
                central=self.transform(starts,n,taper,[actual_f])[0]
                save_windows(path,start_sample=starts,coefficients=central,**result)
            summary=dict(method='music',music_snapshots=self.c['music_snapshots'],snapshot_half_bandwidth_hz=self.c['music_time_bandwidth']*self.fs/n,frequency_hz=nominal_f,actual_frequency_hz=actual_f,total_windows=len(starts),quality_windows=int(result['quality'].sum()),cache_reused=reused,
                median_velocity_m_s=float(np.median(result['velocity_m_s'][result['quality']])) if result['quality'].any() else None,
                reasons={REASONS[k]:int((result['reason']==k).sum()) for k in np.unique(result['reason'])})
            log.update(summary);self.summaries.append(summary)
            for e in examples:
                if len(self.examples)<3:self.examples.append(dict(f=actual_f,**e))
        rows,prs=joint_estimate(central,self.xy,starts,n,actual_f,result,self.c,self.amp)
        for row in rows:self.joint.append(dict(frequency_hz=nominal_f,actual_frequency_hz=actual_f,**row))
        for row in prs:self.joint_pairs.append(dict(frequency_hz=nominal_f,actual_frequency_hz=actual_f,pair=self.names[row.pop('i')]+'--'+self.names[row.pop('j')],**row))
        save_windows(self.dest/(token+'.npz'),start_sample=starts,frequency_hz=nominal_f,actual_frequency_hz=actual_f,
            window_samples=n,slide_samples=step,fs_hz=self.fs,start_utc=str(self.origin),**result)
        self.cache[key]=path
        return dict(start_sample=starts,coefficients=central,**result)

    def groups(self, azimuth):
        width=self.c['direction_group_width_deg']
        return np.floor((np.nan_to_num(azimuth)+width/2)/width).astype(int)%round(360/width)

    def transform(self,starts,n,taper,freqs,snapshots=False):
        from .phase_stack import goertzel_windows
        from scipy.signal.windows import dpss
        tapers=[taper]
        if snapshots:
            if n<=2*self.c['music_time_bandwidth']:raise ValueError('MUSIC窗口太短，无法构造DPSS快照')
            tapers=dpss(n,self.c['music_time_bandwidth'],Kmax=self.c['music_snapshots'],sym=False)
        coefficients=np.empty((len(freqs),len(self.names),len(starts),len(tapers)),complex)
        for si,name in enumerate(self.names):
            first,last=self.limits[name];good=(starts>=max(0,first+self.stable))&(starts+n<=min(len(self.arrays[name]),last))
            for fi,f in enumerate(freqs):
                for ti,wt in enumerate(tapers):
                    for lo in range(0,len(starts),4096):
                        ix=starts[lo:lo+4096];a=int(ix[0]);b=int(ix[-1])+n
                        coefficients[fi,si,lo:lo+len(ix),ti]=goertzel_windows(self.arrays[name][a:b],ix-a,n,float(f),self.fs,wt)
                coefficients[fi,si,~good]=np.nan+1j*np.nan
        return coefficients if snapshots else coefficients[...,0]

    def finish(self, frequencies, scan):
        from .phase_bands import assign_bands
        columns=['frequency_hz','actual_frequency_hz','direction_group','phase_velocity_m_s','accepted','independent_windows']
        q=csv(self.dest/'joint_dispersion.csv',self.joint,columns=None if self.joint else columns)
        observed=q[q.frequency_hz.round(9).isin(np.round(frequencies,9))] if len(q) else q
        bands=assign_bands(observed,frequencies,scan,group_key='direction_group') if len(observed) else pd.DataFrame(columns=['band_id','accepted'])
        csv(self.dest/'joint_continuous_bands.csv',bands)
        csv(self.dest/'joint_pair_phases.csv',self.joint_pairs,columns=None if self.joint_pairs else ['pair','frequency_hz','final_phase_deg'])
        flat=[dict(**{k:v for k,v in r.items() if k!='reasons'},reason=name,windows=count) for r in self.summaries for name,count in r['reasons'].items()]
        csv(self.dest/'frequency_diagnostics.csv',flat,columns=None if flat else ['frequency_hz','reason','windows'])
        write_json(self.dest/'direction_complete.json',dict(method='music',music_module_sha256=sha(Path(__file__).with_name('music.py')),settings=self.c,frequencies=len(self.summaries),
            raw_joint_points=len(q),accepted_continuous_points=int(bands.accepted.sum()) if len(bands) else 0,
            module_sha256=sha(__file__),direction_convention='propagation: north=0 clockwise',
            confidence_note='High-score grid envelope and temporal spread are not calibrated confidence intervals.'))
        self.plot(q)
        import html,base64
        doc='<html><meta charset="utf-8"><style>body{font:15px/1.6 sans-serif;max-width:1200px;margin:30px auto}img{max-width:100%}td,th{padding:5px;border-bottom:1px solid #ccc}</style><h1>MUSIC阵列方向与平面波联合估计</h1>'
        doc+='<p>方向估计使用曲线筛选之前的共同UTC窗口。北为0°、顺时针；传播方向与来波方向相差180°。高分候选范围是分辨率诊断，不是统计置信区间。联合速度不使用CCF源相位或输入曲线；整周数在通过质量检查的MUSIC候选附近解析，仍需检查空间混叠与平面波假设。</p>'
        doc+='<h2>实际参数</h2><pre>'+html.escape(json.dumps(self.c,ensure_ascii=False,indent=2))+'</pre>'
        doc+='<h2>逐频剔除原因</h2>'+pd.DataFrame(flat).to_html(index=False)
        doc+='<h2>联合频散与连续带宽检查</h2>'+bands.to_html(index=False)
        for name in ['direction_overview.png','beam_maps.png']:
            if (self.dest/name).exists():doc+='<img src="data:image/png;base64,'+base64.b64encode((self.dest/name).read_bytes()).decode()+'">'
        doc+='<p><a href="../phase_report.html">相位筛选报告</a> · <a href="../phase_trace.json">全阶段追溯</a></p></html>'
        (self.dest/'report.html').write_text(doc,encoding='utf-8')
        self.trace.note('方向汇总',f'{len(self.summaries)}个频率；联合原始点 {len(q)}；连续频带保留点 {int(bands.accepted.sum()) if len(bands) else 0}',output=str(self.dest/'report.html'))

    def plot(self,q):
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(2,1,figsize=(11,7),constrained_layout=True)
        for path in sorted(self.dest.glob('*.npz')):
            with np.load(path) as z:
                pick=np.arange(0,len(z['start_sample']),max(1,len(z['start_sample'])//150))
                good=z['quality'][pick];f=float(z['frequency_hz'])
                axes[0].scatter(np.full(good.sum(),f),z['azimuth_deg'][pick][good],s=6,alpha=.3,color='#177c8c')
        axes[0].set(xlabel='Frequency (Hz)',ylabel='Propagation azimuth (deg)',ylim=(0,360),title='Quality-qualified windows before curve screening')
        if len(q):
            for g,r in q.groupby('direction_group'):
                r=r.sort_values('frequency_hz');axes[1].plot(r.frequency_hz,r.phase_velocity_m_s,'.-',label=f'direction group {g}')
        axes[1].set(xlabel='Frequency (Hz)',ylabel='Joint horizontal phase velocity (m/s)');
        if len(q):axes[1].legend()
        fig.savefig(self.dest/'direction_overview.png',dpi=140);plt.close(fig)
        if self.examples:
            fig,axes=plt.subplots(2,len(self.examples),figsize=(5*len(self.examples),8),squeeze=False,constrained_layout=True)
            for i,e in enumerate(self.examples):
                ax=axes[0,i]
                p=ax.pcolormesh(e['azimuth'],e['velocity'],e['power'],shading='auto',vmin=0,vmax=1)
                ax.set(xlabel='Propagation azimuth (deg)',ylabel='Velocity (m/s)',title=f"{e['f']:.4g} Hz | window {e['window_index']}");fig.colorbar(p,ax=ax,label='Normalised MUSIC pseudospectrum')
                vectors,_,_=grid(self.c); coords=e['coordinates'];offset=vectors-e['reference_slowness']
                response=sum(abs(np.exp(2j*np.pi*f*(offset@coords.T)).sum(axis=1))**2/len(coords)**2 for f in e['frequencies'])/len(e['frequencies'])
                ax=axes[1,i];p=ax.pcolormesh(e['azimuth'],e['velocity'],response.reshape(e['power'].shape),shading='auto',vmin=0,vmax=1)
                ax.set(xlabel='Propagation azimuth (deg)',ylabel='Velocity (m/s)',title='Array response around peak candidate');fig.colorbar(p,ax=ax)
            fig.savefig(self.dest/'beam_maps.png',dpi=140);plt.close(fig)
