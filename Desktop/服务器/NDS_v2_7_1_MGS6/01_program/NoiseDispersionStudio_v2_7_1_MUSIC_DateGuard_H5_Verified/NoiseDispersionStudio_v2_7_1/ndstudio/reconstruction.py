"""Prior-guided low-frequency measurements using signed wave-vector projections."""
from pathlib import Path
import hashlib,json,shutil
import numpy as np
import pandas as pd
from scipy import signal
from dispersion.common import csv,sha,write_json

PROJECTION_NOTE='若最终要输出沿波矢的真实相速度，还必须一致地处理台距投影，不能直接把台对视速度当作真实波速。'
DEFAULTS=dict(enabled=False,lower_frequency_hz=.2,step_hz=.05,minimum_direction_frequencies=2,null_trials=40,null_p_max=.05)

def settings(value=None):
    value=value or {}
    if set(value)-set(DEFAULTS):raise ValueError('未知低频重建参数: '+str(sorted(set(value)-set(DEFAULTS))))
    c=DEFAULTS|value
    if not isinstance(c['enabled'],bool):raise ValueError('reconstruction.enabled must be boolean')
    for k in ('lower_frequency_hz','step_hz','null_p_max'):
        if not np.isfinite(c[k]) or c[k]<=0:raise ValueError(k+' must be finite and positive')
    for k in ('minimum_direction_frequencies','null_trials'):
        if int(c[k])!=c[k] or c[k]<1:raise ValueError(k+' must be a positive integer')
        c[k]=int(c[k])
    if c['minimum_direction_frequencies']<2 or c['null_p_max']>.05 or c['null_trials']<19:
        raise ValueError('方向至少跨2个低频点；null_trials>=19，null_p_max<=0.05')
    return c

def low_targets(observed,options,config):
    out=[]
    for cid,g in observed.groupby('curve_id',sort=True):
        g=g.sort_values('frequency_hz').copy()
        if len(g)<2 or np.any(np.diff(g.phase_velocity_m_s)>1e-8):continue
        low=float(options['lower_frequency_hz']);step=float(options['step_hz'])
        if low<=config['bandpass_hz'][0] or low<=config['direction_half_width']:
            raise ValueError('低频重建下限必须高于预处理通带下限及MUSIC频带半宽；先修改预处理参数')
        use=g.head(min(4,len(g)));slope=float(np.clip(np.polyfit(np.log(use.frequency_hz),np.log(use.phase_velocity_m_s),1)[0],-3,0))
        count=max(0,int(np.floor((g.frequency_hz.iloc[0]-low+1e-9)/step)))
        for f in g.frequency_hz.iloc[0]-step*np.arange(count,0,-1):
            r=g.iloc[0].to_dict();v=float(np.clip(g.phase_velocity_m_s.iloc[0]*(f/g.frequency_hz.iloc[0])**slope,*config['velocity_bounds']))
            r.update(frequency_hz=float(round(f,9)),phase_velocity_m_s=v,is_model_prediction=True,is_prior_prediction=True,background_discriminated=False);out.append(r)
        g['is_prior_prediction']=False;out.extend(g.to_dict('records'))
    if not out:raise ValueError('没有至少两个频点且速度不随频率增加的输入曲线')
    return pd.DataFrame(out)

def align_cross(raw_cross,frequency,projected_distance,reference_distance,prior_true_velocity):
    rotation=-2*np.pi*frequency*(reference_distance-np.asarray(projected_distance))/prior_true_velocity
    return raw_cross*np.exp(1j*rotation),rotation

def phase_speed(frequency,phase,reference_distance):return -2*np.pi*np.asarray(frequency)*reference_distance/np.asarray(phase)
def derived_group(f,v):return 1/np.gradient(np.asarray(f)/np.asarray(v),np.asarray(f),edge_order=1)

def shifted_null(ga,gb,starts,n,eligible,projection,reference,f,v,center,lo,hi,options,amp):
    from .direction import independent_indices
    from .phase_screening import screen_cross
    ids=independent_indices(starts,n,eligible)
    if len(ids)<3:return 1.,0.,0
    def statistic(a,b):
        cross,_=align_cross(b*np.conj(a),f,projection[ids],reference,v)
        m=screen_cross(cross,center,lo,hi,'curve',amp,False)
        return float(m['accepted'].sum()/max(1,m['amplitude_valid'].sum()))
    actual=statistic(ga[ids],gb[ids]);rng=np.random.default_rng(9075)
    offsets=rng.permutation(np.arange(1,len(ids)))[:options['null_trials']]
    null=[statistic(ga[ids],np.roll(gb[ids],int(k))) for k in offsets]
    return float((1+sum(x>=actual-1e-12 for x in null))/(1+len(null))),actual,len(null)

def run(root,dest,targets,session,xy,config,options,repeat_cycles,slide_cycles,taper,amp,trace,meta):
    from .phase_records import phase_window_geometry
    from .phase_stack import save_windows,stack_stats
    from .phase_screening import screen_cross
    from .direction import pair_gate,independent_indices,angular_distance
    from .phase_bands import assign_bands
    from .phase_velocity import lift_band
    observed=targets[~targets.is_prior_prediction.astype(bool)].copy()
    for name in ('report.html','comparison.png','comparison_complete.json'):(dest/'paired_comparison'/name).unlink(missing_ok=True)
    targets=low_targets(observed,options,dict(bandpass_hz=config['bandpass_hz'],direction_half_width=session.c['frequency_half_width_hz'],velocity_bounds=config['phase']['velocity_bounds_m_s']))
    csv(dest/'selected_phase_targets.csv',targets);csv(dest/'screening_prior_predictions.csv',targets[targets.is_prior_prediction])
    frequencies=np.sort(targets.frequency_hz.unique());fs=session.fs
    scan=dict(config['scan'],step_hz=max(options['step_hz'],config['scan']['step_hz']))
    windows=dest/'windows';windows.mkdir(exist_ok=True)
    for p in windows.glob('*.npz'):p.unlink()
    def fetch(pair,f):
        geometry=phase_window_geometry(float(f),fs,repeat_cycles,slide_cycles);n=geometry['window_samples'];step=geometry['slide_samples']
        wt=np.ones(n) if taper=='none' else signal.windows.hann(n,sym=True);d=session.get(float(f),float(f),n,step,wt)
        a,b=pair.split('--');delta=xy.loc[b,['x_m','y_m']].to_numpy(float)-xy.loc[a,['x_m','y_m']].to_numpy(float)
        gate=pair_gate(d,*delta,session.c)
        first=max(session.limits[a][0],session.limits[b][0])+session.stable;last=min(session.limits[a][1],session.limits[b][1])
        gate['direction_valid'] &= (d['start_sample']>=first)&(d['start_sample']+n<=last)
        return d,gate,d['coefficients'][session.names.index(a)],d['coefficients'][session.names.index(b)],n,step,delta
    references=[];rows=[];before=[]
    for (pair,branch),g in targets.groupby(['pair','branch'],sort=True):
        for f,q in g.groupby('frequency_hz',sort=True):
            d,gate,ga,gb,n,step,delta=fetch(pair,f);starts=d['start_sample'];raw=gb*np.conj(ga)
            from .music_front import frequency_reference
            ref=frequency_reference(d,starts,n,delta,session.c,branch,angle_gate=True)
            if ref and not np.allclose(q.distance_m,ref['distance_m'],rtol=.001,atol=.01):raise ValueError('CCF台对距离与MUSIC坐标不一致: '+pair)
            references.append(dict(pair=pair,branch=branch,frequency_hz=float(f),confirmed=bool(ref),**({k:v for k,v in ref.items() if k not in ('indices','eligible')} if ref else {})))
            for _,r in q.iterrows():
                v=(float(r.phase_velocity_m_s) if r.get('projection_applied',False) else r.phase_velocity_m_s*abs(ref['projected_distance_m'])/r.distance_m) if ref else np.nan
                reference=ref['projected_distance_m'] if ref else np.nan;center=-2*np.pi*f*reference/v if ref else 0.
                eligible=gate['direction_valid'].copy()
                if ref:eligible &= (gate['projected_distance_m']*reference>0)&(angular_distance(d['azimuth_deg'],ref['azimuth_deg'])<=session.c['direction_group_width_deg']+1e-9)
                else:eligible[:]=False
                aligned,rotation=align_cross(raw,f,gate['projected_distance_m'],reference,v) if ref else (raw,np.zeros(len(raw)))
                lo=np.radians(r.phase_offset_low_deg);hi=np.radians(r.phase_offset_high_deg)
                if not 0<hi-lo<np.pi:raise ValueError('投影重建门限宽度须小于180°且大于0°')
                masks=screen_cross(aligned,center,lo,hi,'curve',amp,False,eligible);accepted=masks['accepted'];independent=len(independent_indices(starts,n,accepted))
                p,stat,trials=shifted_null(ga,gb,starts,n,eligible,gate['projected_distance_m'],reference,f,v,center,lo,hi,options,amp) if ref else (1.,0.,0)
                qualified=bool(ref and accepted.sum()>=session.c['minimum_accepted_windows'] and independent>=session.c['minimum_independent_windows'] and p<=options['null_p_max'])
                token=hashlib.sha256(f'{r.curve_id}|{f:.9f}'.encode()).hexdigest()[:16]
                save_windows(windows/(token+'.npz'),start_sample=starts,cross=aligned,raw_cross=raw,alignment_rotation_rad=rotation,window_projected_distance_m=gate['projected_distance_m'],reference_projected_distance_m=reference,prior_true_velocity_m_s=v,propagation_azimuth_deg=d['azimuth_deg'],confirmed_azimuth_deg=ref['azimuth_deg'] if ref else np.nan,**masks,**gate,direction_eligible=eligible,direction_applied=True,frequency_hz=f,actual_frequency_hz=f,window_samples=n,slide_samples=step,fs_hz=fs,start_utc=str(session.origin),pair=pair,predicted_phase_rad=center,phase_offset_low_rad=lo,phase_offset_high_rad=hi,screening_mode='curve',split_branches=False,amplitude_rejection=amp)
                row=r.to_dict();row['ccf_curve_id']=r.get('parent_curve_id',r.curve_id);row.update(gate_id=token,windows_file='windows/'+token+'.npz',accepted=qualified,accepted_windows=int(accepted.sum()),valid_windows=int(masks['valid'].sum()),total_windows=len(starts),direction_independent_windows=independent,selected_direction_group=ref['selected_direction_group'] if ref else -1,projected_distance_m=reference,confirmed_azimuth_deg=ref['azimuth_deg'] if ref else np.nan,predicted_unwrapped_phase_rad=center,predicted_phase_deg=float(np.degrees(np.angle(np.exp(1j*center)))),prior_true_velocity_m_s=v,background_null_p=p,background_null_trials=trials,background_observed_fraction=stat,window_duration_s=n/fs,rejection_reason=';'.join(reason for failed,reason in [(not bool(ref),'MUSIC方向未确认'),(int(accepted.sum())<session.c['minimum_accepted_windows'],'接受窗不足'),(independent<session.c['minimum_independent_windows'],'独立窗不足'),(trials<19,'时移背景样本不足19组'),(p>options['null_p_max'],'时移背景检验未通过')] if failed))
                row.update({'before_'+k:value for k,value in stack_stats(aligned,masks['valid']).items()});row.update({'after_'+k:value for k,value in stack_stats(aligned,accepted).items()});rows.append(row)
                if ref and not r.is_prior_prediction:
                    b=r.to_dict();b['is_model_prediction']=bool(r.get('origin_ccf_is_model_prediction',r.is_model_prediction));b['conditional_comparison_only']=r.get('target_role','')=='apparent_screening_prior_only';b.update(parent_curve_id=r.curve_id,phase_velocity_m_s=v,velocity_type='wave_vector_phase_velocity',projection_applied=True,projected_distance_m=reference,confirmed_azimuth_deg=ref['azimuth_deg']);before.append(b)
                trace.note('低频实测与投影叠加',f'{r.curve_id} · {f:g} Hz · {int(accepted.sum())}窗 · 独立{independent} · null p={p:g} · 合格={qualified}')
    csv(dest/'direction_reference.csv',references)
    frame=pd.DataFrame(rows);csv(dest/'phase_stack_summary.csv',frame);csv(dest/'phase_gate_summary.csv',frame)
    gate_targets=frame.copy();gate_targets['origin_source_phase_rad']=gate_targets.source_phase_rad;gate_targets['origin_constant_cycle']=gate_targets.constant_cycle
    gate_targets['source_phase_rad']=0.;gate_targets['phase_velocity_m_s']=gate_targets.prior_true_velocity_m_s;gate_targets['velocity_type']='wave_vector_phase_velocity';gate_targets['projection_applied']=True;gate_targets['is_model_prediction']=True;gate_targets['target_role']='screening_prior_not_observation';csv(dest/'selected_phase_targets.csv',gate_targets)
    checked=pd.concat([assign_bands(g,np.sort(g.frequency_hz.unique()),scan) for _,g in frame.groupby('curve_id',sort=False)]);results=[];catalog=[]
    for band,g in checked[checked.accepted].groupby('band_id',sort=True):
        g=g.sort_values('frequency_hz').copy()
        try:phase,shift=lift_band(np.radians(g.after_raw_phase_deg),g.predicted_unwrapped_phase_rad.to_numpy())
        except ValueError:checked.loc[g.index,['accepted','band_reason']]=[False,'全频带不能由同一整数解缠'];continue
        velocity=phase_speed(g.frequency_hz,phase,g.projected_distance_m.to_numpy())
        if not np.isfinite(velocity).all() or np.any(velocity<=0) or np.any(np.diff(velocity)>1e-8):checked.loc[g.index,['accepted','band_reason']]=[False,'速度非正/不连续/随频率增加'];continue
        g['parent_curve_id']=g.curve_id;g['curve_id']=band;g['origin_source_phase_rad']=g.source_phase_rad;g['origin_constant_cycle']=g.constant_cycle
        g['origin_CCF_unwrapped_phase_rad']=g.measured_unwrapped_phase_rad;g['measured_unwrapped_phase_rad']=phase
        for col in ('velocity_min_m_s','velocity_max_m_s','relative_phase_q025_rad','relative_phase_q975_rad','prediction_reconstruction_error_rad'):
            if col in g:g['origin_'+col]=g[col];g[col]=np.nan
        g['source_phase_rad']=0.;g['source_phase_verified']=False;g['raw_unwrapped_phase_rad']=phase;g['raw_phase_cycle']=np.rint((phase-np.radians(g.after_raw_phase_deg))/(2*np.pi)).astype(int)
        g['phase_velocity_m_s']=velocity;g['group_velocity_m_s']=derived_group(g.frequency_hz,velocity);g['apparent_pair_velocity_m_s']=velocity*g.distance_m/abs(g.projected_distance_m)
        g['is_low_frequency_reconstruction']=g.is_prior_prediction;g['is_prior_prediction']=False;g['is_model_prediction']=False;g['prior_used']=True;g['background_discriminated']=True;g['group_velocity_is_derived']=True;g['projection_applied']=True
        g['velocity_type']='wave_vector_phase_velocity';g['measurement_method']='direction_projected_goertzel';g['direction_method']='music';g['direction_quality_verified']=True;g['unique_physical_dispersion_verified']=False
        g['distance_over_wavelength']=abs(g.projected_distance_m)*g.frequency_hz/g.phase_velocity_m_s;g['cycle_interpretation']='origin_CCF_hypothesis; raw_phase_cycle is separately measured relative to prior'
        if not np.isfinite(g.group_velocity_m_s).all() or (g.group_velocity_m_s<=0).any():checked.loc[g.index,['accepted','band_reason']]=[False,'派生群速度非正'];continue
        results.append(g);catalog.append(dict(curve_id=band,segment_id=g.segment_id.iloc[0],cohort=g.cohort.iloc[0],pair_index=int(g.pair_index.iloc[0]),pair=g.pair.iloc[0],branch=g.branch.iloc[0],constant_cycle=int(g.constant_cycle.iloc[0]),source_phase_rad=0.,fmin_hz=float(g.frequency_hz.min()),fmax_hz=float(g.frequency_hz.max()),points=len(g),background_discriminated=True,compatible_other_pairs=0,parent_curve_id=g.parent_curve_id.iloc[0],direction_method='music',measurement_method='direction_projected_goertzel',direction_quality_verified=True))
    final=pd.concat(results,ignore_index=True) if results else pd.DataFrame(columns=list(dict.fromkeys(list(frame.columns)+['parent_curve_id','distance_over_wavelength','measurement_method'])))
    csv(dest/'inversion_curves.csv',final);csv(dest/'inversion_catalog.csv',catalog,columns=None if catalog else ['curve_id','segment_id','cohort','pair_index','pair','branch','constant_cycle','source_phase_rad','fmin_hz','fmax_hz','points','background_discriminated','compatible_other_pairs','parent_curve_id','measurement_method','direction_quality_verified'])
    csv(dest/'phase_continuous_bands.csv',checked);csv(dest/'qualified_phase_summary.csv',checked[checked.accepted])
    before=pd.DataFrame(before,columns=list(dict.fromkeys(list(observed.columns)+['parent_curve_id','projection_applied','projected_distance_m','confirmed_azimuth_deg'])))
    for _,g in before.groupby('curve_id'):
        g=g.sort_values('frequency_hz');before.loc[g.index,'group_velocity_m_s']=derived_group(g.frequency_hz,g.phase_velocity_m_s)
    before['group_velocity_is_derived']=True;csv(dest/'before_projected_curves.csv',before);session.finish(frequencies,scan)
    marker=dict(complete=True,direction_method='music',music_module_sha256=sha(Path(__file__).with_name('music.py')),reconstruction=options,reconstruction_module_sha256=sha(__file__),module_sha256=sha(Path(__file__).with_name('phase_stack.py')),direction_sha256=sha(Path(__file__).with_name('direction.py')),phase_records_sha256=sha(Path(__file__).with_name('phase_records.py')),phase_screening_sha256=sha(Path(__file__).with_name('phase_screening.py')),settings_sha256=sha(root/'config.json'),phase_settings_sha256=sha(root/'phase_settings.json') if (root/'phase_settings.json').exists() else None,aligned_h5_sha256=sha(dest/'aligned_waveforms.h5'),station_geometry_sha256=sha(dest/'direction/station_geometry.csv'),stations_csv_sha256=sha(root/'stations.csv'),direction=session.c,frequency_method='exact',screening_mode='curve',repeat_cycles=repeat_cycles,slide_cycles=slide_cycles,taper=taper,amplitude_rejection=amp,split_phase_branches=False,projection_applied=True,velocity_type='wave_vector_phase_velocity',raw_cross_source_phase_rad=0.,monotonic_only=True,prior=None,shared_utc_windows=True,stabilization_s=session.stable/fs,response_corrected=True,valid_window_evaluations=int(frame.valid_windows.sum()),accepted_window_evaluations=int(frame.accepted_windows.sum()),inversion_handoff_points=len(final),low_frequency_measured_points=int(final.get('is_low_frequency_reconstruction',pd.Series(dtype=bool)).sum()),selected_curves=int(targets.curve_id.nunique()),unique_pair_frequency_gates=len(frame),station_metadata=meta,warnings=['Prior-guided gates do not independently identify the absolute cycle or Rayleigh mode.','Derived group velocities are not independent observations.'])
    write_json(dest/'phase_stack_complete.json',marker);write_json(dest/'phase_progress.json',dict(complete=True,processed_gates=len(frame)))
    report='<meta charset="utf-8"><h1>低频实测重建与沿波矢相速度</h1><p>'+PROJECTION_NOTE+'</p><p>逐频点MUSIC确认方向；每频率逐窗旋转到该频点的投影距离后筛选和叠加。c = −2π f ℓref / φ。先验预测不能直接反演，原CCF源相位及周期只保留为先验来源。</p>'
    report+=f'<p>合格实测点 {len(final)}；其中新增低频 {marker["low_frequency_measured_points"]}。</p><p><a href="direction/report.html">方向报告</a> · <a href="velocity_comparison/report.html">频散对比</a> · <a href="paired_comparison/report.html">结构对比（反演后）</a></p>'
    report+=pd.DataFrame(references).to_html(index=False)+checked[['curve_id','frequency_hz','accepted','band_reason','background_null_p']].to_html(index=False)
    (dest/'phase_report.html').write_text(report,encoding='utf-8');dispersion_report(dest,dest/'velocity_comparison');trace.note('投影重建完成',f'{len(final)}实测点；低频{marker["low_frequency_measured_points"]}点');return marker

def dispersion_report(folder,output):
    folder=Path(folder);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    before=pd.read_csv(folder/'before_projected_curves.csv');after=pd.read_csv(folder/'inversion_curves.csv');fig,ax=plt.subplots(figsize=(9,5))
    for label,table,style in [('Before (projected)',before,'--'),('After (measured)',after,'o-')]:
        for i,(_,g) in enumerate(table.groupby('curve_id')):
            g=g.sort_values('frequency_hz');ax.plot(g.frequency_hz,g.phase_velocity_m_s,style,ms=3,label=label if i==0 else None)
    ax.set(xlabel='Frequency / Hz',ylabel='Wave-vector phase velocity / m/s');ax.grid(alpha=.25)
    if ax.get_legend_handles_labels()[0]:ax.legend()
    fig.tight_layout();fig.savefig(output/'dispersion.png',dpi=150);plt.close(fig)
    (output/'report.html').write_text('<meta charset="utf-8"><h1>沿波矢相速度：修正前后</h1><p>'+PROJECTION_NOTE+'</p><p>虚线为原曲线经同一方向投影；实线为合格实测叠加相位反算。先验预测不绘为实测。</p><img src="dispersion.png" style="max-width:100%">',encoding='utf-8')
    return dict(projection_applied=True,measured_points=len(after),output=str(output/'report.html'))

def comparison_tables(before,after,scan):
    from .phase_bands import assign_bands
    before=before.drop(columns=['parent_curve_id'],errors='ignore').rename(columns={'curve_id':'original_curve_id'}).copy()
    joined=after[['curve_id','parent_curve_id','frequency_hz']].merge(before,left_on=['parent_curve_id','frequency_hz'],right_on=['original_curve_id','frequency_hz'],validate='many_to_one');joined['accepted']=True
    checked=pd.concat([assign_bands(g,np.sort(after.loc[after.curve_id==cid,'frequency_hz'].unique()),scan) for cid,g in joined.groupby('curve_id',sort=False)]) if len(joined) else assign_bands(joined,[],scan)
    parts=[]
    for band,g in checked[checked.accepted].groupby('band_id',sort=False):
        g=g.sort_values('frequency_hz').copy();old=g.curve_id.iloc[0];g['curve_id']=band
        other=after[(after.curve_id==old)&after.frequency_hz.isin(g.frequency_hz)].sort_values('frequency_hz').copy();other['curve_id']=band
        for table in (g,other):table['group_velocity_m_s']=derived_group(table.frequency_hz,table.phase_velocity_m_s);table['group_velocity_is_derived']=True
        parts.append((g,other))
    return (pd.concat([x[0] for x in parts],ignore_index=True),pd.concat([x[1] for x in parts],ignore_index=True)) if parts else (joined.iloc[:0],after.iloc[:0])

def paired_compare(project,selection=None,folder=None):
    from dispersion.common import load_config
    from .real_compare import paired_inversion
    from .inversion import selection_digest
    from .candidates import read_curves
    root=Path(project).resolve();phase=root/'output/phase_guided';marker=json.loads((phase/'phase_stack_complete.json').read_text(encoding='utf-8-sig'))
    if not marker.get('reconstruction',{}).get('enabled'):raise ValueError('配对投影反演需要先完成低频投影重建')
    if marker['reconstruction_module_sha256']!=sha(__file__) or marker['settings_sha256']!=sha(root/'config.json'):raise ValueError('低频重建已过期，须重跑相位阶段')
    after,_=read_curves(root);before=pd.read_csv(phase/'before_projected_curves.csv');scan=load_config(root/'config.json')['scan'];selection=selection or json.loads((root/'selection.json').read_text(encoding='utf-8-sig'))
    keep=np.zeros(len(after),bool)
    for case in selection['cases']:
        for band in case['bands']:keep |= after.curve_id.eq(band['curve_id']) & after.frequency_hz.between(band.get('fmin_hz',-np.inf),band.get('fmax_hz',np.inf))
    after=after[keep].copy();dest=Path(folder or phase/'paired_comparison');dest.mkdir(parents=True,exist_ok=True)
    for name in ('comparison_complete.json','report.html','comparison.png','DARE_Full_Comparison.html'):(dest/name).unlink(missing_ok=True)
    b,a=comparison_tables(before,after,scan);csv(dest/'before_common_curves.csv',b);csv(dest/'after_common_curves.csv',a);csv(dest/'after_extended_curves.csv',after)
    cases=[]
    for case in selection['cases']:
        bands=[dict(band,curve_id=cid) for band in case['bands'] for cid in a.curve_id.unique() if cid.startswith(band['curve_id']+'__B')]
        if all(any(x['curve_id'].startswith(v['curve_id']+'__B') for x in bands) for v in case['bands']):cases.append(dict(case,bands=bands))
    common=dict(selection,cases=cases);write_json(dest/'common_selection.json',common)
    status=dict(complete=True,phase_completion_sha256=sha(phase/'phase_stack_complete.json'),selection_sha256=selection_digest(selection),common_points=len(a),extended_points=len(after),same_frequency_and_band_comparison=True,projection_applied_to_both=True,extrapolated_predictions_in_inversion=False,common_cases=len(cases),extended_cases=len(selection['cases']))
    if cases:
        for name,table in [('before_common',b),('after_common',a)]:status[name]=paired_inversion(table,common,dest/name,scan)
    else:
        for name in ('before_common','after_common'):
            if (dest/name).exists():shutil.rmtree(dest/name)
        status['common_inversion_skipped']='共同频带不足，未伪造配对结构'
    if selection['cases'] and len(after):status['after_extended']=paired_inversion(after,selection,dest/'after_extended',scan)
    elif (dest/'after_extended').exists():shutil.rmtree(dest/'after_extended')
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,5))
    for label,table,style in [('Before (common)',b,'--'),('After (common)',a,'-'),('After (extended)',after,':')]:
        for i,(_,g) in enumerate(table.groupby('curve_id')):axes[0].plot(g.frequency_hz,g.phase_velocity_m_s,style,label=label if i==0 else None)
    axes[0].set(xlabel='Frequency / Hz',ylabel='Wave-vector phase velocity / m/s');axes[0].grid(alpha=.25);plotted=False
    for label,name in [('Before (common)','before_common'),('After (common)','after_common'),('After (extended)','after_extended')]:
        path=dest/name/'structure_layers.csv'
        if not path.exists():continue
        layers=pd.read_csv(path)
        for i,(_,g) in enumerate(layers.groupby('model_id')):
            depth=[];vel=[]
            for r in g.itertuples():depth.extend([r.top_depth_m,r.bottom_depth_m if np.isfinite(r.bottom_depth_m) else max(100.,r.top_depth_m*1.3)]);vel.extend([r.vs_m_s,r.vs_m_s])
            axes[1].plot(vel,depth,alpha=.35,label=label if i==0 else None);plotted=True
    axes[1].set(xlabel='Vs / m/s',ylabel='Depth / m');axes[1].invert_yaxis();axes[1].grid(alpha=.25)
    for ax in axes:
        if ax.get_legend_handles_labels()[0]:ax.legend(fontsize=8)
    if not plotted:axes[1].text(.5,.5,'No accepted structural models',transform=axes[1].transAxes,ha='center')
    fig.tight_layout();fig.savefig(dest/'comparison.png',dpi=150);plt.close(fig)
    page='<meta charset="utf-8"><h1>修正前后频散与结构反演</h1><p>'+PROJECTION_NOTE+'</p><p>共同频段采用相同频点、模态、层数、边界、种子和求解器。扩展反演另用完整实测频段；预测点不进入反演，结构仍是条件模型。</p><img src="comparison.png" style="max-width:100%"><p>'
    for name in ('before_common','after_common','after_extended'):
        if (dest/name/'structure_layers.csv').exists():page+=f'<a href="{name}/structure_layers.csv">{name}结构表</a> · '
    page+='</p>';(dest/'report.html').write_text(page,encoding='utf-8');write_json(dest/'comparison_complete.json',status);return status
