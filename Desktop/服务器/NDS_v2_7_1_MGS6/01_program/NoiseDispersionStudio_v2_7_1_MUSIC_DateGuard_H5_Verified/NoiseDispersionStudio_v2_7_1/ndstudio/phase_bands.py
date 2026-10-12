"""Recheck contiguous measured bands after window rejection; no gap filling."""
import numpy as np
import pandas as pd
from dispersion.common import csv, write_json, sha


def assign_bands(frame, frequencies, scan, group_key='curve_id'):
    q=frame.copy();q['band_id']='';q['point_accepted']=q.accepted.astype(bool);q['band_reason']='';q['band_points']=0;q['bandwidth_hz']=0.
    expected=np.sort(np.unique(frequencies)); index={round(float(f),9):i for i,f in enumerate(expected)}
    maximum_gap=1.5*float(scan['step_hz'])+1e-9
    for key,g in q.groupby(group_key,sort=False):
        g=g.sort_values('frequency_hz'); groups=[];current=[];previous=None
        for ix,row in g.iterrows():
            fi=index[round(float(row.frequency_hz),9)]
            if not row.accepted:
                if current:groups.append(current);current=[]
                previous=None;continue
            direction=row.get('selected_direction_group',-1)
            smooth_direction='confirmed_azimuth_deg' in row and np.isfinite(row.confirmed_azimuth_deg)
            changed=direction!=previous[2] if previous else False
            if smooth_direction and previous and len(previous)>3:
                changed=abs((row.confirmed_azimuth_deg-previous[3]+180)%360-180)>scan.get('max_frequency_direction_change_deg',45.)
            if current and (fi!=previous[0]+1 or row.frequency_hz-previous[1]>maximum_gap or changed):groups.append(current);current=[]
            current.append(ix);previous=(fi,float(row.frequency_hz),direction,row.get('confirmed_azimuth_deg',np.nan))
        if current:groups.append(current)
        for number,ids in enumerate(groups,1):
            width=float(q.loc[ids,'frequency_hz'].max()-q.loc[ids,'frequency_hz'].min())
            passed=len(ids)>=scan['min_points'] and width>=scan['min_band_hz']-1e-9
            q.loc[ids,'band_id']=str(key)+'__B'+str(number);q.loc[ids,'band_points']=len(ids);q.loc[ids,'bandwidth_hz']=width
            q.loc[ids,'accepted']=passed;q.loc[ids,'band_reason']='accepted' if passed else '连续点数或最小连续带宽不足'
    q.loc[~q.point_accepted,'band_reason']='窗口/质量检查未通过'
    return q


def handoff(dest, targets, summaries, scan, direction, trace,monotonic_only=False):
    """Preserve CCF interpretation, inherited source/cycle, and explicit band gaps."""
    from pathlib import Path
    dest=Path(dest); rows=summaries[['curve_id','frequency_hz','after_count','valid_windows','windows_file','gate_id']].copy()
    extra=['direction_independent_windows','selected_direction_group']
    for k in extra:
        if k in summaries:rows[k]=summaries[k]
    rows['accepted']=rows.after_count>0
    if direction['mode'] in ('filter','joint'):
        rows['accepted'] &= rows.after_count>=direction['minimum_accepted_windows']
        rows['accepted'] &= rows.direction_independent_windows>=direction['minimum_independent_windows']
    frequencies=np.sort(targets.frequency_hz.unique()); checked=assign_bands(rows,frequencies,scan)
    csv(dest/'phase_continuous_bands.csv',checked)
    ready=checked[checked.accepted].copy()
    for band,g in checked[checked.band_id!=''].groupby('band_id',sort=False):
        trace.note('连续频带',f'{band} · {g.frequency_hz.min():g}–{g.frequency_hz.max():g} Hz · {len(g)}点 · 带宽{g.bandwidth_hz.iloc[0]:g} Hz · {g.band_reason.iloc[0]}')
    trace.note('连续频带参数',f'min_band_hz={scan["min_band_hz"]:g}；min_points={scan["min_points"]}；合格曲线频点={len(ready)}',output=str(dest/'phase_continuous_bands.csv'))
    csv(dest/'qualified_phase_summary.csv',summaries.merge(ready[['curve_id','frequency_hz','band_id']],on=['curve_id','frequency_hz']))
    # Inherit physical and background status, rather than upgrading them after selection.
    columns=list(targets.columns)+['parent_curve_id','band_id']
    result=pd.DataFrame(columns=columns)
    if len(ready):
        result=targets.merge(ready[['curve_id','frequency_hz','band_id']],on=['curve_id','frequency_hz'])
        from .phase_velocity import calculate
        # Each kept band has its own single integer alignment; never unwrap across rejected gaps.
        selected=summaries.merge(ready[['curve_id','frequency_hz','band_id']],on=['curve_id','frequency_hz'])
        selected['curve_id']=selected.band_id;result['parent_curve_id']=result.curve_id;result['curve_id']=result.band_id
        try:
            details,_=calculate(selected,result)
            result['phase_velocity_m_s']=details.set_index(['curve_id','frequency_hz']).loc[pd.MultiIndex.from_frame(result[['curve_id','frequency_hz']]),'raw_velocity_m_s'].to_numpy()
            result['wavenumber_rad_m']=2*np.pi*result.frequency_hz/result.phase_velocity_m_s
            result['distance_over_wavelength']=result.distance_m*result.frequency_hz/result.phase_velocity_m_s
        except ValueError as exc:
            trace.note('反演交接未生成',str(exc));result=result.iloc[:0]
        for k,default in dict(cohort='phase',pair_index=0,group_velocity_m_s=np.nan,instrument_response_corrected=True,
            background_discriminated=False,source_phase_verified=False,is_model_prediction=False).items():
            if k not in result:result[k]=default
        if 'projection_applied' in result and result.projection_applied.eq(True).all():
            from .reconstruction import derived_group
            result['velocity_type']='wave_vector_phase_velocity';result['direction_method']='music'
            result['apparent_pair_velocity_m_s']=result.phase_velocity_m_s*result.distance_m/abs(result.projected_distance_m)
            result['distance_over_wavelength']=abs(result.projected_distance_m)*result.frequency_hz/result.phase_velocity_m_s
            for _,g in result.groupby('curve_id'):
                g=g.sort_values('frequency_hz');result.loc[g.index,'group_velocity_m_s']=derived_group(g.frequency_hz,g.phase_velocity_m_s)
            result['group_velocity_is_derived']=True
        result['conditional_interpretation']=True;result['direction_constrained']=direction['mode'] in ('filter','joint')
        if monotonic_only and len(result):
            from .curve_prior import monotone_curves
            result,audit=monotone_curves(result);csv(dest/'corrected_curve_trend_screen.csv',audit)
            trace.note('修正后升频增速筛除',f'{len(audit)}条连续曲线 → 保留{int(audit.accepted.sum())}条。')
    csv(dest/'inversion_curves.csv',result)
    if monotonic_only:
        selected_keys=result[['parent_curve_id','frequency_hz','band_id']].rename(columns={'parent_curve_id':'curve_id'}) if len(result) else ready.iloc[:0][['curve_id','frequency_hz','band_id']]
        csv(dest/'qualified_phase_summary.csv',summaries.merge(selected_keys,on=['curve_id','frequency_hz']))
    catalog=[]
    for cid,g in result.groupby('curve_id',sort=False):
        catalog.append(dict(curve_id=cid,segment_id=g.segment_id.iloc[0],cohort=g.cohort.iloc[0],pair_index=int(g.pair_index.iloc[0]),
            pair=g.pair.iloc[0],branch=g.branch.iloc[0],constant_cycle=int(g.constant_cycle.iloc[0]),source_phase_rad=float(g.source_phase_rad.iloc[0]),
            fmin_hz=float(g.frequency_hz.min()),fmax_hz=float(g.frequency_hz.max()),points=len(g),
            background_discriminated=bool(g.background_discriminated.all()),compatible_other_pairs=0,parent_curve_id=g.parent_curve_id.iloc[0]))
    csv(dest/'inversion_catalog.csv',catalog,columns=None if catalog else ['curve_id','segment_id','cohort','pair_index','pair','branch','constant_cycle','source_phase_rad','fmin_hz','fmax_hz','points','background_discriminated','compatible_other_pairs','parent_curve_id'])
    trace.note('反演交接',f'{len(catalog)}条合格连续曲线；{len(result)}个频点。沿用原始源相位、周期及背景检验状态。',output=str(dest/'inversion_curves.csv'))
    return len(result)
