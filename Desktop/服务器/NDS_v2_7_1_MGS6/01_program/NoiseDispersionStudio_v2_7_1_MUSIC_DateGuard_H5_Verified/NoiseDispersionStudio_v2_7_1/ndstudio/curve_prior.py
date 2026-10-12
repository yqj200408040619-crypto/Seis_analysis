"""Explicit monotone curve selection and bounded screening-only extrapolation."""
import numpy as np
import pandas as pd
from dispersion.common import csv


def monotone_curves(frame, velocity='phase_velocity_m_s', tolerance=1e-8):
    rows=[]; keep=[]
    for cid,g in frame.groupby('curve_id',sort=False):
        g=g.sort_values('frequency_hz'); dv=np.diff(g[velocity].to_numpy(float))
        accepted=len(g)>=2 and np.isfinite(g[velocity]).all() and bool((dv<=tolerance).all())
        rows.append(dict(curve_id=cid,points=len(g),increasing_steps=int((dv>tolerance).sum()),
            largest_increase_m_s=float(dv.max()) if len(dv) else None,accepted=accepted,
            reason='accepted' if accepted else 'velocity_increases_with_frequency_or_invalid'))
        if accepted:keep.append(cid)
    return frame[frame.curve_id.isin(keep)].copy(),pd.DataFrame(rows,columns=['curve_id','points','increasing_steps','largest_increase_m_s','accepted','reason'])


def screening_prior(targets,dest,scan,options=None,monotonic_only=False,trace=None):
    options=options or {}; unknown=set(options)-{'enabled','extension_hz','maximum_log_slope'}
    if unknown:raise ValueError('未知外延参数: '+str(sorted(unknown)))
    q=targets.copy(); q['is_prior_prediction']=False
    if monotonic_only:
        q,audit=monotone_curves(q);csv(dest/'input_curve_trend_screen.csv',audit)
        if trace:trace.note('升频增速曲线筛除',f'{len(audit)}条候选 → 保留{int(audit.accepted.sum())}条；逐相邻频点要求非增速。')
    if q.empty:raise ValueError('没有符合升频非增速要求的输入曲线')
    extension=float(options.get('extension_hz',.1));limit=float(options.get('maximum_log_slope',3.))
    if not np.isfinite([extension,limit]).all() or extension<0 or limit<=0:raise ValueError('外延范围或斜率限制无效')
    predictions=[]
    if options.get('enabled',False) and extension>0:
        grid=np.round(np.arange(scan['min_hz'],scan['max_hz']+scan['step_hz']/2,scan['step_hz']),9)
        for cid,g in q.groupby('curve_id',sort=False):
            g=g.sort_values('frequency_hz');f=g.frequency_hz.to_numpy();v=g.phase_velocity_m_s.to_numpy()
            # Local log-linear endpoint fits, clipped to a decreasing power law.
            for left,endpoint in [(True,g.iloc[0]),(False,g.iloc[-1])]:
                use=g.iloc[:min(5,len(g))] if left else g.iloc[-min(5,len(g)):]
                slope=float(np.clip(np.polyfit(np.log(use.frequency_hz),np.log(use.phase_velocity_m_s),1)[0],-limit,0.))
                extra=grid[(grid>=f[0]-extension-1e-9)&(grid<f[0]-1e-9)] if left else grid[(grid>f[-1]+1e-9)&(grid<=f[-1]+extension+1e-9)]
                for freq in extra:
                    row=endpoint.to_dict();velocity=float(endpoint.phase_velocity_m_s*(freq/endpoint.frequency_hz)**slope)
                    row.update(frequency_hz=float(freq),phase_velocity_m_s=velocity,is_prior_prediction=True,
                        is_model_prediction=True,background_discriminated=False,velocity_min_m_s=np.nan,velocity_max_m_s=np.nan,
                        prior_method='bounded_endpoint_power_law',prior_log_slope=slope,
                        measured_unwrapped_phase_rad=np.nan,relative_phase_q025_rad=np.nan,relative_phase_q975_rad=np.nan,
                        prediction_reconstruction_error_rad=np.nan)
                    # Different cycle aliases share observed phase gates, but
                    # their separate endpoint fits need not predict the same
                    # phase outside the observed band. Never collapse those.
                    row['segment_id']=str(row.get('segment_id',cid))+'__prior_'+cid
                    row['phase_equivalence_id']=row['segment_id']
                    sign=1. if row['branch']=='positive' else -1.
                    theta=sign*(row['source_phase_rad']-2*np.pi*freq*row['distance_m']/velocity+2*np.pi*row['constant_cycle'])
                    row.update(predicted_unwrapped_phase_rad=theta,predicted_phase_deg=float(np.degrees(np.angle(np.exp(1j*theta)))))
                    predictions.append(row)
    if predictions:q=pd.concat([q,pd.DataFrame(predictions)],ignore_index=True)
    csv(dest/'screening_prior_predictions.csv',q[q.is_prior_prediction])
    if trace:trace.note('相位筛选先验外延',f'{len(predictions)}个外延先验点；均标记为预测，禁止流入观测频散及结构反演。',extension_hz=extension)
    return q
