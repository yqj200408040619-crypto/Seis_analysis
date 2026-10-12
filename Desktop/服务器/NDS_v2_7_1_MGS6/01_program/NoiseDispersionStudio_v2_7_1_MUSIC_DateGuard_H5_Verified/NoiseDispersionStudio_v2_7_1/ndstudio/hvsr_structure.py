"""Conditional shallow constraints and opt-in Rayleigh-only ellipticity search.

Quarter-wave estimates constrain integrated S-wave travel time, not a unique
velocity profile. Ellipticity search is a one-layer exploratory approximation,
never silently equated with mixed-wavefield ambient-noise HVSR.
"""
import numpy as np
import pandas as pd
from .hvsr import write_json


def ellipticity_forward(f,thickness_m,vs_m_s,bedrock_vs_m_s,c):
    from disba import Ellipticity
    # disba uses km, km/s, g/cm3, and periods in increasing order.
    solver=Ellipticity(np.array([thickness_m/1000,0.]),
                       np.array([vs_m_s*c['vp_vs'],bedrock_vs_m_s*c['bedrock_vp_vs']])/1000,
                       np.array([vs_m_s,bedrock_vs_m_s])/1000,
                       np.array([c['density_g_cm3'],c['bedrock_density_g_cm3']]))
    periods=1/np.asarray(f)[::-1]; curve=solver(periods,mode=0)
    if len(curve.period)!=len(periods) or not np.allclose(curve.period,periods): raise ValueError('Ellipticity solver incomplete period coverage')
    hv=np.abs(curve.ellipticity[::-1])
    if c['horizontal'] in ('rms','geometric'): hv/=np.sqrt(2)
    if not np.isfinite(hv).all() or np.any(hv<=1e-8): raise ValueError('Ellipticity undefined')
    return hv


def ellipticity_search(f,mean,std,c,output):
    rng=np.random.default_rng(int(c['ellipticity_seed']))
    # Work on at most 48 log-spaced frequencies; never manufacture degrees of freedom by oversampling.
    ix=np.unique(np.round(np.linspace(0,len(f)-1,min(48,len(f)))).astype(int));ff=f[ix]
    vv=mean[ix];sigma=np.maximum(np.log(std[ix]),.15)
    rows=[];curves=[];failures=0
    for k in range(int(c['ellipticity_trials'])):
        vs=float(np.exp(rng.uniform(np.log(c['vs_min_m_s']),np.log(c['vs_max_m_s']))))
        h=float(np.exp(rng.uniform(np.log(c['thickness_min_m']),np.log(c['thickness_max_m']))))
        try:
            pred=ellipticity_forward(ff,h,vs,c['bedrock_vs_m_s'],c)
            residual=np.log(pred/vv)
            score=float(np.sqrt(np.mean((residual/sigma)**2)))
            raw=float(np.sqrt(np.mean(residual**2)))
            rows.append(dict(trial=k,thickness_m=h,vs_m_s=vs,bedrock_vs_m_s=c['bedrock_vs_m_s'],weighted_log_rms=score,log_rms=raw))
            curves.append(pred)
        except Exception as exc:
            # Numerical trial failures are counted. Dependency/installation errors must not be hidden.
            if isinstance(exc,(ImportError,ModuleNotFoundError)): raise
            failures+=1
        if (k+1)%40==0: print(f"Rayleigh ellipticity: {k+1}/{c['ellipticity_trials']} trials, {failures} solver failures",flush=True)
    if not rows: raise ValueError('Rayleigh ellipticity: all model trials failed')
    order=np.argsort([r['weighted_log_rms'] for r in rows]); rows=[rows[i] for i in order];curves=np.asarray(curves)[order]
    best=rows[0]['weighted_log_rms'];keep=[i for i,r in enumerate(rows) if r['weighted_log_rms']<=max(1.,best*1.25)]
    pd.DataFrame(rows).to_csv(output/'ellipticity_trials.csv',index=False,encoding='utf-8-sig')
    fit=pd.DataFrame({'frequency_hz':ff,'observed_hv':vv,'best_rayleigh_hv':curves[0],
                      'retained_min_hv':curves[keep].min(axis=0),'retained_max_hv':curves[keep].max(axis=0)})
    fit.to_csv(output/'ellipticity_fit.csv',index=False,encoding='utf-8-sig')
    result=dict(method='conditional fundamental-mode Rayleigh ellipticity; one finite layer over fixed halfspace',
                assumption='isotropic Rayleigh-dominated surface wavefield; Love/body waves neglected',
                geometric_mean_approximation=c['horizontal']=='geometric',
                model_trials=int(c['ellipticity_trials']),successful_trials=len(rows),solver_failures=failures,
                best_model=rows[0],retained_models=[rows[i] for i in keep],
                acceptable_fit=bool(best<=1. and rows[0]['log_rms']<=.35),
                uncertainty='trial ensemble only, not posterior confidence; increase trials/check priors',
                fixed_prior=dict(bedrock_vs_m_s=c['bedrock_vs_m_s'],vp_vs=c['vp_vs'],bedrock_vp_vs=c['bedrock_vp_vs'],
                                 density_g_cm3=c['density_g_cm3'],bedrock_density_g_cm3=c['bedrock_density_g_cm3']))
    write_json(output/'ellipticity_models.json',result)
    return result


def analyze_structure(f,mean,std,peak,c,output):
    result=dict(method='quarter-wave conditional travel-time constraint',status='no_reliable_peak',
                explanation='平坦、边缘或不可靠峰不自动生成结构。H/V单独不能唯一确定Vs、层厚或Vs30。')
    if not peak or not peak['clear_peak']:
        write_json(output/'shallow_structure.json',result);return result
    # SESAME permits five of six clarity checks. A one-sided shoulder may
    # pass that score without establishing a fundamental resonance. Keep the
    # score intact but require observed decay on both sides for structure.
    if not (peak.get('left_drop',False) and peak.get('right_drop',False)):
        result.update(status='incomplete_peak_shape',candidate_f0_hz=peak['f0_hz'],
                      explanation='所选峰虽通过常规清晰度评分，但未观测到两侧充分下降，不能自动认定为基频；暂停层厚估计，须扩大频段并复核窗长、响应和噪声来源。')
        write_json(output/'shallow_structure.json',result);return result
    f0=peak['f0_hz'];sigma=peak['sigma_f_hz'] or 0.
    # Observed-window scatter band, explicitly NOT an uncertainty interval on the mean.
    flo=max(float(f[0]),f0-1.96*sigma);fhi=min(float(f[-1]),f0+1.96*sigma)
    flo=max(flo,f0*.5);fhi=max(fhi,f0)
    v=np.linspace(c['vs_min_m_s'],c['vs_max_m_s'],60)
    table=pd.DataFrame(dict(vs_m_s=v,depth_at_f0_m=v/(4*f0),depth_lower_m=v/(4*fhi),depth_upper_m=v/(4*flo),
                            assumed_bedrock_vs_m_s=c['bedrock_vs_m_s']))
    table.to_csv(output/'depth_vs_tradeoff.csv',index=False,encoding='utf-8-sig')
    result.update(status='conditional',f0_hz=f0,period_s=1/f0,
                  s_wave_one_way_travel_time_s=1/(4*f0),frequency_scatter_band_hz=[flo,fhi],
                  vs_prior_m_s=[c['vs_min_m_s'],c['vs_max_m_s']],
                  depth_at_f0_range_m=[float(v[0]/(4*f0)),float(v[-1]/(4*f0))],
                  depth_with_window_scatter_range_m=[float(v[0]/(4*fhi)),float(v[-1]/(4*flo))],
                  assumed_bedrock_vs_m_s=c['bedrock_vs_m_s'],
                  formula='f0 ~= 1/(4*sum(h_i/Vs_i)); uniform equivalent layer H ~= Vs/(4*f0)',
                  explanation='一维高阻抗反差、所选峰为基频的条件下，输出等效覆盖层厚度范围；半空间Vs为输入先验。未自动输出Vs30。')
    if c['rayleigh_assumption']:
        if c['bedrock_vs_m_s']<=c['vs_max_m_s']: raise ValueError('Rayleigh单覆盖层模型要求半空间Vs大于覆盖层Vs先验上限')
        result['ellipticity']=ellipticity_search(f,mean,std,c,output)
    write_json(output/'shallow_structure.json',result)
    return result
