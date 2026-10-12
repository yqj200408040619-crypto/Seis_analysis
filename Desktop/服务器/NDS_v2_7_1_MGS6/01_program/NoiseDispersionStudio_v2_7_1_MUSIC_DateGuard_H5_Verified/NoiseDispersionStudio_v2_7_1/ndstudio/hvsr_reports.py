"""Offline plots and self-contained HTML report for single-station HVSR."""
from pathlib import Path
import base64
import html
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def report(output,table,curves,starts,summary):
    output=Path(output);figures=output/'figures';figures.mkdir(exist_ok=True)
    f=table.frequency_hz.to_numpy();mean=table.hv.to_numpy();std=table.hv_std_factor.to_numpy()
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(10,5.4),constrained_layout=True)
    for y in curves[::max(1,len(curves)//70)]: ax.semilogx(f,y,color='#bac8d4',alpha=.3,lw=.6)
    ax.fill_between(f,mean/std,mean*std,color='#4287b9',alpha=.18,label='Window scatter: +/- 1 log sigma')
    ax.semilogx(f,mean,color='#155784',lw=2.2,label='Geometric mean H/V')
    ax.semilogx(f,table.n_over_v,'--',color='#d48736',alpha=.7,lw=1,label='N/V')
    ax.semilogx(f,table.e_over_v,':',color='#429376',alpha=.8,lw=1.4,label='E/V')
    for p in summary['peaks'][:6]:
        ax.plot(p['f0_hz'],p['amplitude'],'o',color='#c74d42' if p['selected'] else '#7d699c')
        if p['selected']: ax.annotate(f"{p['f0_hz']:.3f} Hz\nclarity {p['clarity_pass_count']}/6",(p['f0_hz'],p['amplitude']),xytext=(10,12),textcoords='offset points')
    ax.axhline(1,color='#999',lw=.6);ax.set(xlabel='Frequency (Hz)',ylabel='H/V amplitude ratio',title=f"{summary['instrument_id']} | {summary['independent_windows']} independent windows")
    ax.set_ylim(0,min(max(float((mean*std).max())*1.12,3),max(25,float(mean.max())*1.4)))
    ax.grid(alpha=.2);ax.legend(fontsize=9);fig.savefig(figures/'hvsr_curve.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(10,7),constrained_layout=True)
    for key,name,color in [('n_amplitude','N','#d48736'),('e_amplitude','E','#429376'),('z_amplitude','Z','#155784')]: axes[0].loglog(f,table[key],label=name,color=color)
    axes[0].set(xlabel='Frequency (Hz)',ylabel='Smoothed Fourier amplitude',title='Component spectra (inspect peaks and vertical minima)');axes[0].legend();axes[0].grid(alpha=.2)
    # Row positions preserve irregular UTC gaps; draw individual time-window cells.
    hours=(starts-starts[0])/3600;edges=np.geomspace(f[0]/np.sqrt(f[1]/f[0]),f[-1]*np.sqrt(f[-1]/f[-2]),len(f)+1)
    for i,h in enumerate(hours):
        mesh=axes[1].pcolormesh(edges,[h,h+summary['settings']['window_s']/3600],np.log10(curves[i:i+1]),vmin=-.3,vmax=max(.6,float(np.log10(curves).max())),cmap='viridis',shading='flat',rasterized=True)
    fig.colorbar(mesh,ax=axes[1],label='log10(H/V)')
    axes[1].set_xscale('log');axes[1].set(xlabel='Frequency (Hz)',ylabel='Hours from first accepted window',title='log10(H/V) per window; blank rows = unavailable/rejected')
    fig.savefig(figures/'hvsr_spectra_time.png',dpi=140);plt.close(fig)
    names=['hvsr_curve.png','hvsr_spectra_time.png']
    structure=summary['structure']
    if structure['status']=='conditional':
        d=pd.read_csv(output/'depth_vs_tradeoff.csv');fig,axes=plt.subplots(1,2,figsize=(11,5.7),constrained_layout=True)
        axes[0].fill_between(d.vs_m_s,d.depth_lower_m,d.depth_upper_m,color='#4287b9',alpha=.2,label='Observed-window frequency scatter')
        axes[0].plot(d.vs_m_s,d.depth_at_f0_m,color='#155784',lw=2,label='H = Vs / (4 f0)')
        axes[0].set(xlabel='Assumed equivalent Vs (m/s)',ylabel='Equivalent layer thickness (m)',title='Depth-Vs trade-off');axes[0].legend(fontsize=9);axes[0].grid(alpha=.2)
        bottom=max(float(d.depth_upper_m.max())*1.3,30)
        for i in np.linspace(0,len(d)-1,10).astype(int):
            row=d.iloc[i];vs=row.vs_m_s;h=row.depth_at_f0_m
            axes[1].plot([vs,vs,structure['assumed_bedrock_vs_m_s'],structure['assumed_bedrock_vs_m_s']],[0,h,h,bottom],color='#4287b9',alpha=.35)
        axes[1].set(xlabel='Vs (m/s)',ylabel='Depth (m)',title='Conditional profiles: halfspace Vs is a prior');axes[1].invert_yaxis();axes[1].grid(alpha=.2)
        fig.suptitle('Equivalent cover over halfspace; not a unique layered inversion',fontsize=12);fig.savefig(figures/'shallow_structure.png',dpi=160);plt.close(fig);names.append('shallow_structure.png')
    ell=structure.get('ellipticity')
    if ell:
        fit=pd.read_csv(output/'ellipticity_fit.csv');fig,axes=plt.subplots(1,2,figsize=(11,5.7),constrained_layout=True)
        axes[0].loglog(fit.frequency_hz,fit.observed_hv,'k',label='Observed H/V');axes[0].loglog(fit.frequency_hz,fit.best_rayleigh_hv,color='#c74d42',label='Best Rayleigh-only trial')
        axes[0].fill_between(fit.frequency_hz,fit.retained_min_hv,fit.retained_max_hv,color='#c74d42',alpha=.15);axes[0].set(xlabel='Frequency (Hz)',ylabel='H/V',title='Conditional ellipticity fit');axes[0].legend(fontsize=9);axes[0].grid(alpha=.2)
        bottom=max(m['thickness_m'] for m in ell['retained_models'])*1.3
        for m in ell['retained_models'][:40]: axes[1].plot([m['vs_m_s'],m['vs_m_s'],m['bedrock_vs_m_s'],m['bedrock_vs_m_s']],[0,m['thickness_m'],m['thickness_m'],bottom],color='#c74d42',alpha=.15)
        axes[1].set(xlabel='Vs (m/s)',ylabel='Depth (m)',title=f"Exploratory ensemble; fit accepted={ell['acceptable_fit']}");axes[1].invert_yaxis();axes[1].grid(alpha=.2)
        fig.savefig(figures/'rayleigh_ellipticity.png',dpi=150);plt.close(fig);names.append('rayleigh_ellipticity.png')
    images=''.join('<img alt="'+name+'" src="data:image/png;base64,'+base64.b64encode((figures/name).read_bytes()).decode()+'">' for name in names)
    peak=summary['selected_peak']; ptext=f"{peak['f0_hz']:.4g} Hz，可靠性：{peak['reliability_pass']}，清晰度：{peak['clarity_pass_count']}/6" if peak else '无内部局部峰'
    if structure['status']=='conditional':
        lo,hi=structure['depth_at_f0_range_m'];stext=f"Vs 先验 {structure['vs_prior_m_s']} m/s → 等效层厚 {lo:.1f}–{hi:.1f} m。半空间 Vs={structure['assumed_bedrock_vs_m_s']:g} m/s 为输入先验。"
    else: stext=structure['explanation']
    content=f'''<!doctype html><html lang="zh"><meta charset="utf-8"><title>单台 H/V 与浅层约束</title>
<style>body{{max-width:1120px;margin:32px auto;padding:0 22px;font:16px/1.7 system-ui;color:#19384a}}h1{{font-size:28px}}img{{width:100%;margin:16px 0}}table{{border-collapse:collapse;font-size:13px}}td,th{{border:1px solid #d9e2e9;padding:6px;white-space:nowrap}}.tablewrap{{overflow-x:auto;width:100%}}.note{{background:#edf4f7;padding:18px;border-radius:8px}}pre{{white-space:pre-wrap}}</style>
<h1>单台 H/V 与浅层结构约束</h1><p>{html.escape(summary['instrument_id'])} · 有效窗 {summary['accepted_windows']} / {summary['windows']} · 不重叠统计窗 {summary['independent_windows']}</p>
<p>所选峰：{ptext}</p><div class="note">{html.escape(stext)}<br>H/V 约束共振周期与 S 波走时，单独不能唯一分离 Vs 与层厚；多峰不自动解释为多层。离散度带是时窗间变化，并非均值置信区间。</div>
{images}<h2>峰值逐项检查（可横向滚动）</h2><div class="tablewrap">{pd.DataFrame(summary['peaks']).to_html(index=False,escape=True)}</div>
<h2>解释与复核</h2><ul>{''.join('<li>'+html.escape(w)+'</li>' for w in summary['warnings'])}</ul>
<p>可选 Rayleigh 拟合使用基阶、各向同性 Rayleigh 主导假设，忽略 Love 波和体波；拟合不佳不构成有效结构结果。所有固定先验保存于结果 JSON。</p>
<p>质量判据依据 <a href="https://sesame.geopsy.org/Papers/HV_User_Guidelines.pdf">SESAME H/V Guidelines</a>：3 项可靠性、6 项清晰度；可靠性全通过且清晰度至少5/6为清晰峰评分条件；自动层厚估计还要求峰两侧均观测到充分下降。该额外条件是程序的解释保护，不改变SESAME评分。窄峰来源、仪器噪声与基频认定仍需复核。</p>
<details><summary>本次参数</summary><pre>{html.escape(json.dumps(summary['settings'],ensure_ascii=False,indent=2))}</pre></details></html>'''
    (output/'hvsr_report.html').write_text(content,encoding='utf-8')
