"""Invert final accepted complex-stack phases on inherited constant cycle branches.

This is a conditional reconstruction: the dispersion curve supplied the phase gate.
It cannot establish the absolute cycle, source phase, or an independent dispersion.
"""
from pathlib import Path
import json
import hashlib
import base64
import html
import numpy as np
import pandas as pd


def wrap(a):
    return np.angle(np.exp(1j * a))


def lift_band(wrapped_outward, anchor):
    """Continuous unwrap followed by ONE integer alignment for the whole band."""
    unwrapped = np.unwrap(np.asarray(wrapped_outward, dtype=float))
    shift = int(np.rint(np.median((anchor - unwrapped) / (2*np.pi))))
    lifted = unwrapped + 2*np.pi*shift
    if np.any(np.abs(lifted - anchor) >= np.pi):
        raise ValueError('Band cannot be aligned by one global integer: inspect gaps/phase jumps')
    return lifted, shift


def calculate(summary, targets):
    keys = ['curve_id', 'frequency_hz']
    if summary.duplicated(keys).any() or targets.duplicated(keys).any():
        raise ValueError('Duplicate curve-frequency key')
    extras = [c for c in targets if c not in summary or c in keys]
    q = summary.merge(targets[extras], on=keys, validate='one_to_one', how='outer', indicator=True)
    if not q['_merge'].eq('both').all():
        raise ValueError('Stack/target frequency coverage differs')
    q = q.drop(columns='_merge')
    if (q.after_count <= 0).any():
        raise ValueError('Empty stack: final phase unavailable; do not invert')
    parts = []
    for cid, g in q.groupby('curve_id', sort=True):
        g = g.sort_values('frequency_hz').copy()
        for c in ['constant_cycle','source_phase_rad','distance_m','branch']:
            if g[c].nunique() != 1:
                raise ValueError(f'{cid}: {c} must be constant across a curve band')
        sign = {'positive':1, 'negative':-1}[g.branch.iloc[0]]
        anchor = sign * g.predicted_unwrapped_phase_rad.to_numpy()
        g['source_phase_deg'] = np.degrees(g.source_phase_rad)
        g['branch_sign'] = sign
        for method in ['raw', 'unit']:
            phase = np.radians(g[f'after_{method}_phase_deg'].to_numpy())
            if not np.isfinite(phase).all():
                raise ValueError(f'{cid}: non-finite {method} phase')
            outward, shift = lift_band(sign * phase, anchor)
            residual = np.degrees(wrap(phase-np.radians(g.predicted_phase_deg)))
            inside = (residual >= g.phase_offset_low_deg.to_numpy()-1e-8) & (residual <= g.phase_offset_high_deg.to_numpy()+1e-8)
            # A cone narrower than pi is convex under positive weights: its sum stays inside.
            if np.all(g.phase_offset_high_deg-g.phase_offset_low_deg < 180) and not inside.all():
                raise ValueError(f'{cid}: final {method} phase outside selected phase arc')
            propagation = g.source_phase_rad.to_numpy()-outward+2*np.pi*g.constant_cycle.to_numpy()
            if np.any(propagation <= 0):
                raise ValueError(f'{cid}: non-positive propagation phase on inherited cycle')
            g[f'final_{method}_wrapped_phase_deg'] = np.degrees(phase)
            g[f'final_{method}_outward_unwrapped_phase_deg'] = np.degrees(outward)
            g[f'{method}_global_unwrap_shift'] = shift
            g[f'{method}_propagation_phase_deg'] = np.degrees(propagation)
            distance=abs(g.projected_distance_m) if 'projection_applied' in g and g.projection_applied.eq(True).all() else g.distance_m
            g[f'{method}_velocity_m_s'] = 2*np.pi*g.frequency_hz*distance/propagation
            g[f'{method}_velocity_relative_change'] = g[f'{method}_velocity_m_s']/g.phase_velocity_m_s-1
            g[f'{method}_gate_residual_deg'] = residual
            g[f'{method}_final_phase_inside_gate'] = inside
        parts.append(g)
    details = pd.concat(parts, ignore_index=True)
    trends = []
    for cid, g in details.groupby('curve_id', sort=True):
        for method in ['raw', 'unit']:
            f = g.frequency_hz.to_numpy()
            v = g[f'{method}_velocity_m_s'].to_numpy()
            dv = np.diff(v)
            trends.append(dict(curve_id=cid, segment_id=g.segment_id.iloc[0], pair=g.pair.iloc[0], branch=g.branch.iloc[0],
                constant_cycle=int(g.constant_cycle.iloc[0]), method=method, points=len(g),
                frequency_min_hz=f[0], frequency_max_hz=f[-1], velocity_at_low_frequency_m_s=v[0],
                velocity_at_high_frequency_m_s=v[-1], endpoint_relative_change=v[-1]/v[0]-1,
                linear_slope_m_s_per_hz=float(np.polyfit(f, v, 1)[0]),
                spearman_rho=float(pd.Series(f).rank().corr(pd.Series(v).rank())),
                decreasing_adjacent_fraction=float(np.mean(dv < -1e-8)), increasing_steps=int(np.sum(dv > 1e-8)),
                strictly_decreasing=bool(np.all(dv < -1e-8)),
                maximum_abs_relative_change_from_input=float(np.max(np.abs(g[f'{method}_velocity_relative_change'])))))
    return details, pd.DataFrame(trends)


def verify_saved_windows(folder, summary):
    """Recompute BOTH final stack phases directly from all saved accepted windows."""
    maximum = 0.0
    count = 0
    for _, row in summary.drop_duplicates('gate_id').iterrows():
        with np.load(folder/row.windows_file) as a:
            c = a['cross'][a['accepted']]
            assert len(c) == int(row.after_count)
            sums = {'raw':np.sum(c), 'unit':np.sum(c/np.abs(c))}
            assert np.allclose(sums['raw'],complex(row.after_raw_real,row.after_raw_imag),rtol=1e-10,atol=1e-25)
            for method, total in sums.items():
                delta = abs(float(np.degrees(wrap(np.angle(total)-np.radians(row[f'after_{method}_phase_deg'])))))
                maximum = max(maximum,delta)
                assert delta < 1e-8, (row.gate_id,method,delta)
            count += 1
    return {'unique_gate_files_checked':count, 'maximum_final_phase_error_deg':maximum}


def figures(details, output, font_path=None, title='曲线相位筛选', focus_curve_ids=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if font_path and Path(font_path).exists():
        from matplotlib import font_manager
        font_manager.fontManager.addfont(str(font_path))
        plt.rcParams['font.family'] = font_manager.FontProperties(fname=str(font_path)).get_name()
    plt.rcParams['axes.unicode_minus'] = False
    curves = list(details.groupby('curve_id',sort=True))
    rows = (len(curves)+2)//2
    fig, axes = plt.subplots(rows,2,figsize=(13,3.6*rows),layout='constrained',squeeze=False)
    for ax, (cid,g) in zip(axes.flat,curves):
        f = g.frequency_hz.to_numpy()
        ax.fill_between(f,g.velocity_min_m_s.to_numpy(),g.velocity_max_m_s.to_numpy(),color='#cbd5e1',alpha=.7,label='原曲线给定的速度范围')
        ax.plot(f,g.phase_velocity_m_s,'--',color='#64748b',lw=1.2,label='筛选输入曲线')
        ax.plot(f,g.raw_velocity_m_s,color='#b45309',lw=1.4,label='原始互谱叠加后反算')
        ax.plot(f,g.unit_velocity_m_s,color='#2563eb',lw=1.4,label='单位相位叠加后反算')
        ax.set_title(f'{cid}\n{g.pair.iloc[0]}  {g.branch.iloc[0]}  d={g.distance_m.iloc[0]:.1f} m',fontsize=10)
        ax.set_xlabel('频率 f (Hz)'); ax.set_ylabel('条件相速度 (m/s)'); ax.grid(alpha=.2)
    for ax in list(axes.flat)[len(curves):]:ax.axis('off')
    handles,labels = axes.flat[0].get_legend_handles_labels()
    axes.flat[-1].legend(handles,labels,loc='upper left',frameon=False,fontsize=11)
    axes.flat[-1].text(.03,.55,'每条曲线固定 N，沿用原曲线源相位。\n相位门限来自输入频散曲线。\n反算结果用于检查筛选后的叠加，\n不能独立证明真实频散趋势。',transform=axes.flat[-1].transAxes,va='top',fontsize=12)
    fig.suptitle(title+'：最终叠加相位与原候选 f–v 曲线',fontsize=16)
    fig.savefig(output/'fv_curves.png',dpi=160); fig.savefig(output/'fv_curves.pdf'); plt.close(fig)
    if focus_curve_ids:
        selected=[(cid,g) for cid,g in curves if cid in focus_curve_ids]
        if set(focus_curve_ids)-set(details.curve_id):raise ValueError('Unknown focus curve IDs')
        fig,axes=plt.subplots(2,len(selected),figsize=(6.5*len(selected),7),layout='constrained',squeeze=False)
        for col,(cid,g) in enumerate(selected):
            ax=axes[0,col];f=g.frequency_hz.to_numpy()
            ax.fill_between(f,g.velocity_min_m_s.to_numpy(),g.velocity_max_m_s.to_numpy(),color='#cbd5e1',alpha=.7)
            ax.plot(f,g.phase_velocity_m_s,'--',color='#64748b',label='原候选曲线')
            ax.plot(f,g.raw_velocity_m_s,'o-',color='#b45309',ms=3,label='幅值加权叠加')
            ax.plot(f,g.unit_velocity_m_s,'s-',color='#2563eb',ms=3,label='等权相位叠加')
            ax.set(title=cid,ylabel='条件相速度 (m/s)');ax.legend(fontsize=9);ax.grid(alpha=.2)
            ax=axes[1,col]
            ax.plot(f,100*g.raw_velocity_relative_change,'o-',color='#b45309',ms=3,label='幅值加权')
            ax.plot(f,100*g.unit_velocity_relative_change,'s-',color='#2563eb',ms=3,label='等权相位')
            ax.axhline(0,color='#64748b',ls='--');ax.set(xlabel='频率 f (Hz)',ylabel='相对原曲线变化 (%)');ax.grid(alpha=.2)
        fig.suptitle(title+'：指定周期分支的叠加对比',fontsize=15)
        fig.savefig(output/'fv_focus.png',dpi=160);plt.close(fig)
    segments=list(details.drop_duplicates(['segment_id','frequency_hz']).groupby('segment_id',sort=True))
    cols=min(2,len(segments));rows=(len(segments)+cols-1)//cols
    fig, axes = plt.subplots(rows,cols,figsize=(6.5*cols,4*rows),layout='constrained',squeeze=False)
    for ax,(sid,g) in zip(axes.flat,segments):
        for r in g.itertuples():
            lo=r.predicted_phase_deg+r.phase_offset_low_deg; hi=r.predicted_phase_deg+r.phase_offset_high_deg
            for k in [-1,0,1]:
                a,b=max(-180,lo+360*k),min(180,hi+360*k)
                if a<=b: ax.vlines(r.frequency_hz,a,b,color='#94a3b8',alpha=.65,lw=2)
        ax.scatter(g.frequency_hz,g.predicted_phase_deg,s=5,color='#64748b',label='门限中心')
        ax.scatter(g.frequency_hz,g.final_raw_wrapped_phase_deg,s=12,color='#b45309',label='原始互谱叠加相位')
        ax.scatter(g.frequency_hz,g.final_unit_wrapped_phase_deg,s=8,color='#2563eb',label='单位相位叠加相位')
        ax.set_title(f'{sid}  {g.pair.iloc[0]}  {g.branch.iloc[0]}',fontsize=10)
        ax.set_ylim(-185,185); ax.set_yticks([-180,-90,0,90,180]);ax.set_xlabel('频率 f (Hz)');ax.set_ylabel('最终互谱相位 (°)');ax.grid(alpha=.2)
    axes.flat[0].legend(fontsize=8,loc='best')
    for ax in list(axes.flat)[len(segments):]:ax.set_visible(False)
    fig.suptitle(title+'：实际互谱相位与接受区间',fontsize=15)
    fig.savefig(output/'final_phases.png',dpi=160);plt.close(fig)


def clear_comparison(output):
    output=Path(output)
    for name in ('final_phase_velocity.csv','velocity_trends.csv','curve_comparison.csv',
                 'verification.json','workbook_data.json','fv_curves.png','fv_curves.pdf',
                 'fv_focus.png','final_phases.png','report.html'):
        (output/name).unlink(missing_ok=True)


def run_phase_velocity(folder, output, verify_windows=True, focus_curve_ids=None, qualified=False):
    marker=Path(folder)/'phase_stack_complete.json'
    if marker.exists() and json.loads(marker.read_text(encoding='utf-8-sig')).get('reconstruction',{}).get('enabled'):
        from .phase_audit import audit
        from .reconstruction import dispersion_report
        if verify_windows:audit(folder)
        return dispersion_report(folder,output)
    folder,output=Path(folder),Path(output)
    output.mkdir(parents=True,exist_ok=True)
    clear_comparison(output)
    summary=pd.read_csv(folder/'phase_stack_summary.csv')
    targets=pd.read_csv(folder/'selected_phase_targets.csv')
    if qualified or ((summary.after_count<=0).any() and (folder/'qualified_phase_summary.csv').exists()):
        kept=pd.read_csv(folder/'qualified_phase_summary.csv')
        targets=targets.merge(kept[['curve_id','frequency_hz','band_id']],on=['curve_id','frequency_hz'])
        summary=kept
        if summary.empty:
            (output/'report.html').write_text('<html><meta charset="utf-8"><h1>没有通过连续频带检查的叠加曲线</h1><p>请查看上级 phase_continuous_bands.csv、方向剔除原因及实际参数。没有生成速度或结构结果。</p></html>',encoding='utf-8')
            return {'qualified_curve_frequency_rows':0,'no_velocity_invented':True}
        summary['parent_curve_id']=summary.curve_id;targets['parent_curve_id']=targets.curve_id
        summary['curve_id']=summary.band_id;targets['curve_id']=targets.band_id
    details,trends=calculate(summary,targets)
    metadata_path=folder/'phase_stack_complete.json'
    metadata=json.loads(metadata_path.read_text(encoding='utf-8')) if metadata_path.exists() else {}
    title=metadata.get('project_name','曲线相位筛选')
    audit=verify_saved_windows(folder,summary) if verify_windows else {'window_verification':'not requested'}
    for method in ['raw','unit']:
        distance=abs(details.projected_distance_m) if 'projection_applied' in details and details.projection_applied.eq(True).all() else details.distance_m
        restored=details.source_phase_rad-2*np.pi*details.frequency_hz*distance/details[f'{method}_velocity_m_s']+2*np.pi*details.constant_cycle
        err=np.max(np.abs(np.degrees(restored)-details[f'final_{method}_outward_unwrapped_phase_deg']))
        if err>=1e-8:raise ValueError('Phase/velocity inverse check failed')
        audit[f'{method}_inverse_roundtrip_error_deg']=float(err)
    sources=sorted(float(x) for x in details.source_phase_deg.unique())
    audit.update(curve_frequency_rows=len(details),unique_segment_frequencies=len(details.drop_duplicates(['segment_id','frequency_hz'])),
        original_constant_cycles_preserved=True,source_phase_deg=sources[0] if len(sources)==1 else sources,
        sources={name:hashlib.sha256((folder/name).read_bytes()).hexdigest() for name in ['phase_stack_summary.csv','selected_phase_targets.csv']},
        module_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        start_utc=metadata.get('start_utc'),end_utc_exclusive=metadata.get('end_utc_exclusive'),
        conditional_comparison_not_independent_validation=True)
    details.to_csv(output/'final_phase_velocity.csv',index=False)
    trends.to_csv(output/'velocity_trends.csv',index=False)
    comparison=details[['curve_id','segment_id','pair','branch','frequency_hz','constant_cycle','source_phase_deg',
        'distance_m','phase_velocity_m_s','raw_velocity_m_s','unit_velocity_m_s','raw_velocity_relative_change',
        'unit_velocity_relative_change','final_raw_wrapped_phase_deg','final_unit_wrapped_phase_deg',
        'raw_propagation_phase_deg','unit_propagation_phase_deg','accepted_windows','valid_windows']]
    comparison.to_csv(output/'curve_comparison.csv',index=False)
    (output/'verification.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    payload={'details':json.loads(details.to_json(orient='records')), 'trends':json.loads(trends.to_json(orient='records'))}
    (output/'workbook_data.json').write_text(json.dumps(payload,ensure_ascii=False),encoding='utf-8')
    font=Path(__file__).resolve().parent.parent/'assets/NotoSansSC.otf'
    figures(details,output,font,title,focus_curve_ids)
    doc='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>叠加相位与原候选曲线对比</title><style>body{max-width:1250px;margin:40px auto;padding:0 20px;font:16px/1.6 system-ui;color:#172033}img{max-width:100%}table{border-collapse:collapse;font-size:12px}td,th{padding:6px;border-bottom:1px solid #ddd}th{background:#f1f5f9}</style>'
    doc+='<h1>'+html.escape(title)+'：Goertzel叠加与原候选曲线</h1>'
    doc+=f'<p>UTC {html.escape(str(metadata.get("start_utc","未记录")))} 至 {html.escape(str(metadata.get("end_utc_exclusive","未记录")))}。{details.curve_id.nunique()}条曲线、{len(details)}个曲线频点，对应{audit["unique_segment_frequencies"]}个独立片段频点。范围来源{html.escape(str(metadata.get("range_source","未记录")))}，额外余量/容差{metadata.get("tolerance_deg","未记录")}°；窗口{metadata.get("repeat_cycles","未记录")}周期，滑动{metadata.get("slide_cycles","未记录")}周期。</p>'
    if metadata.get('frequency_method')=='matlab_bin':doc+='<p>本次与附件相同，Goertzel使用最近DFT频点。曲线频率为frequency_hz，实际变换频率为actual_frequency_hz；两者差值可在完整CSV中核查。相位与速度对应关系为条件近似；需要严格逐频对应时，选择曲线精确频率模式重算。</p>'
    if metadata.get('split_phase_branches'):doc+='<p>正、负互谱相位分别累加；此对比使用曲线预测相位同号分支。两分支相位及复数和均保存在上级phase_branch_summary.csv。</p>'
    doc+='<p>互谱C=Gtarget·conj(Greference)。幅值加权叠加相位arg(ΣC)，等权相位叠加arg(ΣC/|C|)。负时延分支先反号，沿升频方向展开，整段仅加一次360°整数偏移对齐原片段，每条曲线保留原N。传播相位D=源相位−展开相位+360°N，速度v=360°·f·d/D。源相位沿用原曲线，并非由本次叠加确定。</p>'
    doc+='<p>灰色虚线为原候选曲线，灰带为原片段相位范围对应的条件速度范围；棕色为幅值加权叠加反算，蓝色为等权相位叠加反算。筛选使用了原曲线的相位范围，因此曲线接近不能独立证明原频散正确，也不能消除整数周期歧义。滑动窗高度重叠，窗口数不是独立样本数。</p><h2>逐曲线趋势</h2>'
    view=trends[['curve_id','method','frequency_min_hz','frequency_max_hz','velocity_at_low_frequency_m_s','velocity_at_high_frequency_m_s','linear_slope_m_s_per_hz','decreasing_adjacent_fraction','increasing_steps','maximum_abs_relative_change_from_input']]
    doc+=view.to_html(index=False,float_format=lambda x:f'{x:.5f}')
    image_names=['fv_curves.png','final_phases.png']
    if (output/'fv_focus.png').exists():image_names.insert(0,'fv_focus.png')
    for name in image_names:
        encoded=base64.b64encode((output/name).read_bytes()).decode('ascii')
        doc+='<img src="data:image/png;base64,'+encoded+'">'
    doc+='<p>curve_comparison.csv：逐频原速度、两种叠加速度、变化幅度、最终相位、传播相位、台距及保留窗数。final_phase_velocity.csv保留完整字段，velocity_trends.csv记录趋势，verification.json记录所有逐窗复数求和和速度回代检查。</p></html>'
    (output/'report.html').write_text(doc,encoding='utf-8')
    return audit
