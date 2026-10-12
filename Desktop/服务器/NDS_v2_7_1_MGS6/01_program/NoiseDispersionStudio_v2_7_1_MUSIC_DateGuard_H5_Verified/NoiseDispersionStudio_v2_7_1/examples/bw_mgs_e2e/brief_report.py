"""Short, standalone report embedding figures produced by Studio itself."""
import base64
import csv
import html
import json
from pathlib import Path

def make_report(workdir):
    root=Path(workdir);out=root/'project/output'
    def read(path,default=None):
        return json.loads(path.read_text(encoding='utf-8-sig')) if path.exists() else default
    selection=read(root/'array_selection.json')
    prep=read(out/'prepare_complete.json',{})
    scan=read(out/'scan_complete.json',{})
    summary=read(out/'studio_summary.json',{})
    modelinfo=read(out/'conditional_models.json',{'models':[]})
    geometry=selection['geometry']
    text='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>BW.MGS 频散与结构反演</title>
<style>body{font:16px/1.7 system-ui,sans-serif;color:#243746;max-width:1150px;margin:36px auto;padding:0 24px}h1,h2{color:#164c60}img{width:100%}table{border-collapse:collapse;width:100%;font-size:14px}td,th{border:1px solid #d7e0e7;padding:8px}th{background:#eef5f8}.note{padding:15px;background:#edf5fa}</style>
<h1>BW.MGS01–MGS05：频散与结构反演</h1>'''
    text+=f"<p>LMU FDSN；{html.escape(selection['start_utc'])} 至 {html.escape(selection['end_utc'])}（UTC）。五站原生 HHZ 200 Hz；台对距离 {geometry['min_pair_distance_m']/1000:.2f}–{geometry['max_pair_distance_m']/1000:.2f} km。</p>"
    text+=f"<p>全部数据下载、响应校正、互相关、多峰FTAN、相位测量、条件频散及Rayleigh反演均由 Noise Dispersion Studio 的现有模块完成。有效窗 {prep.get('accepted_windows',0)} / {prep.get('scheduled_windows',0)}；通过波形验证的片段 {scan.get('segments',0)}；条件曲线 {summary.get('candidate_curves',0)}；兼容结构 {summary.get('conditional_models',0)}。</p>"
    text+='<p>本次最低连续带宽由0.10 Hz调整为0.05 Hz，逐点相位稳定性门槛保持原值；该调整使用了本次数据，验证标记为内部复用。默认距离/波长门槛3得到14个兼容模型；图中补充展示门槛2的探索性对照，保留N=2和N=3两种周期分支。</p>'
    figures=[('频散曲线：实测条件分支与结构正演','11_selected_dispersion.png'),('Vs–深度结构：全部兼容模型','12_vs_structure.png')]
    for title,name in figures:
        text+='<h2>'+title+'</h2>'
        path=out/'figures'/name
        if path.exists():
            text+='<img alt="'+title+'" src="data:image/png;base64,'+base64.b64encode(path.read_bytes()).decode()+'">'
        else:
            text+='<p>当前数据与检验条件下无通过结果。完整拒绝原因保存在项目诊断表中。</p>'
    models=modelinfo['models']
    if models:
        cases={}
        for m in models:
            if m['case'] not in cases or m['phase_weighted_rms_fraction']<cases[m['case']]['phase_weighted_rms_fraction']:cases[m['case']]=m
        text+='<h2>每个周期分支的拟合实例</h2><p>下表仅列每组相速度RMS最小的实例，便于查看。图中和CSV保留全部兼容模型；拟合较小不代表结构较可信。</p><table><tr><th>模型</th><th>周期分支</th><th>有限层厚度（m）</th><th>各层与半空间 Vs（m/s）</th><th>相速度 RMS 相对误差</th></tr>'
        for m in cases.values():
            cycle=m['curve_ids'][0].rsplit('_N',1)[-1]
            text+='<tr>'+''.join('<td>'+html.escape(v)+'</td>' for v in [m['model_id'],'N='+cycle,', '.join(f'{x:.1f}' for x in m['thickness_m']),', '.join(f'{x:.1f}' for x in m['vs_m_s']),f"{100*m['phase_weighted_rms_fraction']:.2f}%"] )+'</tr>'
        text+='</table>'
        with (out/'selected_observed_dispersion.csv').open(encoding='utf-8-sig',newline='') as handle:
            rows=list(csv.DictReader(handle))
        branches={}
        for row in rows:branches.setdefault(row['curve_id'],{})[float(row['frequency_hz'])]=float(row['phase_velocity_m_s'])
        text+='<h2>本频段的速度变化</h2><table><tr><th>周期分支</th><th>0.270 Hz相速度（m/s）</th><th>0.340 Hz相速度（m/s）</th><th>随频率的变化</th></tr>'
        for cid,values in branches.items():
            low,high=values[min(values)],values[max(values)]
            text+=f'<tr><td>N={html.escape(cid.rsplit("_N",1)[-1])}</td><td>{low:.2f}</td><td>{high:.2f}</td><td>{"下降" if high<low else "上升"}</td></tr>'
        text+='</table><p>本次一个台对的窄频段同时支持不同整数周期解释；不能仅以“低频速度较高”来判定唯一分支或地层。深度分辨率尚未建立。</p>'
    diagnostic=out/'scan_band_diagnostics.csv'
    if diagnostic.exists():
        with diagnostic.open(encoding='utf-8-sig',newline='') as handle:bands=list(csv.DictReader(handle))
        text+='<h2>扫描后的连续稳健区间</h2><table><tr><th>台对 / 时段</th><th>频率范围（Hz）</th><th>点数</th><th>带宽（Hz）</th><th>结果</th></tr>'
        for row in bands[:8]:
            values=[row['pair']+' / '+row['cohort'],f"{float(row['fmin_hz']):.3f}–{float(row['fmax_hz']):.3f}",row['points'],f"{float(row['bandwidth_hz']):.3f}",row['reason']]
            text+='<tr>'+''.join('<td>'+html.escape(value)+'</td>' for value in values)+'</tr>'
        text+='</table><p>逐频点检验见point_diagnostics.csv；连续区间见scan_band_diagnostics.csv；候选门窗的主要限制与失败计数见scan_candidate_diagnostics.csv。失败计数可重叠，候选排序不代表通过概率。</p>'
    text+='<p class="note">频散速度和 Vs 剖面保留源相位、整数周期、Rayleigh 模态及共同1D结构假设；这些是条件兼容结果。反演使用程序内置 disba 求解路径；本环境未安装可选的独立 pysurf96 求解器，因此本次结果未做第二求解器复验。频散图中的包络群速度不等于相速度或 Vs。</p></html>'
    destination=root/'BW_MGS_Brief_Report.html';destination.write_text(text,encoding='utf-8')
    return destination

if __name__=='__main__':
    import sys
    print(make_report(sys.argv[1]))
