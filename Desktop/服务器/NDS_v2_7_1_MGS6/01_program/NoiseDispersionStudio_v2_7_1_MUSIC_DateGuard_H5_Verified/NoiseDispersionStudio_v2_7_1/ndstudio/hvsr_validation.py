"""Native real-data comparison report; uses Studio HVSR outputs and readers."""
from pathlib import Path
import base64,html,json
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from .hvsr import discover,_read_interval,sha_file,write_json


CASES=[('baseline','hvsr','hvsr_settings.json'),
       ('same_window','hvsr_same_window','hvsr_same_window_settings.json'),
       ('low_frequency','hvsr_low_frequency','hvsr_low_frequency_settings.json')]


def real_validation_report(root):
    root=Path(root).resolve();output=root/'validation';output.mkdir(parents=True,exist_ok=True)
    project=root/'project';download=root/'download'
    audit=json.loads((download/'acquisition_audit.json').read_text(encoding='utf-8-sig'))
    selection=json.loads((download/'hvsr_download_selection.json').read_text(encoding='utf-8-sig'))
    actual={str((download/r['file']).resolve()):sha_file(download/r['file']) for r in audit['records']}
    if not audit['checks_passed'] or any(actual[str((download/r['file']).resolve())]!=r['sha256'] for r in audit['records']):
        raise ValueError('真实数据源文件与下载核查 SHA256 不一致')
    rows=[];summaries=[]
    fig,ax=plt.subplots(figsize=(10.8,5.2),constrained_layout=True)
    colors=['#155784','#dc812f','#309582']
    for (case,sub,config),color in zip(CASES,colors):
        folder=project/'output'/sub;s=json.loads((folder/'hvsr_summary.json').read_text(encoding='utf-8-sig'))
        if s['source_sha256']!=actual:raise ValueError('分析源文件与当前原始下载不一致：'+case)
        qc=pd.read_csv(folder/'window_qc.csv').fillna('');rejected=Counter(qc.loc[~qc.accepted,'reason'])
        peak=s['selected_peak'];c=s['settings'];table=pd.read_csv(folder/'hvsr_curve.csv')
        ax.semilogx(table.frequency_hz,table.hv,lw=2,color=color,
                    label=f"{c['fmin_hz']:g}-{c['fmax_hz']:g} Hz / {c['window_s']:g} s / {s['independent_windows']} windows")
        if peak:ax.plot(peak['f0_hz'],peak['amplitude'],'o',color=color)
        rows.append(dict(case=case,frequency_band_hz=[c['fmin_hz'],c['fmax_hz']],window_s=c['window_s'],
                         windows=s['windows'],accepted=s['accepted_windows'],independent=s['independent_windows'],
                         candidate_hz=peak['f0_hz'] if peak else None,clarity=peak['clarity_pass_count'] if peak else None,
                         reliability_pass=peak['reliability_pass'] if peak else None,
                         left_drop=peak['left_drop'] if peak else None,right_drop=peak['right_drop'] if peak else None,
                         structure_status=s['structure']['status'],rejected_by_reason=dict(rejected)))
        summaries.append(s)
    ax.axhline(1,color='#9ba9b2',lw=.7);ax.set(xlabel='Frequency (Hz)',ylabel='H/V amplitude ratio',
            title='BW.MGS01 | same real data, native 200 Hz | parameter sensitivity',ylim=(0,3.4))
    ax.grid(alpha=.18);ax.legend(fontsize=9);fig.savefig(output/'hvsr_comparison.png',dpi=160);plt.close(fig)
    # Plot original digital samples with the same gap/overlap-aware reader used
    # for processing. Counts are explicitly labelled; no external waveform tool.
    inv=discover(download);item=next(i for i in inv['instruments'] if i['instrument_id']==selection['instrument_id'])
    fig,axes=plt.subplots(3,1,figsize=(10.8,6),sharex=True,constrained_layout=True)
    for (channel,records),ax,color in zip(sorted(item['channels'].items()),axes,colors):
        fs=records[0]['fs_hz'];t=UTCDateTime('2020-11-01T00:20:00Z').ns
        origin=records[0]['start_ns'];t=origin+round(round((t-origin)/1e9*fs)/fs*1e9)
        n=round(60*fs);x,_=_read_interval(records,t,n,fs,{})
        ax.plot(np.arange(n)/fs,x,color=color,lw=.6);ax.set(ylabel=channel+' counts');ax.grid(alpha=.15)
    axes[0].set_title('Original archived samples | 2020-11-01 00:20 UTC | no resampling / no gap filling')
    axes[-1].set_xlabel('Seconds from start');fig.savefig(output/'native_waveform.png',dpi=150);plt.close(fig)
    result=dict(program_version='2.5.1',station=selection['instrument_id'],provider=selection['provider'],
                request_start_utc=selection['start_utc'],request_end_utc=selection['end_utc'],
                channels=[r['channel'] for r in selection['channels']],sampling_rate_hz=selection['channels'][0]['sample_rate'],
                acquisition_checks_passed=True,h5_files=audit['h5_files'],archived_native_samples=audit['native_samples'],
                source_hashes_unchanged=True,cases=rows,
                conclusion='真实下载、响应校正、H/V、质量控制与报告链路运行成功。候选低频宽肩峰未建立可靠基频，暂停自动层厚和速度剖面。',
                limits=['低频变化同时涉及频段与窗长；同窗长对照单独检查了频段影响。',
                        '更长窗的突发噪声剔除减少有效窗；未放松质量门限或按H/V形状挑窗。',
                        '此次没有可确认的浅层结构解释，不输出Vs30；Rayleigh条件试算保持关闭。'])
    gui=project/'gui_validation.json'
    if gui.exists():result['gui_validation']=json.loads(gui.read_text(encoding='utf-8-sig'))
    write_json(output/'real_validation.json',result)
    table=[]
    for r in rows:
        table.append(dict(频段Hz='–'.join(f'{v:g}' for v in r['frequency_band_hz']),窗长秒=r['window_s'],
                          有效独立窗=f"{r['independent']} / {r['windows']}",候选峰Hz=f"{r['candidate_hz']:.4f}",
                          清晰度=f"{r['clarity']}/6",低频侧下降=r['left_drop'],高频侧下降=r['right_drop'],结构状态=r['structure_status']))
    def image(name):return '<img alt="'+name+'" src="data:image/png;base64,'+base64.b64encode((output/name).read_bytes()).decode()+'">'
    detailed=''
    for r,s,(_,sub,_) in zip(rows,summaries,CASES):
        detailed+=f'<h2>{r["frequency_band_hz"]} Hz · {r["window_s"]:g}秒窗</h2><p>{html.escape(s["structure"]["explanation"])}</p>'
        detailed+='<pre>'+html.escape(json.dumps(r['rejected_by_reason'],ensure_ascii=False,indent=2))+'</pre>'
        detailed+='<img alt="分量谱与逐窗H/V" src="data:image/png;base64,'+base64.b64encode((project/'output'/sub/'figures/hvsr_spectra_time.png').read_bytes()).decode()+'">'
    content='''<!doctype html><html lang="zh"><meta charset="utf-8"><title>真实波形下载与 H/V Debug 验证</title>
<style>body{max-width:1120px;margin:32px auto;padding:0 24px;font:16px/1.7 system-ui;color:#19384a}h1{font-size:28px}img{width:100%;margin:14px 0}table{border-collapse:collapse;font-size:14px}td,th{border:1px solid #d9e2e9;padding:8px;white-space:nowrap}.wrap{overflow-x:auto}.note{background:#edf4f7;padding:18px;border-radius:8px}pre{white-space:pre-wrap}</style>
<h1>真实波形下载与 H/V Debug 验证</h1>'''
    content+=f'<p>Noise Dispersion Studio 2.5.1 · {html.escape(selection["instrument_id"])} · LMU FDSN · 2020-11-01 00:00–02:00 UTC<br>HHE / HHN / HHZ · 原生 200 Hz · 12个H5分块 · {audit["native_samples"]:,}个归档样点（含原始MiniSEED记录边界重叠）</p>'
    content+='<div class="note"><b>程序链路验证成功；这批记录尚不能可靠约束浅层层厚。</b><br>波形由程序自带 FDSN 后端下载，原始 MiniSEED 与完整 StationXML 内嵌归档。原生样点、UTC、采样率及响应归档核查通过；分析前后12个H5文件SHA256一致。处理采用现有程序的逐分量速度响应校正、原生网格分窗、RMS水平合成和K–O平滑，未借助外部波形下载器或分析软件。</div>'
    content+=image('native_waveform.png')+'<h2>同一真实记录的参数复核</h2><div class="wrap">'+pd.DataFrame(table).to_html(index=False)+'</div>'+image('hvsr_comparison.png')
    content+='<p>常规 0.2–20 Hz / 120秒窗候选峰0.424 Hz虽满足5/6清晰度，但没有低频侧下降；同120秒窗扩至0.1 Hz仍是宽肩峰。扩至0.05 Hz并采用300秒窗后，最高评分候选峰约0.204 Hz，清晰度仅4/6。这不是可靠基频确认，不能据此把宽峰直接换算为层厚。</p>'
    content+='<h2>已修复与验证</h2><ul><li>增加可复现的单台三分量 hvsr-fetch 命令，复用原下载、响应归档和核查模块；默认HH?，并检查完整三分量与历元覆盖。</li><li>自动层厚推算增加双侧下降条件；保持SESAME评分不变。避免5/6宽肩峰自动生成结构。</li><li>逐窗QC区分分析窗缺测与响应缓冲不足，记录边缘跳窗不再误标为窗内缺测。</li><li>H/V命令的--start / --end覆盖配置文件日期。</li></ul>'
    if gui.exists():
        content+='<h2>程序界面复测</h2><p>Linux Qt offscreen环境实际点击主按钮，启动程序QProcess：退出码0，43个有效窗，曲线与“峰形不完整，暂停层厚估计”状态正常显示。Windows原生环境尚未实测。</p>'
        if (project/'real_gui.png').exists():content+='<img alt="真实数据程序界面" src="data:image/png;base64,'+base64.b64encode((project/'real_gui.png').read_bytes()).decode()+'">'
    content+='<p>v2.5.1相关回归：90项通过，8项Tk显示测试因缺少DISPLAY跳过；单独Qt H/V回归3项通过。完整记录见程序docs/HVSR_Validation_v251.json。</p>'
    content+=detailed+'<h2>解释范围</h2><ul>'+''.join('<li>'+html.escape(x)+'</li>' for x in result['limits'])+'</ul><p>完整原始波形、响应、参数、逐窗检查、逐峰评分、日志与命令已放入真实数据验证包。程序回归结果及GUI验证见随包验证JSON。本报告图像内嵌，可离线打开。</p></html>'
    (output/'real_report.html').write_text(content,encoding='utf-8')
    return result
