"""Expose contiguous robust runs and rejection bottlenecks after every scan."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from dispersion.common import csv,sha,write_json

RUN_COLUMNS=['rank','gate_id','cohort','pair','branch','fmin_hz','fmax_hz','points','bandwidth_hz',
             'required_points','required_bandwidth_hz','passes_minimum_points',
             'passes_minimum_bandwidth','all_point_checks_pass','accepted_continuous_band','reason']
GATE_COLUMNS=['rank','gate_id','cohort','pair','branch','fmin_hz','fmax_hz','evaluated_points',
              'robust_points','robust_fraction','longest_run_fmin_hz','longest_run_fmax_hz',
              'longest_run_points','longest_run_bandwidth_hz','main_bottlenecks']

def _bool(values):
    return values.map(lambda x:str(x).strip().lower() in ('true','1','yes')).to_numpy(bool)

def _runs(q,step):
    good=_bool(q.robust_point);f=q.frequency_hz.to_numpy(float);current=[];result=[]
    for i,ok in enumerate(good):
        if current and (not ok or abs(f[i]-f[current[-1]]-step)>max(1e-8,step*1e-5)):
            result.append(current);current=[]
        if ok:current.append(i)
    if current:result.append(current)
    return result

def write_scan_diagnostics(c,top_k=8):
    if not isinstance(top_k,int) or top_k<1:raise ValueError('Diagnostic top_k must be a positive integer')
    out=Path(c['output']);out.mkdir(parents=True,exist_ok=True)
    point_file=out/'point_diagnostics.csv';gate_file=out/'candidate_gates.csv'
    if not point_file.exists():raise ValueError('扫描逐频点诊断表不存在，请先执行scan')
    point=pd.read_csv(point_file);gates=pd.read_csv(gate_file) if gate_file.exists() else pd.DataFrame()
    provenance=out/'cache/provenance.json';scan_c=c;threshold_source='supplied_config'
    if provenance.exists():
        scan_c=json.loads(provenance.read_text(encoding='utf-8-sig'))['config'];threshold_source='preparation_provenance'
    marker=out/'scan_complete.json'
    if marker.exists():
        completed=json.loads(marker.read_text(encoding='utf-8-sig'))
        if 'config' in completed:
            for name,key in [('point_diagnostics.csv','point_diagnostics_sha256'),('candidate_gates.csv','candidate_gates_sha256')]:
                if completed.get(key) and sha(out/name)!=completed[key]:raise ValueError('诊断表与已完成扫描不一致，请重新扫描：'+name)
            scan_c=completed['config'];threshold_source='scan_completion'
    matches=scan_c==c;s=scan_c['scan'];v=scan_c['validation'];rows=[];gate_rows=[]
    metadata=gates.set_index('gate_id').to_dict('index') if len(gates) else {}
    checks=[('maximum_route_phase_change_rad','处理路线相位',v['phase_route_max_rad'],'max'),
            ('maximum_split_phase_difference_rad','时段相位',v['phase_split_max_rad'],'max'),
            ('gate_change_rad','门窗扰动相位',v['gate_phase_max_rad'],'max'),
            ('maximum_split_arrival_fraction','时段到时',v['arrival_fraction_max'],'max'),
            ('maximum_phase_vs_packet_delay_fraction','相位时延/波包到时',v['phase_delay_fraction_max'],'max'),
            ('narrow_arrival_fraction','窄频宽到时',v['arrival_fraction_max'],'max'),
            ('narrow_development_fraction','窄频宽开发到时',v['arrival_fraction_max'],'max'),
            ('block_phase_fraction','分块波包与相位联合比例',v['block_fraction_min'],'min'),
            ('block_packet_fraction','分块波包比例',v['block_fraction_min'],'min'),
            ('maximum_leaveout_phase_rad','留一块相位',v['leaveout_phase_max_rad'],'max'),
            ('maximum_leaveout_arrival_fraction','留一块到时',v['leaveout_arrival_fraction_max'],'max')]
    if len(point):
        for gid,q in point.groupby('gate_id',sort=True):
            q=q.sort_values('frequency_hz').reset_index(drop=True);m=metadata.get(gid,{})
            info={k:m.get(k,q[k].iloc[0] if k in q else '') for k in ('cohort','pair','branch')}
            runs=_runs(q,s['step_hz']);longest=max(runs,key=lambda ix:(q.frequency_hz.iloc[ix[-1]]-q.frequency_hz.iloc[ix[0]],len(ix)),default=[])
            for ix in runs:
                lo,hi=float(q.frequency_hz.iloc[ix[0]]),float(q.frequency_hz.iloc[ix[-1]])
                width=hi-lo;enough=len(ix)>=s['min_points'];wide=width>=s['min_band_hz']-1e-8
                reasons=[]
                if not enough:reasons.append(f"点数不足({len(ix)}<{s['min_points']})")
                if not wide:reasons.append(f"带宽不足({width:.3f}<{s['min_band_hz']:.3f} Hz)")
                rows.append(dict(gate_id=gid,**info,fmin_hz=lo,fmax_hz=hi,points=len(ix),bandwidth_hz=width,
                    required_points=s['min_points'],required_bandwidth_hz=s['min_band_hz'],passes_minimum_points=enough,
                    passes_minimum_bandwidth=wide,all_point_checks_pass=True,accepted_continuous_band=enough and wide,
                    reason='通过连续频段条件' if enough and wide else '；'.join(reasons)))
            counts={};metric_failures={}
            for column,label,limit,kind in checks:
                if column not in q:continue
                x=pd.to_numeric(q[column],errors='coerce').to_numpy(float)
                failed=(~np.isfinite(x))|((x>limit) if kind=='max' else (x<limit))
                counts[label]=int(failed.sum());metric_failures['failed_'+column]=int(failed.sum())
            for column,label in [('required_routes_pass','处理路线综合检验（含SNR/边界）'),('phase_sampling_pass','相位采样密度')]:
                if column in q:counts[label]=int((~_bool(q[column])).sum())
            bottlenecks=sorted(((n,k) for k,n in counts.items() if n),reverse=True)
            row=dict(gate_id=gid,**info,fmin_hz=float(q.frequency_hz.min()),fmax_hz=float(q.frequency_hz.max()),
                evaluated_points=len(q),robust_points=int(_bool(q.robust_point).sum()),robust_fraction=float(_bool(q.robust_point).mean()),
                longest_run_fmin_hz=float(q.frequency_hz.iloc[longest[0]]) if longest else None,
                longest_run_fmax_hz=float(q.frequency_hz.iloc[longest[-1]]) if longest else None,
                longest_run_points=len(longest),longest_run_bandwidth_hz=float(q.frequency_hz.iloc[longest[-1]]-q.frequency_hz.iloc[longest[0]]) if longest else 0.,
                main_bottlenecks='；'.join(f'{label}:{n}/{len(q)}点' for n,label in bottlenecks[:3]))
            row.update(metric_failures)
            gate_rows.append(row)
    rows.sort(key=lambda x:(-x['bandwidth_hz'],-x['points'],x['gate_id'],x['fmin_hz']))
    gate_rows.sort(key=lambda x:(-x['longest_run_bandwidth_hz'],-x['robust_fraction'],-x['evaluated_points'],x['gate_id']))
    for rank,row in enumerate(rows,1):row['rank']=rank
    for rank,row in enumerate(gate_rows,1):row['rank']=rank
    csv(out/'scan_band_diagnostics.csv',rows,RUN_COLUMNS)
    gate_columns=GATE_COLUMNS+sorted(set().union(*(set(x) for x in gate_rows))-set(GATE_COLUMNS)) if gate_rows else GATE_COLUMNS
    csv(out/'scan_candidate_diagnostics.csv',gate_rows,gate_columns)
    summary={'robust_contiguous_runs':len(rows),'accepted_continuous_runs':sum(x['accepted_continuous_band'] for x in rows),
        'candidate_gates':len(gate_rows),'threshold_source':threshold_source,'thresholds_match_current_config':matches,
        'minimum_points':s['min_points'],'minimum_bandwidth_hz':s['min_band_hz'],'frequency_step_hz':s['step_hz'],
        'rank_is_not_a_probability':True,'failure_counts_can_overlap':True,'block_phase_metric_also_requires_packet_pass':True,
        'source_sha256':{p.name:sha(p) for p in (point_file,gate_file) if p.exists()}}
    write_json(out/'scan_diagnostics_summary.json',summary)
    print(f"scan diagnostics: 连续稳健区间{len(rows)}个，通过频段条件{summary['accepted_continuous_runs']}个；要求≥{s['min_points']}点、≥{s['min_band_hz']:.3f} Hz",flush=True)
    if not matches:print('scan diagnostics: 当前配置与生成诊断表时的参数不同，摘要采用历史参数；修改后须重新运行相应阶段。',flush=True)
    for row in rows[:top_k]:
        print(f"  {row['gate_id']} | {row['pair']} {row['cohort']} {row['branch']} | {row['fmin_hz']:.3f}–{row['fmax_hz']:.3f} Hz | {row['points']}点，带宽{row['bandwidth_hz']:.3f} Hz | {row['reason']}",flush=True)
    print('scan diagnostics: 较好的候选门窗与主要限制（失败计数可重叠）：',flush=True)
    for row in gate_rows[:top_k]:
        print(f"  {row['gate_id']} | {row['pair']} {row['cohort']} {row['branch']} | {row['fmin_hz']:.3f}–{row['fmax_hz']:.3f} Hz | 稳健{row['robust_points']}/{row['evaluated_points']}点；最长连续{row['longest_run_points']}点、{row['longest_run_bandwidth_hz']:.3f} Hz | {row['main_bottlenecks'] or '无逐点失败项'}",flush=True)
    if not summary['accepted_continuous_runs'] and (out/'cohort_counts.csv').exists():
        cohorts=pd.read_csv(out/'cohort_counts.csv')
        print('scan diagnostics: 固定时段的基础可用性：',flush=True)
        for row in cohorts.head(top_k).itertuples():
            print(f"  {row.cohort} | 开发{row.train_windows}窗、验证{row.validation_windows}窗（各需≥{v['min_windows_each_split']}）；可用块{row.eligible_blocks}（需≥{v['min_eligible_blocks']}）",flush=True)
        if not gate_rows:print('scan diagnostics: 尚无候选门窗；训练波包频点与SNR可查看local_packet_peaks.csv。',flush=True)
    print('scan diagnostics: 完整表见',out/'scan_band_diagnostics.csv','和',out/'scan_candidate_diagnostics.csv',flush=True)
    return summary
