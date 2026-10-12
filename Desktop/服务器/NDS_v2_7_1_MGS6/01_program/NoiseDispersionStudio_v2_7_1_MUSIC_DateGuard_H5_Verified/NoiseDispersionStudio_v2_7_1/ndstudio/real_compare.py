"""Resume the program's native acquisition/processing stages and compare matching bands."""
from pathlib import Path
import argparse,json,time,hashlib
import numpy as np
import pandas as pd
from dispersion.common import write_json,sha,load_config,ensure_ready
from .acquisition import audit_acquisition
from .pipeline import execute


def preprocess(root,expected_stations=None):
    project=root/'project'; state=root/'real_compare_checkpoint.json'
    acquisition=json.loads((root/'download/acquisition_status.json').read_text(encoding='utf-8-sig'))
    if not acquisition.get('complete'):raise ValueError('下载尚未完成；请先用程序fetch续传')
    old=json.loads(state.read_text(encoding='utf-8-sig')) if state.exists() else {'stages':[]}
    markers={'acquisition-audit':root/'download/acquisition_audit.json',
        'import':project/'hdf_import_audit.json',**{s:project/'output'/f'{s}_complete.json' for s in ['prepare','scan','controls','families']}}
    def check_stations():
        if expected_stations is not None:
            t=pd.read_csv(project/'stations.csv')
            if len(t)!=len(expected_stations) or {s.split('.')[1] for s in t.seed_id}!=set(expected_stations):raise ValueError('导入的物理台站与六站名单不一致')
    downstream_invalid=False
    for stage,path in markers.items():
        reusable=not downstream_invalid and any(r['stage']==stage and r['status']=='complete' for r in old['stages']) and path.exists()
        if reusable:
            try:
                if stage=='import':
                    records=json.loads(path.read_text(encoding='utf-8-sig'))['records']
                    if {Path(r['file']).resolve() for r in records}!={p.resolve() for p in (root/'download').rglob('*.h5')}:
                        raise ValueError('下载文件集合已改变')
                    for r in records:
                        if sha(r['file'])!=r['source_sha256'] or sha(project/'canonical'/r['canonical_file'])!=r['canonical_sha256']:
                            raise ValueError('原始或导入波形已改变')
                elif stage in ('prepare','scan','controls','families'):
                    ensure_ready(load_config(project/'config.json'),'prepare_complete.json' if stage=='prepare' else 'scan_complete.json')
            except (ValueError,KeyError,OSError) as e:
                print('【检查点失效】'+stage+' · '+str(e),flush=True);reusable=False
        if reusable:
            if stage=='import':check_stations()
            print('【检查点复用】'+stage,flush=True);continue
        downstream_invalid=True
        print('【真实全流程】'+stage,flush=True);start=time.time()
        try:
            if stage=='acquisition-audit':audit_acquisition(root/'download')
            else:execute(project,stage)
            if stage=='import':check_stations()
        except Exception as e:
            old['failed_stage']=stage;old['error']=str(e);write_json(state,old);raise
        old['stages'].append(dict(stage=stage,status='complete',elapsed_s=time.time()-start))
        old.update(complete=False,project=str(project),input_folder=str(root/'download'))
        old.pop('failed_stage',None);old.pop('error',None);write_json(state,old)
    return old


def completed_phase(folder,options):
    """A completion marker is reusable only for the actual requested calculation."""
    from .direction import settings
    marker=folder/'phase_stack_complete.json'
    if not marker.exists():return False
    recorded=json.loads(marker.read_text(encoding='utf-8-sig'))
    expected={k:options[k] for k in ('repeat_cycles','slide_cycles','taper','frequency_method','monotonic_only')}
    expected.update(direction=settings(options.get('direction')),prior=options.get('prior'))
    if options.get('reconstruction') or recorded.get('reconstruction'):
        from .reconstruction import settings as reconstruction_settings
        expected['reconstruction']=reconstruction_settings(options.get('reconstruction'))
    if any(recorded.get(k)!=v for k,v in expected.items()):
        raise ValueError('相位检查点参数与当前请求不同；请使用新的相位输出目录')
    if recorded.get('module_sha256')!=sha(Path(__file__).with_name('phase_stack.py')) or recorded.get('direction_sha256')!=sha(Path(__file__).with_name('direction.py')):
        return False
    if recorded.get('music_module_sha256')!=sha(Path(__file__).with_name('music.py')):return False
    if recorded.get('aligned_h5_sha256')!=sha(folder/'aligned_waveforms.h5'):
        raise ValueError('相位检查点波形已改变')
    if recorded.get('reconstruction_module_sha256') and recorded['reconstruction_module_sha256']!=sha(Path(__file__).with_name('reconstruction.py')):return False
    return True


def paired_inversion(table,selection,dest,scan):
    from .inversion import invert_tables,code_digest,selection_digest
    marker=dest/'inversion_complete.json'
    signature=hashlib.sha256(table.to_csv(index=False,float_format='%.17g').encode()).hexdigest()
    expected=dict(input_curves_sha256=signature,selection_sha256=selection_digest(selection),studio_code_sha256=code_digest(),
        minimum_points=scan['min_points'],minimum_band_hz=scan['min_band_hz'])
    previous=json.loads(marker.read_text(encoding='utf-8-sig')) if marker.exists() else {}
    if any(previous.get(k)!=v for k,v in expected.items()):
        invert_tables(table,selection,dest,scan['min_points'],scan['min_band_hz'])
        previous=json.loads(marker.read_text(encoding='utf-8-sig'));previous.update(expected);write_json(marker,previous)
    return previous


def stage_statistics(root):
    out=root/'project/output';q=pd.read_csv(out/'conditional_dispersion.csv')
    from .curve_prior import monotone_curves
    decreasing,_=monotone_curves(q)
    result=dict(scan_band_hz=[load_config(root/'project/config.json')['scan'][k] for k in ('min_hz','max_hz')],
        conditional_curves=int(q.curve_id.nunique()),nonincreasing_curves=int(decreasing.curve_id.nunique()),
        observed_phase_segments=json.loads((out/'scan_complete.json').read_text(encoding='utf-8-sig'))['segments'],
        measured_frequency_rows=json.loads((out/'scan_complete.json').read_text(encoding='utf-8-sig'))['measured_points'])
    result['preparation']=json.loads((out/'prepare_complete.json').read_text(encoding='utf-8-sig'))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workdir',type=Path,required=True)
    p.add_argument('--compare',action='store_true');p.add_argument('--control-only',action='store_true');p.add_argument('--fine-cycles',type=float,default=200.)
    p.add_argument('--low-frequency-hz',type=float,default=.2);p.add_argument('--legacy-phase',action='store_true')
    a=p.parse_args();root=a.workdir.resolve();preprocess(root)
    if a.compare:compare(root,a.fine_cycles,not a.legacy_phase,a.low_frequency_hz)
    if a.control_only:phase_control(root,a.fine_cycles)


def select_inputs(root):
    from .candidates import auto_selection
    project=root/'project';out=project/'output';spec=json.loads((project/'project.json').read_text(encoding='utf-8-sig'))
    frozen=root/'before_selection.json';saved=root/'before_selected_curves.csv';signature=root/'comparison_input_signature.json'
    current=dict(config_sha256=sha(project/'config.json'),curves_sha256=sha(out/'conditional_dispersion.csv'))
    if frozen.exists() and saved.exists():
        if signature.exists() and json.loads(signature.read_text(encoding='utf-8-sig'))!=current:
            raise ValueError('成对对比的原始配置/曲线已改变；请使用新比较目录，避免混用历史频段')
        if not signature.exists():
            original=pd.read_csv(out/'conditional_dispersion.csv');previous=pd.read_csv(saved)
            merged=previous.merge(original,on=['curve_id','frequency_hz'],suffixes=('_saved','_current'),validate='one_to_one')
            if len(merged)!=len(previous) or not np.allclose(merged.phase_velocity_m_s_saved,merged.phase_velocity_m_s_current,rtol=1e-12):
                raise ValueError('历史对比输入与当前原始候选不同')
            write_json(signature,current)
        selection=json.loads(frozen.read_text(encoding='utf-8-sig'));ids=sorted({b['curve_id'] for case in selection['cases'] for b in case['bands']})
        print('【冻结对比输入复用】',len(ids),'条曲线',flush=True);return ids,selection
    spec['exclude_increasing_velocity']=True;write_json(project/'project.json',spec)
    write_json(project/'inversion_settings.json',dict(require_independent_solver=True))
    selection=auto_selection(project,max_cases=12,modes=(0,),minimum_wavelengths=3)
    catalog=pd.read_csv(out/'curve_catalog.csv').set_index('curve_id')
    # Use one representative validated band (and preserve separate cycle
    # hypotheses), so the before/after comparison has a stable common scope.
    preferred=[case for case in selection['cases'] if all(
        catalog.loc[b['curve_id'],'fmax_hz']-catalog.loc[b['curve_id'],'fmin_hz']<=1.+1e-9 for b in case['bands'])]
    selection['cases']=(preferred or selection['cases'])[:3]
    if not selection['cases']:raise ValueError('没有同时通过背景检验、三倍波长和升频非增速的候选')
    write_json(root/'before_selection.json',selection)
    write_json(signature,current)
    ids=sorted({b['curve_id'] for case in selection['cases'] for b in case['bands']})
    curves=pd.read_csv(out/'conditional_dispersion.csv');chosen=curves[curves.curve_id.isin(ids)]
    chosen.to_csv(root/'before_selected_curves.csv',index=False)
    print('【同组对比输入】',len(ids),'曲线；',len(chosen),'频点；',chosen.frequency_hz.nunique(),'独立频率',flush=True)
    return ids,selection


def compare(root, fine_cycles=200., reconstruct_low=True, low_frequency_hz=.2):
    from .phase_stack import run_phase_stack
    from .inversion import invert_tables
    from .curve_prior import monotone_curves
    project=root/'project';out=project/'output';ids,selection=select_inputs(root)
    options=dict(curve_ids=ids,repeat_cycles=20.,slide_cycles=20.,tolerance_deg=15.,frequency_method='exact',
        taper='hann',monotonic_only=True,prior=dict(enabled=True,extension_hz=.1),
        direction=dict(mode='filter',max_pair_angle_deg=45.))
    coarse=out/'phase_coarse'
    if not (coarse/'aligned_waveforms.h5').exists():run_phase_stack(project,output=coarse,**options)
    options.update(repeat_cycles=fine_cycles,slide_cycles=fine_cycles,
        direction=dict(mode='filter',max_pair_angle_deg=45.,direction_step_deg=1.,velocity_points=256,
            frequency_half_width_hz=.025,frequency_bins=3))
    if reconstruct_low:options.update(prior=None,reconstruction=dict(enabled=True,lower_frequency_hz=low_frequency_hz,step_hz=.05))
    write_json(project/'phase_settings.json',options)
    phase=out/'phase_guided'
    if reconstruct_low and (phase/'phase_stack_complete.json').exists():
        if not json.loads((phase/'phase_stack_complete.json').read_text(encoding='utf-8-sig')).get('reconstruction',{}).get('enabled'):
            print('【检查点更新】旧相位结果未做投影与低频实测；复用校正波形重算相位。',flush=True);(phase/'phase_stack_complete.json').unlink()
    if not completed_phase(phase,options):
        run_phase_stack(project,aligned_h5=coarse/'aligned_waveforms.h5',**options)
    after=pd.read_csv(phase/'inversion_curves.csv');before=pd.read_csv(out/'conditional_dispersion.csv')
    final=root/'comparison';final.mkdir(exist_ok=True)
    if reconstruct_low:
        from .reconstruction import paired_compare
        cases=[]
        for case in selection['cases']:
            bands=[dict(b,curve_id=cid) for b in case['bands'] for cid in after.loc[after.parent_curve_id==b['curve_id'],'curve_id'].unique()]
            if all(any(after.loc[after.curve_id==x['curve_id'],'parent_curve_id'].eq(b['curve_id']).any() for x in bands) for b in case['bands']):cases.append(dict(case,bands=bands))
        summary=paired_compare(project,dict(selection,cases=cases),final);summary.update(fine_cycles=fine_cycles,stages=stage_statistics(root));write_json(final/'comparison_summary.json',summary)
        (final/'DARE_Full_Comparison.html').write_text((final/'report.html').read_text(),encoding='utf-8')
        state=json.loads((root/'real_compare_checkpoint.json').read_text(encoding='utf-8-sig'));state.update(complete=True,comparison=summary);write_json(root/'real_compare_checkpoint.json',state);return
    after.to_csv(final/'after_curves.csv',index=False)
    summary=dict(final_curves=int(after.curve_id.nunique()),final_points=len(after),fine_cycles=fine_cycles,
        same_frequency_and_band_comparison=bool(len(after)),extrapolated_points_used_in_inversion=False,
        all_final_curves_nonincreasing=True,source_phase_and_cycles_conditional=True,
        mode_assumed=0,unique_structure_identified=False,stages=stage_statistics(root))
    if len(after):
        matching=before.merge(after[['parent_curve_id','frequency_hz','curve_id']].rename(columns={'curve_id':'matched_id','parent_curve_id':'curve_id'}),on=['curve_id','frequency_hz'])
        matching['parent_curve_id']=matching.curve_id;matching['curve_id']=matching.pop('matched_id')
        matching.to_csv(final/'before_matched_curves.csv',index=False)
        cases=[]
        for case in selection['cases']:
            bands=[]
            for b in case['bands']:
                for cid in after.loc[after.parent_curve_id==b['curve_id'],'curve_id'].unique():bands.append({'curve_id':cid,'mode':b['mode']})
            if bands:cases.append(dict(case,bands=bands))
        paired=dict(selection,cases=cases);write_json(final/'paired_selection.json',paired)
        for name,table in [('before',matching),('after',after)]:
            summary[name+'_inversion']=paired_inversion(table,paired,final/name,load_config(project/'config.json')['scan'])
    else:
        matching=before.iloc[:0];matching.to_csv(final/'before_matched_curves.csv',index=False)
        # Preserve the original selected inversion as a diagnostic reference.
        chosen=before[before.curve_id.isin(ids)]
        summary['before_diagnostic_inversion']=paired_inversion(chosen,selection,final/'before_diagnostic',load_config(project/'config.json')['scan'])
        summary['after_inversion_skipped']='没有通过方向、连续频带与趋势检查的观测；禁止外延先验替代实测反演'
    write_json(final/'comparison_summary.json',summary)
    comparison_report(root,matching,after,summary)
    state=json.loads((root/'real_compare_checkpoint.json').read_text(encoding='utf-8-sig'));state.update(complete=True,comparison=summary)
    write_json(root/'real_compare_checkpoint.json',state)


def phase_control(root,cycles=200.):
    """Separate module control; cannot replace the direction-qualified product."""
    from .phase_stack import run_phase_stack
    from .inversion import invert_tables
    project=root/'project';selection=json.loads((root/'before_selection.json').read_text(encoding='utf-8-sig'))
    ids=sorted({b['curve_id'] for c in selection['cases'] for b in c['bands']})
    phase=project/'output/phase_control'
    options=dict(repeat_cycles=cycles,slide_cycles=cycles,taper='hann',frequency_method='exact',monotonic_only=True,direction={'mode':'off'})
    if not completed_phase(phase,options):
        run_phase_stack(project,output=phase,curve_ids=ids,repeat_cycles=cycles,slide_cycles=cycles,taper='hann',
            frequency_method='exact',monotonic_only=True,direction={'mode':'off'},
            aligned_h5=project/'output/phase_coarse/aligned_waveforms.h5')
    after=pd.read_csv(phase/'inversion_curves.csv');before=pd.read_csv(root/'before_selected_curves.csv')
    final=root/'phase_control_comparison';final.mkdir(exist_ok=True)
    matching=before.merge(after[['parent_curve_id','frequency_hz','curve_id']].rename(columns={'curve_id':'matched_id','parent_curve_id':'curve_id'}),on=['curve_id','frequency_hz'])
    matching['parent_curve_id']=matching.curve_id;matching['curve_id']=matching.pop('matched_id')
    after.to_csv(final/'after_curves.csv',index=False);matching.to_csv(final/'before_matched_curves.csv',index=False)
    cases=[]
    for case in selection['cases']:
        bands=[{'curve_id':cid,'mode':b['mode']} for b in case['bands'] for cid in after.loc[after.parent_curve_id==b['curve_id'],'curve_id'].unique()]
        if bands:cases.append(dict(case,bands=bands))
    paired=dict(selection,cases=cases);write_json(final/'paired_selection.json',paired)
    summary=dict(diagnostic_only=True,passed_final_direction_screen=False,final_curves=int(after.curve_id.nunique()),
        final_points=len(after),cycles=cycles,extrapolated_points_used_in_inversion=False,
        same_frequency_and_band_comparison=True,unique_structure_identified=False)
    if len(after):
        delta=after[['curve_id','frequency_hz','phase_velocity_m_s']].merge(matching[['curve_id','frequency_hz','phase_velocity_m_s']],on=['curve_id','frequency_hz'],suffixes=('_after','_before'),validate='one_to_one')
        summary['maximum_velocity_change_percent']=float(abs(delta.phase_velocity_m_s_after/delta.phase_velocity_m_s_before-1).max()*100)
    if len(after):
        for name,table in [('before',matching),('after',after)]:
            summary[name+'_inversion']=paired_inversion(table,paired,final/name,load_config(project/'config.json')['scan'])
    write_json(final/'comparison_summary.json',summary);comparison_report(root,matching,after,summary,final,True)


def comparison_report(root,before,after,summary,final=None,diagnostic=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import base64,html
    final=Path(final or root/'comparison');fig,axes=plt.subplots(1,2,figsize=(13,5),constrained_layout=True)
    if len(after):
        for cid,g in after.groupby('curve_id'):
            raw=before[before.curve_id==cid];axes[0].plot(raw.frequency_hz,raw.phase_velocity_m_s,'--',alpha=.6)
            axes[0].plot(g.frequency_hz,g.phase_velocity_m_s,'o-',ms=3,label=cid)
        axes[0].legend(fontsize=7)
        colors={'before':'#64748b','after':'#0d9488'}
        for name in ['before','after']:
            p=final/name/'conditional_models.json'
            if not p.exists():continue
            models=json.loads(p.read_text(encoding='utf-8-sig'))['models']
            for m in models:
                depths=np.r_[0,np.cumsum(m['thickness_m'])];depths=np.r_[depths,max(800.,depths[-1]*1.25)]
                axes[1].step(np.r_[m['vs_m_s'],m['vs_m_s'][-1]],depths,where='post',color=colors[name],alpha=.35)
        axes[1].plot([],[],color=colors['before'],label='Before');axes[1].plot([],[],color=colors['after'],label='After');axes[1].legend()
    else:
        q=pd.read_csv(root/'before_selected_curves.csv')
        for cid,g in q.groupby('curve_id'):axes[0].plot(g.frequency_hz,g.phase_velocity_m_s,label=cid)
        axes[0].text(.02,.98,'No qualified corrected band',transform=axes[0].transAxes,va='top',color='#b91c1c')
        p=final/'before_diagnostic/conditional_models.json'
        if p.exists():
            for m in json.loads(p.read_text(encoding='utf-8-sig'))['models']:
                depths=np.r_[0,np.cumsum(m['thickness_m'])];depths=np.r_[depths,max(800.,depths[-1]*1.25)]
                axes[1].step(np.r_[m['vs_m_s'],m['vs_m_s'][-1]],depths,where='post',color='#64748b',alpha=.25)
        axes[1].text(.02,.98,'Before: conditional models\nAfter: no qualified measured band',ha='left',va='top',transform=axes[1].transAxes,color='#b91c1c')
    axes[0].set(xlabel='Frequency (Hz)',ylabel='Conditional phase velocity (m/s)',title='Matched dispersion: dashed before, solid after' if len(after) else 'Original candidates; corrected bands rejected');axes[0].grid(alpha=.2)
    axes[1].set(xlabel='Vs (m/s)',ylabel='Depth (m)',title='Same bands, modes, seed and inversion settings' if len(after) else 'Original conditional models; no corrected inversion');axes[1].invert_yaxis();axes[1].grid(alpha=.2)
    if diagnostic:fig.suptitle('PHASE MODULE CONTROL ONLY - direction screen not passed',color='#b91c1c')
    fig.savefig(final/'Dispersion_Structure_Comparison.png',dpi=160);plt.close(fig)
    doc='<html lang="zh-CN"><meta charset="utf-8"><title>DARE全流程对比</title><style>body{font:16px/1.65 sans-serif;max-width:1150px;margin:35px auto}img{width:100%}pre{white-space:pre-wrap}table{border-collapse:collapse;font-size:13px}td,th{padding:6px;border:1px solid #ddd}</style><h1>DARE六台站96小时：频散与结构对比</h1>'
    doc+='<p>程序原生GFZ下载：2020-10-12至2020-10-16（UTC），A00、B02、C00、C01、C02、D02；台距0.93–3.83 km，平均最近邻间距1.43 km。96个H5下载分块审计通过。相位模块实际使用全部校正记录，长度96小时21.3秒。预处理、扫描、错时对照、周期枚举、Goertzel相位筛选、方向约束和条件结构反演均由本程序执行。</p>'
    doc+='<p>只保留升频非增速曲线。端点外延0.1 Hz只作为筛选先验，预测点不参与最终观测或反演。方向夹角上限45°；曲线源相位、整数周期及基阶Rayleigh模态仍为条件假设，结构不是唯一解。</p>'
    if diagnostic:doc+='<p style="color:#b91c1c;font-weight:bold">此页仅为相位模块控制试验，关闭方向门限；不是通过最终方向筛选的频散或结构结果。用于在相同频点和反演设置下检查相位筛选模块的作用。</p>'
    elif not len(after):
        doc+='<p style="color:#b91c1c;font-weight:bold">严格筛选最终保留0条连续修正频散，修正后结构反演跳过。下面灰色结构仅是原始曲线的条件相容模型，不是已验证的真实地层结构。</p>'
    if not diagnostic:
        doc+='<p>冻结输入为同一台对C02–D02的三条候选，2.45–3.15 Hz，各15个频点；源相位45°，整周数N=4或5。P00031与P00032来自不同测量门限，不是独立阵列证据。200周期Hann窗，1°×256慢度方向网格，邻频三点波束功率平均；相位仍在中心精确频率计算。最少3个独立且相位接受的窗口，最少6点且连续带宽0.25 Hz，所有质量阈值保留。</p>'
        raw=pd.read_csv(root/'before_selected_curves.csv');wavelength=raw.phase_velocity_m_s/raw.frequency_hz
        doc+=f'<p>由候选速度推算，输入条件波长约{wavelength.min()/1000:.2f}–{wavelength.max()/1000:.2f} km，明显短于最近邻间距。这提示空间混叠风险；加密数值扫描网格不能恢复缺失的空间采样。不能仅凭台距满足下载约束认定阵列适合该频段。</p>'
    doc+='<img src="data:image/png;base64,'+base64.b64encode((final/'Dispersion_Structure_Comparison.png').read_bytes()).decode()+'">'
    doc+='<pre>'+html.escape(json.dumps(summary,ensure_ascii=False,indent=2))+'</pre>'
    diag=root/'project/output/phase_guided/direction/frequency_diagnostics.csv'
    if diag.exists():doc+='<h2>方向检查原因</h2>'+pd.read_csv(diag).groupby('reason').windows.sum().reset_index().to_html(index=False)
    supplement=root/'supplementary_scan.json'
    if not diagnostic and supplement.exists():
        doc+='<h2>低频补充核验</h2><p>同一批原始H5重新执行导入、预处理、扫描、错时对照及候选枚举。扫描0.2–1.5 Hz，步长0.025 Hz；质量标准、最少6点、最小0.25 Hz连续带宽均保留。该频段没有通过验证的测量段，因此相位修正和结构反演没有观测输入。</p><pre>'+html.escape(supplement.read_text())+'</pre>'
    control=root/'phase_control_comparison'
    if not diagnostic and (control/'Dispersion_Structure_Comparison.png').exists():
        doc+='<h2>相位模块控制试验：未通过最终方向筛选</h2><p style="color:#b91c1c">下面仅比较相位筛选模块作用，关闭方向门限。该控制试验不能替代上面的严格筛选结果，曲线接近也不能独立证明真实频散。</p>'
        doc+='<img src="data:image/png;base64,'+base64.b64encode((control/'Dispersion_Structure_Comparison.png').read_bytes()).decode()+'">'
        doc+='<pre>'+html.escape((control/'comparison_summary.json').read_text())+'</pre>'
    doc+='</html>';(final/'DARE_Full_Comparison.html').write_text(doc,encoding='utf-8')


if __name__=='__main__':main()
