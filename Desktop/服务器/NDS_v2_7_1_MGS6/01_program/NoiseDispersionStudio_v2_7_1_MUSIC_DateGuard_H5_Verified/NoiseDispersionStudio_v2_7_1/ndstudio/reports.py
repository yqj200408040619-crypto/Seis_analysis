"""Portable offline report and figure generation for conditional structure groups."""
from pathlib import Path
import json,base64,html
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def studio_report(project):
 root=Path(project);out=root/'output';out.mkdir(exist_ok=True);figdir=out/'figures';figdir.mkdir(exist_ok=True)
 spec=json.loads((root/'project.json').read_text(encoding='utf-8-sig'));info={'project':spec.get('name',root.name),'conditional_only':True};models=[];pics=[]
 def table(df):return '<div class="table">'+df.head(40).to_html(index=False,border=0)+'</div>' if len(df) else '<p>无通过结果。</p>'
 from .inversion import inversion_current
 complete=inversion_current(root)
 if complete and (out/'conditional_models.json').exists():models=json.loads((out/'conditional_models.json').read_text(encoding='utf-8-sig'))['models']
 if not models:
  for name in ('11_selected_dispersion.png','12_vs_structure.png'):(figdir/name).unlink(missing_ok=True)
 if models:
  fig,axes=plt.subplots(1,2,figsize=(12,5));maxdepth=max(50,max(sum(m['thickness_m']) for m in models)*1.3)
  obs=pd.read_csv(out/'selected_observed_dispersion.csv')
  branch_color={cid:f'C{i%10}' for i,cid in enumerate(sorted(obs.curve_id.unique()))}
  branch_cycle={cid:int(q.constant_cycle.iloc[0]) for cid,q in obs.groupby('curve_id')}
  labelled=set()
  for m in models[:80]:
   depth=np.r_[0,np.cumsum(m['thickness_m']),maxdepth];xx=[];yy=[]
   for i,v in enumerate(m['vs_m_s']):xx.extend([v,v]);yy.extend([depth[i],depth[i+1]])
   cid=m['curve_ids'][0];label=f'N={branch_cycle[cid]}' if len(m['curve_ids'])==1 else m['case']
   axes[0].plot(xx,yy,alpha=.45,lw=1,color=branch_color[cid],label=label if label not in labelled else None)
   labelled.add(label)
  axes[0].invert_yaxis();axes[0].set(xlabel='Conditional Vs (m/s)',ylabel='Model depth (m)',title=f'{len(models)} compatible models; not a posterior')
  if len(labelled)<=12:axes[0].legend(fontsize=8)
  for cid,q in obs.groupby('curve_id'):
   one=q.drop_duplicates('frequency_hz').sort_values('frequency_hz')
   label=f'N={int(one.constant_cycle.iloc[0])}, {one.cohort.iloc[0]}, p{int(one.pair_index.iloc[0])}'
   axes[1].plot(one.frequency_hz,one.phase_velocity_m_s,'o',ms=3,alpha=.7,label=label,color=branch_color[cid])
  for mid,q in obs.groupby('model_id'):
   for cid,b in q.groupby('curve_id'):
    b=b.sort_values('frequency_hz');axes[1].plot(b.frequency_hz,b.predicted_phase_velocity_m_s,alpha=.3,lw=.8,color=branch_color[cid])
  axes[1].set(xlabel='Frequency (Hz)',ylabel='Phase speed (m/s)',title='Selected observed points / conditional forward curves')
  if len(branch_color)<=12:axes[1].legend(fontsize=8)
  fig.tight_layout();path=figdir/'10_structure_models.png';fig.savefig(path,dpi=150,bbox_inches='tight')
  fig.canvas.draw();renderer=fig.canvas.get_renderer()
  for ax,name in zip(axes,['12_vs_structure.png','11_selected_dispersion.png']):
   extent=ax.get_tightbbox(renderer).transformed(fig.dpi_scale_trans.inverted()).expanded(1.035,1.04)
   fig.savefig(figdir/name,dpi=150,bbox_inches=extent)
  plt.close(fig)
 if (figdir/'10_structure_models.png').exists() and models:pics.append(figdir/'10_structure_models.png')
 prep=json.loads((out/'prepare_complete.json').read_text(encoding='utf-8-sig')) if (out/'prepare_complete.json').exists() else {}
 review=pd.read_csv(out/'candidate_review.csv') if (out/'candidate_review.csv').exists() else pd.DataFrame()
 attempts=pd.read_csv(out/'inversion_attempts.csv') if (out/'inversion_attempts.csv').exists() else pd.DataFrame()
 summary={'project':spec.get('name',root.name),'accepted_windows':prep.get('accepted_windows'),
  'candidate_curves':len(review),'conditional_models':len(models),'inversion_complete':complete,'unique_physical_curve_or_structure_verified':False,
  'independent_solver_passed_models':sum(bool(m['independent_solver_pass']) for m in models)}
 from dispersion.common import write_json
 write_json(out/'studio_summary.json',summary)
 css='body{font:16px/1.75 system-ui,"Microsoft YaHei",sans-serif;color:#21364b;max-width:1200px;margin:32px auto;padding:0 24px}h1,h2{color:#174963}.notice{background:#edf5fa;border-left:5px solid #3980a0;padding:18px}.table{overflow:auto}table{border-collapse:collapse;font-size:12px}th,td{border:1px solid #dbe4ec;padding:7px}img{max-width:100%}code{background:#edf1f4;padding:2px 5px}'
 text=f'<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>频散与条件结构报告</title><style>{css}</style><h1>{html.escape(spec.get("name",root.name))}</h1><div class="notice">候选曲线{len(review)}条；通过设置条件的结构{len(models)}个。曲线周期、源相位、模态及共同1D路径假设均保留，不宣称唯一地下结构或地层置信区间。</div>'
 config=json.loads((root/'config.json').read_text(encoding='utf-8-sig'));trim=config.get('discard_first_s',600)
 text+=f'<h2>输入和处理</h2><p>HDF5原始数据以只读方式打开，按明确的起始UTC和原生采样率转换，不使用文件修改时间作为授时。缺少必需字段时阻止分析，映射可保存复用。统一对齐后去共同起点前{trim:g}秒；短孤立尖峰在端点可靠时局部修复，长异常/缺测局部掩码。高维或无分量顺序的数据需要明确映射。通用核心提供多路互相关、多峰FTAN、固定开发/验证窗、门窗扰动、分块、留一块和错时背景对照。</p>'
 text+='<p><a href="workflow_trace.json">全流程阶段摘要</a> · <a href="workflow_processing.log">完整处理日志</a> · <a href="workflow_settings_used.json">实际参数快照</a></p>'
 if (out/'geothermal/report.html').exists():
  text+='<p><a href="geothermal/report.html">地热地质融合与空间层析成果</a> · <a href="geothermal/tomography_paths.csv">层析射线路径交接表</a>。地热模块不会把共同1D条件结构平铺成二维Vs。</p>'
 if (out/'geothermal/tomography2d/tomography_2d_report.html').exists():
  text+='<p><a href="geothermal/tomography2d/tomography_2d_report.html">二维相速度层析报告</a> · <a href="geothermal/tomography2d/local_dispersion.csv">局部频散交接表</a>。二维相速度层析仍不是Vs深度结构。</p>'
 if (out/'music_front/report.html').exists():text+='<p><a href="music_front/report.html">前段MUSIC逐频点方向与投影速度</a>（前段不使用45°台对门槛）。</p>'
 if (out/'phase_guided/phase_stack_complete.json').exists():
  text+='<p><a href="phase_guided/phase_report.html">逐频相位筛选与叠加</a> · <a href="phase_guided/phase_continuous_bands.csv">剔除后的连续频带</a>。启用方向筛选时，反演候选从通过连续频带检查的叠加曲线交接，源相位、周期及背景检验状态仍沿用原曲线。</p>'
  if (out/'phase_guided/direction/report.html').exists():text+='<p><a href="phase_guided/direction/report.html">阵列方向、45°筛选与独立平面波联合估计</a>。联合MUSIC结果与CCF条件曲线分别展示；未经模态与背景确认，不自动当作已验证的结构输入。</p>'
 skipped=json.loads((out/'inversion_skipped.json').read_text(encoding='utf-8-sig')) if (out/'inversion_skipped.json').exists() else None
 if skipped:text+='<p>自动反演已跳过：'+html.escape(skipped['reason'])+'。这是有效的无通过结果，不表示处理仍在运行。</p>'
 elif not complete:text+='<p>反演尚未完成。检查点中的中间模型不计入本报告的通过结果。</p>'
 if spec.get('analysis_mode')=='curves_only':text+='<p>本项目从已生成的候选曲线CSV导入，没有在本项目重新测量波形。输入表的质量声明来自其生产流程，不能解释为本次新波形验证。</p>'
 else:text+='<p>完整波形图件和频散检验表见同目录 <a href="report.html">核心处理报告</a>；实际参数及文件指纹在缓存provenance.json中。</p>'
 text+='<h2>候选与组选择</h2><p>候选按背景通过、频带覆盖和其他路径兼容性排序。自动组合须来自同一时间组，组内曲线两两兼容；不同台对的传播方位按正/负支路与实际几何换算。整数N固定于每个连续段，不要求不同台对编号相同，不逐频点改N，不合并全天/夜间的重叠观测。自动排序不是概率；仍保留未进入反演的完整候选。</p>'+table(review)
 text+='<h2>反演设置与限制</h2><p>每个明确的曲线/模态组合分别做多起点、1–4个有限层加半空间的Rayleigh条件反演。固定Vp/Vs和密度可编辑；层厚与Vs有明确搜索范围。相速度和实测群速度参与拟合，每条频段等总权；稀疏选点仅避免扫描密度主导，不把相关频点当作独立样本。接受前检查全段相速度/群速度误差、优化收敛与数值根；安装pysurf96时再做独立求解器比较。未安装时明确标记未独立复验；GUI可要求它为必需，此时缺少即停止。</p><p>Vs剖面是所选周期、源相位、模态和1D结构假设下的兼容模型；模型数量、包络或拟合误差均不是正式后验概率和置信区间。有限频段、固定Vp/Vs、层数先验和横向变化限制深度分辨率；接近边界的参数单独标记。未得到兼容模型时保留全部失败原因，不自动放宽门槛。</p>'
 display=pd.DataFrame(models)
 if len(display):text+=table(display[['model_id','case','finite_layers','thickness_m','vs_m_s','phase_weighted_rms_fraction','maximum_group_median_fraction','independent_solver_pass','parameters_near_search_bound']])
 text+='<h2>优化与拒绝记录</h2>'+table(attempts)
 for path in pics:text+='<figure><img src="data:image/png;base64,'+base64.b64encode(path.read_bytes()).decode()+'"></figure>'
 text+='<h2>输出分级</h2><p>conditional_dispersion.csv为实测相位条件解释；selected_observed_dispersion.csv为已选择的观测及对应预测列；model_predictions.csv只存模型预测，缺根为空，不能填进观测空档。conditional_models.json保留每组结构与假设，inversion_attempts.csv记录失败。修改输入、配置、核心代码或选择后必须重新运行相应阶段。</p><h2>方法和依赖</h2><p><a href="https://doi.org/10.1111/j.1365-246X.2007.03374.x">Bensen等(2007)</a>环境噪声互相关/FTAN；<a href="https://docs.h5py.org/en/stable/">h5py官方文档</a>；<a href="https://keurfonluu.github.io/disba/">disba官方文档</a>；<a href="https://doc.qt.io/qtforpython-6/">Qt for Python</a>。具体筛选阈值为本项目的探索性设置。</p></html>'
 (out/'studio_report.html').write_text(text,encoding='utf-8')
 return out/'studio_report.html'

