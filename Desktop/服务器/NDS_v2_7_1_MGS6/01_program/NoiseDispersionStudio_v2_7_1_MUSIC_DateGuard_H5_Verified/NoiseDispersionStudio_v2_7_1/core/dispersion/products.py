"""Offline Chinese report, scientific plots and fail-closed output checks."""
from pathlib import Path
import base64,html,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .common import ensure_ready,csv,write_json

def read_optional(path,columns=()):
    if not Path(path).exists(): return pd.DataFrame(columns=columns)
    return pd.read_csv(path)

def report(c):
    out,cache,meta=ensure_ready(c,"families_complete.json")
    stations=pd.read_csv(out/"stations.csv");pairs=pd.read_csv(out/"pairs.csv")
    cat=pd.read_csv(out/"measured_segment_catalog.csv");phase=pd.read_csv(out/"measured_phase_segments.csv")
    quality=pd.read_csv(out/"window_quality.csv");counts=pd.read_csv(out/"cohort_counts.csv")
    curves=pd.read_csv(out/"conditional_dispersion.csv");family=pd.read_csv(out/"curve_catalog.csv")
    control=read_optional(out/"timeshift_summary.csv");diag=pd.read_csv(out/"point_diagnostics.csv")
    figdir=out/"figures";figdir.mkdir(exist_ok=True);images=[]
    def save(fig,name):
        path=figdir/name;fig.savefig(path,dpi=150,bbox_inches="tight");plt.close(fig)
        from PIL import Image
        with Image.open(path) as im: im.load()
        images.append(path)
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    axes[0].scatter(stations.x_m,stations.y_m)
    for r in stations.itertuples(): axes[0].annotate(r.station,(r.x_m,r.y_m),xytext=(4,4),textcoords="offset points")
    axes[0].set(xlabel="East / local x (m)",ylabel="North / local y (m)",title="Actual station geometry");axes[0].axis("equal")
    axes[1].plot(quality.timestamp-quality.timestamp.min(),quality.clean_fraction,".")
    axes[1].axhline(c["min_clean_fraction"],color="C3",ls="--")
    axes[1].set(xlabel="Seconds since first scheduled window",ylabel="Effective common clean fraction",ylim=(-.03,1.05),title="Local masking and window retention")
    save(fig,"01_geometry_quality.png")
    fig,ax=plt.subplots(figsize=(10,max(3,len(counts)*.5)))
    y=np.arange(len(counts));ax.barh(y-.18,counts.train_windows,height=.35,label="Development");ax.barh(y+.18,counts.validation_windows,height=.35,label="Internal validation")
    ax.set_yticks(y,counts.cohort);ax.set(xlabel="Windows",title="Fixed cohorts: windows are shared across frequency");ax.legend();save(fig,"02_cohort_windows.png")
    fig,ax=plt.subplots(figsize=(11,max(3,min(16,len(cat)*.24))))
    for i,r in enumerate(cat.itertuples()):
        ax.plot([r.fmin_hz,r.fmax_hz],[i,i],lw=3)
    ax.set_yticks(np.arange(len(cat)),cat.segment_id+" "+cat.cohort+" p"+cat.pair_index.astype(str)+" "+cat.branch,fontsize=6)
    ax.set(xlabel="Frequency (Hz)",ylabel="Validated phase segments",title="Measured coverage; gaps remain unmeasured")
    if cat.empty: ax.text(.5,.5,"No band passed all waveform checks",transform=ax.transAxes,ha="center")
    save(fig,"03_measured_coverage.png")
    if len(cat) and len(curves):
        shown=cat.sort_values(["points","segment_id"],ascending=[False,True]).head(4)
        fig,axes=plt.subplots(2,2,figsize=(11,8));axes=axes.ravel()
        for ax,(_,r) in zip(axes,shown.iterrows()):
            q=curves[curves.segment_id==r.segment_id];ids=q.curve_id.unique()
            take=ids[np.unique(np.linspace(0,len(ids)-1,min(50,len(ids))).astype(int))]
            for cid in take:
                part=q[q.curve_id==cid];ax.plot(part.frequency_hz,part.phase_velocity_m_s,color="C0",alpha=.15,lw=.7)
            meas=phase[phase.segment_id==r.segment_id];ax.plot(meas.frequency_hz,meas.group_velocity_m_s,color="C3",ls=":",label="Envelope group speed")
            ax.set(title=f"{r.segment_id}: {len(ids)} cycle/source families",xlabel="Frequency (Hz)",ylabel="Conditional c (m/s)");ax.legend(fontsize=7)
        for ax in axes[len(shown):]: ax.axis("off")
        fig.suptitle("Representative segments; every family is in the CSV. No unique mode assumed.")
        fig.tight_layout();save(fig,"04_conditional_families.png")
    def table(df,limit=40):
        if df.empty: return "<p>无对应结果。</p>"
        return df.head(limit).to_html(index=False,border=0,escape=True)+("<p>完整表格见CSV文件。</p>" if len(df)>limit else "")
    summary=json.loads((out/"families_complete.json").read_text(encoding='utf-8-sig'));prep=json.loads((out/"prepare_complete.json").read_text(encoding='utf-8-sig'))
    modelblock="<p>未运行结构解释，或没有模型通过当前条件。</p>"
    model_current=(out/"inversion_complete.json").exists()
    if model_current:
        marker=json.loads((out/"inversion_complete.json").read_text(encoding='utf-8-sig'))
        if 'selection_sha256' in marker:
            try:
                from ndstudio.inversion import inversion_current
                model_current=inversion_current(out.parent)
            except ImportError:model_current=False
    if model_current and (out/"conditional_models.json").exists():
        modelinfo=json.loads((out/"conditional_models.json").read_text(encoding='utf-8-sig'));models=modelinfo["models"]
        modelblock="<p>当前条件兼容结构实例："+str(len(models))+"。这些不是唯一结构或地层置信区间；合成平面波示例中的结构仅检验数值运行，不能当作真实地下结构。</p>"
        if models:
            display=pd.DataFrame(models)[["model_id","case","thickness_m","vs_m_s","phase_weighted_rms_fraction","maximum_group_median_fraction","independent_solver_pass"]]
            modelblock+=table(display)
            fig,ax=plt.subplots(figsize=(7,5))
            maxdepth=max(50.,max(sum(m["thickness_m"]) for m in models)*1.4)
            for m in models[:20]:
                h=m["thickness_m"];vs=m["vs_m_s"];depth=np.r_[0.,np.cumsum(h),maxdepth]
                xx=[];yy=[]
                for i,v in enumerate(vs): xx += [v,v];yy += [depth[i],depth[i+1]]
                ax.plot(xx,yy,alpha=.6,label=m["model_id"])
            ax.invert_yaxis();ax.set(xlabel="Conditional Vs (m/s)",ylabel="Model depth (m)",title="Conditional structure examples; depth sensitivity not established")
            ax.legend(fontsize=7);save(fig,"05_optional_conditional_models.png")
    passes=int(control.background_threshold_exceeded.sum()) if len(control) else 0
    rejected={}
    for key in ("required_routes_pass","phase_sampling_pass"):
        if key in diag: rejected[key]=int((~diag[key].astype(bool)).sum())
    steps=[
        ("输入与几何", "读取任意台数台站、真实坐标、通道和采集UTC；保存文件指纹。", "防止把距离、通道或旧缓存误用于新数据。"),
        ("共同授时与裁剪", f"按UTC重采样到同一网格，统一删除共同记录前{c['discard_first_s']:g}秒；显式记录用户给定的钟差。", "避免逐台从首样点截取造成相位错位。"),
        ("局部异常处理", "只修复极短、端点可靠的孤立尖峰；较长异常、平线和缺测保留局部掩码。", "保留更多窗内有效片段，避免伪造长段波形。"),
        ("响应与单台预处理", "按NSLC和完整采集历元核对StationXML响应，统一物理量，带通滤波。", "区分仪器相位与传播相位；记录标称方向及未核实的现场极性。"),
        ("多路互相关", "运行绝对均值归一化、谱白化及one-bit挑战；在有效连续片段上计算线性互相关。", "压制强瞬态和窄带源，检查结果是否依赖特定处理参数。"),
        ("固定时段和状态", "预先划分开发/验证块；状态分类器仅用开发窗拟合，整段频带共用窗集。", "降低逐频点挑窗造成的选择偏差。"),
        ("多峰FTAN与追踪", "保留局部峰，按到时连续性追踪多个波包；固定门窗仅从开发集学习。", "减少次峰遗漏，允许多个传播状态并存。"),
        ("波形稳健性", "独立时段重复、参数对照、门窗扩缩/整体移动、频宽改变、分块和留一块验证。", "检验相位是否重复，排除峰混合和不稳定门窗。"),
        ("错时与重复性诊断", "台站一致错时、重新应用冻结分类器；整块重采样并整段对齐2π。", "检查背景分离和时间稳定性；诊断统计不当作校准p值或置信区间。"),
        ("条件频散与跨路径比较", "每个连续段枚举常数周期N和源相位；比较不同台对的重叠频段。", "保留多解；用额外路径约束，而不把缺测点补成观测。"),
        ("可选结构解释", "用户明确指定周期/模态组合，独立多起点反演并用第二套求解器核对。", "检验假设下的结构兼容性；模型预测单独保存。")]
    rows="".join(f"<tr><td>{html.escape(a)}</td><td>{html.escape(b)}</td><td>{html.escape(d)}</td></tr>" for a,b,d in steps)
    pics="".join(f'<figure><img src="data:image/png;base64,{base64.b64encode(p.read_bytes()).decode()}"><figcaption>{p.name}</figcaption></figure>' for p in images)
    text=f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>传统噪声频散提取报告</title>
<style>body{{font:16px/1.7 system-ui,"Microsoft YaHei",sans-serif;max-width:1180px;margin:35px auto;padding:0 22px;color:#1f2937}}h1,h2{{color:#123d52}}table{{border-collapse:collapse;width:100%;font-size:13px;display:block;overflow:auto}}th,td{{border:1px solid #d7e1e8;padding:7px;vertical-align:top}}th{{background:#edf4f7}}figure{{margin:24px 0}}img{{max-width:100%}}.notice{{padding:16px;background:#eff6ff;border-left:5px solid #3280a0}}code{{background:#f1f5f9;padding:2px 4px}}</style>
<h1>传统环境噪声频散提取：本次数据</h1>
<div class="notice">共同有效窗 {prep["accepted_windows"]} / {prep["scheduled_windows"]}；通过波形验证的相位片段 {len(cat)}；错时背景分离通过 {passes} 个片段；条件周期/源相位曲线 {summary["conditional_curves"]} 条。重叠状态与片段不是独立观测条数。这些输出不表示已确认唯一Rayleigh频散或地下结构。</div>
<p>开发/验证分界：{html.escape(meta["split_utc"])}；验证状态：{html.escape(c["validation"]["status"])}。初始共同裁剪 {c["discard_first_s"]} 秒。0.025Hz等扫描步长仅增加扫描密度，实际分辨能力受记录长度、FTAN频宽、几何和信噪比约束。</p>
<h2>完整处理流程与意义</h2><table><tr><th>步骤</th><th>处理</th><th>意义</th></tr>{rows}</table>
<h2>结果分级</h2><p><code>measured_phase_segments.csv</code>保存实测相位及包络群速度；<code>conditional_dispersion.csv</code>保存 c=2πfr/(φ−θ+2πN) 的条件解释，每个连续段固定一个N；<code>model_predictions.csv</code>只保存可选反演预测。源相位由配置给出（常用π/4），方向性噪声可能破坏这个假设。周期、极性、源相位与模态不通过优化相干度自动确定。</p>
<p>没有通过全部验证时，程序保留诊断表并输出空曲线表。不得为了得到曲线而自动降低阈值。缺失错时对照记为空值，不能记作零噪声。bootstrap范围是重复性描述，频点相关，不能据此给出正式地层置信区间。</p>
<h2>测得相位覆盖</h2>{table(cat)}<h2>错时对照</h2>{table(control)}
<h2>可选结构解释</h2>{modelblock}
<h2>仪器与授时审计</h2>{table(pd.read_csv(out/"record_audit.csv"))}
<h2>图件</h2>{pics}<p>候选图只显示最长的最多4个片段，每段最多50条曲线用于阅读；完整多解见CSV。图中群速度不等于相速度，也不能直接当作Vs。</p>
<h2>复用时需要修改</h2><p>更换stations.csv、原始波形与采集日期匹配的响应文件；按采样率、台间距和预期频段设置带通、到时搜索、晚到噪声区和分析窗。新的验证数据须在调整参数前划分；本次验证类型以配置中的validation.status为准，内部时段重复不等于未触碰的外部验证。线阵也可输入本互相关流程，但本程序不输出F–K二维传播方向，也不把阵列投影速度自动当作真实相速度。</p>
<h2>方法来源</h2><p><a href="https://doi.org/10.1111/j.1365-246X.2007.03374.x">Bensen et al. (2007)</a>：环境噪声预处理、互相关、FTAN及稳定性质量控制；<a href="https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.remove_response.html">ObsPy响应文档</a>；<a href="https://docs.fdsn.org/projects/stationxml/en/latest/reference.html">FDSN StationXML方向/响应定义</a>。当前阈值及多解筛选为本项目的探索性设置，不能宣称为通用统计标准。</p></html>'''
    (out/"report.html").write_text(text,encoding="utf-8")

def verify(c):
    out,cache,meta=ensure_ready(c,"families_complete.json");ph=pd.read_csv(out/"measured_phase_segments.csv");cat=pd.read_csv(out/"measured_segment_catalog.csv")
    curves=pd.read_csv(out/"conditional_dispersion.csv");point=pd.read_csv(out/"point_diagnostics.csv");pairs=pd.read_csv(out/"pairs.csv").set_index("pair_index")
    z=np.load(cache/"window_routes.npz");shape=z["ram5"].shape
    assert shape[1]==len(pairs) and shape[-1]==2*round(c["max_lag_s"]*c["sample_rate_hz"])+1
    assert len(z["training"])==shape[0] and z["training"].any() and (~z["training"]).any()
    for r in c["routes"]: assert z[r["name"]].shape==shape and np.isfinite(z[r["name"]]).all()
    if len(ph):
        join=ph.merge(point[["gate_id","frequency_hz","robust_point"]],on=["gate_id","frequency_hz"],how="left",validate="many_to_one")
        assert join.robust_point.eq(True).all()
        assert np.isfinite(ph[["measured_unwrapped_phase_rad","arrival_s","group_velocity_m_s"]]).all().all()
        assert not ph.calibrated_confidence_interval.any() and not ph.absolute_cycle_resolved.any() and not ph.rayleigh_mode_verified.any()
        for sid,q in ph.groupby("segment_id"):
            q=q.sort_values("frequency_hz");assert len(q)==int(cat.set_index("segment_id").loc[sid,"points"])
            assert np.allclose(np.diff(q.frequency_hz),c["scan"]["step_hz"],atol=1e-8)
            assert q.cohort.nunique()==1 and q.pair_index.nunique()==1
    if len(curves):
        assert not curves.is_model_prediction.any() and curves.conditional_interpretation.all()
        assert not curves.absolute_cycle_resolved.any() and not curves.rayleigh_mode_verified.any()
        merged=curves.merge(ph[["segment_id","frequency_hz","measured_unwrapped_phase_rad"]],on=["segment_id","frequency_hz"],validate="many_to_one")
        r=merged.pair_index.map(pairs.distance_m.to_dict()).to_numpy();kr=merged.source_phase_rad-merged.measured_unwrapped_phase_rad+2*np.pi*merged.constant_cycle
        expected=2*np.pi*merged.frequency_hz*r/kr
        assert np.allclose(expected,merged.phase_velocity_m_s,rtol=1e-10)
        assert np.isfinite(curves.phase_velocity_m_s).all() and (curves.phase_velocity_m_s>0).all()
        for _,q in curves.groupby("curve_id"):
            assert q.constant_cycle.nunique()==q.cohort.nunique()==q.source_phase_rad.nunique()==1
    if (out/"timeshift_summary.csv").exists():
        control=pd.read_csv(out/"timeshift_summary.csv")
        if len(control):
            assert control.loc[control.valid_trials==0,"null_q95"].isna().all()
            assert not control.calibrated_p_value.any()
            assert (control.loc[control.background_threshold_exceeded,"valid_trials"]>=20).all()
    if (out/"model_predictions.csv").exists():
        predicted=pd.read_csv(out/"model_predictions.csv");assert not predicted.is_observation.any()
        if len(predicted): assert predicted.loc[~predicted.root_available,"phase_velocity_m_s"].isna().all()
    result=dict(checks_passed=True,windows=shape[0],pairs=len(pairs),measured_segments=len(cat),conditional_curves=curves.curve_id.nunique(),
                empty_curves_are_valid_outcome=True,verified_unique_physical_dispersion=False,raw_files_unchanged_and_fingerprint_checked=True)
    write_json(out/"verification.json",result);print("verify:",json.dumps(result),flush=True)
    return result
