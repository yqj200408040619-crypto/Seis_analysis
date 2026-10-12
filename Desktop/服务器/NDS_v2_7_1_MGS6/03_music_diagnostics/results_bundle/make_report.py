#!/usr/bin/env python3
"""Build the Chinese report and a compact reproducibility bundle."""
from __future__ import annotations

import hashlib
import json
import sys
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent / "project"
WORKSPACE = HERE.parents[1]
SOURCE = WORKSPACE / "source_2023" / "NoiseDispersionStudio_v2_7_1"
sys.path[:0] = [str(SOURCE), str(SOURCE / "core")]
from dispersion.common import sha
from ndstudio.music_front import current as music_front_current


def read_csv(name):
    return pd.read_csv(HERE / name)


def fmt(x, digits=1):
    if pd.isna(x):
        return "—"
    return f"{float(x):,.{digits}f}"


def md_table(frame, columns, headers=None, digits=1):
    headers = headers or columns
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"]*len(headers)) + "|"]
    for row in frame[columns].itertuples(index=False, name=None):
        vals = []
        for column, value in zip(columns, row):
            if pd.isna(value):
                vals.append("—")
            elif isinstance(value, (float, np.floating)):
                if column.endswith("_hz"):
                    vals.append(fmt(value, 3))
                elif column == "min_coherence":
                    vals.append(fmt(value, 3))
                elif column.endswith("_fraction") or column.endswith("_rate") or column in {"primary_share", "bootstrap_positive_fraction"}:
                    vals.append(fmt(100 * value, 2) + "%")
                else:
                    vals.append(fmt(value, digits))
            else:
                vals.append(str(value))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def figure_bootstrap():
    d = read_csv("low_frequency_joint_bootstrap_summary.csv")
    f = d.frequency_hz.to_numpy(float)
    y = d.original_velocity_m_s.to_numpy(float)
    lo = y - d.velocity_q025_m_s.to_numpy(float)
    hi = d.velocity_q975_m_s.to_numpy(float) - y
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    ax.errorbar(f, y, yerr=np.vstack([lo, hi]), fmt="o", ms=6, capsize=4,
                color="#205493", ecolor="#6c8ebf", lw=1.8)
    ax.set(xlabel="Frequency (Hz)", ylabel="Measured phase velocity (m/s)",
           title="0.08–0.13 Hz: frequency-specific estimates and block-bootstrap intervals")
    ax.grid(alpha=.25)
    fig.text(.12, .015, "Windows are selected separately by frequency; bars are marginal 95% intervals. "
             "Cross-frequency paired differences are analyzed separately; no curve is fitted.", fontsize=8)
    fig.tight_layout(rect=(0, .045, 1, .94))
    path = HERE / "low_frequency_joint_bootstrap.png"
    fig.savefig(path, dpi=170)
    plt.close(fig)
    return path


def figure_threshold():
    d = read_csv("coherence_threshold_tradeoff.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
    for ax, f in zip(axes, sorted(d.frequency_hz.unique())):
        g = d[d.frequency_hz == f].sort_values("min_coherence")
        ax.plot(g.min_coherence, g.end_to_end_pure_noise_false_accept_rate,
                "o-", label="Pure-noise false accepts", color="#b22222")
        ax.plot(g.min_coherence, g.end_to_end_single_wave_snr0_recall,
                "s-", label="Single-wave recall (0 dB)", color="#176b3a")
        ax.set(title=f"{f:.2f} Hz", xlabel="Minimum coherence")
        ax.set_ylim(-.05, 1.05)
        ax.grid(alpha=.25)
        for t in [.98, .999]:
            hit = g[np.isclose(g.min_coherence, t)]
            if len(hit):
                ax.annotate(f"{t:g}", (t, hit.end_to_end_single_wave_snr0_recall.iloc[0]),
                            xytext=(3, 4), textcoords="offset points", fontsize=7)
    axes[0].set_ylabel("Rate")
    axes[0].legend(fontsize=8, loc="lower left")
    fig.suptitle("Coherence cutoff sensitivity under known-answer simulation")
    fig.tight_layout()
    path = HERE / "coherence_cutoff_tradeoff.png"
    fig.savefig(path, dpi=170)
    plt.close(fig)
    return path


def figure_direction():
    d = read_csv("real_multipeak_frequency_summary.csv").sort_values("frequency_hz")
    fig, ax = plt.subplots(figsize=(7.8, 4.4))
    sizes = 45 + 120*d.primary_share.to_numpy(float)
    ax.scatter(d.frequency_hz, d.median_azimuth_deg, s=sizes, color="#7851a9", alpha=.9)
    for row in d.itertuples():
        ax.annotate(f"{row.primary_windows}/{row.primary_windows_total}",
                    (row.frequency_hz, row.median_azimuth_deg), xytext=(4, 5),
                    textcoords="offset points", fontsize=8)
    ax.set(xlabel="Frequency (Hz)", ylabel="Dominant sampled primary direction (°)",
           title="Low-frequency primary peaks across five sampled windows per frequency")
    ax.grid(alpha=.25)
    ax.text(.01, .02, "Point size and label show the dominant peak share; no direction curve was fitted.",
            transform=ax.transAxes, fontsize=8, va="bottom")
    fig.tight_layout()
    path = HERE / "low_frequency_peak_tracking.png"
    fig.savefig(path, dpi=170)
    plt.close(fig)
    return path


def build_report():
    sim = read_csv("simulation_summary.csv")
    noise = sim[sim.case == "pure_noise"]
    wave = sim[(sim.case == "single_wave") & sim.snr_db.isin([-5, 0, 10, 20])]
    gain = sim[sim.case == "single_wave_gain_imbalance"]
    dual = sim[sim.case == "two_independent_narrowband_waves"]
    coherent = sim[sim.case == "two_coherent_tones"]
    ab = read_csv("nw_snapshot_ablation_summary.csv")
    ab_noise = ab[ab.case == "pure_noise"]
    leave = read_csv("real_leave_one_station_out.csv")
    leave_low = leave[leave.label.str.startswith("0.100")].copy()
    leave_low = leave_low[leave_low.dropped_station != "none"]
    grid = read_csv("real_grid_refinement.csv")
    boot = read_csv("low_frequency_joint_bootstrap_summary.csv")
    diffs = read_csv("low_frequency_adjacent_differences.csv")
    paired_diffs = read_csv("low_frequency_common_window_block_deltas.csv")
    common_summary = read_csv("low_frequency_common_window_summary.csv")
    rejection = read_csv("rejection_condition_summary.csv")
    reasons = read_csv("rejection_reason_reproduction.csv")
    lower_sim = read_csv("coherence_lower_floor_simulation_summary.csv")
    real_floor = read_csv("real_coherence_floor_counterfactual.csv")
    front_floor = read_csv("coherence_front_counterfactual_summary.csv")
    front_refs = read_csv("coherence_front_counterfactual_references.csv")
    hold = read_csv("holdout_date_validation.csv")
    hold_base = hold[(hold.NW == 3.5) & (hold.snapshots == 6)]
    h01 = hold_base[(np.isclose(hold_base.frequency_hz,.10)) & (hold_base.min_coherence == .70)]
    h045 = hold_base[(np.isclose(hold_base.frequency_hz,.45)) & (hold_base.min_coherence == .70)]
    peaks = read_csv("real_multipeak_frequency_summary.csv")
    bootdraws = read_csv("low_frequency_joint_bootstrap_draws.csv")
    complete = bootdraws[bootdraws.all_six_valid]

    marker_music = json.loads((PROJECT / "output/music_front/music_complete.json").read_text())
    marker_phase = json.loads((PROJECT / "output/phase_guided/phase_stack_complete.json").read_text())
    integrity = dict(baseline_music_current=bool(music_front_current(PROJECT)),
        baseline_music_marker_sha256=sha(PROJECT/"output/music_front/music_complete.json"),
        baseline_phase_marker_sha256=sha(PROJECT/"output/phase_guided/phase_stack_complete.json"),
        aligned_waveform_sha256=sha(PROJECT/"output/music_front/aligned_waveforms.h5"),
        music_settings_sha256=sha(PROJECT/"output/music_front/music_settings_used.json"),
        phase_inversion_handoff_points=marker_phase.get("inversion_handoff_points"),
        phase_accepted_window_evaluations=marker_phase.get("accepted_window_evaluations"),
        phase_valid_window_evaluations=marker_phase.get("valid_window_evaluations"))
    (HERE / "baseline_integrity.json").write_text(json.dumps(integrity,ensure_ascii=False,indent=2),encoding="utf-8")

    fig1 = figure_bootstrap()
    fig2 = figure_threshold()
    fig3 = figure_direction()

    real = md_table(leave_low,
        ["dropped_station","accepted","windows","paired_median_azimuth_shift_deg","paired_median_velocity_change_fraction"],
        ["Dropped station","Accepted windows / 32","N","Median paired azimuth shift (°)","Median paired velocity change (%)"])
    real += "\n\n" + md_table(grid,
        ["frequency_hz","windows","coarse_accepted","fine_accepted","paired_q1_windows","paired_median_azimuth_difference_deg","paired_median_velocity_change_fraction"],
        ["f (Hz)","N","2° grid accepted","1° grid accepted","Paired single-source windows","Median direction difference (°)","Median velocity change (%)"])

    low = boot.copy()
    low = low[["frequency_hz","original_azimuth_deg","azimuth_q025_deg","azimuth_q975_deg",
               "original_projected_distance_m","projection_q025_m","projection_q975_m",
               "original_velocity_m_s","velocity_q025_m_s","velocity_q975_m_s"]]
    low_text = md_table(low,
        list(low.columns),
        ["f (Hz)","Observed azimuth","Azimuth 2.5%","Azimuth 97.5%","Observed projection (m)",
         "Projection 2.5%","Projection 97.5%","Observed velocity (m/s)","Velocity 2.5%","Velocity 97.5%"])
    paired_text = md_table(paired_diffs,
        ["frequency_low_hz","frequency_high_hz","paired_centers","paired_six_hour_blocks",
         "bootstrap_estimable","observed_delta_median_m_s","delta_q025_m_s","delta_median_m_s",
         "delta_q975_m_s","positive_fraction"],
        ["Step (Hz)","to (Hz)","Paired UTC centers","6-hour blocks","Interval estimable",
         "Observed paired Δv median (m/s)","Bootstrap Δv 2.5%","Bootstrap Δv median",
         "Bootstrap Δv 97.5%","Bootstrap Δv>0 (%)"])

    gate_totals = rejection.sum(numeric_only=True)
    gate_frame = pd.DataFrame([
        {"condition":"Window-frequency evaluations", "count":int(gate_totals.windows)},
        {"condition":"Not estimable: station/geometry", "count":int(gate_totals.not_estimable)},
        {"condition":"q=0 among estimable", "count":int(gate_totals.q0)},
        {"condition":"q=1", "count":int(gate_totals.q1)},
        {"condition":"q≥2", "count":int(gate_totals.q2plus)},
        {"condition":"q=1 fails coherence", "count":int(gate_totals.low_coherence_q1)},
        {"condition":"q=1 fails direction uncertainty", "count":int(gate_totals.uncertain_q1)},
        {"condition":"q=1 fails secondary-peak ratio", "count":int(gate_totals.ambiguous_q1)},
        {"condition":"q=1 fails grid resolution", "count":int(gate_totals.grid_resolution_q1)},
        {"condition":"q=1 fails velocity boundary", "count":int(gate_totals.boundary_q1)},
        {"condition":"q=1 fails phase residual", "count":int(gate_totals.phase_residual_fail_q1)},
        {"condition":"q=1 passes every direction gate", "count":int(rejection.accepted.sum())},
    ])
    gate_text = md_table(gate_frame,["condition","count"],["Condition","Window-frequency rows"])
    floor_text = md_table(real_floor,
        ["min_coherence","counterfactual_accepted_windows","n_frequencies_with_any",
         "longest_contiguous_points","longest_band_min_hz","longest_band_max_hz"],
        ["Coherence min","Passing q=1 windows","Frequencies with windows","Longest raw frequency run (points)","Run start (Hz)","Run end (Hz)"])
    front_text = md_table(front_floor,
        ["min_coherence","accepted_frequency_references","corrected_points","candidate_bands"],
        ["Coherence min","Accepted pair-frequency references / 50","Corrected CCF points","Continuous candidate bands"])
    lower_wave = lower_sim[lower_sim.case == "single_wave_snr0"]
    lower_noise = lower_sim[lower_sim.case == "pure_noise"]
    lower_sim_text = md_table(lower_wave,
        ["frequency_hz","min_coherence","trials","accepted","acceptance_rate"],
        ["f (Hz)","Coherence min","0 dB single-wave trials","Accepted","Recall"])
    lower_noise_text = md_table(lower_noise,
        ["frequency_hz","min_coherence","trials","accepted","acceptance_rate"],
        ["f (Hz)","Coherence min","Pure-noise trials","False single-window accepts","Rate"])
    common_table = common_summary.copy()
    common_table = common_table[["frequency_hz","common_utc_centers","q1_direction_accepted","pair_direction_usable",
                                 "k0_velocity_count","k0_velocity_median_m_s"]]
    common_text = md_table(common_table, list(common_table.columns),
        ["f (Hz)","Common UTC centers","q=1 direction-gate windows","Pair-direction usable","k=0 phase-speed windows",
         "k=0 median speed (m/s)"])
    accepted_refs = front_refs[front_refs.accepted].copy()
    refs_text = md_table(accepted_refs,
        ["min_coherence","frequency_hz","independent_windows","azimuth_deg","projected_distance_m"],
        ["Coherence min","f (Hz)","Independent windows","Reference azimuth (°)","Projection (m)"])

    sim_table = []
    for row in wave.itertuples():
        sim_table.append(dict(frequency_hz=row.frequency_hz, snr_db=row.snr_db,
            accepted=f"{row.accepted}/{row.trials}", median_azimuth_error_deg=row.median_azimuth_error_deg,
            p90_azimuth_error_deg=row.p90_azimuth_error_deg,
            median_abs_velocity_error_fraction=row.median_abs_velocity_error_fraction))
    sim_frame = pd.DataFrame(sim_table)
    sim_text = md_table(sim_frame,
        ["frequency_hz","snr_db","accepted","median_azimuth_error_deg","p90_azimuth_error_deg","median_abs_velocity_error_fraction"],
        ["f (Hz)","SNR (dB)","Accepted","Median azimuth error (°)","P90 azimuth error (°)","Median |velocity error| (%)"])
    noise_text = md_table(noise,
        ["frequency_hz","trials","q0","q1","q2plus","accepted"],
        ["f (Hz)","Pure-noise trials","q=0","q=1","q≥2","Accepted"])
    gain_text = md_table(gain,
        ["frequency_hz","trials","accepted","median_azimuth_error_deg","median_abs_velocity_error_fraction"],
        ["f (Hz)","Trials","Accepted","Median azimuth error (°)","Median |velocity error| (%)"])
    dual_text = md_table(dual,
        ["frequency_hz","trials","q1","q2plus","accepted"],
        ["f (Hz)","Independent dual-wave trials","q=1","q≥2","Accepted as one usable wave"])
    coh_text = md_table(coherent,
        ["frequency_hz","trials","q1","q2plus","accepted"],
        ["f (Hz)","Coherent same-tone trials","q=1","q≥2","Accepted"])
    ab_noise_text = md_table(ab_noise,
        ["frequency_hz","NW","snapshots","trials","q0","q1","q2plus","accepted"],
        ["f (Hz)","NW","K","Noise trials","q=0","q=1","q≥2","Accepted"])
    htable = pd.concat([h01,h045]).copy()
    htable["split"] = htable.split.str.replace("train_Nov01_03","Train Nov1–3",regex=False).str.replace("holdout_Nov04","Holdout Nov4",regex=False)
    hold_text = md_table(htable,
        ["frequency_hz","split","min_coherence","windows","accepted","q1","q2plus","median_accepted_azimuth_deg","median_accepted_velocity_m_s"],
        ["f (Hz)","Split","Coherence min","N","Accepted","q=1","q≥2","Median accepted azimuth (°)","Median accepted velocity (m/s)"])
    peak_text = md_table(peaks,
        ["frequency_hz","median_azimuth_deg","primary_windows","primary_windows_total","primary_share"],
        ["f (Hz)","Dominant peak azimuth (°)","Dominant windows","Sampled windows","Share (%)"])

    report = f"""# MGS 六台 MUSIC 筛选验证报告

**数据时段：** 2023-11-01 至 2023-11-05 00:00 UTC（四个完整 UTC 日）  
**台站：** BW_MGS01、BW_MGS02、BW_MGS03、BW_MGS04、BW_MGS05、BW_MGSBH  
**处理基线：** 200 周期窗，NW=3.5、DPSS 快照 K=6；方向网格 2°、速度网格 128 点；自动 MDL 信号数；最低相干度 0.70。

## 结论

当前结果**不支持直接放宽 MUSIC 门槛**。已知答案单波仿真说明估计器能恢复理想单波，但不是实际多分支噪声场的充分验证。实际数据在 0.70 门槛下没有可用的前段参考；把相干度降到 0.40 或 0.50 后虽能恢复少量逐频参考，方向分支散乱且没有连续候选带。较低门槛的纯噪声仿真也出现了单窗误收。

此前按各频率分别挑窗计算出的 0.09→0.10 Hz 正增速，并未在相同 UTC 中心的配对重算中复现：观测配对差中位数约 **−5 m/s**，时间块 bootstrap 中位数约 −2 m/s、95% 区间 **−489 至 +673 m/s**。六个低频点只有 1/96 个共同 UTC 中心具备完整的单源、台对方向和相位速度。因此当前不做真实结构反演；对这六点也不能用先前的非配对正增速判定其不满足单调条件。

## 仿真恢复与误接受

仿真使用本次真实台站坐标，合成波形经过生产版 DPSS 窗、快照构造和 MUSIC 估计。单波测试每频率每 SNR 8 次，另有台站振幅相差约 28 dB 的测试每频率 12 次。下表是方向/速度误差；SNR 是时域信号对噪声功率比。

{sim_text}

200 个独立纯噪声试验（每频率 100 个）没有一个被最终接受；自动信号数为 q=0 的比例为 94%，其余 q=1 或 q≥2 样本也被后续质量门槛拒绝。此处观测误接受率为 0/100/频率；样本量有限，不能把它解释成误接受概率为零。

{noise_text}

相干度从 0.70 扫描至 0.98 时，已知答案单波在 0 dB SNR 下仍保持 24/24 与 23/24 的接受率，纯噪声保持 0/40 与 0/40。更高阈值会逐渐损失单波：模拟种子留出选择的阈值为 0.999（0.10 Hz）和 0.98（0.45 Hz），其留出单波保留率为 90%。但这些阈值迁移到真实数据后会拒绝 11 月 4 日 0.10 Hz 的全部窗口，因此没有应用到正式流程。

![相干度阈值敏感性](coherence_cutoff_tradeoff.png)

{gain_text}

独立窄带双波的 48 次测试中，46 次识别为 q≥2；仅 1 次低频样本被当作 q=1 接受，0.45 Hz 没有双波样本被接受。两个完全相干的同频正弦波会退化为 rank-1 的空间快照，20 次低频测试均被识别为 q=1；这类信号本身不能靠当前协方差快照区分两个波源，需与独立双波结果分开解释。

{dual_text}

{coh_text}

## NW、快照数和纯噪声特征值

单因素对照固定 NW=3.5 改 K=4/6/7，再固定 K=6 改 NW=3.5/4.0。每频率、SNR、参数组合的已知答案单波测试为 10 次。各组合都通过单波质量门槛；方向误差中位数范围约 0.009–0.116°，速度相对误差中位数均低于 0.30%。纯噪声对照每组合 12 次，最终接受数仍为零；0.45 Hz 下 K=7 的 12 个噪声窗均判为 q=0，K=6 有 9/12 判零源。

{ab_noise_text}

因此快照数对纯噪声信号数的影响比 NW 对单波方向精度的影响更明显。当前仿真规模支持保留 NW=3.5、K=6 作为基线；尚不足以证明它在真实背景噪声波场中最优。

## 真实数据：台站、网格和多峰

0.10 Hz 逐站剔除对同一批 q=1 窗的方向中位变化为 0.17–1.46°，相速度变化中位数为 −0.65% 至 +0.86%。接受窗数从全六站的 18/32 变化到剔除单站后的 18–23/32。0.45 Hz 的全六站与所有逐站剔除对照均为 0/32 接受。网格由 2°/128 点加密到 1°/256 点后，0.10 Hz 两种网格各接受 4/12 窗；共同 q=1 窗方向与速度估计数值相同。

{real}

每个低频点抽取 5 个有效窗口查看谱峰。占优峰方向在 0.08 Hz 的 4/5 个窗口中约为 108°；0.09–0.12 Hz 的占比降到 2/5；0.13 Hz 主峰分散，最高方向组仅占 1/5。该采样只能用于观察分支稳定性，不能独立证明哪个峰是物理传播分支。

{peak_text}

![低频方向主峰跟踪](low_frequency_peak_tracking.png)

## 低频六点：边际区间与同中心配对结果

此前的 499 次时间块重采样分别汇总各频率的可用窗口，频率之间没有要求同一个 UTC 中心同时通过单源与台对门槛。因此这些逐频区间只能描述各频率的边际变化，不能把相邻频率差的区间当作同一窗口配对差，也不能据此断言速度随频率上升。原始测量点仍保留在表中；没有拟合或平滑曲线。

{low_text}

现补充同一 UTC 小时中心、同一六小时块同步重采样的相邻频率配对差。只有同一中心两频率均通过 q=1 方向、台对方向及 k=0 相位速度检查时才计算 Δv；少于 10 个配对中心或 4 个时间块时不报告区间。

下表中的 k=0 速度由未解缠的主值互谱相位计算，只作同时间中心的配对诊断，不是正式频散反演速度。

{common_text}

{paired_text}

0.09→0.10 Hz 的配对区间跨过零，不能确认增速；0.12→0.13 Hz 只有 9 个配对中心，区间不可估。六频点完整的 q=1 配对速度只出现在 1/96 个共同 UTC 中心，所以这批数据不能建立足够稳健的同分支六点曲线。前一版报告把非配对边际重采样误解释为明确正增速，现已更正。

![低频联合重采样速度区间](low_frequency_joint_bootstrap.png)

## 留出日核验

11 月 1–3 日作为训练期，11 月 4 日作为未参与参数选择的留出期；每个分组最多均匀抽取 16 个有效窗口。当前基线 NW=3.5、K=6、相干度 0.70 在 0.10 Hz 的训练期接受 7/16，留出期接受 15/16；q=1 窗方向中位数由约 105.3°变为 95.1°。0.45 Hz 两期均没有通过窗。相干度 0.98 在真实 0.10 Hz 两期均拒绝全部单源窗，说明按理想单波仿真选出的高阈值不能直接迁移到现场数据。

{hold_text}

## MUSIC 门槛与相干度反事实复核

前段 44 个频率合计 40,130 个窗频点。逐窗重算后，最终拒绝标签与正式基线逐项完全一致（最大差值 0）。以下门槛条件互相重叠，计数不能相加：

{gate_text}

2,798 个 q=1 窗中，360 个只失败一项门槛，2,438 个同时失败至少两项，没有一个通过全部当前方向质量门槛。相干度失败很常见，但网格分辨率、方向不确定性和次峰歧义也广泛存在，不能把零接受简单归因于单个阈值。另有 2,811 个 q≥2 窗的生产诊断字段把相位残差记为 0°；该残差实际未评估，应记为未评估，不能当作通过证据。多源窗仍被单独拒绝，不进入候选引用。

保持其余 q=1 方向门槛不变、仅下调相干度下限时，单窗层面可见：

{lower_sim_text}

{lower_noise_text}

仿真每组 100 个纯噪声窗和 24 个 0 dB 单波窗，样本量有限；0.40–0.50 的纯噪声误收提示降低门槛会增加风险。真实数据的单窗计数与生产逐频引用、连续带结果如下：

{floor_text}

再按生产版逐频方向分组、独立窗抽取、台对传播方向、投影和连续频带规则重建 50 组现有 pair-branch-frequency 引用：

{front_text}

0.40 只得到 4 个逐频引用，0.50 得到 3 个，0.60 及以上为 0；五种设置都没有形成任何修正 CCF 点或连续频散带。少量保留下来的引用方向也彼此跳变：

{refs_text}

因此原始单窗数量在 0.40–0.50 下看似连成 40 个频点，不能等同于同一 pair、同一传播分支通过生产引用与连续性检查。**本轮没有调整正式相干度阈值，也没有把反事实点送入相位反演。**

## 对全流程和报告质量的影响

原始全流程基线仍为完整状态，验证签名匹配；已完成 120,927 个相位窗频点评估、150 个逐窗通过，但连续带质量规则未产生反演交接点（0 点），所以没有真实结构模型。本轮全部测试结果写入独立目录，没有改写基线输入、MUSIC 输出或相位输出。`baseline_integrity.json` 记录基线完整性和哈希。

审计时曾发现第一版“纯噪声”测试仍叠加了单波信号。该批纯噪声结果已作废；本报告和结果包只使用显式零信号的修正版并完整重跑。仿真采用独立站点白噪声与有限窄带波场，不能覆盖真实非平稳噪声场；SNR 试验每格 8 次，NW/K 单波每格 10 次、噪声每格 12 次。真实方向峰跟踪每频率抽取 5 个窗口。因此结论用于检验当前门槛和估计器风险，不替代更广的独立数据集验证。

## 复核文件

- `run_validation.py`：可复跑的模拟、真实数据与联合重采样脚本。
- `run_next_diagnostics.py`、`run_coherence_front_rebuild.py`：门槛重算、共同 UTC 中心配对与前段反事实重建。
- `simulation_trials.csv`、`nw_snapshot_ablation_trials.csv`：已知答案仿真逐次结果及特征值。
- `real_leave_one_station_out.csv`、`real_grid_refinement.csv`：真实数据敏感性对照。
- `real_multipeak_observations.csv`、`real_multipeak_tracks.csv`：频率内局部谱峰与方向组汇总。
- `low_frequency_joint_bootstrap_draws.csv`、`low_frequency_adjacent_differences.csv`：同步重采样结果。
- `low_frequency_common_window_results.csv`、`low_frequency_common_window_block_deltas.csv`：同 UTC 窗与配对速度差。
- `rejection_condition_windows.csv`、`rejection_reason_reproduction.csv`：门槛重叠统计与正式拒绝标签复核。
- `coherence_lower_floor_simulation_trials.csv`、`coherence_front_counterfactual_summary.csv`：低相干度模拟与前段带宽反事实。
- `holdout_date_validation.csv`、`coherence_threshold_holdout.csv`：真实留出日与模拟种子留出核验。
- `baseline_integrity.json`：基线签名与输入哈希。
"""
    (HERE / "report.md").write_text(report, encoding="utf-8")

    files = sorted([p for p in HERE.iterdir() if p.is_file() and p.name not in
        {"results_bundle.zip", "results_manifest.json"} and
        (p.suffix in {".csv", ".npz", ".png", ".md", ".py", ".json"}) and
        p.name != "__init__.py"])
    manifest = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    (HERE / "results_manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    bundle_members = [p for p in files if p.name not in {"results_manifest.json"}]
    bundle_members.append(HERE / "results_manifest.json")
    with zipfile.ZipFile(HERE / "results_bundle.zip", "w", compression=zipfile.ZIP_DEFLATED, compresslevel=8) as z:
        for p in bundle_members:
            z.write(p, arcname=p.name)
    print("report", HERE / "report.md")
    print("bundle", HERE / "results_bundle.zip")
    print("bundle files", len(bundle_members), "bytes", (HERE/"results_bundle.zip").stat().st_size)


if __name__ == "__main__":
    build_report()
