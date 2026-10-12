# 2.4.1：BW.MGS 真实数据流程完善

- `array-select --station MGS01,MGS02,MGS03,MGS04,MGS05` 支持明确台站名单，名单与前缀、地理范围取交集，随后才进行阵列选择，避免相近前缀台站替代指定台站。
- 经原始 MiniSEED 和 StationXML 核验的 `h5data.ndstudio/1.0` 归档，在重叠样点完全一致时容许最多 2 微秒的数据块时间残差。这一容差仅用于重复样点识别；保留数据的原始 UTC 不被平移。相位分析仍按每条保留记录的真实时间对齐。
- 普通 HDF5 输入不继承这个 MiniSEED 容差；不同采样网格、超过容差的残差或任何重叠样点冲突仍会停止导入。审计表新增容差、最大时间残差及原始时间是否平移的记录。
- 导入重试会先核对现有 canonical 记录的完整样点、采样率和起始 UTC，相同记录直接复用，避免失败后重复压缩已完成的数据。
- 新增 `examples/bw_mgs_e2e/run_e2e.py`，仅编排本程序已有下载、审计、预处理、互相关、FTAN、相位频散、背景检验、条件反演和报告模块；新增简要报告生成器。
- CLI 增加 `--minimum-wavelengths`，用于显式指定自动反演选择的距离/波长下限。默认仍为3；本次另外使用2做探索性对照，所有逐点相位验收与背景检验条件保持不变，选择文件记录实际下限。
- 原有结构报告增加按分支着色的频散图，并分别导出 `11_selected_dispersion.png` 和 `12_vs_structure.png`。重新分析、未完成或无兼容结构时清除对应旧图，避免误展示上一轮结果。
- 每次 `scan` 完成后自动打印连续稳健区间和较好的候选门窗，包括频率范围、点数、带宽与未通过原因；默认各打印8个，`--diagnostic-top-k` 可调整。未产生合格频段时也会打印，不会自动修改阈值。
- 新增 `scan_band_diagnostics.csv`（所有逐点检查通过的连续区间，包含因带宽/点数不足而未收录的区间）和 `scan_candidate_diagnostics.csv`（候选门窗及主要限制、各项失败点数）。原始逐频点诊断仍在 `point_diagnostics.csv`。分块相位指标包含波包条件，打印名称明确为“分块波包与相位联合比例”，防止误当成纯相位失败。

只重新查看已有扫描的诊断，无须重算波形：

```bash
python run.py scan-diagnostics --project BW_MGS_run/project --diagnostic-top-k 10
```

摘要采用生成诊断表时的缓存参数。当前配置与历史参数不同时会提示；不会因修改了配置就把旧诊断重新标记成通过结果。

## 完整流程

联网重做：

```bash
python run.py doctor
python examples/bw_mgs_e2e/run_e2e.py --workdir BW_MGS_run
python examples/bw_mgs_e2e/brief_report.py BW_MGS_run
```

复用本次原始 HDF5 归档，先将归档解压至 `BW_MGS_run/download/`，并把选择文件放到 `BW_MGS_run/array_selection.json`：

```bash
python examples/bw_mgs_e2e/run_e2e.py --workdir BW_MGS_run --offline
python examples/bw_mgs_e2e/brief_report.py BW_MGS_run
```

本次五站 HHZ 数据时段为 2020-11-01 至 2020-11-05（UTC）；原始采样率 200 Hz。按照真实 5.92–18.40 km 台对距离设置 0.15–1.5 Hz 扫描、10 Hz 分析采样率、±100 s 互相关时延及 0.5–45 s 到时搜索。反演有限层厚度范围 10–10000 m、Vs 范围 200–5000 m/s，沿用明确固定的 Vp/Vs 和密度假设。

本环境未安装可选独立求解器 pysurf96，本次仅使用程序原有 disba 反演路径，并明确标注未作第二求解器复验。程序不以调低验收条件、补造频散点或模型来保证产生结构图。

## 本次测试结果与重现

有效共同窗110/287；相位片段1个（MGS01–MGS04，0.270–0.340 Hz，15点）；40次错时背景对照全部有效，观测中位SNR为29.45，对照95%分位2.39；条件周期曲线9条。

首次连续带宽下限0.10 Hz排除了15个已逐点通过的频点；改为0.05 Hz后提取该片段，所有相位门槛原样保留。因为使用本次诊断调整了参数，项目验证状态记为 `internal_reuse`。

默认距离/波长下限3得到14个兼容结构。最终额外对照使用下限2，共23个兼容结构，保留N=2、N=3两种周期解释；所有曲线、失败记录及默认下限3的原始结果均随结果包保留。

最终结果的离线重现命令：

```bash
python examples/bw_mgs_e2e/run_e2e.py --workdir BW_MGS_run --offline --minimum-wavelengths 2
python examples/bw_mgs_e2e/brief_report.py BW_MGS_run
```

最低带宽0.05 Hz只是连续观测的保留条件，扫描0.005 Hz步长不等于真实频率分辨率。N=2分支相速度随频率下降，N=3分支相速度随频率略升；不能用期望的变化方向自动选择唯一分支，也不能把窄频段的兼容Vs模型当成唯一地下结构。
