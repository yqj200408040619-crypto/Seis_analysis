# v2.7.1：统一 MUSIC 逐频点方向与相速度

前段频散分析、后段 Goertzel 筛窗调用同一个 MUSIC 实现；旧相位叠加定向和基于台对方位推测传播方向的关联门限已移除。相位联合速度计算只拟合 MUSIC 方向上的标量慢度，不再用二维相位直线拟合重新定向。

## 计算

1. 使用响应校正、统一 UTC、保留缺口的全台阵波形。每频率按 T=N/f、hop=S/f 遍历完整记录，MUSIC 模式补上最后一个完整窗；缺口不能拼接成连续数据。
2. 在同一中心频率，使用 DPSS 正交窗生成多组复数谱快照。构造 Hermitian 互谱协方差、特征值分解、信号/噪声子空间。
3. 在二维水平慢度网格上搜索 MUSIC 伪谱，再最小化噪声子空间投影进行细化。北=0°、顺时针为传播方向；来波方位另相差180°。
4. 方向质量检查包括有效台站数、二维几何、信号数、相干性、残差、搜索边界、分离峰及高分候选范围。近直线阵列和无可识别信号不输出可用于修正的方向。自动 MDL 还要求次信号达到可配置的最小能量比例，避免有限带宽泄漏误判多波源；多波源窗不混成一个校正方向。
5. 每个台对、每个频点分别确定方向 θ(f) 和带符号投影 ℓ(f)=Δx sinθ+Δy cosθ。前段 c(f)=c_app(f)|ℓ(f)|/d。保留原视速度、MUSIC 独立估计速度、投影和来源字段，不用 MUSIC 速度直接替换条件周期假设。
6. Goertzel 使用同一实现重新定向，保留默认45°台对夹角、方向分辨率和独立窗要求。每窗按它自己的投影距离旋转到该频点的参考投影，再筛窗、叠加和反算速度。频率间不共用一个平均方向或一个固定投影距离。

此处速度是沿水平波矢的相速度；表面波解释仍需台阵与波场满足相应假设。高分候选范围和时间分位数是分辨率诊断，不是已标定置信区间。

## 两阶段规则

| 部分 | 方向算法 | 45°台对门槛 | 不合格频点 |
|---|---|---|---|
| CCF/FTAN候选的前段校正 | MUSIC | 不使用 | 在审计表中保留，合格连续曲线中剔除 |
| Goertzel与低频实测重建 | MUSIC | 默认使用 | 不入最终实测曲线和反演 |

前段接近垂直入射、投影过小或不确定范围跨零时也不输出稳定速度；这不是45°筛窗。

低频外延仅提供筛窗先验。新频点必须有实际波形、MUSIC方向、足够独立窗、时移背景检验和连续频带；预测值不能直接生成频散或结构。修正后随频率增加的曲线继续排除。方向随频率平滑变化可连续通过；方向分组的15°分箱边界不会单独截断低频曲线。

## 参数和结果

- 方向参数保存在项目 `phase_settings.json` 的 `direction` 中。`method` 只允许 `music`；缺少该字段的旧设置自动补为 MUSIC。
- `music_snapshots=6`、`music_time_bandwidth=3.5`：DPSS 快照数和 NW；实际谱半宽 NW/T 会逐频打印。
- `music_front_cycles=20`：前段谱窗周期数；六站大孔径验证设200，后段也使用200周期以收窄DPSS实际带宽。
- `music_sources=0`：自动信号数；`music_max_sources=2`；`music_min_source_fraction=0.1`；`music_diagonal_loading=1e-5`。可在高级 JSON 调整。
- `velocity_min_m_s`、`velocity_max_m_s` 是 MUSIC 搜索范围，独立于用于到时几何设置的2000m/s。
- 前段：`output/music_front/frequency_directions.csv`、`dispersion.csv`、`point_audit.csv`、`report.html`。
- 后段：`output/phase_guided/direction/`、`direction_reference.csv`（低频模式逐频参考）、`windows/`、`inversion_curves.csv`、`phase_report.html`。
- 窗档案保存原始互谱、逐窗方向/投影、旋转角、参考投影、接受掩码；`phase-audit` 重算并核查。
- 完成标记和方向缓存纳入 MUSIC 源码哈希；旧算法缓存不能复用为 MUSIC 结果。

## 使用与更新

新项目和主流程默认使用 MUSIC。完整分析会在候选生成后自动运行前段 MUSIC，再交接给后段。已有项目可先运行：

```text
python run.py music --project 项目目录
python run.py phase-stack --project 项目目录 --phase-settings 项目目录/phase_settings.json
```

六站流程：

```text
python -m ndstudio.mgs_e2e --workdir MGS6_run
```

只在下载已经完整且已核查时使用 `--reuse-download`。中断后同目录续行；结果改变会自动废止相应缓存。前段连续带未通过时，只允许原CCF候选作为后段低频筛窗先验；只有重新实测的合格点可流入反演。全流程会记录科学门槛未通过的原因并跳过自动反演。

关闭正在运行的 Studio 和下载窗口，将补丁解压覆盖程序根目录，重新打开 `run_windows.bat`。无需重启 Windows；项目设置不会因单纯覆盖程序而自动改写。旧频散/相位/反演要重新运行相应阶段，不能把旧图片当作新算法结果。

固定40秒阻断已删除。自动到时范围为 ceil(max(已配置到时上限,1.1×最大台距/2000))，噪声区、最大 lag 和连续窗要求随之向外调整；2000m/s只用于这一步，后续速度可以超过2000m/s。

参考：Goldstein & Archuleta (1987), Array Analysis of Seismic Signals, GRL 14(1), 13–16，https://archuleta.faculty.geol.ucsb.edu/library/pdf/1987-PG-RJA-GL014i001p00013.pdf 。
