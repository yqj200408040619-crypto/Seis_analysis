# v2.5.4：方向约束、连续处理与人工调参

第6页现在把整段波形、逐频方向估计、相位筛选、分组叠加、连续频带检查连在一起。一键分析也执行这个相位阶段，再将合格的连续频带交给候选审查和条件结构反演。界面日志同步打印实际参数、输入范围、窗口数量、淘汰原因、最终相位和输出位置。

## 时间与相位的含义

默认处理全部已导入记录。开始/结束UTC只描述数据覆盖范围；高级区间选项仅供对照试验，不决定单窗长度。每个频点的目标窗长为 `N/f` 秒，滑动为 `M/f` 秒，取整到采样点；只使用完整有效窗。默认共同起点顺延600秒，可在“起点稳定时间”调整。

开启方向功能时，全阵列使用同一UTC滑窗网格，每个台对只使用它自己的共同有效窗。台站缺测不通过拼接消除，也不把所有台站的记录先裁成最短共同区间。关闭方向时保留附件的台对局部滑窗逻辑。

互谱为 `G目标 × conj(G参考)`。正、负相位分别累计，最终相位为复数和的辐角；不通过反号或共轭强行合并。最近DFT频点与精确曲线频率均可选，实际计算频率写入报告。

## 方向怎么参与

1. 用所有有效台站在同一窗口的逐频复系数扫描二维慢度，得到传播方向与速度。北为0°，顺时针增加；记录的是传播方向。
2. 邻近频点分别计算波束功率再平均。输入频散曲线只提供待处理频点，不为F-K设置目标速度；搜索上下限来自方向参数。
3. 检查台站数、二维几何、相干度、扫描网格、边界峰、多峰歧义及空间相位残差。
4. 检查台对连线与波矢的轴向夹角，默认超过45°排除。反向传播保留有符号投影并单独分组。
5. 默认同时要求近峰候选范围不越过45°，并检查该范围引起的余弦修正误差预算（高级JSON默认15%）。近峰候选范围是分辨率诊断，不是校准后的统计置信区间。
6. 在曲线相位筛选之前选定主方向组，再执行幅值检查和相位门限。不同方向不直接混叠成一个平均角度。

方向模式 `filter` 与 `joint` 均执行方向筛窗并输出独立平面波联合估计；`diagnostic` 只输出方向诊断，保留附件式相位筛窗；`off` 关闭方向模块。界面默认 `filter`；未指定方向参数的旧CLI/API调用默认 `off`，保证旧脚本兼容。近直线阵列不能唯一辨识二维方向与速度，宜先用诊断模式查看几何限制。

`direction/joint_dispersion.csv` 是不依赖输入CCF源相位与周期的独立平面波估计。它使用有符号空间相位、F-K支持的整数周期和多台对联合拟合，按方向组对齐后叠加，并用不重叠窗口计数。报告中的时间分位范围表示记录内变化，不是速度置信区间。

CCF叠加频散沿用原有源相位、周期、时延分支、群速度与背景检验状态。独立平面波结果单独显示；没有自动指定Rayleigh模态或将它解释成唯一浅层结构。

## 连续带宽与旧参数问题

第2页“最小连续带宽”对应 `config.json → scan.min_band_hz`，不是第6页F-K频带半宽。高级JSON与可视控件采用三方合并：只有从加载值实际改动的控件会覆盖JSON，保存后同步显示。计算开始打印实际参数，并保存本次参数快照。

叠加完成后再次按当前 `min_band_hz` 和 `min_points` 检查连续频带。淘汰频点、频率缺口、方向组改变均会分段，不跨缺口展开相位。每段单独计算叠加速度和反演交接曲线。无合格频带时明确显示原因并跳过自动反演。

改动方向、相位或带宽设置后应重跑相应阶段。新生成的校正波形缓存把校正参数与扫描参数分开；仅改变带宽无需重新响应校正。旧版没有校正签名的H5仍要求配置一致，升级时先从原始记录重建一次。参数或输入改变会使历史交接结果失效。

只修改45°阈值、角度误差预算或方向分组时，可以复用逐频变换与F-K缓存，重新生成掩码、叠加及频带结果。缓存包含波形、几何、UTC、计算参数与源码指纹；日志会打印 `cache_reused`。

## 从哪里追溯

| 输出（相对于项目目录） | 用途 |
|---|---|
| `output/workflow_processing.log` | 一键流程的各阶段参数、结果与失败信息 |
| `output/workflow_trace.json` | 阶段状态、耗时、结果与运行编号 |
| `output/workflow_settings_used.json` | 一键流程实际参数快照 |
| `output/phase_guided/phase_processing.log` | 每频点方向原因、每台对筛窗数量与最终相位 |
| `output/phase_guided/phase_trace.json` | 相位阶段结构化追溯 |
| `output/phase_guided/phase_settings_used.json` | 相位阶段实际参数与方向设置 |
| `output/phase_guided/phase_gate_summary.csv` | 逐台对/频点接受数、方向组与淘汰原因 |
| `output/phase_guided/windows/*.npz` | 原互谱、UTC起点、逐项掩码、夹角与候选范围 |
| `output/phase_guided/direction/report.html` | 方位、速度、波束图和阵列响应图 |
| `output/phase_guided/direction/frequency_diagnostics.csv` | 方向阶段各原因的窗口数 |
| `output/phase_guided/phase_continuous_bands.csv` | 当前带宽和连续点数的接受/淘汰结果 |
| `output/phase_guided/inversion_curves.csv` | 合格的条件CCF叠加频散交接表 |

第6页提供“查看逐阶段日志”“查看方向与联合频散”“查看连续频带诊断”按钮。阶段失败不保留本次完成标记；界面会提示参数改变后的历史结果需要重跑。日志与追溯文件保存当前运行，保留人工调参历史时请复制本次参数、日志和诊断表。

## 按原因调参

| 看到的原因 | 优先检查 |
|---|---|
| 阵列孔径与频率要求更细扫描网格 | 减小方向步长、增加慢度点数；搜索速度范围应有独立依据 |
| 有效台站不足 / 几何退化 | 台站UTC、缺测、坐标和二维覆盖 |
| 存在近等强的分离峰 | 阵列响应与空间混叠；不要仅凭最大峰强制指定方向 |
| 最优速度在搜索边界 | 有依据地扩大速度范围，检查空间混叠 |
| 夹角或候选范围超限 | 改用更接近传播轴的台对；核对方向分辨率 |
| 相干度或空间相位残差不合格 | 窗长、频点、缺测及多方向波场 |
| 方向合格但曲线门限不合格 | 源相位、周期、时延分支、输入曲线及相位容差 |
| 连续点数或最小连续带宽不足 | 先看淘汰频点与方向组变化，再有依据地调整带宽 |

大量重叠窗口不等于独立样本；最少不重叠窗口默认3个。45°是控制误差放大的工程阈值，没有在45°处发生数学奇异。

## 可复现实验与启动

Windows完整解压后运行 `install_windows.bat`，再运行 `run_windows.bat`。本次验证环境为Linux/Qt离屏，未做Windows原生实测。

在程序目录运行：

```bash
python examples/create_direction_example.py --directory examples/direction_plane/project
python run.py gui --project examples/direction_plane/project
python run.py phase-stack --project examples/direction_plane/project --phase-settings examples/direction_plane/project/phase_settings.json --aligned-h5 examples/direction_plane/project/aligned.h5
python run.py phase-audit --folder examples/direction_plane/project/output/phase_guided
```

示例是7台站、600秒、20Hz、已知30°方向与500m/s的受控平面波。它从预设的条件CCF目标开始，用于验证方向、共享窗口、45°筛选与交接；不是实际观测，也不是完整扫描或浅层结构真值。CCF源相位/周期别名的速度与独立500m/s结果不同，报告保留这一区别。

完整回归：175通过、12跳过。真实DARE已下载校正波形的整段重放和细网格补充核对见 `Direction_Validation_v254.html`，本轮没有重新向服务商下载。用户MATLAB方向提取附件原样保存在 `matlab/reference/traditional_fk_mseed_direction.m`，Python实现修正传播符号并加入分辨率与歧义检查；MATLAB运行时未验证。
