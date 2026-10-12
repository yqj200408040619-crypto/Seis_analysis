> 此目录是底层兼容计算核心。新版 HDF5 / GUI / 分组反演请从上级目录 `run.py` 启动，以根目录 README.md 为准。

# 可复用传统环境噪声频散程序 v1.0

这是此次0901分析流程的通用版本：从原始波形开始，完成UTC对齐、共同去前600秒、局部异常处理、响应校正、多路线性互相关、多峰FTAN追踪、固定窗验证、错时背景对照和条件频散输出。**没有0901专用台站、日期、固定窗数、旧状态分类器或旧地下模型依赖。**

程序默认保留多解。数据不能支持可靠片段时输出空曲线表和拒绝原因；它不会靠降低阈值或插值空档制造一条完整曲线。它实现的是本次四台0901的传统互相关/FTAN路线；早期13台、100m线阵的F–K路线是不同的测量方法。

## 这次0901提取走过的步骤

| 步骤 | 做了什么 | 意义与需要检查的内容 |
|---|---|---|
| 1. 输入和几何审计 | 核对波形、UTC、采样率、通道、台站坐标、台对距离和数据指纹 | 相速度依赖真实距离；防止不同采集批次、分量和旧缓存混用。 |
| 2. 共同时间轴 | 各台按真实UTC重采样到20Hz，统一去共同记录前600s，保护滤波边界 | 保证同一分析窗来自同一物理时段，而非各台首样点后的同一序号。授时正确仍依赖记录时间及外部授时信息。 |
| 3. 数据清理 | 检查缺测、平线、异常尖峰；建立共同有效样点 | 短孤立尖峰可局部修复，较长异常局部掩码；不因一个尖峰直接丢掉整个600s窗。 |
| 4. 仪器响应/极性核对 | 用采集历元匹配的StationXML校正响应，核对标称传感器方向与增益 | 传播相位不能混入仪器相位；标称同极性不等于现场极性已经校准，也不能单独消除绝对周期。 |
| 5. 单台降噪和互相关 | 零相位带通、运行绝对均值归一化、谱白化；对各台对计算正/负时延线性互相关 | 压制强瞬态和窄带源，累积共同传播信号。RAM5/20、不同白化频宽和one-bit用于比较处理依赖。 |
| 6. 宽频诊断、固定窗分组 | 先看各频段窗数和谱结构，再按全天、UTC昼夜与开发集学习的状态分组 | 频段由证据决定；同一支路的整段频带共用时间窗，避免逐频点“挑好数据”。 |
| 7. 多峰FTAN | 在频率—到时图中保留多个局部峰并追踪连续波包 | 次峰可能属于真实传播状态；不能永远只取最强峰。到时给群速度，尚不是相速度。 |
| 8. 固定门窗与波形验证 | 开发时段学习门窗；验证时段重复，并检查参数、FTAN频宽及门窗扩缩/整体移动 | 拒绝混合波包、门窗敏感相位和只在一个时段出现的峰。增加扫描频点不能自动提高物理分辨率。 |
| 9. 分块和错时对照 | 6h分块、留一块、整块bootstrap；40次台站一致错时，重新应用冻结分类器 | 检查时间稳定性和背景分离。缺少对照记NaN；结果不是校准p值或正式置信区间。 |
| 10. 周期/源相位多解 | 对每个连续实测段枚举整段固定的2π周期N和源相位，再比较不同台对的重叠频段 | 相对相位可重复，但绝对相速度仍可能多解；不逐点换N，不跨空档展开相位。 |
| 11. 可选结构兼容性 | 在明确的周期和Rayleigh模态假设下做多起点分层反演，并用独立求解器核对 | 缩小可兼容解释；保留多模型。小拟合残差不意味着唯一结构或可信地层界面。 |

0901已有结果：285个共同分析窗；35个通过波形验证的相位片段，片段和状态存在重叠，不能视为35条独立曲线。主低频段约1.65–2.75Hz，辅助低频段2.325–2.65Hz；夜间慢波包高频段扩大到4.30–5.00Hz。**主支路2.75–4.30Hz仍缺少实测连接。**

响应核对和高频扩展后，原探索性范围内保留32个结构兼容实例；更严格的群速度优先集仍有4种周期/模态解释。它们是源相位、周期、模态、Vp/Vs和密度等假设下的结果。这些数字属于旧0901结果，不是新程序对任何数据的固定输出数目。新程序重新发现与验证全部候选，不能承诺复现相同的候选编号或数量。

## 安装

推荐Python 3.12，Linux/macOS/Windows均可执行波形和相位流程：

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

可选结构反演另需：

```bash
python -m pip install -r requirements-inversion.txt
```

`pysurf96`含Fortran扩展，部分平台需要编译环境或WSL。缺少它不影响波形/频散候选提取；默认结构验证会明确失败，不把缺少第二求解器标记为通过。把`require_independent_solver`改为false只能得到未独立验证的探索性模型。

## 新数据最少需要什么

1. 每台一个明确的垂直分量通道，以及带UTC和采样率的波形。MiniSEED/SAC等由ObsPy读取；也支持下面定义的NPZ。
2. 台站真实坐标，选经纬度或本地平面坐标之一。
3. 采集日期匹配的StationXML完整响应，或明确声明输入已经响应校正，并给出统一物理单位。
4. 一份JSON配置。默认值只适合类似0901的分析尺度，换台间距或频段后必须调整相关参数。

```bash
python run.py init --directory my_project
```

随后编辑`my_project/stations.csv`和`my_project/config.json`。相对文件路径均以配置所在目录为基准；`response_xml`以stations.csv所在目录为基准；波形`file`以`data_root`为基准。支持波形通配符及用分号分隔多个文件，但`seed_id`必须是一个精确NSLC，不能用通配符。

台站表例子（坐标只是格式示例，必须换成真实值）：

```csv
station,file,seed_id,x_m,y_m,response_xml,precorrected,units,polarity,polarity_verified,clock_offset_s
A,A.mseed,XX.A..HHZ,0,0,response.xml,false,counts,1,false,0
B,B.mseed,XX.B..HHZ,1000,0,response.xml,false,counts,1,false,0
```

| 字段 | 含义 |
|---|---|
| station | 唯一标签，字母/数字/下划线/连字符；不必等于SEED台名 |
| seed_id | `NET.STA.LOC.CHA`，无位置码时保留两个点，如`CH.SBAJ2..HGZ` |
| latitude, longitude | WGS84经纬度；程序计算大地距离和本地坐标 |
| x_m, y_m | 没有经纬度时使用本地米制坐标 |
| response_xml | StationXML文件；响应必须唯一且覆盖完整输入历元 |
| precorrected, units | 已校正输入填true；ACC填`m/s^2`，VEL填`m/s`，DISP填`m` |
| polarity | 仅应用明确声明的+1/-1，不用更好拟合来自动翻转台站 |
| polarity_verified | 是否有独立现场极性证据；仅标称方向不能自动填true |
| clock_offset_s | 记录时间减真实时间；程序按`真实UTC=记录UTC−该值`修正，0表示未声明修正 |

NPZ必须包含一维`data`、标量`fs_hz`、ISO UTC字符串`start_utc`；可选同长度布尔`valid`。无效样点既可用`valid=false`也可用NaN。NPZ必须由外部记录保证绝对时间和物理单位真实，程序不会猜采集时间。

同一台站输入内采样率改变、互相矛盾的重叠记录、响应历元切换，程序会停下并说明原因，需要先按实际历元整理；不会用单个响应强行覆盖不同仪器。真正缺测不插值成波形。

## 完整运行和分步运行

```bash
python run.py all --config my_project/config.json
```

结果写到配置的`output`目录，先打开`report.html`。Windows也可把配置文件路径作为参数传给`run_windows.bat`。

```bash
python run.py prepare  --config my_project/config.json
python run.py scan     --config my_project/config.json
python run.py controls --config my_project/config.json
python run.py families --config my_project/config.json
python run.py report   --config my_project/config.json
python run.py verify   --config my_project/config.json
```

改变参数、代码或原始输入后重新执行`prepare`或`all`；配置、代码和原始文件指纹改变时缓存会失效。分步命令不会静默使用不匹配的缓存。重新测量会清除旧的可选反演产物，须再运行`invert`，以防把旧结构接到新数据上。读取已有结果报告无需计算。数据大时可按阶段执行，程序不修改原始波形。

## 哪些参数需要按新数据调整

| 参数 | 默认/作用 | 调整依据 |
|---|---|---|
| sample_rate_hz | 20Hz，共同输出采样率 | 至少覆盖目标带宽且低于各台原生Nyquist；采用抗混叠重采样，不通过插值创造高频。 |
| discard_first_s | 600s，只对共同记录去除一次 | 本项目要求；仪器进入采集时段不要逐台按序号截掉。 |
| bandpass_hz, response.pre_filt_hz | 0.1–7.142857Hz；响应预滤四角频率包围带通 | 目标频带须处于共同仪器有效频带内，远离Nyquist及响应极弱区。 |
| window_s, min_clean_segment_s | 600s窗；有效连续片段至少150s | 连续片段必须长于最大时延，窗内有效样点比例默认≥0.7；未满足者输出拒绝记录。 |
| spikes.* | 高阈值、极短修复上限0.04s；较长异常掩码 | 应先对波形、修复日志和敏感性结果检查。不能保证自动识别每一种异常；高阈值是保守起点，不是把所有高幅信号判为坏数据。 |
| max_lag_s, scan.arrival_* | 互相关±60s；波包搜索0.2–25s | 根据真实台间距和群速度范围设定；到时、晚到噪声区与时延边缘必须分开。 |
| scan.noise_lag_s | 35–55s，包络晚到噪声RMS | 必须检查该区是否也有真实到达；更短台间距需相应缩小或重新设置。 |
| scan.min/max/step_hz | 1.5–5.0Hz，扫描0.025Hz | 先宽频诊断再锁定候选；更密扫描不增加有效独立频点。相位跳跃采样不足的点会拒绝。 |
| scan.ftan_width/narrow_width | 高斯σ=max(0.04,0.15f)；窄频宽0.08f | 时频分辨率有取舍，必须比较到时及相位稳定性。 |
| validation.block_s, split_utc | 6h块，开发/验证各一半 | 按数据持续时间和源状态相关时间设置，分界必须落在块边界；分析前冻结。 |
| validation.status | 新数据默认内部时间顺序验证 | 对已反复分析的0901设`internal_reuse`；外部验证必须由独立数据的实际使用情况支持。 |
| cohorts.states / feature_band_hz | 开发集PCA+k-means，默认2种状态 | 0关闭状态分组；只在开发数据拟合，验证和错时数据均只应用冻结分类器。分类不是物理模态识别。 |
| routes | RAM/白化多路必需对照，one-bit为挑战 | 1Hz陷波仅在0901示例中作为必需对照，新的数据按观测到的干扰决定，不能照搬干扰频率。 |
| controls.trials | 默认40次 | 候选至少20次有效错时对照且观测SNR超过背景95%分位才标记背景分离；这仍不是校准显著性。 |
| phase.* | π/4源相位；N=0..80；c=150..3000m/s | 明确来源相位、速度范围和有限搜索范围。绝对周期不是自动确定；真解可能不在设定范围内。 |

波形验证默认：两时段包络SNR≥10；处理路线相位差≤0.2rad、两时段相位差≤0.35rad、到时变化≤10%、相位导数时延与包络到时差≤20%；至少4块、75%块支持；门窗相位变化≤0.35rad；留一块相位变化≤0.2rad。连续片段至少8点且宽≥0.25Hz。这些是本项目的探索性筛选条件，不能理解为通用“提取标准”。

不要直接把某一个coherence门槛当作最终条件。平方相干通常在0–1之间；原讨论的2.5以上指标是峰比，并不等于coherence。真实传播信号还需要重复性、波包、相位和几何证据。

## 三个层级的输出

| 文件 | 保存什么 | 怎么用 |
|---|---|---|
| record_audit.csv / local_quality_events.csv | UTC、钟差声明、响应、标称方向、修复/掩码事件 | 先检查基础数据与处理是否正确。 |
| window_quality.csv / window_schedule.csv | 每窗有效比例、有效片段、实际参与的共同窗 | 评估丢数据原因；修复只能恢复局部，不无限放宽窗质量。 |
| cohort_counts.csv / cohort_assignments.csv | 各固定群的开发/验证窗数 | 先看哪些群有足够的独立时间块。 |
| local_packet_peaks.csv / candidate_gates.csv | 开发峰、追踪学习的固定门窗 | 检查选的是否为同一波包。 |
| point_diagnostics.csv / route_diagnostics.csv | 包括拒绝点的全部验证指标 | 查清具体失败环节，再决定延长观测或检查波包混合。 |
| measured_phase_segments.csv | 相对展开相位、到时、包络群速度、重复性范围 | **第一层：实测量。**没有绝对N或模态确认。 |
| timeshift_summary.csv | 背景分位、有效对照数和分离标记 | 缺少控制数据是NaN，不能当作零噪声；只作诊断。 |
| conditional_dispersion.csv / curve_catalog.csv | 整段固定N、源相位假设下的候选c(f)及kr、波长比 | **第二层：条件频散。**所有点有实测相位，未夹入预测。先看背景/响应标记和跨台对兼容链接。 |
| cross_pair_compatible_curves.csv | 不同台对的重叠频段速度一致组合 | 约束备选而不是证明同模态；共享台站的台对也不完全独立。 |
| segment_overlap_connections.csv | 同台对、同窗群、同分支的实测重叠连接及整段整数对齐 | 没有实测重叠不自动拼接；不能跨缺测区补点。 |
| conditional_models.json | 可选分层结构兼容实例、假设及独立求解器状态 | **第三层：模型解释。**仅表示在假设下兼容，不是唯一Vs或地层置信区间。 |
| inversion_attempts.csv | 各次优化的误差、数值失败及拒绝原因 | 区分没通过物理兼容条件、漏根和优化未收敛；不用零填补漏根。 |
| selected_observed_dispersion.csv / model_predictions.csv | 分开的观测曲线和模型预测；预测缺根留NaN | 绘图时把空档预测画成虚线，不放进实测数据或拟合权重。 |
| report.html / figures / verification.json | 离线报告、图和程序一致性校验 | `checks_passed`表示数值/文件逻辑检查通过，不表示已证明唯一物理频散。 |

相位约定：互相关为conj(X_i)X_j，j延迟对应正时延；负分支反转为正时间后再分析。门窗傅里叶相位为θ=arg∫CCF(t)exp(−i2πft)dt。在假设θ=φ−kr+2πN下，

```text
c(f) = 2π f r / [φ − θ(f) + 2πN]
U(f) = r / 包络到时
```

每个连续片段的N为常数。不同片段的局部展开基准可以相差整数2π，不能机械要求编号相同。π/4源相位仅在适当格林函数远场约定下合理；强方向性、近场、体波及多模态混合可能使其不成立。相位导数帮助检查时延，不自动消除N。只用Z分量和数值匹配不能证明Rayleigh模态。

本地平面坐标可以输入13台线阵，但这是互相关流程；**它不实现F–K反演，不从线阵估计唯一二维方向，也不把沿阵列投影速度自动当作真实相速度。**短台间距下，到时可能小于波包分辨能力，即使空间Nyquist没有违反也未必有可用频散。

## 可选：从候选频散估计结构

先阅读曲线目录和跨台对链接，为每个备选明确选择整段`curve_id`和假设模态。根据`examples/selection_template.json`建立文件：

```json
{
  "settings": {"layers": [1, 2], "starts": 12, "vp_vs_ratio": 2.0,
               "density_g_cm3": 2.0, "require_independent_solver": true},
  "cases": [
    {"name": "cycle_mode_case_1", "bands": [
      {"curve_id": "把这里换成curve_catalog中的实际ID", "mode": 0}
    ]}
  ]
}
```

```bash
python run.py invert --config my_project/config.json --selection my_project/selection.json
python run.py report --config my_project/config.json
```

可用`fmin_hz/fmax_hz`限定某个片段中的观测子区间；同台对的重叠观测不能重复计权。`mode=0`基阶、`mode=1`第一高阶；这是用户声明的假设，不能由同一波包的峰次序推断。

反演从新随机起点开始，不读取旧0901模型。每个选中观测带等权、带内点平均分配权重；默认相速度加权RMS≤3%、最大点差≤5%，每带群速度中位误差≤15%；独立surf96相速度差≤0.5m/s，超过60个频点时自动分批求出全部实根，漏模态根、非物理根不以零补齐。1–4个有限层可配置，厚度和Vs界限、Vp/Vs、密度均是模型假设。有限多起点搜索可能漏解，未找到模型不是地质无解证明；多个模式和周期组合需分别运行。

## 自带测试与示例

```bash
python -m unittest discover -s tests -v
python run.py demo --directory demo
python run.py demo --directory demo_noise --field noise
```

`demo`包含4台已知沿线传播方向的色散随机场：c(f)=600−20f m/s，原生40Hz、各台起始网格相差0/0.007/0.013/0.021s，加入单样点尖峰、2s平线与缺测，按真实UTC处理到20Hz。其源相位0是合成设计的已知量，不能照搬到真实噪声。测试要求至少两个台对恢复已知频散，且保留包含局部异常的部分有效窗。它不是Rayleigh分层地球模拟，不能用来验证地下结构恢复。

`demo_noise`是独立随机噪声，要求没有通过全部验证的频散候选。可选的disba/pysurf96测试核对均匀半空间解析Rayleigh根、分层基阶/第一高阶以及单位/模式索引；缺少可选库时明确skip。

`validation/`提供此次运行的测试记录。0901回归核对范围是35个旧已验证固定窗片段、1042个相位点的FTAN与门窗相位复算，不是重新处理两天250Hz原始数据；记录里明确写了`raw_0901_end_to_end_replay=false`。新数据使用全部流程必须实际运行自己的原始波形。

`examples/config_0901.json`、`stations_0901.csv`和对应响应快照方便重跑0901。把原始4个mseed放到程序根目录的`data_0901/0901/`；其他数据不得使用这份响应或该日期分界。

## 仍然提取不出来时

1. 先检查时轴、单位、响应、距离和异常处理，查看保留窗数；不要直接降相干阈值。
2. 若开发/验证SNR不足，优先增加观测时间、有效台对与独立时间块；短时各向同性随机场也可能被正确拒绝。
3. 若门窗敏感、到时和相位导数不一致，检查波包混合、局部近场源、多模态和方向性噪声；继续增加扫描点不能解决物理混合。
4. 若相位稳定但周期/模态多解，增加不同距离的路径、另一日未使用记录、三分量/极化和独立地质信息。几何和仪器极性是优先约束。
5. 缺测频带保持空白；模型跨空档预测可以帮助设计补测，但不是新观测。

## 方法与软件来源

- Bensen et al. (2007), Processing seismic ambient noise data to obtain reliable broad-band surface wave dispersion measurements: https://doi.org/10.1111/j.1365-246X.2007.03374.x
- ObsPy仪器响应文档：https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.remove_response.html
- FDSN StationXML定义：https://docs.fdsn.org/projects/stationxml/en/latest/reference.html
- disba官方实现：https://github.com/keurfonluu/disba
- pysurf96官方实现：https://github.com/miili/pysurf96

FTAN多峰追踪、固定门窗和分块验证沿用本项目阶段4–7的算法思想，已移除旧数据专用依赖。本包保存完整源码、配置、示例及测试；报告无需联网。
