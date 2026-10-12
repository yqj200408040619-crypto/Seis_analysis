# 2.5.6：低频实测重建、台距投影与成对结构反演

若最终要输出沿波矢的真实相速度，还必须一致地处理台距投影，不能直接把台对视速度当作真实波速。

完整流程：程序下载与完整性检查 → 导入和预处理 → CCF/FTAN扫描、背景检验和候选频散 → 排除升频增速曲线 → 向低频外延速度先验 → 独立阵列方向确认 → 逐频Goertzel相差及方向筛窗 → 投影相位对齐、复数叠加 → 合格实测频散 → 修正前后结构反演。下载、预处理、扫描与反演求解器继续使用原模块。

**45°台对角度限制只用于后面的方向/Goertzel筛窗，不对前面的下载、预处理或扫描新增这个限制。**

## 低频预测不是测量

用非增速候选的低频端点对数斜率外延，并限制斜率和项目速度边界。低频下限可以低于原扫描频带，但必须高于预处理通带下限及F-K邻频半宽。记录不足完整窗口、方向不明或质量不合格时不制造低频曲线；频带内缺口不会直接填为观测。

`screening_prior_predictions.csv` 保存外延预测；`selected_phase_targets.csv` 保存投影后的原始平面波筛选先验，标记 `is_model_prediction=True`、`target_role=screening_prior_not_observation`。只有真实波形的合格叠加相位才进入 `inversion_curves.csv`，其中 `is_model_prediction=False`、`is_prior_prediction=False`；新增低频实测另标记 `is_low_frequency_reconstruction=True`。

方向先由阵列F-K独立估计：同一传播方向组至少跨两个低频点，每点有足够不重叠窗口。此后才用曲线速度筛选相差。不能在每个频点重新挑一个最匹配曲线相位的方向；各频点的独立方向估计用于排除偏离已确认方向、反向传播、夹角越界、方向候选范围越界或角度误差超预算的窗口。

## 投影及相位约定

Δx是台对向量，p_w是每窗波矢单位方向，p_ref是已确认参考方向，d是原CCF所用台距。

| 量 | 定义 |
| --- | --- |
| 每窗有符号投影 | ℓ_w = Δx·p_w |
| 参考有符号投影 | ℓ_ref = Δx·p_ref |
| 沿波矢速度先验 | c_prior = c_app·\|ℓ_ref\|/d |
| 原始互谱 | C_w = G_目标·conj(G_参考) |
| 平面波预测相位 | −2πfℓ_w/c_prior |
| 共同参考互谱 | C_ref,w = C_w·exp[−i2πf(ℓ_ref−ℓ_w)/c_prior] |
| 叠加后相速度 | c = −2πfℓ_ref/φ_stack |

先把每窗互谱旋转到同一投影距离，再按预测相位筛选和复数叠加。不能叠加不同投影的相位后，再乘一个平均角度修正因子。保留投影符号区分正反传播；不能先取绝对值再混合。输出真实相速度 `phase_velocity_m_s` 和台对视速度 `apparent_pair_velocity_m_s`，并标记 `velocity_type=wave_vector_phase_velocity`。

原CCF源相位（例如π/4）及整数周期保存为 `origin_source_phase_rad`、`origin_constant_cycle`，用于解释原候选速度。公共源相位在瞬时单一平面波的两台互谱中抵消，原始平面波模型源相位取0；不能直接照搬CCF源相位。兼容表的 `constant_cycle` 保留原候选假设，实际解缠整周另存 `raw_phase_cycle`，二者含义不同。原CCF相位和bootstrap范围移至 `origin_*` 字段，不充当新相位或新速度的置信区间。

原始互谱、每窗投影、参考投影、旋转角、方向和掩码保存于窗口NPZ。`phase-audit` 从校正波形重算Goertzel互谱，再由台站坐标重算投影、相位旋转及筛选掩码。

新模式按物理传播方向分组，允许同一方向的包裹相位跨越±π，不以相位正负符号剔除该窗口；旧模式继续使用旧正负相位分支。

## 质量和反演

逐点检查接受窗口及不重叠窗口；对已固定方向集合做无自配对的循环错移背景检验，最多40个不同错移，实际次数另存。p≤0.05才能进入下一步。随后重查最少点数和最小连续带宽；解缠失败、速度非正、派生群速度非正或任何相邻升频增速的整条频带排除。不同曲线的频率网格不相互打断连续频带。

旧CCF自动选择保留台距/波长门槛。新投影Goertzel测量采用阵列方向、独立窗口、角度误差和错移背景标准，并用明确方法标记，避免原CCF门槛阻断全部低频重建。群速度由d(f/c)/df派生，标记 `group_velocity_is_derived=True`，不是新的独立群速度观测。

修正前曲线也投影到同一已确认方向，保证比较的是同一种速度。配对目录下分别输出：

| 目录 | 反演输入 |
| --- | --- |
| `before_common/` | 原曲线投影后的共同合格频点 |
| `after_common/` | 修正后的完全相同共同频点 |
| `after_extended/` | 修正后的完整实测频段，含新低频点 |

共同频段两组使用相同模态、层数、边界、种子和求解器；扩展频段结果单列。共同频段不足或拟合未通过时不生成合格结构；空结果清除旧结构，避免混用上轮结果。先验筛选、绝对整周、平面波和Rayleigh模态仍是条件假设，不能声称唯一真实结构。

## 操作和输出

第6页勾选“向低频外延先验，确认方向后用实测相位重建”，设置低频下限和步长；勾选时同步切换至精确Goertzel、曲线门限和方向筛窗。旧项目默认保持旧模式，需主动启用。处理后在候选页选曲线，再点“按当前选择反演修正前后曲线”；一键分析启用自动反演时也会执行配对。

```bash
python run.py phase-stack --project PROJECT --reconstruct-low --low-frequency-hz 0.2 --goertzel-frequency exact --taper hann --auto-invert
python run.py phase-audit --folder PROJECT/output/phase_guided
python run.py phase-compare --project PROJECT
```

可在 `phase_settings.json` 保存：

```json
{"frequency_method":"exact","screening_mode":"curve","direction":{"mode":"filter"},"reconstruction":{"enabled":true,"lower_frequency_hz":0.2,"step_hz":0.05}}
```

原生真实数据续跑 `ndstudio.real_compare.compare()` 默认新模式；复用下载、导入、预处理和扫描检查点，旧相位结果不作为投影重建完成证据。显式 `reconstruct_low=False` 或CLI `--legacy-phase` 可调用旧诊断路径。

主要报告：`phase_report.html`、`velocity_comparison/report.html`、`paired_comparison/report.html`；三组结构各存 `structure_layers.csv`。更新相位、曲线选择或反演参数后不能将旧对比报告当作本次结果。

新增路径在Python原生流程实现；旧MATLAB入口未宣称与此新投影路径等价，未做Windows原生或MATLAB运行时验证。本次不重跑DARE下载及完整现场流程；定向合成检查及具体结果见 `Validation_v256.md`。
