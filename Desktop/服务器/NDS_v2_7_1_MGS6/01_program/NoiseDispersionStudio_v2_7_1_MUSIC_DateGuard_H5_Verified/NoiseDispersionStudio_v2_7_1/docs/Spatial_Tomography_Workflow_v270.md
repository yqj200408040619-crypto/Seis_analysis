# 2.7.0 二维相速度层析与局部频散交接

## 1. 本阶段在完整流程中的位置

```text
背景噪声观测
  -> 预处理 / 扫描 / 条件频散
  -> 方向确认 + Goertzel + 低频实测修正
  -> tomography_paths.csv
  -> 逐频二维相速度层析 c(x,y,f)
  -> ray coverage / resolution / checkerboard
  -> local_dispersion.csv
  -> 下一阶段：逐网格 Vs(z) 反演
  -> Vs(x,y,z)
  -> 地质 / 断裂 / 钻井 / 温泉 / MT / 重力 / 热流融合
```

本版完成到 `local_dispersion.csv`。它**没有**把二维相速度直接叫做 Vs，也没有把已有共同 1-D 结构平铺为空间模型。

## 2. 输入数据和分支唯一性

求解器读取 `output/geothermal/tomography_paths.csv`。默认要求同一台对、同一频率只有一条观测。如果候选页仍保留多个整数周期/源相位竞争分支，二维层析会停止并要求先在 `selection.json` 中确定要检验的条件分支。

可以明确设置 `tomo_duplicate_policy=median`，但这只是用户指定的合并策略，程序不会把它写成“唯一真实频散”。

默认 `tomo_quality_mode=default_pass`，保留当前频散交接中的背景检验等质量标志；`strict_direction` 还要求方向质量通过；`all_valid` 仅做数值有效性检查，主要用于诊断。

## 3. 方向投影数据如何进入直线路径层析

2.5.6 的 Goertzel 新链路输出沿已确认波矢的真实相速度：

```text
velocity_type = wave_vector_phase_velocity
```

它不能不经处理就当成台站 A-B 直线路径速度。本版按相位时延保持一致：

```text
Δt_phase = |projected_distance| / c_wavevector
c_pair_equivalent = distance_AB / Δt_phase
```

若已有 `apparent_pair_velocity_m_s`，优先直接使用；否则由投影距离和真实波矢速度计算。该等效量只用于直线路径方程。默认还要求波矢轴与台对轴夹角不大于 `tomo_max_direction_angle_deg=15°`；强斜入射数据不会被硬塞进直线路径层析。

这与前面的 45° Goertzel 窗口筛选不是一回事：45°是相位重建阶段的窗口质量上限；二维直线路径层析默认采用更严格的近共线门限，因为两种模型的物理假设不同。

## 4. 网格和正演

每个频率独立建立规则矩形网格，默认：

```text
tomo_nx = 12
tomo_ny = 12
tomo_margin_fraction = 0.05
```

经纬度项目先用以阵列中心为原点的局部方位等距投影（AEQD）转换到米制坐标。每条台站线段与每个网格单元的交长度通过精确线段-网格边界分割计算，不用简单“采样点计数”近似。

正演方程：

```text
t = G s
```

其中 `G` 为各路径在单元内的长度（km），`s` 为慢度（s/km）。路径矩阵行会缩放到项目记录的实际 `distance_m`，消除局部投影与大地距离之间的微小数值差。

## 5. 反演

以每频率观测速度的中位慢度作为参考：

```text
s = s_ref + δs
```

求解：

```text
min || W(G δs - δt) ||²
  + λ_d² ||δs||²
  + λ_s² ||L δs||²
```

- `λ_d`：阻尼，对应 `tomo_damping`；
- `λ_s`：相邻网格一阶平滑，对应 `tomo_smoothing`；
- `L`：水平/垂直相邻单元差分矩阵；
- 如果输入含速度标准差，转换为时延标准差做相对权重；否则等权。

用户填写的是相对正则系数，程序按设计矩阵正列范数的中位数自动转换为实际系数，避免网格单位改变导致数量级失控。

## 6. 必须随结果查看的诊断

每个频率输出：

- `travel_time_rms_s`：时延拟合 RMS；
- `ray_count`：穿过单元的路径数；
- `path_length_m`：单元内累计路径长度；
- `resolution`：当前阻尼/平滑下模型分辨率矩阵对角元；
- `phase_velocity_std_m_s`：线性化、正则化条件下的误差尺度估计；
- checkerboard 相关系数与振幅恢复比。

`resolution` 和标准差不是概率意义的置信区间；checkerboard 也不是“通过/失败”的万能判据。解释时应结合射线方向覆盖、路径数量、频率覆盖和正则参数敏感性。

## 7. 输出

目录：

```text
output/geothermal/tomography2d/
```

主要文件：

```text
phase_velocity_tomography.csv   # c(x,y,f)
local_dispersion.csv            # 每个网格单元的 c(f)，供下一阶段 Vs(z)
data_fit.csv                    # 每条路径实测/预测时延及残差
frequency_catalog.csv           # 每频率路径数、RMS、覆盖、分辨率、checkerboard
checkerboard_summary.csv
checkerboard/checkerboard_*.csv
grid_definition.json
tomography_2d_complete.json
tomography_2d_report.html
figures/phase_velocity_*.png/.pdf
```

## 8. GUI 与 CLI

GUI 第 8 页：

1. 先在候选/反演流程中确定要检验的条件分支；
2. 点“导出层析射线路径”；
3. 设置网格、质量模式、最大方向夹角、阻尼/平滑、最少路径；
4. 点“运行二维相速度层析”；
5. 打开二维层析报告，检查 coverage / resolution / checkerboard；
6. `local_dispersion.csv` 作为下一阶段节点 Vs 反演输入。

CLI：

```bash
python run.py tomography-export --project PROJECT --selected-only
python run.py tomography-invert --project PROJECT
python run.py tomography-example --directory examples/generated_tomography
```

## 9. 科学限制

- 直线路径是假设，不适合强横向梯度导致明显弯曲射线的情况；后续可扩展 bending/FMM 路径更新。
- 一个台站对/频率的多个周期别名不能当成独立观测重复增加权重。
- 平面波方向测量与台间路径积分不是同一个观测算子，因此强斜入射数据默认排除。
- 当前输出是二维相速度频率切片，不是三维 Vs。
- 下一阶段逐网格 Vs(z) 反演仍必须保留模态、频带、频散不确定性和深度分辨率限制。
