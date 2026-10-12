# Noise Dispersion Studio 2.6.0：地热解释与空间层析接口

本版以 2.5.6 为基线。原下载、HDF5核查、预处理、扫描、错时背景、候选频散、方向/Goertzel、低频实测重建、条件结构反演全部保留；地热能力作为下游第 8 页加入，不反向改变成熟频散处理。

## 1. 全流程

背景噪声观测 → 频散提取 → 方向/低频实测修正 → 条件 Vs 结构 → 层析路径交接 → 空间 Vs 输入/层析结果 → 区域地质融合 → 地热异常解释 → 后续靶区评价。

重要边界：当前共同 1-D 条件结构不会被平铺成二维或三维 Vs。二维/三维 Vs 必须来自真实空间层析或明确外部空间模型。

## 2. 第 8 页输入

- 空间 Vs：CSV/TXT，或安装可选依赖后使用 NetCDF/GRD。
- 区域地质：Shapefile / GeoJSON / GeoPackage（需要 GeoPandas）。
- 断裂：Shapefile / GeoJSON / GeoPackage。
- 钻井：CSV 或矢量点。
- 温泉：CSV 或矢量点。
- Studio 自身 stations.csv 自动作为面波台站图层。

标准空间 Vs 长表：

```text
longitude,latitude,depth_m,vs_m_s
```

可选质量字段：

```text
vs_std_m_s,resolution,ray_count
```

也接受 `depth_km`、`vs_km_s`。若速度列只写 `vs`，必须在第 8 页明确单位，程序不猜。

## 3. 不确定性/分辨率门限

可设置：

- 最大相对 Vs 标准差：`vs_std_m_s / vs_m_s`
- 最小 `resolution`
- 最少 `ray_count`

缺少被要求的质量字段时，对应点不会被默认视为通过。被屏蔽的数据保留在 `spatial_vs_standardized.csv`，并记录 `quality_pass` 与 `quality_reasons`。

## 4. 层析路径交接

“导出层析射线路径”生成：

```text
output/geothermal/tomography_paths.csv
output/geothermal/tomography_paths_manifest.json
```

它直接读取当前有效频散交接源：若 2.5.6 的方向/Goertzel重建结果有效，则继承 `wave_vector_phase_velocity`、投影距离和方向质量字段；否则保留原条件频散及固定周期/源相位假设。

核心字段包括：

- curve_id / pair / branch / frequency_hz / period_s
- phase_velocity_m_s / group_velocity_m_s
- station_a / station_b
- 两端经纬度或本地 x/y
- distance_m / pair_azimuth_deg
- projected_distance_m / confirmed_azimuth_deg（若有）
- measurement_method / velocity_type
- background_discriminated / direction_quality_verified
- unique_physical_dispersion_verified
- tomography_handoff_status

此表是下一阶段二维相速度层析的输入协议，不代表每条条件曲线已经自动成为唯一物理频散。

## 5. 多深度 Vs 成果

“生成地热成果”将空间 Vs 按设置深度提取，并输出：

```text
output/geothermal/vs_000500m.csv
output/geothermal/figures/vs_000500m.png
output/geothermal/figures/vs_000500m.pdf
...
output/geothermal/report.html
```

规则网格优先使用 PyGMT（backend=auto 且运行时可用）；没有 GMT/PyGMT 时自动使用 Matplotlib。对不规则空间点，Matplotlib仅按原始点显示，不做无约束插值。

## 6. Windows 安装

基础 Studio：

```bat
install_windows.bat
```

CSV + Matplotlib 地热功能使用基础环境即可。若需要 Shapefile/GeoJSON、NetCDF/GRD 和 PyGMT：

```bat
install_geothermal_windows.bat
```

PyGMT 还需要系统中可用的 GMT 6.5+。如果 GMT 不可用，第 8 页仍可选择 Matplotlib 后端。

## 7. CLI

```bat
.venv\Scripts\python.exe run.py tomography-export --project "C:\path\project"
.venv\Scripts\python.exe run.py geothermal-validate --project "C:\path\project"
.venv\Scripts\python.exe run.py geothermal --project "C:\path\project"
```

仅导出当前 selection.json 中曲线：

```bat
.venv\Scripts\python.exe run.py tomography-export --project "C:\path\project" --selected-only
```

合成空间 Vs 示例：

```bat
run_geothermal_example_windows.bat
```

## 8. 下一阶段

2.6.0 已经固定“频散路径 → 空间层析 → 深度 Vs → 地质解释”的接口。下一步适合在该接口上增加：

1. 直线路径二维相速度层析（按频率/周期）；
2. 路径覆盖、阻尼和分辨率诊断；
3. 每个网格点的局部频散曲线及 Vs(z) 反演；
4. 生成标准 `longitude, latitude, depth_m, vs_m_s` 三维长表；
5. 再引入 MT、重力、热流和钻井测温，建立多证据地热靶区评价。

这些步骤必须保留路径覆盖和不确定性，不用“低 Vs = 地热”的单一规则代替地质解释。
