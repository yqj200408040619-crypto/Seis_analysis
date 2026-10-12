# 2.6.0 变更摘要

基线：完整 Noise Dispersion Studio 2.5.6。

## 合并

- 主窗口新增第 8 页：地热与空间层析。
- 地热配置统一保存为 `geothermal_settings.json`。
- 主 Studio 报告自动链接地热报告与层析交接表。
- 新增 `tomography-export`、`geothermal-validate`、`geothermal`、`geothermal-example` CLI。
- Windows 增加可选地热依赖安装脚本，不改变基础 Studio 安装方式。

## 空间层析交接

- 导出 `tomography_paths.csv`，包含频率、速度、台站端点、距离、方位、投影距离与科学质量状态。
- 显式 `pair` 标签优先匹配当前台站几何；历史 `pair_index` 仅作审计，避免错误射线映射。
- 有效的方向/Goertzel `inversion_curves.csv` 优先进入层析交接。

## 地热成果

- 空间 Vs 标准协议：`longitude, latitude, depth_m, vs_m_s`。
- 可选质量字段：`vs_std_m_s, resolution, ray_count`。
- CSV/TXT 直接支持；NetCDF/GRD 通过可选 xarray/netCDF4 支持。
- Matplotlib 基础后端；PyGMT/GMT 可选后端。
- 多深度 Vs 切片与地质、断裂、台站、钻井、温泉叠加。
- 输出独立 HTML 报告并把主要图件同步到 Studio `output/figures`。

## 科学限制

- 不把共同 1-D 条件结构转换成伪 2-D/3-D Vs。
- 不对不规则空间 Vs 做无约束插值。
- 不把低 Vs 单独等同于地热靶区。
- 层析路径交接保留周期/方向/背景检验等状态，后续求解器不得静默丢弃。
