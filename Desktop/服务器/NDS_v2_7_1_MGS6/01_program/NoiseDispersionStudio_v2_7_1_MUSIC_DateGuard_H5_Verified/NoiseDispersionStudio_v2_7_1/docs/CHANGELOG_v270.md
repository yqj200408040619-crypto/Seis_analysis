# 2.7.0 变更摘要

基线：2.6.0 地热整合版；原 2.5.6 成熟下载、预处理、扫描、方向/Goertzel、低频重建和条件结构反演不改。

## 新增二维相速度层析

- 新增 `ndstudio/tomography.py`。
- 新增 CLI：`tomography-invert`、`tomography-example`。
- 第 8 页新增直线路径二维层析设置和“运行二维相速度层析”按钮。
- 经度纬度自动转局部 AEQD 米制坐标；本地 x/y 也可直接处理。
- 精确计算直线在每个矩形单元中的路径长度。
- 逐频阻尼 + 一阶平滑慢度反演。
- 输出 travel-time residual、ray_count、累计 path length、resolution diagonal、速度标准差尺度。
- 每频率可自动运行 checkerboard 恢复测试。
- 自动形成 `local_dispersion.csv`，作为下一阶段网格节点 Vs(z) 反演的统一输入。

## 方向数据的物理处理

- `wave_vector_phase_velocity` 不直接当作 A-B 路径速度。
- 优先使用 `apparent_pair_velocity_m_s`；否则由 `projected_distance_m` 与真实波矢速度转换成台对等效相位时延速度。
- 默认仅接受波矢轴/台对轴夹角 <=15° 的投影测量进入直线路径反演。
- 该 15° 门限属于二维层析观测算子约束，不替代 2.5.6 Goertzel 阶段的 45°筛窗。

## 分支和质量保护

- 同一台对/频率存在多个竞争候选时默认报错，不把周期别名重复计权。
- 可显式使用 `median` 合并策略，但报告中保留条件解释警告。
- 选定曲线存在时，层析求解自动优先使用 selection.json 对应的路径交接。

## 地热报告

- 主 Studio 报告链接二维相速度层析报告和 `local_dispersion.csv`。
- 地热报告会显示已完成的二维相速度层析状态，但仍明确 c(x,y,f) != Vs(x,y,z)。
