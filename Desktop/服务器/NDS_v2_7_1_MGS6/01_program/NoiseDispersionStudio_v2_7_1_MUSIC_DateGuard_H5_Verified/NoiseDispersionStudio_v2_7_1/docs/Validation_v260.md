# Noise Dispersion Studio 2.6.0 — 地热/空间层析合并验证

## 基线

本版直接从用户提供的完整 `Studio_v2_5_6_Program` 源码合并，不以此前独立 Geothermal Mapper 包替换核心流程。2.5.6 的下载、H5 数据审计、预处理、扫描、错时背景、候选族、方向/Goertzel、低频实测重建、结构反演、H/V 和报告模块均保留。

## 新增链路

```text
现有 Studio 频散结果
        ↓
当前有效曲线选择
（方向约束 inversion_curves 优先）
        ↓
tomography_paths.csv
        ↓
二维相速度层析【下一求解阶段接口】
        ↓
longitude / latitude / depth / Vs
        ↓
不确定性 / resolution / ray_count 遮罩
        ↓
地质 + 断裂 + 钻井 + 温泉
        ↓
多深度 PNG/PDF/CSV + HTML 报告
```

## 本轮自动测试

执行：

```text
PYTHONPATH=.:core pytest -q tests/test_geothermal.py
```

结果：`3 passed`。

覆盖：

1. 空间 Vs 的 km/s、km 深度统一换算到 m/s、m，并验证 Vs 标准差、resolution、ray_count 三类质量门限。
2. 层析射线导出时优先依据显式 `pair` 标签重新匹配当前台站几何；即使历史 `pair_index` 已过期，也不会把频散挂到错误台对上，同时保留 `source_pair_index` 与 `geometry_pair_index` 供审计。
3. 自包含空间 Vs 示例用 Matplotlib 生成两个深度切片，PNG/PDF/CSV、标准化空间表和 HTML 报告全部落盘。

另外使用 2.5.6 自带 `examples/direction_plane/project` 定向验证层析交接：程序优先读取当前有效的 `output/phase_guided/inversion_curves.csv`，并发现示例中的历史 `pair_index` 与当前 station ordering 不一致；修正后 `S0--S2`、`S1--S3` 均按显式台对名称匹配到正确几何。这一检查避免了静默的射线路径错配。

## 当前环境限制

本容器缺少 ObsPy 与 PySide6，因此不能在这里重新执行全部 2.5.6 原始波形和 Qt GUI 回归。新增 Python 文件和被修改入口已通过 `py_compile`；地热核心定向测试可在当前环境完成。Windows GUI、GMT DLL 与 PyGMT 仍应在 Windows 目标机上运行一次 `doctor`/示例进行平台验证。

## 科学边界

- 1-D 条件结构只作为 1-D 结果；不生成伪空间 Vs。
- `tomography_paths.csv` 是带科学状态字段的交接表，不代表每条条件候选都已经被证明为唯一物理频散。
- 低频方向/Goertzel结果保留方向、错时背景、投影距离等字段，后续层析必须继续传播这些质量信息。
- 空间 Vs 切片只显示提供的空间反演结果；程序不为不规则输入做无约束插值。
- 低 Vs 异常不能单独作为“可打井”结论，后续需和断裂、热流、MT、重力及钻井温度等独立证据联合解释。
