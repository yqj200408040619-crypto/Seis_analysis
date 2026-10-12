# 2.7.0 二维相速度层析定向验证

本轮只验证新增/修改链路，不重复消耗时间重跑 2.5.6 已成熟的下载、预处理、扫描和反演历史回归。

## 1. 自动测试

执行：

```text
pytest -q tests/test_geothermal.py tests/test_tomography2d.py
```

结果：

```text
8 passed
```

覆盖：

- 线段穿越规则网格时总路径长度守恒；
- 方向投影真实波矢速度按相位时延转换为台对等效速度；
- 超过二维层析最大方向夹角的数据被拒绝；
- 同台对/同频率竞争分支默认拒绝重复计权；
- `selection.json` 的标准 `bands[].curve_id` 可按 case 正确提取；
- 多频直线路径慢度反演；
- resolution diagonal、ray count、局部频散和 checkerboard 输出；
- 原 2.6.0 地热空间 Vs 输出回归仍通过。

## 2. 自包含合成二维层析

命令：

```text
python run.py tomography-example --directory <workdir>
```

合成设置：8 个边界台站，20 条长路径/频率，4×4 网格，3 个频率 0.8 / 1.2 / 1.8 Hz。观测由已知 checker 型慢度场正演得到。

结果：

| 频率 Hz | 路径 | 时延 RMS s | 覆盖格 | 中位 resolution | checker 相关 | 振幅恢复比 |
|---:|---:|---:|---:|---:|---:|---:|
| 0.8 | 20 | 0.02101 | 16/16 | 0.9638 | 0.99994 | 0.98951 |
| 1.2 | 20 | 0.02358 | 16/16 | 0.9638 | 0.99994 | 0.98951 |
| 1.8 | 20 | 0.02685 | 16/16 | 0.9638 | 0.99994 | 0.98951 |

这是受控合成数值验证，不代表真实地质分辨率。

## 3. 与 Studio 真实结构的交接保护

使用包内 `examples/direction_plane/project` 的实际 Studio 相位输出复制测试：

- `tomography-export` 成功识别方向约束后的 `inversion_curves.csv`；
- 导出 2 条台站路径、2 条曲线、12 个频点；
- 每个频率只有 2 条路径；
- 默认要求每频率至少 6 条路径，因此 1.3–1.8 Hz 六个频率全部明确跳过；
- 最终不生成二维层析结果，并报“没有频率达到最少路径数”。

这验证了程序不会因为用户点击“运行二维层析”就在明显欠定的两条射线数据上制造彩色速度图。

## 4. 原地热输出回归

`geothermal-example` 仍成功生成 500 / 1000 / 2000 / 4000 m 四个空间 Vs 切片及 PNG/PDF/CSV/HTML。没有候选频散的纯空间 Vs 示例只记录“层析射线路径未导出”提示，不影响地热制图。

## 5. 尚未宣称完成的内容

- 本轮未用新的现场大阵列重跑完整下载到二维层析；包内 direction_plane 的射线数故意不足，正确拒绝。
- 当前为 straight-ray；没有实现弯曲射线/FMM 更新。
- 当前输出是 `c(x,y,f)` 与 `local_dispersion.csv`，尚未完成网格节点 `Vs(z)` 批量反演。
- 没有把 checkerboard、resolution 或线性化标准差解释为概率置信区间。
- Windows 原生 GUI 仍建议用户目标机运行 `run_tomography_example_windows.bat` 做最终平台确认。
