# 2.6.0 空间层析数据协议

## A. 频散路径输入：tomography_paths.csv

每一行是一条台站路径在一个频率上的观测/条件频散点。至少应使用：

```text
pair_index,station_a,station_b,frequency_hz,phase_velocity_m_s,distance_m
```

地理层析还需要：

```text
longitude_a,latitude_a,longitude_b,latitude_b
```

2.5.6 低频方向重建产生的沿波矢真实相速度会额外保留：

```text
projected_distance_m,confirmed_azimuth_deg,velocity_type,measurement_method,direction_quality_verified
```

任何二维层析求解器都应读取并保留 `tomography_handoff_status`、背景检验和条件解释字段，不应静默丢掉周期/方向假设。

## B. 空间 Vs 输出

Studio 2.6.0 地热模块接受标准长表：

```text
longitude,latitude,depth_m,vs_m_s
```

建议二维/三维层析程序同时输出：

```text
vs_std_m_s,resolution,ray_count
```

其中 `resolution` 的定义必须由求解器自己的报告明确说明；Studio 不把任意数值自动解释成概率或置信度。

## C. 科学约束

- 1-D 条件反演不能生成二维空间信息。
- 深度切片只显示输入空间模型在该深度的值，不增加深度分辨率。
- 频散路径覆盖不足的区域必须通过 ray_count/resolution/其他诊断进行屏蔽或明确标记。
- 后续局部 Vs 反演须保留频率覆盖、模态、绝对周期和方向约束来源。
