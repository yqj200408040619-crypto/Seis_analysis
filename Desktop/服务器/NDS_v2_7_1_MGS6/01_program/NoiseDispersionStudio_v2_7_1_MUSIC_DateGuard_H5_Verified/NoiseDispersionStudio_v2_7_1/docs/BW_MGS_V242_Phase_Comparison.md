# 2.4.2：Goertzel相位筛选叠加与原候选曲线对比

第6页“曲线相位筛选与叠加”完成后，自动计算最终相位、继承原曲线整数周期并反算速度，保存到 `output/phase_guided/velocity_comparison/`。点击“打开叠加曲线对比”查看；“重新生成叠加曲线对比”只重读已经保存的逐窗结果，不重算波形。

图中灰色虚线为原候选速度，棕色为幅值加权复数叠加 `arg(sum(C))` 的反算速度，蓝色为等权相位叠加 `arg(sum(C/abs(C)))` 的反算速度。灰带来自原片段bootstrap相位范围，属于条件范围，不是单窗校准置信区间。

已修复旧输出写死DARE名称、日期、频点数量和源相位的问题；曲线面板数随实际数量变化，本次N=1至N=9全部展示。窗口归档增加临时写入、关闭与CRC校验后提交；失败时保留旧文件。重新处理前清除已知旧对比输出，空叠加不补造速度。

## 本次MGS参数与核对

- UTC：2020-11-01T00:00:00Z至2020-11-05T00:00:00Z，96小时。
- 台对：BW_MGS01--BW_MGS04，距离18404.892 m。
- 频段0.270–0.340 Hz，15个频点，9个整数周期解释，135行曲线频点。
- 每窗20周期、滑动1周期、矩形窗none。
- 接受范围来自原片段bootstrap相位范围，额外余量0°；幅值及突变只标记。
- 保留59561 / 1584557个窗口频率评估；逐频保留率3.27%–4.26%。频率与滑动窗口相互依赖，这不是独立样本计数。
- 15份逐窗记录全部复算两种最终叠加；45个抽样窗口互谱从已对齐波形重新计算，所有相位门限掩码逐项核对通过。
- 速度回代误差小于1e-8度；原候选曲线与原结构文件经SHA256比对保持不变。
- 118项测试通过，12项跳过（GUI显示环境与可选独立求解器）；另外用Qt offscreen核对新增按钮。

相位门限由原曲线给出，最终速度接近输入曲线受这个规则约束；不能据此独立证明原频散正确或选择唯一N。原源相位45°沿用，本次没有再次进行结构反演。

## 复现

结果包包含响应校正并按UTC对齐的 `aligned_waveforms.h5`，可完全使用本程序离线重做逐频Goertzel筛选。将程序与BW_MGS_Results.zip解压到同一父目录，在程序根目录执行：

```bash
python run.py phase-stack --project ../BW_MGS_run/project --start 2020-11-01T00:00:00Z --end 2020-11-05T00:00:00Z --phase-range-source measured --phase-tolerance-deg 0 --repeat-cycles 20 --slide-cycles 1 --taper none --aligned-h5 ../BW_MGS_run/project/output/phase_guided/aligned_waveforms.h5
```

上述命令自动输出全部9条曲线。另生成N=2、N=3放大对比（全部曲线仍保留在报告中）：

```bash
python run.py phase-velocity --folder ../BW_MGS_run/project/output/phase_guided --output ../BW_MGS_run/project/output/phase_guided/velocity_comparison --curve-ids P00003_b1_S0_N2,P00003_b1_S0_N3
python run.py phase-audit --folder ../BW_MGS_run/project/output/phase_guided
```

不用已对齐H5时，先按原程序完整导入原始H5数据，再去掉 `--aligned-h5` 即由程序自行校正响应与对齐。只查看报告无需重跑数据。

## 输出

- `curve_comparison.csv`：逐频原速度、两种叠加速度、相对变化、最终相位、传播相位、距离与窗口数。
- `final_phase_velocity.csv`：以上结果的完整字段。
- `velocity_trends.csv`：逐曲线、逐叠加方式的速度趋势。
- `fv_curves.png/pdf`：全部候选的曲线对比；`fv_focus.png`：指定分支及百分比偏移。
- `final_phases.png`：实际最终互谱相位及接受范围。
- `report.html`：嵌入图片的完整对比报告。
- `verification.json`和上级`phase_audit.json`：复算及来源核对。
