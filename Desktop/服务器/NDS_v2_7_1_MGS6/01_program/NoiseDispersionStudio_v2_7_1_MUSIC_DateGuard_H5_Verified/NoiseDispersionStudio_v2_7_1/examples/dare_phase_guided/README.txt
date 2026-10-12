在程序目录运行 run_dare_phase_test_windows.bat 即可离线重放；或运行：
.venv\Scripts\python.exe run.py phase-stack --project examples\dare_phase_guided\project --phase-settings examples\dare_phase_guided\settings.json --aligned-h5 examples\dare_phase_guided\aligned_waveforms.h5 --output examples\dare_phase_guided\replay
本示例含7条DARE真实候选与30分钟校正波形，264个去重逐频任务。
默认自动使用此H5的全部30分钟，按附件共同起点顺延600秒后遍历剩余20分钟。
这里没有四天原始波形，不能以此证明多日原始下载已全部重跑。
若要检查附件原门限，复制settings.json并把screening_mode改为matlab，使用另一个输出目录。
若希望不跳过起点，把stabilization_s改为0；严格曲线频点计算改frequency_method为exact。
两个相位分支完整统计见phase_branch_summary.csv。
