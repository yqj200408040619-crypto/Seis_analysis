NDS v2.7.1 日期目录防混用补丁（H5 实际波形校验版）

将压缩包内 ndstudio/date_guard.py、ndstudio/mgs_e2e.py 覆盖到程序根目录相同路径。tests/test_date_guard.py 是回归测试，可选放入根目录 tests/ 后运行 pytest。

保护逻辑会在下载前和续跑时校验工作目录请求时段、下载清单、H5 内原始请求时段/分块/台站号，以及实际波形起止时间；检查已对齐波形缓存是否属于本次时段。MiniSEED 记录边界允许最多 60 秒合理补齐。检测不符时停止并提示使用新工作目录。

Windows：关闭 Studio，将文件覆盖到程序根目录，再重新打开 Studio。运行四天验证脚本时按 README 指引运行 run_mgs_20231101_1105_windows.bat。
验证：专项 8 项通过；程序回归总计 214 项通过、8 项跳过、0 项失败。
