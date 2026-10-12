Noise Dispersion Studio 2.5.6 / 2.5.6.2
Windows JSON Invalid \\escape 热修复补丁

修复目标
========
报错示例：
json.decoder.JSONDecodeError: Invalid \\escape: line 5 column 32 (...)

该错误发生在 import 已完成后读取 hdf_import_audit.json / 追溯摘要时。
Windows 路径若以单反斜杠写入旧 JSON，例如：
C:\Users\...（文件里实际只有单个反斜杠）
会被 JSON 解析器解释为非法转义。

补丁内容
========
1. ndstudio\trace.py
   - 正常 JSON 仍按原方式读取。
   - 仅当 JSON 解析失败时，尝试修复 Windows 路径反斜杠。
   - 对生成型 JSON 自动保留 *.invalid_escape_backup 备份后重写为标准 JSON。
   - 追溯/摘要 JSON 解析失败不再让已经完成的科学处理阶段整体失败。

2. tools\repair_json_escape.py
   - 可手动扫描并修复现有项目中的旧 JSON。

3. repair_project_json.bat
   - Windows 下交互式调用上述修复工具。

安装方法
========
A. 先完全关闭 Noise Dispersion Studio 主窗口。
B. 将本补丁中的 ndstudio 文件夹复制到原程序根目录：
   NoiseDispersionStudio_v2_5_6\
   选择“替换目标中的文件”。
C. tools 和 repair_project_json.bat 也可一并复制到程序根目录。
D. 重新启动 run_windows.bat（或你平时启动 Studio 的方式）。

不需要
======
- 不需要重启 Windows。
- 不需要重新运行 install_windows.bat。
- 不需要重新下载 HDF5。
- 不需要删除 canonical、output 或项目 config.json。

出错项目的推荐恢复顺序
====================
1. 重启 Studio。
2. 打开原来的同一个项目。
3. 先运行“仅导入与核查”一次，确认 import 正常结束。
4. 再运行“一键分析”。

如果旧 hdf_import_audit.json 仍报错：
- 双击程序根目录的 repair_project_json.bat；
- 粘贴项目目录，例如：
  C:\Users\...\NoiseDispersionStudio_v2_5_6\downloads_dispersion_project
- 修复完成后重新打开 Studio，再跑“仅导入与核查”。

说明
====
这不是分析参数修改补丁，不会改变频散、Goertzel、反演或速度参数。
它只修复 Windows JSON 路径转义和追溯摘要容错。
