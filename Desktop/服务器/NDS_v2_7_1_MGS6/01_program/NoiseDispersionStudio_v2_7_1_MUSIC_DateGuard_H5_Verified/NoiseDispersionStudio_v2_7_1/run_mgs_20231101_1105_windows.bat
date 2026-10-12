@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
 echo Run install_windows.bat first.
 pause
 exit /b 1
)
.venv\Scripts\python.exe -m ndstudio.mgs_e2e --workdir "%~dp0MGS6_20231101_1105" --start 2023-11-01 --end 2023-11-05
pause
