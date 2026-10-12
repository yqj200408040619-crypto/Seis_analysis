@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
 echo Run install_windows.bat first.
 pause
 exit /b 1
)
.venv\Scripts\python.exe -m ndstudio.mgs_e2e --workdir "%~dp0MGS6_run" %*
pause
