@echo off
setlocal
cd /d "%~dp0"
where py >nul 2>nul
if errorlevel 1 (
  echo Install 64-bit Python 3.12 from python.org, including the Python launcher.
  pause
  exit /b 1
)
if not exist ".venv\Scripts\python.exe" (
  py -3.12 -m venv .venv
  if errorlevel 1 goto failed
)
.venv\Scripts\python.exe -m pip install --upgrade pip
if errorlevel 1 goto failed
.venv\Scripts\python.exe -m pip install -r requirements.txt
if errorlevel 1 goto failed
.venv\Scripts\python.exe run.py doctor
if errorlevel 1 goto failed
echo Installation complete. Open run_windows.bat.
pause
exit /b 0
:failed
echo Installation failed. Keep this window's error and retry after resolving it.
pause
exit /b 1
