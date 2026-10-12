@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  echo Run install_windows.bat first.
  pause
  exit /b 1
)
echo Installing optional geospatial Python packages...
.venv\Scripts\python.exe -m pip install -r requirements-geothermal.txt
if errorlevel 1 goto failed
echo.
where gmt >nul 2>nul
if errorlevel 1 (
  echo NOTE: GMT executable was not found in PATH.
  echo The integrated geothermal page can still use the Matplotlib backend.
  echo For the PyGMT backend, install GMT 6.5+ separately and make it visible to PyGMT.
) else (
  gmt --version
)
.venv\Scripts\python.exe run.py doctor
echo Optional geothermal Python stack installation finished.
pause
exit /b 0
:failed
echo Optional geothermal installation failed. Keep this window's error.
pause
exit /b 1
