@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  echo Run install_windows.bat first.
  pause
  exit /b 1
)
if exist "examples\generated_geothermal" rmdir /s /q "examples\generated_geothermal"
.venv\Scripts\python.exe run.py geothermal-example --directory "examples\generated_geothermal"
if errorlevel 1 pause
