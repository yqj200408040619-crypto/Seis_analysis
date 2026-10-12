@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  echo Run install_windows.bat first.
  pause
  exit /b 1
)
if exist "examples\generated_tomography" rmdir /s /q "examples\generated_tomography"
.venv\Scripts\python.exe run.py tomography-example --directory "examples\generated_tomography"
if errorlevel 1 pause
