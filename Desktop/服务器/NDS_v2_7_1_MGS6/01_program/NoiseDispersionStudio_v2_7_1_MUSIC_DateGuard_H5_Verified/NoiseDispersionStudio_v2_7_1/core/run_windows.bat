@echo off
cd /d "%~dp0"
if "%~1"=="" (
  echo Usage: run_windows.bat path_to_config.json
  pause
  exit /b 1
)
python run.py all --config "%~1"
pause
