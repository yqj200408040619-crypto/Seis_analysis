@echo off
cd /d "%~dp0"
set "PYTHONPATH=%~dp0"
set "PYTHONUTF8=1"
set "PYTHONIOENCODING=utf-8"
py -3.14 -X utf8 run_extended_direction.py --workdir runs/extended_la > "runs\extended_la\driver_ext.log" 2>&1
echo EXIT %ERRORLEVEL% >> "runs\extended_la\driver_ext.log"
