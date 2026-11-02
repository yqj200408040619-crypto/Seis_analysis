@echo off
cd /d "%~dp0"
set "PYTHONPATH=%~dp0"
set "PYTHONUTF8=1"
set "PYTHONIOENCODING=utf-8"
set "EXT_START=2023-11-01"
set "EXT_END=2023-11-15"
py -3.14 -X utf8 run_extended_direction.py --workdir runs/extended_la_14d > "runs\extended_la_14d\driver_ext14.log" 2>&1
echo EXIT %ERRORLEVEL% >> "runs\extended_la_14d\driver_ext14.log"
