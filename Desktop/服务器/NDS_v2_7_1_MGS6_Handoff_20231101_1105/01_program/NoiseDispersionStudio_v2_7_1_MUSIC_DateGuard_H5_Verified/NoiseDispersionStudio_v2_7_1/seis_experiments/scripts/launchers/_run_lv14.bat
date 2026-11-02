@echo off
cd /d "%~dp0"
set "PYTHONPATH=%~dp0"
set "PYTHONUTF8=1"
set "PYTHONIOENCODING=utf-8"
py -3.14 -X utf8 -m ndstudio.array_search_e2e --workdir runs/array_search_14d --start 2023-11-01 --end 2023-11-15 --only LongValley_NC > "runs\array_search_14d\driver_lv14.log" 2>&1
echo EXIT %ERRORLEVEL% >> "runs\array_search_14d\driver_lv14.log"
