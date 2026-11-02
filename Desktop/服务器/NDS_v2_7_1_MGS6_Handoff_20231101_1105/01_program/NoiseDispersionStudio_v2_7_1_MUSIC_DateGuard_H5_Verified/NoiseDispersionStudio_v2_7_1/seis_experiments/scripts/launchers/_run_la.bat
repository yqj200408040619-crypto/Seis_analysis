@echo off
cd /d "%~dp0"
set "PYTHONPATH=%~dp0"
set "PYTHONUTF8=1"
set "PYTHONIOENCODING=utf-8"
py -3.14 -X utf8 -m ndstudio.array_search_e2e --workdir runs/array_search_labasin --start 2023-11-01 --end 2023-11-04 --only LABasin_CI > "runs\array_search_labasin\driver_la.log" 2>&1
echo EXIT %ERRORLEVEL% >> "runs\array_search_labasin\driver_la.log"
