@echo off
setlocal
cd /d "%~dp0"
set /p PROJECT=Paste Studio project folder path and press Enter: 
if "%PROJECT%"=="" exit /b 1
if exist ".venv\Scripts\python.exe" (
  ".venv\Scripts\python.exe" tools\repair_json_escape.py "%PROJECT%"
) else (
  py tools\repair_json_escape.py "%PROJECT%"
)
pause
