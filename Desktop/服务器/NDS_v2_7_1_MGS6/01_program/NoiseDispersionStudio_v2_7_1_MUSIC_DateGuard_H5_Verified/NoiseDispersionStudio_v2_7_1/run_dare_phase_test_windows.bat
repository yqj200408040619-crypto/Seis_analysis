@echo off
cd /d "%~dp0"
.venv\Scripts\python.exe run.py phase-stack --project examples\dare_phase_guided\project --phase-settings examples\dare_phase_guided\settings.json --aligned-h5 examples\dare_phase_guided\aligned_waveforms.h5 --output examples\dare_phase_guided\replay
if errorlevel 1 goto done
.venv\Scripts\python.exe run.py phase-audit --folder examples\dare_phase_guided\replay
:done
pause
