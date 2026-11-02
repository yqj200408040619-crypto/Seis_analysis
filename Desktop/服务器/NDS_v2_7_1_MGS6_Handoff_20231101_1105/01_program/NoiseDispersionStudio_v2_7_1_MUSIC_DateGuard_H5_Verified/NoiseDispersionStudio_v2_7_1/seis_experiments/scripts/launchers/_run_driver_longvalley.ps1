$ErrorActionPreference = 'Continue'
$prog = 'C:\Users\游齐家\Desktop\服务器\NDS_v2_7_1_MGS6_Handoff_20231101_1105\01_program\NoiseDispersionStudio_v2_7_1_MUSIC_DateGuard_H5_Verified\NoiseDispersionStudio_v2_7_1'
Set-Location $prog
$env:PYTHONPATH = $prog
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
$log = Join-Path $prog 'runs\array_search\driver_longvalley.log'
"=== launch $(Get-Date -Format o) ===" | Out-File -FilePath $log -Encoding utf8
& py -3.14 -X utf8 -m ndstudio.array_search_e2e --workdir runs/array_search --start 2023-11-01 --end 2023-11-04 --only LongValley_NC *>> $log
"=== exit $LASTEXITCODE at $(Get-Date -Format o) ===" | Out-File -FilePath $log -Append -Encoding utf8
