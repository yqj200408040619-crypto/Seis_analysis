NDS v2.7.1 MUSIC - 2023-11-01 to 2023-11-05 validation launcher

1. Run install_windows.bat once if .venv does not exist.
2. Close any existing Studio/download window.
3. Double-click run_mgs_20231101_1105_windows.bat.
4. The work directory is MGS6_20231101_1105. Do not reuse a 2022 workdir.

The native workflow requests BW.MGS01-MGS05 and BW.MGSBH, HHZ, start inclusive
2023-11-01T00:00:00Z and end exclusive 2023-11-05T00:00:00Z (96 h), then runs:
array selection -> lossless H5 acquisition -> preprocessing/CCF scan -> front MUSIC
per-frequency direction and projected velocity -> review/verify -> Goertzel 200-cycle
low-frequency reconstruction with MUSIC direction -> phase audit -> qualified band selection
-> paired inversion only if a physical continuous corrected band survives.

Safety fix in this patch: --start/--end are explicit, and an existing workdir from a
different time range is rejected so old H5 cannot be mixed into the 2023 run.
