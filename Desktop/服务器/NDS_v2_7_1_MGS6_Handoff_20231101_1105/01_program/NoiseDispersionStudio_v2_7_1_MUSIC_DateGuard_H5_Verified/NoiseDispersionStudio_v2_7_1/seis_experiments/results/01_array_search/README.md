# 01 — Array search & aperture/frequency matching

Live-FDSN probing + array-search driver outputs. Metadata probes are cheap; driver statuses
summarize multi-GB downloads that are not committed.

- `catalog_probe.json`         — Parkfield box has only 2 broadband HHZ stns (preset fails).
- `parkfield_fail_probe.json`  — select_array failure for the preset candidate.
- `multibox_probe.json`        — scan of candidate dense boxes; LongValley & LA basin stand out.
- `geometry_dump.json`         — full pairwise geometry for LongValley & LA basin candidates.
- `longvalley_diag.json`       — compact 7-station subset geometry (aperture 16.2 km).
- `newregion_probe.json`       — slightly-larger-aperture (<30 km) candidates (Anza, LA basin...).
- `extend_probe.json`          — sparse stations around the LA core for aperture extension.
- `array_search_status_*.json` — end-to-end driver stage status for LongValley 3d, LA basin 3d/14d.

Key fact: NC waveforms live at NCEDC, CI waveforms at SCEDC (not IRIS). See ../../docs.
