"""Geothermal interpretation and spatial-tomography handoff.

This module is intentionally downstream of the mature ambient-noise workflow.
It never converts a common 1-D conditional structure into a fake 2-D/3-D
velocity volume. Spatial Vs products must come from a real spatial inversion or
an explicitly supplied spatial model.
"""
from __future__ import annotations

from pathlib import Path
import json, math, shutil, html
from typing import Any, Optional

import numpy as np
import pandas as pd

from dispersion.common import write_json, sha, load_config


DEFAULT_SETTINGS: dict[str, Any] = {
    "tomography_file": "",
    "tomography_variable": "vs_m_s",
    "tomography_vs_unit": "auto",      # auto / m/s / km/s (for ambiguous 'vs')
    "grid_depth_m": None,               # only for a 2-D grid with no depth coordinate
    "slice_depths_m": [500.0, 1000.0, 2000.0, 5000.0],
    "depth_tolerance_m": None,
    "backend": "auto",                 # auto / pygmt / matplotlib
    "cmap": "turbo",
    "geology": "",
    "faults": "",
    "wells": "",
    "hot_springs": "",
    "label_points": True,
    "selected_paths_only": False,
    "max_relative_vs_std": None,
    "min_resolution": None,
    "min_ray_count": None,
    "title": "Geothermal interpretation",
    # Straight-ray 2-D phase-velocity tomography (next-stage solver).
    "tomography_paths_file": "",
    "tomo_selection_case": "",          # optional selection.json case name/index
    "tomo_quality_mode": "default_pass",  # all_valid / default_pass / strict_direction
    "tomo_duplicate_policy": "error",     # error / median
    "tomo_max_direction_angle_deg": 15.0,
    "tomo_nx": 12,
    "tomo_ny": 12,
    "tomo_margin_fraction": 0.05,
    "tomo_damping": 0.15,
    "tomo_smoothing": 1.0,
    "tomo_min_paths_per_frequency": 6,
    "tomo_min_cell_rays": 1,
    "tomo_checkerboard": True,
    "tomo_checkerboard_amplitude_fraction": 0.08,
    "tomo_checkerboard_block_cells": 2,
    "tomo_max_figures": 12,
}


class GeothermalError(RuntimeError):
    pass


def settings_path(project: str | Path) -> Path:
    return Path(project).resolve() / "geothermal_settings.json"


def load_settings(project: str | Path) -> dict[str, Any]:
    path = settings_path(project)
    saved = json.loads(path.read_text(encoding="utf-8-sig")) if path.exists() else {}
    unknown = set(saved) - set(DEFAULT_SETTINGS)
    if unknown:
        raise GeothermalError("未知 geothermal_settings 参数: " + str(sorted(unknown)))
    result = dict(DEFAULT_SETTINGS)
    result.update(saved)
    _validate_settings(result)
    return result


def save_settings(project: str | Path, value: dict[str, Any]) -> Path:
    unknown = set(value) - set(DEFAULT_SETTINGS)
    if unknown:
        raise GeothermalError("未知 geothermal_settings 参数: " + str(sorted(unknown)))
    result = dict(DEFAULT_SETTINGS)
    result.update(value)
    _validate_settings(result)
    path = settings_path(project)
    write_json(path, result)
    return path


def _validate_settings(s: dict[str, Any]) -> None:
    if s["backend"] not in ("auto", "pygmt", "matplotlib"):
        raise GeothermalError("backend 必须为 auto / pygmt / matplotlib")
    if s["tomography_vs_unit"] not in ("auto", "m/s", "km/s"):
        raise GeothermalError("tomography_vs_unit 必须为 auto / m/s / km/s")
    if s.get("tomo_quality_mode") not in ("all_valid", "default_pass", "strict_direction"):
        raise GeothermalError("tomo_quality_mode 必须为 all_valid / default_pass / strict_direction")
    if s.get("tomo_duplicate_policy") not in ("error", "median"):
        raise GeothermalError("tomo_duplicate_policy 必须为 error / median")
    for key in ("tomo_nx", "tomo_ny", "tomo_min_paths_per_frequency", "tomo_min_cell_rays", "tomo_checkerboard_block_cells", "tomo_max_figures"):
        if int(s.get(key, 0)) < 1:
            raise GeothermalError(key + " 必须为正整数")
    for key in ("tomo_margin_fraction", "tomo_damping", "tomo_smoothing", "tomo_max_direction_angle_deg", "tomo_checkerboard_amplitude_fraction"):
        if not np.isfinite(float(s.get(key))) or float(s.get(key)) < 0:
            raise GeothermalError(key + " 必须为有限非负数")
    if float(s.get("tomo_max_direction_angle_deg")) > 90:
        raise GeothermalError("tomo_max_direction_angle_deg 不能超过90度")
    depths = s.get("slice_depths_m") or []
    if not isinstance(depths, list):
        raise GeothermalError("slice_depths_m 必须是数值列表")
    for d in depths:
        if not np.isfinite(float(d)) or float(d) < 0:
            raise GeothermalError("slice_depths_m 必须为有限非负深度")
    for key in ("depth_tolerance_m", "grid_depth_m", "max_relative_vs_std", "min_resolution", "min_ray_count"):
        v = s.get(key)
        if v is not None and (not np.isfinite(float(v)) or float(v) < 0):
            raise GeothermalError(key + " 必须为空或有限非负数")


def resolve_project_path(project: str | Path, value: str | Path | None) -> Optional[Path]:
    if value is None or str(value).strip() == "":
        return None
    p = Path(value)
    if not p.is_absolute():
        p = Path(project).resolve() / p
    return p.resolve()


def _truth(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float)) and value in (0, 1):
        return bool(value)
    return str(value).strip().lower() in ("true", "1", "yes", "y")


def _selected_curve_ids(root: Path, selection_case: str | int | None = None) -> set[str]:
    path = root / "selection.json"
    if not path.exists():
        return set()
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    cases = list(obj.get("cases", []))
    if selection_case not in (None, ""):
        key = str(selection_case).strip()
        chosen = [c for c in cases if str(c.get("name", "")) == key]
        if not chosen and key.isdigit():
            idx = int(key)
            if 0 <= idx < len(cases): chosen = [cases[idx]]
            elif 1 <= idx <= len(cases): chosen = [cases[idx-1]]
        if not chosen:
            raise GeothermalError("selection.json 中找不到 tomo_selection_case=" + key)
        cases = chosen
    ids: set[str] = set()
    for case in cases:
        for cid in case.get("curve_ids", []): ids.add(str(cid))
        for band in case.get("bands", []):
            if band.get("curve_id") is not None: ids.add(str(band["curve_id"]))
    return ids


def _pair_geometry(root: Path):
    """Return canonical station and pair tables using the same geometry as the solver."""
    from dispersion.prepare import station_table

    c = load_config(root / "config.json")
    stations, pairs, _ = station_table(c)
    return stations.reset_index(drop=True), pairs.reset_index(drop=True)


def export_tomography_paths(project: str | Path, selected_only: Optional[bool] = None, selection_case: str | int | None = None) -> dict[str, Any]:
    """Export current measured/conditional dispersion as straight-ray tomography input.

    The export preserves scientific-status columns.  It is a handoff dataset, not
    a claim that every curve is the unique physical dispersion branch.
    """
    root = Path(project).resolve()
    out = root / "output" / "geothermal"
    out.mkdir(parents=True, exist_ok=True)
    s = load_settings(root)
    if selected_only is None:
        selected_only = bool(s["selected_paths_only"])

    from .candidates import read_curves
    curves, _ = read_curves(root)
    curves = curves.copy()
    if "is_model_prediction" in curves:
        curves = curves[~curves["is_model_prediction"].map(_truth)].copy()
    if curves.empty:
        raise GeothermalError("当前没有可导出的实测/条件频散点")

    selected = _selected_curve_ids(root, selection_case)
    if selected_only:
        if not selected:
            raise GeothermalError("selected_paths_only=true，但 selection.json 中没有 curve_ids")
        curves = curves[curves.curve_id.astype(str).isin(selected)].copy()
        if curves.empty:
            raise GeothermalError("selection.json 中的曲线与当前频散交接结果不匹配")

    stations, pairs = _pair_geometry(root)
    has_geo = {"latitude", "longitude"} <= set(stations.columns)
    if has_geo:
        lat = pd.to_numeric(stations["latitude"], errors="coerce")
        lon = pd.to_numeric(stations["longitude"], errors="coerce")
        if lat.isna().any() or lon.isna().any():
            has_geo = False

    geom_rows = []
    for p in pairs.itertuples(index=False):
        a = stations.iloc[int(p.i)]
        b = stations.iloc[int(p.j)]
        row = {
            "geometry_pair_index": int(p.pair_index),
            "pair_geometry": str(p.pair),
            "station_a": str(a.station),
            "station_b": str(b.station),
            "geometry_distance_m": float(p.distance_m),
            "pair_azimuth_deg": float(p.azimuth_deg),
            "x_a_m": float(a.x_m), "y_a_m": float(a.y_m),
            "x_b_m": float(b.x_m), "y_b_m": float(b.y_m),
        }
        if has_geo:
            row.update(
                longitude_a=float(a.longitude), latitude_a=float(a.latitude),
                longitude_b=float(b.longitude), latitude_b=float(b.latitude),
            )
        geom_rows.append(row)
    geom = pd.DataFrame(geom_rows)

    # Resolve geometry by the explicit pair label first.  Historical/example
    # archives can carry pair_index values from an earlier station ordering, so
    # silently trusting pair_index would attach a curve to the wrong ray.
    exact = {str(r["pair_geometry"]): r for _, r in geom.iterrows()}
    resolved = []
    for row in curves.to_dict("records"):
        pair_name = str(row.get("pair", ""))
        g = exact.get(pair_name)
        reversed_pair = False
        if g is None and "--" in pair_name:
            a_name, b_name = pair_name.split("--", 1)
            g = exact.get(b_name + "--" + a_name)
            reversed_pair = g is not None
        if g is None and row.get("pair_index") is not None:
            hit = geom[geom.geometry_pair_index == int(row["pair_index"])]
            if len(hit) == 1:
                g = hit.iloc[0]
        if g is None:
            raise GeothermalError("频散台对无法与当前 stations.csv 几何对应: " + pair_name)
        outrow = dict(row)
        outrow["source_pair_index"] = row.get("pair_index")
        outrow["geometry_pair_index"] = int(g.geometry_pair_index)
        outrow["pair_geometry"] = str(g.pair_geometry)
        if reversed_pair:
            outrow.update(
                station_a=str(g.station_b), station_b=str(g.station_a),
                x_a_m=float(g.x_b_m), y_a_m=float(g.y_b_m),
                x_b_m=float(g.x_a_m), y_b_m=float(g.y_a_m),
                pair_azimuth_deg=float((g.pair_azimuth_deg + 180.0) % 360.0),
            )
            if has_geo:
                outrow.update(longitude_a=float(g.longitude_b), latitude_a=float(g.latitude_b), longitude_b=float(g.longitude_a), latitude_b=float(g.latitude_a))
        else:
            for key in ("station_a","station_b","x_a_m","y_a_m","x_b_m","y_b_m","pair_azimuth_deg"):
                outrow[key] = g[key]
            if has_geo:
                for key in ("longitude_a","latitude_a","longitude_b","latitude_b"):
                    outrow[key] = g[key]
        geomdist = float(g.geometry_distance_m)
        if "distance_m" in row and pd.notna(row["distance_m"]):
            measured = float(row["distance_m"])
            if not np.isclose(measured, geomdist, rtol=1e-3, atol=.05):
                raise GeothermalError("频散表的 distance_m 与台对标签对应的 stations.csv 几何不一致: " + pair_name)
            outrow["distance_m"] = measured
        else:
            outrow["distance_m"] = geomdist
        resolved.append(outrow)
    merged = pd.DataFrame(resolved)

    merged["period_s"] = 1.0 / pd.to_numeric(merged.frequency_hz, errors="coerce")
    merged["coordinate_system"] = "WGS84" if has_geo else "local_project_xy_m"
    if "unique_physical_dispersion_verified" not in merged:
        merged["unique_physical_dispersion_verified"] = False
    if "direction_quality_verified" not in merged:
        merged["direction_quality_verified"] = False
    if "measurement_method" not in merged:
        merged["measurement_method"] = "conditional_ccf_cycle_hypothesis"
    if "velocity_type" not in merged:
        merged["velocity_type"] = "conditional_pair_velocity"
    if "background_discriminated" not in merged:
        merged["background_discriminated"] = False
    if "source_phase_verified" not in merged:
        merged["source_phase_verified"] = False
    if "is_low_frequency_reconstruction" not in merged:
        merged["is_low_frequency_reconstruction"] = False

    merged["tomography_handoff_status"] = np.where(
        merged["measurement_method"].astype(str).eq("direction_projected_goertzel")
        & merged["direction_quality_verified"].map(_truth)
        & merged["background_discriminated"].map(_truth),
        "measured_direction_projected_conditional",
        "conditional_cycle_hypothesis",
    )
    merged["default_quality_pass"] = (
        pd.to_numeric(merged.phase_velocity_m_s, errors="coerce").gt(0)
        & merged["background_discriminated"].map(_truth)
    )

    preferred = [
        "curve_id", "segment_id", "cohort", "pair_index", "source_pair_index", "geometry_pair_index", "pair", "branch",
        "frequency_hz", "period_s", "phase_velocity_m_s", "group_velocity_m_s",
        "station_a", "station_b", "longitude_a", "latitude_a", "longitude_b", "latitude_b",
        "x_a_m", "y_a_m", "x_b_m", "y_b_m", "distance_m", "pair_azimuth_deg",
        "projected_distance_m", "confirmed_azimuth_deg", "velocity_type", "measurement_method",
        "direction_quality_verified", "background_discriminated", "source_phase_verified",
        "unique_physical_dispersion_verified", "is_low_frequency_reconstruction",
        "tomography_handoff_status", "default_quality_pass", "coordinate_system",
    ]
    columns = [c for c in preferred if c in merged.columns] + [c for c in merged.columns if c not in preferred]
    merged = merged[columns].sort_values(["frequency_hz", "pair_index", "curve_id"]).reset_index(drop=True)
    csv_path = out / "tomography_paths.csv"
    merged.to_csv(csv_path, index=False)

    manifest = {
        "complete": True,
        "source": str(getattr(curves, "attrs", {}).get("source_file", "current Studio dispersion handoff")),
        "selected_only": bool(selected_only),
        "selection_case": None if selection_case in (None, "") else str(selection_case),
        "curve_count": int(merged.curve_id.nunique()),
        "pair_count": int(merged.pair_index.nunique()),
        "frequency_count": int(merged.frequency_hz.nunique()),
        "rows": int(len(merged)),
        "coordinate_system": "WGS84" if has_geo else "local_project_xy_m",
        "output": str(csv_path),
        "scientific_guards": [
            "Exported conditional curves are not automatically declared the unique physical dispersion branch.",
            "Direction-projected Goertzel measurements retain their direction/background quality flags.",
            "A tomography solver must propagate branch, cycle, uncertainty, and path-coverage assumptions.",
        ],
    }
    write_json(out / "tomography_paths_manifest.json", manifest)
    print(f"【层析交接】{manifest['pair_count']}条台站路径，{manifest['curve_count']}条曲线，{manifest['rows']}个频点 -> {csv_path}", flush=True)
    if not has_geo:
        print("【层析交接提示】当前项目只有本地x/y，没有经纬度；可做局部坐标层析，但不能直接叠加区域地质图。", flush=True)
    return manifest


_LON = ("longitude", "lon")
_LAT = ("latitude", "lat")
_DEPTH_M = ("depth_m", "z_m")
_DEPTH_KM = ("depth_km", "z_km")
_VS_M = ("vs_m_s", "vs_mps", "shear_velocity_m_s")
_VS_KM = ("vs_km_s", "shear_velocity_km_s")
_STD_M = ("vs_std_m_s", "vs_sigma_m_s", "sigma_vs_m_s")
_STD_KM = ("vs_std_km_s", "vs_sigma_km_s")
_RES = ("resolution", "resolution_value", "checkerboard_resolution")
_RAYS = ("ray_count", "path_count", "rays")


def _find_col(df: pd.DataFrame, names: tuple[str, ...]) -> Optional[str]:
    low = {str(c).strip().lower(): c for c in df.columns}
    return next((low[n] for n in names if n in low), None)


def _normalise_spatial_frame(df: pd.DataFrame, vs_unit: str = "auto", grid_depth_m: Optional[float] = None) -> pd.DataFrame:
    lon = _find_col(df, _LON); lat = _find_col(df, _LAT)
    dm = _find_col(df, _DEPTH_M); dk = _find_col(df, _DEPTH_KM)
    vm = _find_col(df, _VS_M); vk = _find_col(df, _VS_KM)
    generic_vs = _find_col(df, ("vs",))
    if lon is None or lat is None:
        raise GeothermalError("空间Vs表必须包含 longitude/lon 与 latitude/lat")
    if dm is None and dk is None and grid_depth_m is None:
        raise GeothermalError("空间Vs表必须包含 depth_m/depth_km；二维网格需在设置中明确 grid_depth_m")
    if vm is None and vk is None and generic_vs is None:
        raise GeothermalError("空间Vs表必须包含 vs_m_s 或 vs_km_s")
    if generic_vs is not None and vm is None and vk is None and vs_unit == "auto":
        raise GeothermalError("列名只有 'vs' 时单位不明确；请设置 tomography_vs_unit=m/s 或 km/s")

    out = pd.DataFrame()
    out["longitude"] = pd.to_numeric(df[lon], errors="coerce")
    out["latitude"] = pd.to_numeric(df[lat], errors="coerce")
    if dm is not None:
        out["depth_m"] = pd.to_numeric(df[dm], errors="coerce")
    elif dk is not None:
        out["depth_m"] = 1000.0 * pd.to_numeric(df[dk], errors="coerce")
    else:
        out["depth_m"] = float(grid_depth_m)
    if vm is not None:
        out["vs_m_s"] = pd.to_numeric(df[vm], errors="coerce")
    elif vk is not None:
        out["vs_m_s"] = 1000.0 * pd.to_numeric(df[vk], errors="coerce")
    else:
        factor = 1000.0 if vs_unit == "km/s" else 1.0
        out["vs_m_s"] = factor * pd.to_numeric(df[generic_vs], errors="coerce")

    sm = _find_col(df, _STD_M); sk = _find_col(df, _STD_KM)
    if sm is not None: out["vs_std_m_s"] = pd.to_numeric(df[sm], errors="coerce")
    elif sk is not None: out["vs_std_m_s"] = 1000.0 * pd.to_numeric(df[sk], errors="coerce")
    rc = _find_col(df, _RES)
    if rc is not None: out["resolution"] = pd.to_numeric(df[rc], errors="coerce")
    rr = _find_col(df, _RAYS)
    if rr is not None: out["ray_count"] = pd.to_numeric(df[rr], errors="coerce")

    out = out.replace([np.inf, -np.inf], np.nan).dropna(subset=["longitude", "latitude", "depth_m", "vs_m_s"])
    out = out[out.longitude.between(-360, 360) & out.latitude.between(-90, 90) & (out.depth_m >= 0) & (out.vs_m_s > 0)]
    if out.empty:
        raise GeothermalError("空间Vs没有有效经纬度/深度/速度点")
    if out.duplicated(["longitude", "latitude", "depth_m"]).any():
        raise GeothermalError("空间Vs存在重复 longitude/latitude/depth_m 网格点")
    return out.sort_values(["depth_m", "latitude", "longitude"]).reset_index(drop=True)


def _read_netcdf_spatial(path: Path, variable: str, vs_unit: str, grid_depth_m: Optional[float]) -> pd.DataFrame:
    try:
        import xarray as xr
    except Exception as exc:
        raise GeothermalError("读取 NetCDF/GRD 需要地热可选依赖 xarray；可运行 install_geothermal_windows.bat") from exc
    ds = xr.open_dataset(path)
    try:
        if variable in ds.data_vars:
            da = ds[variable]
        elif len(ds.data_vars) == 1:
            da = next(iter(ds.data_vars.values()))
        else:
            aliases = [n for n in ds.data_vars if str(n).lower() in ("vs", "vs_m_s", "vs_km_s", "z")]
            if len(aliases) != 1:
                raise GeothermalError("NetCDF含多个变量，请在 tomography_variable 指定Vs变量")
            da = ds[aliases[0]]
        frame = da.to_dataframe(name=str(da.name or variable)).reset_index()
        # Rename chosen variable to a unit-explicit name when attrs are trustworthy.
        units = str(da.attrs.get("units", "")).lower().replace(" ", "")
        name = str(da.name or variable)
        if name not in frame:
            raise GeothermalError("无法从NetCDF展开Vs变量")
        if units in ("m/s", "m*s-1", "ms-1"):
            frame = frame.rename(columns={name: "vs_m_s"})
        elif units in ("km/s", "km*s-1", "kms-1"):
            frame = frame.rename(columns={name: "vs_km_s"})
        elif name.lower() not in _VS_M + _VS_KM + ("vs",):
            frame = frame.rename(columns={name: "vs"})
        # Coordinate aliases from common GMT/xarray names.
        ren = {}
        for c in frame.columns:
            low = str(c).lower()
            if low in ("x", "longitude", "lon"): ren[c] = "longitude"
            elif low in ("y", "latitude", "lat"): ren[c] = "latitude"
            elif low in ("z", "depth", "depth_m"): ren[c] = "depth_m"
        frame = frame.rename(columns=ren)
        return _normalise_spatial_frame(frame, vs_unit, grid_depth_m)
    finally:
        ds.close()


def normalize_spatial_vs(path: str | Path, *, variable: str = "vs_m_s", vs_unit: str = "auto", grid_depth_m: Optional[float] = None) -> pd.DataFrame:
    path = Path(path).resolve()
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() in (".csv", ".txt"):
        return _normalise_spatial_frame(pd.read_csv(path), vs_unit, grid_depth_m)
    if path.suffix.lower() in (".nc", ".netcdf", ".grd"):
        return _read_netcdf_spatial(path, variable, vs_unit, grid_depth_m)
    raise GeothermalError("空间Vs当前支持 CSV/TXT 或 NetCDF/GRD")


def apply_quality_mask(frame: pd.DataFrame, settings: dict[str, Any]) -> pd.DataFrame:
    q = frame.copy()
    passed = np.ones(len(q), dtype=bool)
    reasons = np.array([""] * len(q), dtype=object)

    def reject(mask, reason):
        nonlocal passed, reasons
        bad = ~mask
        reasons[bad] = np.where(reasons[bad] == "", reason, reasons[bad] + ";" + reason)
        passed &= mask

    max_rel = settings.get("max_relative_vs_std")
    if max_rel is not None:
        if "vs_std_m_s" not in q:
            reject(np.zeros(len(q), dtype=bool), "missing_vs_std")
        else:
            rel = q.vs_std_m_s / q.vs_m_s
            q["relative_vs_std"] = rel
            reject(np.isfinite(rel) & (rel <= float(max_rel)), "vs_std_too_large")
    min_res = settings.get("min_resolution")
    if min_res is not None:
        if "resolution" not in q:
            reject(np.zeros(len(q), dtype=bool), "missing_resolution")
        else:
            reject(np.isfinite(q.resolution) & (q.resolution >= float(min_res)), "resolution_too_low")
    min_rays = settings.get("min_ray_count")
    if min_rays is not None:
        if "ray_count" not in q:
            reject(np.zeros(len(q), dtype=bool), "missing_ray_count")
        else:
            reject(np.isfinite(q.ray_count) & (q.ray_count >= float(min_rays)), "ray_count_too_low")
    q["quality_pass"] = passed
    q["quality_reasons"] = reasons
    return q


def validate_spatial_vs(project: str | Path) -> dict[str, Any]:
    root = Path(project).resolve(); s = load_settings(root)
    source = resolve_project_path(root, s["tomography_file"])
    if source is None:
        raise GeothermalError("尚未设置 tomography_file")
    frame = normalize_spatial_vs(source, variable=s["tomography_variable"], vs_unit=s["tomography_vs_unit"], grid_depth_m=s["grid_depth_m"])
    q = apply_quality_mask(frame, s)
    out = root / "output" / "geothermal"; out.mkdir(parents=True, exist_ok=True)
    standardized = out / "spatial_vs_standardized.csv"; q.to_csv(standardized, index=False)
    depths = sorted(float(x) for x in q.depth_m.unique())
    result = {
        "complete": True,
        "source": str(source), "source_sha256": sha(source),
        "rows": int(len(q)), "quality_pass_rows": int(q.quality_pass.sum()),
        "depth_levels_m": depths,
        "longitude_range": [float(q.longitude.min()), float(q.longitude.max())],
        "latitude_range": [float(q.latitude.min()), float(q.latitude.max())],
        "vs_range_m_s": [float(q.vs_m_s.min()), float(q.vs_m_s.max())],
        "has_vs_std": "vs_std_m_s" in q, "has_resolution": "resolution" in q, "has_ray_count": "ray_count" in q,
        "output": str(standardized),
    }
    write_json(out / "spatial_vs_validation.json", result)
    print(f"【空间Vs检查】{len(q)}点，{len(depths)}个深度层，质量门限保留{int(q.quality_pass.sum())}点", flush=True)
    return result


def _point_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    lon = _find_col(df, _LON); lat = _find_col(df, _LAT)
    if lon is None or lat is None:
        raise GeothermalError(f"点位文件 {path.name} 需要 lon/lat 或 longitude/latitude")
    out = df.copy(); out["longitude"] = pd.to_numeric(df[lon], errors="coerce"); out["latitude"] = pd.to_numeric(df[lat], errors="coerce")
    return out.dropna(subset=["longitude", "latitude"])


def _read_stations(root: Path) -> pd.DataFrame:
    path = root / "stations.csv"
    if not path.exists(): return pd.DataFrame()
    df = pd.read_csv(path)
    if not {"latitude", "longitude"} <= set(df): return pd.DataFrame()
    out = pd.DataFrame({
        "name": df.get("station", pd.Series([f"STA{i+1}" for i in range(len(df))])).astype(str),
        "longitude": pd.to_numeric(df.longitude, errors="coerce"),
        "latitude": pd.to_numeric(df.latitude, errors="coerce"),
    })
    return out.dropna(subset=["longitude", "latitude"])


def _vector_gdf(path: Optional[Path]):
    if path is None: return None
    try:
        import geopandas as gpd
    except Exception as exc:
        raise GeothermalError("读取 Shapefile/GeoJSON/GeoPackage 需要 geopandas；可运行 install_geothermal_windows.bat") from exc
    gdf = gpd.read_file(path)
    if gdf.crs is None:
        raise GeothermalError(f"{path.name} 没有 CRS，不能安全叠加地质图")
    if gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs(4326)
    return gdf


def _regular_grid(slice_df: pd.DataFrame):
    lon = np.sort(slice_df.longitude.unique()); lat = np.sort(slice_df.latitude.unique())
    if len(slice_df) != len(lon) * len(lat):
        return None
    p = slice_df.pivot(index="latitude", columns="longitude", values="vs_m_s").reindex(index=lat, columns=lon)
    if p.isna().any().any(): return None
    return lon, lat, p.to_numpy(float)


def _region_from_inputs(slice_df: Optional[pd.DataFrame], stations: pd.DataFrame, points: list[pd.DataFrame], vectors: list[Any]) -> list[float]:
    xs=[]; ys=[]
    if slice_df is not None and len(slice_df): xs += slice_df.longitude.tolist(); ys += slice_df.latitude.tolist()
    if len(stations): xs += stations.longitude.tolist(); ys += stations.latitude.tolist()
    for p in points:
        if p is not None and len(p): xs += p.longitude.tolist(); ys += p.latitude.tolist()
    for g in vectors:
        if g is not None and len(g):
            b = g.total_bounds; xs += [b[0], b[2]]; ys += [b[1], b[3]]
    if not xs or not ys: raise GeothermalError("没有可确定地图范围的空间数据")
    w,e,s,n=float(min(xs)),float(max(xs)),float(min(ys)),float(max(ys)); dx=max(e-w,1e-6);dy=max(n-s,1e-6)
    return [w-.06*dx,e+.06*dx,s-.06*dy,n+.06*dy]


def _plot_matplotlib(root: Path, slice_df: Optional[pd.DataFrame], actual_depth: Optional[float], settings: dict[str,Any], basename: str, warnings: list[str]) -> dict[str,str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    geology_p=resolve_project_path(root,settings["geology"]); faults_p=resolve_project_path(root,settings["faults"])
    wells_p=resolve_project_path(root,settings["wells"]); springs_p=resolve_project_path(root,settings["hot_springs"])
    geology=_vector_gdf(geology_p) if geology_p else None; faults=_vector_gdf(faults_p) if faults_p else None
    wells=_point_table(wells_p) if wells_p else pd.DataFrame(); springs=_point_table(springs_p) if springs_p else pd.DataFrame(); stations=_read_stations(root)
    region=_region_from_inputs(slice_df,stations,[wells,springs],[geology,faults])

    fig,ax=plt.subplots(figsize=(9,7))
    mappable=None
    if slice_df is not None and len(slice_df):
        regular=_regular_grid(slice_df)
        if regular is not None:
            lon,lat,z=regular; mappable=ax.pcolormesh(lon,lat,z/1000.0,shading="auto",cmap=settings["cmap"])
        else:
            mappable=ax.scatter(slice_df.longitude,slice_df.latitude,c=slice_df.vs_m_s/1000.0,cmap=settings["cmap"],s=28,marker="s")
            warnings.append("空间Vs不是完整规则网格；Matplotlib按原始点显示，不做无约束空间插值。")
        fig.colorbar(mappable,ax=ax,label="Vs (km/s)")
    if geology is not None and len(geology): geology.boundary.plot(ax=ax,linewidth=.7,color="0.25")
    if faults is not None and len(faults): faults.plot(ax=ax,linewidth=1.2,color="crimson")
    if len(stations):
        ax.scatter(stations.longitude,stations.latitude,marker="v",s=34,facecolors="white",edgecolors="black",label="Seismic stations",zorder=5)
        if settings["label_points"]:
            for r in stations.itertuples(): ax.annotate(str(r.name),(r.longitude,r.latitude),xytext=(3,3),textcoords="offset points",fontsize=7)
    if len(wells): ax.scatter(wells.longitude,wells.latitude,marker="^",s=38,facecolors="none",edgecolors="black",label="Wells",zorder=6)
    if len(springs): ax.scatter(springs.longitude,springs.latitude,marker="o",s=32,facecolors="none",edgecolors="darkred",label="Hot springs",zorder=6)
    ax.set_xlim(region[:2]);ax.set_ylim(region[2:]);ax.set_xlabel("Longitude");ax.set_ylabel("Latitude");ax.grid(alpha=.18)
    title=settings.get("title") or "Geothermal interpretation"
    if actual_depth is not None: title += f" · Vs depth {actual_depth/1000:g} km"
    ax.set_title(title)
    if len(stations) or len(wells) or len(springs): ax.legend(fontsize=8,loc="best")
    fig.tight_layout()
    out=root/"output"/"geothermal"/"figures";out.mkdir(parents=True,exist_ok=True)
    png=out/(basename+".png");pdf=out/(basename+".pdf");fig.savefig(png,dpi=250);fig.savefig(pdf);plt.close(fig)
    return {"png":str(png),"pdf":str(pdf),"backend":"matplotlib"}


def _plot_pygmt(root: Path, slice_df: Optional[pd.DataFrame], actual_depth: Optional[float], settings: dict[str,Any], basename: str, warnings: list[str]) -> dict[str,str]:
    try:
        import pygmt
        import xarray as xr
    except Exception as exc:
        raise GeothermalError("PyGMT后端需要 pygmt + GMT + xarray；请运行 install_geothermal_windows.bat 并安装GMT") from exc
    geology_p=resolve_project_path(root,settings["geology"]);faults_p=resolve_project_path(root,settings["faults"])
    wells_p=resolve_project_path(root,settings["wells"]);springs_p=resolve_project_path(root,settings["hot_springs"])
    geology=_vector_gdf(geology_p) if geology_p else None;faults=_vector_gdf(faults_p) if faults_p else None
    wells=_point_table(wells_p) if wells_p else pd.DataFrame();springs=_point_table(springs_p) if springs_p else pd.DataFrame();stations=_read_stations(root)
    region=_region_from_inputs(slice_df,stations,[wells,springs],[geology,faults])
    fig=pygmt.Figure();title=settings.get("title") or "Geothermal interpretation"
    if actual_depth is not None:title+=f" · Vs depth {actual_depth/1000:g} km"
    fig.basemap(region=region,projection="M16c",frame=["af",f'+t"{title}"'])
    if slice_df is not None and len(slice_df):
        regular=_regular_grid(slice_df)
        if regular is None:
            raise GeothermalError("PyGMT后端当前要求Vs切片为完整规则经纬度网格；不规则点请选 matplotlib，程序不会自动插值。")
        lon,lat,z=regular;da=xr.DataArray(z/1000.0,coords={"latitude":lat,"longitude":lon},dims=("latitude","longitude"),name="Vs")
        fig.grdimage(grid=da,cmap=settings["cmap"],region=region,projection="M16c");fig.colorbar(frame='af+l"Vs (km/s)"')
    if geology is not None and len(geology):fig.plot(data=geology,pen="0.5p,gray30")
    if faults is not None and len(faults):fig.plot(data=faults,pen="1.1p,red")
    if len(stations):fig.plot(x=stations.longitude,y=stations.latitude,style="i0.22c",fill="white",pen="0.5p,black",label="Seismic stations")
    if len(wells):fig.plot(x=wells.longitude,y=wells.latitude,style="t0.24c",fill="white",pen="0.5p,black",label="Wells")
    if len(springs):fig.plot(x=springs.longitude,y=springs.latitude,style="c0.20c",fill="red",pen="0.4p,black",label="Hot springs")
    if len(stations) or len(wells) or len(springs):fig.legend(position="JTR+jTR+o0.2c",box="+gwhite+p0.5p")
    out=root/"output"/"geothermal"/"figures";out.mkdir(parents=True,exist_ok=True);png=out/(basename+".png");pdf=out/(basename+".pdf")
    fig.savefig(png,dpi=300);fig.savefig(pdf)
    return {"png":str(png),"pdf":str(pdf),"backend":"pygmt"}


def _choose_backend(settings: dict[str,Any], slice_df: Optional[pd.DataFrame]) -> str:
    b=settings["backend"]
    if b!="auto":return b
    # Prefer PyGMT only when its runtime is actually usable and the Vs grid is regular.
    try:
        import pygmt  # noqa:F401
        ok = slice_df is None or _regular_grid(slice_df) is not None
        if ok:return "pygmt"
    except Exception:pass
    return "matplotlib"


def _pick_depth(frame: pd.DataFrame, requested: float, tolerance: Optional[float]) -> tuple[float,pd.DataFrame]:
    levels=np.sort(frame.depth_m.unique().astype(float));actual=float(levels[np.argmin(abs(levels-requested))])
    if tolerance is None:
        if len(levels)>1:
            dif=np.diff(levels); tolerance=max(1e-6,float(np.min(dif))/3.0)
        else:tolerance=1.0
    if abs(actual-requested)>float(tolerance)+1e-9:
        raise GeothermalError(f"请求深度{requested:g} m附近没有空间Vs层；最近{actual:g} m，容差{float(tolerance):g} m")
    return actual,frame[np.isclose(frame.depth_m,actual,rtol=0,atol=1e-8)].copy()


def _base_map_possible(root: Path, s: dict[str,Any]) -> bool:
    if len(_read_stations(root)):return True
    return any(resolve_project_path(root,s[k]) is not None for k in ("geology","faults","wells","hot_springs"))


def generate_geothermal(project: str | Path) -> dict[str,Any]:
    root=Path(project).resolve();s=load_settings(root);out=root/"output"/"geothermal";figdir=out/"figures";figdir.mkdir(parents=True,exist_ok=True)
    warnings=[];products=[]
    print("【地热成果】开始：Studio结果 → 层析交接/空间Vs → 地质融合",flush=True)

    # Export current pair dispersion when possible. This gives the next tomography
    # stage a stable handoff even if no spatial inversion has been supplied yet.
    try:
        path_manifest=export_tomography_paths(root)
    except Exception as exc:
        path_manifest=None;warnings.append("层析射线路径未导出: "+str(exc));print("【地热成果提示】"+warnings[-1],flush=True)

    spatial_validation=None;frame=None
    source=resolve_project_path(root,s["tomography_file"])
    if source is not None:
        spatial_validation=validate_spatial_vs(root)
        frame=pd.read_csv(spatial_validation["output"])
        frame=frame[frame.quality_pass.map(_truth)].copy()
        if frame.empty:raise GeothermalError("不确定性/分辨率门限后没有剩余空间Vs点")
        for req in s["slice_depths_m"]:
            try:
                actual,sl=_pick_depth(frame,float(req),s["depth_tolerance_m"])
            except GeothermalError as exc:
                warnings.append(str(exc));print("【地热切片跳过】"+str(exc),flush=True);continue
            tag=f"vs_{int(round(actual)):06d}m";csv_path=out/(tag+".csv");sl.to_csv(csv_path,index=False)
            backend=_choose_backend(s,sl)
            plot=_plot_pygmt(root,sl,actual,s,tag,warnings) if backend=="pygmt" else _plot_matplotlib(root,sl,actual,s,tag,warnings)
            products.append({"requested_depth_m":float(req),"actual_depth_m":actual,"csv":str(csv_path),**plot})
            print(f"【地热切片】{actual:g} m · {len(sl)}点 · {plot['backend']}",flush=True)
    elif _base_map_possible(root,s):
        backend=_choose_backend(s,None)
        plot=_plot_pygmt(root,None,None,s,"geological_context",warnings) if backend=="pygmt" else _plot_matplotlib(root,None,None,s,"geological_context",warnings)
        products.append({"requested_depth_m":None,"actual_depth_m":None,**plot})
        warnings.append("未提供空间Vs层析；当前只输出地质/断裂/台站/钻井/温泉背景图，不把1-D结构平铺为2-D速度场。")
    else:
        warnings.append("没有空间Vs或可绘制的经纬度图层；仅保留层析交接表。")

    # Surface geothermal figures in the existing Studio figures page.
    studio_fig=root/"output"/"figures";studio_fig.mkdir(parents=True,exist_ok=True)
    for i,p in enumerate(products[:6],start=1):
        if p.get("png"):
            dest=studio_fig/f"20_geothermal_{i:02d}_{Path(p['png']).name}"
            shutil.copy2(p["png"],dest)

    try:
        from .tomography import tomography_status
        phase_tomography=tomography_status(root)
    except Exception:
        phase_tomography={"complete":False}
    result={
        "complete":True,"settings_sha256":sha(settings_path(root)) if settings_path(root).exists() else None,
        "tomography_paths":path_manifest,"phase_velocity_tomography":phase_tomography,"spatial_vs":spatial_validation,"products":products,
        "warnings":warnings,
        "scientific_guards":[
            "A common 1-D conditional structure is never expanded into a fake spatial Vs volume.",
            "Low-Vs is not interpreted as geothermal resource by itself; geology, faults, fluids, temperature and independent geophysics are required.",
            "Depth slices inherit the uncertainty/resolution of the supplied spatial inversion; plotting depth is not proof of depth resolution.",
        ],
    }
    write_json(out/"geothermal_complete.json",result);_write_report(root,result,s)
    print("【地热成果】完成 -> "+str(out/"report.html"),flush=True)
    return result


def _write_report(root: Path, result: dict[str,Any], settings: dict[str,Any]) -> None:
    out=root/"output"/"geothermal";parts=['<meta charset="utf-8"><style>body{font:16px/1.7 system-ui,"Microsoft YaHei",sans-serif;max-width:1180px;margin:30px auto;padding:0 22px;color:#203548}h1,h2{color:#174963}.warn{background:#fff4e5;border-left:5px solid #d97706;padding:12px}.guard{background:#eef6f8;border-left:5px solid #24748a;padding:12px}img{max-width:100%}table{border-collapse:collapse}td,th{border:1px solid #ccd7df;padding:6px}</style>']
    parts.append('<h1>地热地质融合与空间层析成果</h1>')
    parts.append('<div class="guard"><b>科学边界：</b><ul>'+''.join('<li>'+html.escape(x)+'</li>' for x in result['scientific_guards'])+'</ul></div>')
    if result['warnings']:parts.append('<div class="warn"><b>提示：</b><ul>'+''.join('<li>'+html.escape(x)+'</li>' for x in result['warnings'])+'</ul></div>')
    p=result.get('tomography_paths')
    if p:parts.append(f"<h2>层析射线交接</h2><p>{p['pair_count']}条台站路径，{p['curve_count']}条曲线，{p['rows']}个频点。<a href=\"tomography_paths.csv\">tomography_paths.csv</a></p>")
    v=result.get('spatial_vs')
    if v:parts.append(f"<h2>空间Vs输入</h2><p>{v['rows']}点，{len(v['depth_levels_m'])}个深度层；质量门限保留{v['quality_pass_rows']}点。<a href=\"spatial_vs_standardized.csv\">规范化空间Vs表</a></p>")
    if result['products']:
        parts.append('<h2>成果图</h2>')
        for p in result['products']:
            if p.get('png'):
                rel=Path(p['png']).relative_to(out).as_posix();parts.append(f'<figure><img src="{html.escape(rel)}"><figcaption>深度 {p.get("actual_depth_m")} m · {html.escape(p.get("backend",""))}</figcaption></figure>')
    parts.append('<h2>下一阶段接口</h2><p>tomography_paths.csv 是从当前频散结果到二维相速度层析的稳定交接表；空间Vs规范表采用 longitude, latitude, depth_m, vs_m_s，并可附 vs_std_m_s、resolution、ray_count。后续可在此协议上加入直线路径/弯曲路径层析、局部频散反演、MT/重力/热流与钻井测温综合评价。</p>')
    (out/'report.html').write_text(''.join(parts),encoding='utf-8')


def geothermal_status(project: str | Path) -> dict[str,Any]:
    root=Path(project).resolve();path=root/'output/geothermal/geothermal_complete.json'
    if not path.exists():return {"complete":False,"message":"尚未生成地热成果"}
    return json.loads(path.read_text(encoding='utf-8-sig'))


def create_example(directory: str | Path) -> Path:
    """Create a small self-contained geographic Vs example for the geothermal tab."""
    from .project import new_project
    root=new_project(directory);rng=np.random.default_rng(9075)
    stations=pd.DataFrame([
        ['S01','XX.S01..HHZ',50.00,10.00],['S02','XX.S02..HHZ',50.01,10.08],['S03','XX.S03..HHZ',50.07,10.02],
        ['S04','XX.S04..HHZ',50.08,10.10],['S05','XX.S05..HHZ',50.04,10.14],['S06','XX.S06..HHZ',49.97,10.12]],columns=['station','seed_id','latitude','longitude'])
    stations['file']='unused.npz';stations['precorrected']=True;stations['units']='m/s';stations['polarity']=1;stations.to_csv(root/'stations.csv',index=False)
    lons=np.linspace(9.96,10.16,31);lats=np.linspace(49.95,50.10,25);depths=[500.,1000.,2000.,4000.];rows=[]
    for d in depths:
        for lat in lats:
            for lon in lons:
                anomaly=np.exp(-((lon-10.075)**2/(2*.035**2)+(lat-50.025)**2/(2*.025**2)))
                vs=2600+0.28*d-420*anomaly*np.exp(-d/4500)
                rows.append([lon,lat,d,vs,90+15*rng.random(),.65+.3*rng.random(),4+int(10*rng.random())])
    pd.DataFrame(rows,columns=['longitude','latitude','depth_m','vs_m_s','vs_std_m_s','resolution','ray_count']).to_csv(root/'spatial_vs_demo.csv',index=False)
    pd.DataFrame([['Well-A',10.07,50.02],['Well-B',10.12,50.06]],columns=['name','longitude','latitude']).to_csv(root/'wells_demo.csv',index=False)
    pd.DataFrame([['Spring-A',10.06,50.03]],columns=['name','longitude','latitude']).to_csv(root/'hot_springs_demo.csv',index=False)
    save_settings(root,dict(DEFAULT_SETTINGS,tomography_file='spatial_vs_demo.csv',slice_depths_m=depths,wells='wells_demo.csv',hot_springs='hot_springs_demo.csv',backend='matplotlib',title='Geothermal spatial Vs demo'))
    return root
