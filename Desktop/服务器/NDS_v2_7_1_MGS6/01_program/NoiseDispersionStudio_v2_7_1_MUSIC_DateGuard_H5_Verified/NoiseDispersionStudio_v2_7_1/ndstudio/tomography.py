"""Straight-ray 2-D phase-velocity tomography for Studio geothermal handoff.

The solver is deliberately conservative:
- one observation per station-pair/frequency is required by default;
- direction-projected plane-wave measurements are converted to the equivalent
  pair-distance apparent velocity before straight-ray inversion;
- strongly oblique plane-wave observations are excluded by a configurable
  pair/wave-vector angle threshold;
- regularization, ray coverage, resolution diagonal, residuals and checkerboard
  recovery are exported with every run.

This module produces phase-velocity maps and local dispersion curves.  It does
not itself convert local dispersion into Vs(z); that remains the next inversion
stage.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional
import json, math, html

import numpy as np
import pandas as pd
from scipy import linalg

from dispersion.common import write_json, sha


class TomographyError(RuntimeError):
    pass


@dataclass(frozen=True)
class Grid2D:
    x_edges: np.ndarray
    y_edges: np.ndarray
    x_centers: np.ndarray
    y_centers: np.ndarray
    nx: int
    ny: int

    @property
    def ncell(self) -> int:
        return self.nx * self.ny


def _truth(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float)) and value in (0, 1):
        return bool(value)
    return str(value).strip().lower() in ("true", "1", "yes", "y")


def _angle_axis_deg(a: float, b: float) -> float:
    """Acute angle between two axes in degrees (0..90)."""
    d = abs((float(a) - float(b) + 180.0) % 360.0 - 180.0)
    return float(min(d, 180.0 - d))


def _choose_coordinate_columns(df: pd.DataFrame) -> str:
    geo = {"longitude_a", "latitude_a", "longitude_b", "latitude_b"} <= set(df.columns)
    local = {"x_a_m", "y_a_m", "x_b_m", "y_b_m"} <= set(df.columns)
    if geo:
        return "WGS84"
    if local:
        return "local_project_xy_m"
    raise TomographyError("层析路径缺少两端经纬度，也缺少 x/y 本地坐标")


def _project_paths(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any], Any]:
    """Add x/y endpoints and return an inverse transformer when geographic."""
    mode = _choose_coordinate_columns(df)
    q = df.copy()
    if mode == "local_project_xy_m":
        for c in ("x_a_m", "y_a_m", "x_b_m", "y_b_m"):
            q[c] = pd.to_numeric(q[c], errors="coerce")
        meta = {"coordinate_system": mode, "projection": "input local Cartesian metres"}
        return q, meta, None

    try:
        from pyproj import CRS, Transformer
    except Exception as exc:
        raise TomographyError("经纬度二维层析需要 pyproj（Studio 核心依赖）") from exc

    for c in ("longitude_a", "latitude_a", "longitude_b", "latitude_b"):
        q[c] = pd.to_numeric(q[c], errors="coerce")
    lon0 = float(pd.concat([q.longitude_a, q.longitude_b]).median())
    lat0 = float(pd.concat([q.latitude_a, q.latitude_b]).median())
    crs_local = CRS.from_proj4(
        f"+proj=aeqd +lat_0={lat0:.12f} +lon_0={lon0:.12f} +datum=WGS84 +units=m +no_defs"
    )
    forward = Transformer.from_crs("EPSG:4326", crs_local, always_xy=True)
    inverse = Transformer.from_crs(crs_local, "EPSG:4326", always_xy=True)
    xa, ya = forward.transform(q.longitude_a.to_numpy(), q.latitude_a.to_numpy())
    xb, yb = forward.transform(q.longitude_b.to_numpy(), q.latitude_b.to_numpy())
    q["x_a_m"] = xa; q["y_a_m"] = ya; q["x_b_m"] = xb; q["y_b_m"] = yb
    meta = {
        "coordinate_system": "WGS84",
        "local_projection": crs_local.to_string(),
        "projection_center": {"longitude": lon0, "latitude": lat0},
    }
    return q, meta, inverse


def _prepare_observations(paths: pd.DataFrame, settings: dict[str, Any]) -> tuple[pd.DataFrame, list[str]]:
    q = paths.copy()
    warnings: list[str] = []
    if "frequency_hz" not in q or "phase_velocity_m_s" not in q:
        raise TomographyError("tomography_paths.csv 至少需要 frequency_hz 与 phase_velocity_m_s")
    q["frequency_hz"] = pd.to_numeric(q.frequency_hz, errors="coerce")
    q["phase_velocity_m_s"] = pd.to_numeric(q.phase_velocity_m_s, errors="coerce")
    q["distance_m"] = pd.to_numeric(q.get("distance_m"), errors="coerce")
    q = q.replace([np.inf, -np.inf], np.nan).dropna(subset=["frequency_hz", "phase_velocity_m_s", "distance_m"])
    q = q[(q.frequency_hz > 0) & (q.phase_velocity_m_s > 0) & (q.distance_m > 0)].copy()
    if q.empty:
        raise TomographyError("没有有限、正值的频率/相速度/台距观测")

    mode = str(settings.get("tomo_quality_mode", "default_pass"))
    if mode not in ("all_valid", "default_pass", "strict_direction"):
        raise TomographyError("tomo_quality_mode 必须为 all_valid / default_pass / strict_direction")
    if mode in ("default_pass", "strict_direction") and "default_quality_pass" in q:
        before = len(q); q = q[q.default_quality_pass.map(_truth)].copy()
        if len(q) < before: warnings.append(f"质量标志剔除 {before-len(q)} 个观测点")
    if mode == "strict_direction":
        if "direction_quality_verified" not in q:
            raise TomographyError("strict_direction 要求 tomography_paths 包含 direction_quality_verified")
        before = len(q); q = q[q.direction_quality_verified.map(_truth)].copy()
        if len(q) < before: warnings.append(f"方向质量剔除 {before-len(q)} 个观测点")

    # Build a quantity mathematically consistent with a station-pair straight ray.
    q["tomography_velocity_m_s"] = q.phase_velocity_m_s.astype(float)
    q["tomography_observable_kind"] = "pair_path_phase_velocity"
    q["tomography_direction_angle_deg"] = np.nan
    max_angle = float(settings.get("tomo_max_direction_angle_deg", 15.0))
    projected = pd.Series(False, index=q.index)
    if "velocity_type" in q:
        projected |= q.velocity_type.astype(str).eq("wave_vector_phase_velocity")
    if "measurement_method" in q:
        projected |= q.measurement_method.astype(str).eq("direction_projected_goertzel")

    for idx in q.index[projected]:
        row = q.loc[idx]
        angle = np.nan
        if pd.notna(row.get("confirmed_azimuth_deg")) and pd.notna(row.get("pair_azimuth_deg")):
            angle = _angle_axis_deg(float(row.confirmed_azimuth_deg), float(row.pair_azimuth_deg))
            q.at[idx, "tomography_direction_angle_deg"] = angle
        if np.isfinite(angle) and angle > max_angle:
            q.at[idx, "tomography_velocity_m_s"] = np.nan
            q.at[idx, "tomography_observable_kind"] = "rejected_oblique_plane_wave"
            continue
        app = row.get("apparent_pair_velocity_m_s")
        if app is not None and pd.notna(app) and float(app) > 0:
            q.at[idx, "tomography_velocity_m_s"] = float(app)
            q.at[idx, "tomography_observable_kind"] = "apparent_pair_velocity_from_projected_phase"
        elif row.get("projected_distance_m") is not None and pd.notna(row.get("projected_distance_m")):
            proj = abs(float(row.get("projected_distance_m")))
            if proj <= 0:
                q.at[idx, "tomography_velocity_m_s"] = np.nan
                q.at[idx, "tomography_observable_kind"] = "rejected_zero_projection"
            else:
                # phase delay = projected_distance / true wave-vector velocity;
                # equivalent pair velocity = pair_distance / delay.
                q.at[idx, "tomography_velocity_m_s"] = float(row.distance_m) * float(row.phase_velocity_m_s) / proj
                q.at[idx, "tomography_observable_kind"] = "apparent_pair_equivalent_from_projected_phase"
        else:
            # If angle is almost parallel, true and pair-equivalent values are close.
            if not np.isfinite(angle) or angle > 3.0:
                q.at[idx, "tomography_velocity_m_s"] = np.nan
                q.at[idx, "tomography_observable_kind"] = "rejected_missing_projection_metadata"

    rejected = int(q.tomography_velocity_m_s.isna().sum())
    if rejected:
        warnings.append(f"方向投影/夹角规则剔除 {rejected} 个不能安全作为直线路径数据的观测点")
    q = q.dropna(subset=["tomography_velocity_m_s"]).copy()
    q = q[q.tomography_velocity_m_s > 0]
    if q.empty:
        raise TomographyError("质量与方向投影规则后没有可用于二维直线路径层析的观测")

    # Keep frequency bins deterministic and avoid float serialization noise.
    q["frequency_bin_hz"] = q.frequency_hz.round(8)
    duplicate_keys = ["frequency_bin_hz"]
    if "pair_geometry" in q: duplicate_keys.append("pair_geometry")
    elif "pair" in q: duplicate_keys.append("pair")
    else: duplicate_keys += ["station_a", "station_b"]
    dup = q.duplicated(duplicate_keys, keep=False)
    if dup.any():
        policy = str(settings.get("tomo_duplicate_policy", "error"))
        sample = q.loc[dup, duplicate_keys + ["tomography_velocity_m_s"]].head(12)
        if policy == "error":
            raise TomographyError(
                "同一台对/频率存在多个候选观测，二维层析默认拒绝把竞争分支重复计权。"
                "请先建立包含每个台对唯一曲线的联合 selection case，并在 tomo_selection_case 指定它；或明确将 tomo_duplicate_policy 设为 median。\n"
                + sample.to_string(index=False)
            )
        if policy != "median":
            raise TomographyError("tomo_duplicate_policy 必须为 error 或 median")
        group_cols = duplicate_keys
        numeric = q.select_dtypes(include=[np.number]).columns.tolist()
        first_cols = [c for c in q.columns if c not in numeric and c not in group_cols]
        med = q.groupby(group_cols, as_index=False)[numeric].median(numeric_only=True)
        first = q.groupby(group_cols, as_index=False)[first_cols].first() if first_cols else q[group_cols].drop_duplicates()
        q = first.merge(med, on=group_cols, how="inner")
        warnings.append("竞争分支按明确的 median 策略合并；该策略是用户选择，不代表证明了唯一物理分支。")
    return q.reset_index(drop=True), warnings


def _make_grid(df: pd.DataFrame, nx: int, ny: int, margin_fraction: float) -> Grid2D:
    if nx < 2 or ny < 2:
        raise TomographyError("tomo_nx/tomo_ny 至少为 2")
    xs = np.r_[df.x_a_m.to_numpy(float), df.x_b_m.to_numpy(float)]
    ys = np.r_[df.y_a_m.to_numpy(float), df.y_b_m.to_numpy(float)]
    xmin, xmax = float(np.min(xs)), float(np.max(xs)); ymin, ymax = float(np.min(ys)), float(np.max(ys))
    dx = max(xmax - xmin, 1.0); dy = max(ymax - ymin, 1.0)
    m = max(float(margin_fraction), 0.0)
    xmin -= m * dx; xmax += m * dx; ymin -= m * dy; ymax += m * dy
    x_edges = np.linspace(xmin, xmax, nx + 1); y_edges = np.linspace(ymin, ymax, ny + 1)
    return Grid2D(x_edges, y_edges, (x_edges[:-1]+x_edges[1:])/2, (y_edges[:-1]+y_edges[1:])/2, nx, ny)


def _segment_lengths(x0: float, y0: float, x1: float, y1: float, grid: Grid2D) -> np.ndarray:
    dx = x1 - x0; dy = y1 - y0; total = math.hypot(dx, dy)
    if total <= 0:
        raise TomographyError("发现零长度台站路径")
    ts = [0.0, 1.0]
    if abs(dx) > 1e-12:
        ts += [float((x - x0) / dx) for x in grid.x_edges[1:-1] if 0 < (x - x0) / dx < 1]
    if abs(dy) > 1e-12:
        ts += [float((y - y0) / dy) for y in grid.y_edges[1:-1] if 0 < (y - y0) / dy < 1]
    ts = sorted(set(round(t, 14) for t in ts))
    row = np.zeros(grid.ncell, dtype=float)
    for ta, tb in zip(ts[:-1], ts[1:]):
        if tb <= ta: continue
        tm = (ta + tb) / 2.0
        xm = x0 + tm * dx; ym = y0 + tm * dy
        ix = int(np.searchsorted(grid.x_edges, xm, side="right") - 1)
        iy = int(np.searchsorted(grid.y_edges, ym, side="right") - 1)
        ix = min(max(ix, 0), grid.nx - 1); iy = min(max(iy, 0), grid.ny - 1)
        row[iy * grid.nx + ix] += (tb - ta) * total
    return row


def _design_matrix(df: pd.DataFrame, grid: Grid2D) -> np.ndarray:
    G = np.zeros((len(df), grid.ncell), dtype=float)
    for r, row in enumerate(df.itertuples(index=False)):
        g = _segment_lengths(float(row.x_a_m), float(row.y_a_m), float(row.x_b_m), float(row.y_b_m), grid)
        # Make the path integral use the same station distance as the measurement.
        # This removes tiny AEQD/geodesic inconsistencies without changing cell fractions.
        geom = float(g.sum()); target = float(row.distance_m)
        if geom <= 0: raise TomographyError("路径矩阵出现零总长度")
        G[r] = g * (target / geom) / 1000.0  # km, while slowness is s/km
    return G


def _smoothing_matrix(nx: int, ny: int) -> np.ndarray:
    rows = []
    for iy in range(ny):
        for ix in range(nx):
            j = iy * nx + ix
            if ix + 1 < nx:
                r = np.zeros(nx*ny); r[j] = 1.0; r[j+1] = -1.0; rows.append(r)
            if iy + 1 < ny:
                r = np.zeros(nx*ny); r[j] = 1.0; r[j+nx] = -1.0; rows.append(r)
    return np.asarray(rows, dtype=float) if rows else np.zeros((0, nx*ny), dtype=float)


def _weighted_system(G: np.ndarray, t: np.ndarray, velocity: np.ndarray, frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    # Optional phase-velocity standard deviation.  Otherwise equal weights.
    sigma_t = np.ones(len(frame), dtype=float)
    candidates = ["phase_velocity_std_m_s", "velocity_std_m_s", "phase_velocity_sigma_m_s"]
    found = next((c for c in candidates if c in frame.columns), None)
    if found is not None:
        sv = pd.to_numeric(frame[found], errors="coerce").to_numpy(float)
        dist = frame.distance_m.to_numpy(float)
        good = np.isfinite(sv) & (sv > 0)
        # t=d/v -> sigma_t ~= d/v^2 sigma_v
        sigma_t[good] = dist[good] / (velocity[good]**2) * sv[good]
        median = np.nanmedian(sigma_t[good]) if np.any(good) else 1.0
        sigma_t[~good] = median
        sigma_t = np.maximum(sigma_t, max(np.nanmedian(sigma_t)*0.05, 1e-6))
    else:
        sigma_t[:] = 1.0
    w = 1.0 / sigma_t
    # Normalize weights so regularization scale is not tied to seconds.
    w /= np.median(w)
    return G * w[:, None], t * w, w


def _solve_frequency(frame: pd.DataFrame, grid: Grid2D, settings: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any], Optional[pd.DataFrame]]:
    G = _design_matrix(frame, grid)
    vobs = frame.tomography_velocity_m_s.to_numpy(float)
    t_obs = frame.distance_m.to_numpy(float) / vobs
    s_ref = float(np.median(1000.0 / vobs))  # s/km
    y = t_obs - G.sum(axis=1) * s_ref
    Gw, yw, weights = _weighted_system(G, y, vobs, frame)
    col_norm = np.sqrt(np.sum(Gw**2, axis=0))
    scale = float(np.median(col_norm[col_norm > 0])) if np.any(col_norm > 0) else 1.0
    damping = float(settings.get("tomo_damping", 0.15)) * scale
    smoothing = float(settings.get("tomo_smoothing", 1.0)) * scale
    L = _smoothing_matrix(grid.nx, grid.ny)
    blocks = [Gw]; rhs = [yw]
    if damping > 0:
        blocks.append(damping * np.eye(grid.ncell)); rhs.append(np.zeros(grid.ncell))
    if smoothing > 0 and len(L):
        blocks.append(smoothing * L); rhs.append(np.zeros(len(L)))
    A = np.vstack(blocks); b = np.concatenate(rhs)
    dm, _, rank, singular = np.linalg.lstsq(A, b, rcond=None)
    slowness = s_ref + dm
    velocity = np.where(slowness > 0, 1000.0/slowness, np.nan)
    pred_t = G @ slowness
    residual = t_obs - pred_t
    rms = float(np.sqrt(np.mean(residual**2)))
    wrms = float(np.sqrt(np.mean((residual*weights)**2)))

    GTWG = Gw.T @ Gw
    H = GTWG + (damping**2)*np.eye(grid.ncell)
    if smoothing > 0 and len(L): H += (smoothing**2) * (L.T @ L)
    try:
        Hinv = linalg.pinvh(H, check_finite=False)
    except Exception:
        Hinv = np.linalg.pinv(H)
    Rdiag = np.clip(np.diag(Hinv @ GTWG), 0.0, 1.0)
    dof = max(len(t_obs) - int(rank), 1)
    sigma2 = float(np.sum((residual*weights)**2) / dof)
    s_std = np.sqrt(np.maximum(np.diag(Hinv)*sigma2, 0.0))
    v_std = np.where(slowness > 0, 1000.0/(slowness**2) * s_std, np.nan)
    ray_count = np.sum(G > 0, axis=0).astype(int)
    path_len_m = np.sum(G, axis=0) * 1000.0
    ray_axis = (np.degrees(np.arctan2(frame.y_b_m.to_numpy(float)-frame.y_a_m.to_numpy(float),
                                      frame.x_b_m.to_numpy(float)-frame.x_a_m.to_numpy(float))) % 180.0)
    az_gap = np.full(grid.ncell, np.nan, dtype=float)
    az_bins = np.zeros(grid.ncell, dtype=int)
    for j in range(grid.ncell):
        a = np.sort(ray_axis[G[:, j] > 0])
        if len(a) == 1:
            az_gap[j] = 180.0; az_bins[j] = 1
        elif len(a) > 1:
            gaps = np.diff(np.r_[a, a[0] + 180.0])
            az_gap[j] = float(np.max(gaps))
            az_bins[j] = int(len(np.unique(np.floor(a / 22.5).astype(int))))

    rows = []
    for iy, yc in enumerate(grid.y_centers):
        for ix, xc in enumerate(grid.x_centers):
            j = iy*grid.nx+ix
            rows.append({
                "cell_index": j, "ix": ix, "iy": iy, "x_m": float(xc), "y_m": float(yc),
                "phase_velocity_m_s": float(velocity[j]) if np.isfinite(velocity[j]) else np.nan,
                "phase_velocity_std_m_s": float(v_std[j]) if np.isfinite(v_std[j]) else np.nan,
                "slowness_s_km": float(slowness[j]), "reference_velocity_m_s": 1000.0/s_ref,
                "ray_count": int(ray_count[j]), "path_length_m": float(path_len_m[j]),
                "ray_azimuthal_gap_deg": float(az_gap[j]) if np.isfinite(az_gap[j]) else np.nan,
                "ray_azimuth_bin_count": int(az_bins[j]),
                "resolution": float(Rdiag[j]),
            })
    model = pd.DataFrame(rows)
    min_cell = int(settings.get("tomo_min_cell_rays", 1))
    model["coverage_pass"] = model.ray_count >= min_cell

    fit = frame.copy().reset_index(drop=True)
    fit["observed_travel_time_s"] = t_obs
    fit["predicted_travel_time_s"] = pred_t
    fit["travel_time_residual_s"] = residual
    fit["tomography_weight"] = weights

    checker = None
    checker_summary = None
    if bool(settings.get("tomo_checkerboard", True)):
        amp = float(settings.get("tomo_checkerboard_amplitude_fraction", 0.08))
        block = max(int(settings.get("tomo_checkerboard_block_cells", 2)), 1)
        pattern = np.empty(grid.ncell, dtype=float)
        for iy in range(grid.ny):
            for ix in range(grid.nx):
                pattern[iy*grid.nx+ix] = 1.0 if ((ix//block + iy//block) % 2 == 0) else -1.0
        s_true = s_ref * (1.0 + amp*pattern)
        t_syn = G @ s_true
        y_syn = t_syn - G.sum(axis=1)*s_ref
        _, yw_syn, _ = _weighted_system(G, y_syn, vobs, frame)
        b_syn_parts = [yw_syn]
        if damping > 0: b_syn_parts.append(np.zeros(grid.ncell))
        if smoothing > 0 and len(L): b_syn_parts.append(np.zeros(len(L)))
        dm_syn, *_ = np.linalg.lstsq(A, np.concatenate(b_syn_parts), rcond=None)
        recovered = dm_syn / (s_ref*amp) if amp > 0 else np.zeros_like(dm_syn)
        valid = ray_count > 0
        corr = float(np.corrcoef(pattern[valid], recovered[valid])[0,1]) if valid.sum() >= 3 and np.std(recovered[valid]) > 0 else np.nan
        ratio = float(np.std(recovered[valid]) / np.std(pattern[valid])) if valid.sum() >= 2 else np.nan
        checker = model[["cell_index","ix","iy","x_m","y_m","ray_count","resolution"]].copy()
        checker["checkerboard_true_sign"] = pattern
        checker["checkerboard_recovered_fraction"] = recovered
        checker_summary = {"correlation": corr, "amplitude_recovery_ratio": ratio, "amplitude_fraction": amp, "block_cells": block}

    summary = {
        "paths": int(len(frame)), "reference_velocity_m_s": float(1000.0/s_ref),
        "travel_time_rms_s": rms, "weighted_rms": wrms, "rank_augmented": int(rank),
        "damping_effective": damping, "smoothing_effective": smoothing,
        "covered_cells": int(np.sum(ray_count > 0)), "coverage_pass_cells": int(np.sum(ray_count >= min_cell)),
        "median_resolution": float(np.median(Rdiag[ray_count > 0])) if np.any(ray_count > 0) else 0.0,
        "median_azimuthal_gap_deg": float(np.nanmedian(az_gap[ray_count > 0])) if np.any(ray_count > 0) else np.nan,
        "checkerboard": checker_summary,
    }
    return model, fit, summary, checker


def _add_geographic_centers(frame: pd.DataFrame, inverse: Any) -> pd.DataFrame:
    q = frame.copy()
    if inverse is None:
        return q
    lon, lat = inverse.transform(q.x_m.to_numpy(float), q.y_m.to_numpy(float))
    q["longitude"] = lon; q["latitude"] = lat
    return q


def _plot_velocity_grid(frame: pd.DataFrame, frequency: float, out_png: Path, out_pdf: Path) -> None:
    import matplotlib.pyplot as plt
    pivot = frame.pivot(index="iy", columns="ix", values="phase_velocity_m_s")
    arr = pivot.to_numpy(float) / 1000.0
    mask = frame.pivot(index="iy", columns="ix", values="coverage_pass").to_numpy(bool)
    arr = np.where(mask, arr, np.nan)
    fig, ax = plt.subplots(figsize=(7.2, 6.1))
    im = ax.imshow(arr, origin="lower", aspect="auto")
    cb = fig.colorbar(im, ax=ax); cb.set_label("Phase velocity (km/s)")
    ax.set_title(f"2-D phase-velocity tomography · {frequency:g} Hz")
    ax.set_xlabel("Grid X index"); ax.set_ylabel("Grid Y index")
    fig.tight_layout(); fig.savefig(out_png, dpi=220); fig.savefig(out_pdf); plt.close(fig)


def run_phase_velocity_tomography(project: str | Path, *, paths_file: str | Path | None = None) -> dict[str, Any]:
    """Invert exported station-pair phase velocities frequency by frequency."""
    root = Path(project).resolve()
    from .geothermal import load_settings, resolve_project_path, export_tomography_paths
    settings = load_settings(root)
    out = root / "output" / "geothermal" / "tomography2d"
    figdir = out / "figures"; checkdir = out / "checkerboard"
    out.mkdir(parents=True, exist_ok=True); figdir.mkdir(exist_ok=True); checkdir.mkdir(exist_ok=True)

    if paths_file is None:
        configured = settings.get("tomography_paths_file", "")
        source = resolve_project_path(root, configured) if configured else root / "output" / "geothermal" / "tomography_paths.csv"
        # When using the standard Studio handoff, refresh it if the project has a
        # current candidate/curve source. A standalone imported tomography_paths.csv
        # (such as the synthetic example) remains usable without fabricating families.
        if not configured:
            spec={}
            if (root/"project.json").exists():
                try: spec=json.loads((root/"project.json").read_text(encoding="utf-8-sig"))
                except Exception: spec={}
            can_export=(root/"output/families_complete.json").exists() or spec.get("analysis_mode")=="curves_only"
            if can_export:
                case_key = str(settings.get("tomo_selection_case", "") or "").strip()
                cases=[]
                sel_path=root/"selection.json"
                if sel_path.exists():
                    try: cases=list(json.loads(sel_path.read_text(encoding="utf-8-sig")).get("cases", []))
                    except Exception: cases=[]
                if case_key:
                    manifest=export_tomography_paths(root, selected_only=True, selection_case=case_key)
                    source=Path(manifest["output"])
                elif len(cases)==1:
                    manifest=export_tomography_paths(root, selected_only=True, selection_case=str(cases[0].get("name", "0")))
                    source=Path(manifest["output"])
                else:
                    manifest=export_tomography_paths(root, selected_only=bool(settings.get("selected_paths_only", False)))
                    source=Path(manifest["output"])
            elif not source.exists():
                raise TomographyError("没有 tomography_paths.csv，且项目尚未完成候选频散阶段，无法自动导出层析路径")
        elif not source.exists():
            raise FileNotFoundError(source)
    else:
        source = Path(paths_file).resolve()
    if not source.exists(): raise FileNotFoundError(source)

    raw = pd.read_csv(source)
    obs, warnings = _prepare_observations(raw, settings)
    obs, coord_meta, inverse = _project_paths(obs)
    obs = obs.dropna(subset=["x_a_m","y_a_m","x_b_m","y_b_m"]).copy()
    if obs.empty: raise TomographyError("坐标投影后没有有效路径")

    nx = int(settings.get("tomo_nx", 12)); ny = int(settings.get("tomo_ny", 12))
    grid = _make_grid(obs, nx, ny, float(settings.get("tomo_margin_fraction", 0.05)))
    min_paths = int(settings.get("tomo_min_paths_per_frequency", 6))
    max_fig = int(settings.get("tomo_max_figures", 12))
    all_models=[]; all_fits=[]; catalog=[]; checker_catalog=[]
    frequencies = sorted(float(x) for x in obs.frequency_bin_hz.unique())
    print(f"【二维层析】开始：{len(obs)} 个观测，{len(frequencies)} 个频率，网格 {nx}×{ny}", flush=True)
    plotted = 0
    for freq in frequencies:
        f = obs[np.isclose(obs.frequency_bin_hz, freq, atol=5e-9, rtol=0)].copy()
        if len(f) < min_paths:
            catalog.append({"frequency_hz":freq,"period_s":1/freq,"paths":len(f),"status":"skipped_too_few_paths"})
            print(f"【二维层析跳过】{freq:g} Hz：{len(f)} 条路径 < 最少 {min_paths}", flush=True)
            continue
        model, fit, summary, checker = _solve_frequency(f, grid, settings)
        model = _add_geographic_centers(model, inverse); model["frequency_hz"]=freq; model["period_s"]=1/freq
        fit["frequency_bin_hz"] = freq
        all_models.append(model); all_fits.append(fit)
        cb=summary.pop("checkerboard", None)
        entry={"frequency_hz":freq,"period_s":1/freq,"status":"accepted",**summary}
        if cb:
            entry.update(checkerboard_correlation=cb.get("correlation"),checkerboard_amplitude_recovery_ratio=cb.get("amplitude_recovery_ratio"),checkerboard_amplitude_fraction=cb.get("amplitude_fraction"),checkerboard_block_cells=cb.get("block_cells"))
        catalog.append(entry)
        summary["checkerboard"]=cb
        if checker is not None:
            checker=_add_geographic_centers(checker,inverse);checker["frequency_hz"]=freq;checker.to_csv(checkdir/f"checkerboard_{freq:.6f}Hz.csv",index=False)
            checker_catalog.append({"frequency_hz":freq,**(summary.get("checkerboard") or {})})
        if plotted < max_fig:
            png=figdir/f"phase_velocity_{freq:.6f}Hz.png";pdf=figdir/f"phase_velocity_{freq:.6f}Hz.pdf"
            _plot_velocity_grid(model,freq,png,pdf);plotted += 1
        print(f"【二维层析】{freq:g} Hz · {len(f)}路径 · RMS={summary['travel_time_rms_s']:.4g}s · 覆盖{summary['covered_cells']}/{grid.ncell}格", flush=True)

    if not all_models:
        raise TomographyError("没有频率达到最少路径数，未生成二维层析结果")
    models=pd.concat(all_models,ignore_index=True);fits=pd.concat(all_fits,ignore_index=True)
    models.to_csv(out/"phase_velocity_tomography.csv",index=False)
    fits.to_csv(out/"data_fit.csv",index=False)
    pd.DataFrame(catalog).to_csv(out/"frequency_catalog.csv",index=False)
    if checker_catalog: pd.DataFrame(checker_catalog).to_csv(out/"checkerboard_summary.csv",index=False)

    # Local dispersion is the direct handoff to the next node-wise Vs inversion stage.
    local = models[models.coverage_pass.map(_truth)].copy()
    keep=[c for c in ["cell_index","ix","iy","x_m","y_m","longitude","latitude","frequency_hz","period_s","phase_velocity_m_s","phase_velocity_std_m_s","resolution","ray_count","path_length_m","ray_azimuthal_gap_deg","ray_azimuth_bin_count","reference_velocity_m_s"] if c in local]
    local=local[keep].sort_values(["cell_index","frequency_hz"])
    local.to_csv(out/"local_dispersion.csv",index=False)

    grid_def={
        "nx":nx,"ny":ny,"x_edges_m":grid.x_edges.tolist(),"y_edges_m":grid.y_edges.tolist(),
        **coord_meta,
    }
    write_json(out/"grid_definition.json",grid_def)
    result={
        "complete":True,"source":str(source),"source_sha256":sha(source),"observations_used":int(len(obs)),
        "frequencies_requested":len(frequencies),"frequencies_inverted":int(sum(x.get("status")=="accepted" for x in catalog)),
        "grid":{"nx":nx,"ny":ny,"cells":grid.ncell},"coordinate_system":coord_meta,
        "outputs":{
            "phase_velocity_tomography":str(out/"phase_velocity_tomography.csv"),
            "local_dispersion":str(out/"local_dispersion.csv"),"data_fit":str(out/"data_fit.csv"),
            "frequency_catalog":str(out/"frequency_catalog.csv"),"grid_definition":str(out/"grid_definition.json")},
        "warnings":warnings,
        "scientific_guards":[
            "Competing branches for the same pair/frequency are rejected by default rather than counted as independent rays.",
            "Direction-projected plane-wave speeds are converted to a station-pair-equivalent observable; strongly oblique directions are excluded by the configured angle limit.",
            "Straight-ray phase-velocity tomography is an approximation; resolution and checkerboard recovery must be inspected before interpreting anomalies.",
            "The output is phase velocity c(x,y,f), not Vs(x,y,z). local_dispersion.csv is the handoff to node-wise depth inversion.",
        ],
    }
    write_json(out/"tomography_2d_complete.json",result)
    _write_tomography_report(out,result,pd.DataFrame(catalog),checker_catalog)
    print("【二维层析】完成 -> "+str(out/"tomography_2d_report.html"),flush=True)
    return result


def _write_tomography_report(out: Path, result: dict[str,Any], catalog: pd.DataFrame, checker: list[dict[str,Any]]) -> None:
    parts=['<meta charset="utf-8"><style>body{font:16px/1.7 system-ui,"Microsoft YaHei",sans-serif;max-width:1180px;margin:30px auto;padding:0 22px;color:#203548}h1,h2{color:#174963}.guard{background:#eef6f8;border-left:5px solid #24748a;padding:12px}.warn{background:#fff4e5;border-left:5px solid #d97706;padding:12px}table{border-collapse:collapse}td,th{border:1px solid #ccd7df;padding:6px}img{max-width:100%}</style>']
    parts.append('<h1>二维相速度层析</h1>')
    parts.append('<div class="guard"><b>科学边界</b><ul>'+''.join('<li>'+html.escape(x)+'</li>' for x in result['scientific_guards'])+'</ul></div>')
    if result['warnings']:parts.append('<div class="warn"><ul>'+''.join('<li>'+html.escape(x)+'</li>' for x in result['warnings'])+'</ul></div>')
    parts.append(f"<p>输入观测 {result['observations_used']}；反演频率 {result['frequencies_inverted']}/{result['frequencies_requested']}；网格 {result['grid']['nx']}×{result['grid']['ny']}。</p>")
    parts.append('<h2>频率诊断</h2>'+catalog.to_html(index=False,float_format=lambda x:f"{x:.5g}"))
    if checker:parts.append('<h2>Checkerboard</h2>'+pd.DataFrame(checker).to_html(index=False,float_format=lambda x:f"{x:.5g}"))
    figs=sorted((out/'figures').glob('*.png'))[:12]
    if figs:
        parts.append('<h2>相速度切片</h2>')
        for f in figs:parts.append(f'<figure><img src="figures/{html.escape(f.name)}"><figcaption>{html.escape(f.stem)}</figcaption></figure>')
    parts.append('<h2>下一阶段</h2><p><code>local_dispersion.csv</code> 已按网格单元组织 c(f)。下一阶段应对满足覆盖和分辨率门限的网格单元进行局部 Rayleigh 频散 Vs(z) 反演，并将不确定性继续传播到三维 Vs。</p>')
    (out/'tomography_2d_report.html').write_text(''.join(parts),encoding='utf-8')


def tomography_status(project: str | Path) -> dict[str,Any]:
    p=Path(project).resolve()/"output/geothermal/tomography2d/tomography_2d_complete.json"
    if not p.exists():return {"complete":False,"message":"尚未运行二维相速度层析"}
    return json.loads(p.read_text(encoding="utf-8-sig"))


def create_tomography_example(directory: str | Path) -> Path:
    """Create a deterministic WGS84 straight-ray tomography example project."""
    from .project import new_project
    from .geothermal import save_settings, DEFAULT_SETTINGS
    root = new_project(directory)
    try:
        from pyproj import CRS, Transformer
    except Exception as exc:
        raise TomographyError("示例需要 pyproj") from exc
    lon0, lat0 = 10.05, 50.02
    local = CRS.from_proj4(f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +datum=WGS84 +units=m +no_defs")
    inv = Transformer.from_crs(local, "EPSG:4326", always_xy=True)
    pts = {
        'A':(0.,0.),'B':(2000.,0.),'C':(4000.,0.),'D':(4000.,2000.),
        'E':(4000.,4000.),'F':(2000.,4000.),'G':(0.,4000.),'H':(0.,2000.),
    }
    station_rows=[]
    for name,(x,y) in pts.items():
        lon,lat=inv.transform(x,y)
        station_rows.append(dict(station=name,seed_id=f'XX.{name}..HHZ',longitude=lon,latitude=lat,file='unused.npz',precorrected=True,units='m/s',polarity=1))
    pd.DataFrame(station_rows).to_csv(root/'stations.csv',index=False)
    grid=Grid2D(np.linspace(0,4000,5),np.linspace(0,4000,5),np.array([500,1500,2500,3500]),np.array([500,1500,2500,3500]),4,4)
    pattern=np.array([1 if ((ix//2+iy//2)%2==0) else -1 for iy in range(4) for ix in range(4)],float)
    rows=[]; names=list(pts)
    for freq,base_v in [(0.8,2300.),(1.2,2050.),(1.8,1800.)]:
        s0=1000/base_v; s=s0*(1+0.055*pattern)
        for i,a in enumerate(names):
            for b in names[i+1:]:
                x0,y0=pts[a];x1,y1=pts[b];d=math.hypot(x1-x0,y1-y0)
                if d<2500:continue
                G=_segment_lengths(x0,y0,x1,y1,grid)/1000.0;t=float(G@s);v=d/t
                lona,lata=inv.transform(x0,y0);lonb,latb=inv.transform(x1,y1)
                rows.append(dict(curve_id=f'{freq:g}_{a}_{b}',pair=f'{a}--{b}',frequency_hz=freq,phase_velocity_m_s=v,distance_m=d,
                    longitude_a=lona,latitude_a=lata,longitude_b=lonb,latitude_b=latb,
                    default_quality_pass=True,background_discriminated=True,velocity_type='conditional_pair_velocity',measurement_method='synthetic_straight_ray'))
    out=root/'output/geothermal';out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(out/'tomography_paths.csv',index=False)
    save_settings(root,dict(DEFAULT_SETTINGS,tomo_nx=4,tomo_ny=4,tomo_margin_fraction=0.0,tomo_damping=0.02,tomo_smoothing=0.08,tomo_min_paths_per_frequency=6,tomo_checkerboard=True,tomo_checkerboard_block_cells=2,tomo_max_figures=3,title='2-D phase-velocity tomography demo'))
    return root
