from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from ndstudio.geothermal import (
    DEFAULT_SETTINGS,
    apply_quality_mask,
    create_example,
    export_tomography_paths,
    generate_geothermal,
    normalize_spatial_vs,
    save_settings,
)
from ndstudio.project import new_project
from ndstudio.candidates import import_curve_csv


def test_normalize_spatial_csv_and_quality_mask(tmp_path: Path):
    source = tmp_path / "spatial.csv"
    pd.DataFrame(
        [
            [10.0, 50.0, 1.0, 2.50, 0.10, 0.80, 8],
            [10.1, 50.0, 1.0, 2.70, 0.50, 0.40, 2],
        ],
        columns=["lon", "lat", "depth_km", "vs_km_s", "vs_std_km_s", "resolution", "ray_count"],
    ).to_csv(source, index=False)

    frame = normalize_spatial_vs(source)
    assert frame.depth_m.tolist() == [1000.0, 1000.0]
    assert frame.vs_m_s.tolist() == [2500.0, 2700.0]
    assert frame.vs_std_m_s.tolist() == [100.0, 500.0]

    settings = dict(DEFAULT_SETTINGS)
    settings.update(max_relative_vs_std=0.10, min_resolution=0.60, min_ray_count=5)
    masked = apply_quality_mask(frame, settings)
    assert masked.quality_pass.tolist() == [True, False]
    assert "vs_std_too_large" in masked.loc[1, "quality_reasons"]
    assert "resolution" in masked.loc[1, "quality_reasons"]
    assert "ray_count" in masked.loc[1, "quality_reasons"]


def _curves_frame(distance_by_pair: dict[str, float]) -> pd.DataFrame:
    rows = []
    for curve_id, pair, wrong_source_index in [
        ("C01", "S0--S2", 0),  # pair label is authoritative; pair_index intentionally stale
        ("C02", "S1--S2", 1),
    ]:
        for f, v in [(0.5, 1800.0), (0.7, 1650.0), (0.9, 1500.0)]:
            rows.append(
                dict(
                    curve_id=curve_id,
                    segment_id="seg_" + curve_id,
                    cohort="all",
                    pair_index=wrong_source_index,
                    pair=pair,
                    branch="causal",
                    frequency_hz=f,
                    phase_velocity_m_s=v,
                    group_velocity_m_s=v * 0.95,
                    constant_cycle=1,
                    source_phase_rad=0.0,
                    distance_over_wavelength=4.0,
                    instrument_response_corrected=True,
                    background_discriminated=True,
                    source_phase_verified=False,
                    is_model_prediction=False,
                    distance_m=distance_by_pair[pair],
                )
            )
    return pd.DataFrame(rows)


def test_tomography_export_resolves_pair_label_before_stale_pair_index(tmp_path: Path):
    root = new_project(tmp_path / "project")
    # Deliberately choose coordinates where S0--S2 is not pair_index 0 in the
    # canonical combinations ordering.  That reproduces archived/stale indices.
    stations = pd.DataFrame(
        [
            ["S0", "XX.S0..HHZ", 50.00, 10.00],
            ["S1", "XX.S1..HHZ", 50.00, 10.05],
            ["S2", "XX.S2..HHZ", 50.06, 10.00],
        ],
        columns=["station", "seed_id", "latitude", "longitude"],
    )
    stations["file"] = "unused.npz"
    stations["precorrected"] = True
    stations["units"] = "m/s"
    stations["polarity"] = 1
    stations.to_csv(root / "stations.csv", index=False)

    # Get canonical distances so import/export distance auditing remains strict.
    from dispersion.common import load_config
    from dispersion.prepare import station_table

    _, pairs, _ = station_table(load_config(root / "config.json"))
    distances = dict(zip(pairs.pair, pairs.distance_m))
    source = tmp_path / "curves.csv"
    _curves_frame(distances).to_csv(source, index=False)
    import_curve_csv(root, source)

    manifest = export_tomography_paths(root)
    result = pd.read_csv(manifest["output"])
    s0s2 = result[result.pair == "S0--S2"]
    assert len(s0s2) == 3
    assert set(s0s2.geometry_pair_index) == {1}  # canonical S0--S2 index
    assert set(s0s2.source_pair_index) == {0}    # stale source index preserved for audit
    assert set(s0s2.station_a) == {"S0"}
    assert set(s0s2.station_b) == {"S2"}
    assert set(s0s2.coordinate_system) == {"WGS84"}


def test_generate_geothermal_example_with_matplotlib(tmp_path: Path):
    root = create_example(tmp_path / "example")
    # Reduce figures for test speed but preserve multi-depth behavior.
    settings = json.loads((root / "geothermal_settings.json").read_text(encoding="utf-8"))
    settings["slice_depths_m"] = [500.0, 2000.0]
    settings["backend"] = "matplotlib"
    save_settings(root, settings)

    result = generate_geothermal(root)
    assert result["complete"] is True
    assert len(result["products"]) == 2
    for product in result["products"]:
        assert Path(product["png"]).is_file()
        assert Path(product["pdf"]).is_file()
        assert Path(product["csv"]).is_file()
    assert (root / "output" / "geothermal" / "report.html").is_file()
    assert (root / "output" / "geothermal" / "spatial_vs_standardized.csv").is_file()
