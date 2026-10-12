"""Lossless FDSN archive and Noise Dispersion Studio 2.2 interchange contract.

No filtering, response removal, resampling, merging, or gap filling here.
Raw MiniSEED and returned StationXML bytes are embedded in each HDF5 archive.
"""
from pathlib import Path
from io import BytesIO
import hashlib
import json
import os
from collections.abc import Mapping

import h5py
import numpy as np
from obspy import UTCDateTime, read, read_inventory
from pyproj import Geod

SCHEMA = "h5data.ndstudio/1.0"
XML_DATASET = "/metadata/stationxml"
GEOD = Geod(ellps="WGS84")


def utc9(value):
    return str(UTCDateTime(ns=UTCDateTime(value).ns, precision=9))


def json_value(value):
    if isinstance(value, UTCDateTime):
        return utc9(value)
    if isinstance(value, Mapping):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_value(value.tolist())
    if isinstance(value, np.generic):
        return json_value(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def geometry(rows):
    """WGS84 surface distances in m, one point per NET.STA, bounded memory.

    Primary spacing is mean nearest-neighbour distance. Pair mean is separate;
    neither depends on alphabetic station ordering or counts channels twice.
    Conflicting historical coordinates are excluded and reported.
    """
    coords = {}
    stations = set()
    conflicts = set()
    for r in rows:
        key = (r.network, r.station)
        stations.add(key)
        lat, lon = r.latitude, r.longitude
        if lat is None or lon is None or not np.isfinite([lat, lon]).all():
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            continue
        if key in coords and not np.allclose(coords[key][0:2], (lat, lon), rtol=0, atol=1e-6):
            conflicts.add(key)
        else:
            coords[key] = (lat, lon, r.elevation)
    points = [(net, sta, *coords[(net, sta)]) for net, sta in sorted(coords) if (net, sta) not in conflicts]
    n = len(points)
    result = dict(stations=len(stations), stations_with_coordinates=n,
                  missing_coordinates=len(stations - set(coords)),
                  coordinate_conflicts=[".".join(k) for k in sorted(conflicts)],
                  mean_nearest_neighbor_m=None, mean_pair_distance_m=None,
                  min_pair_distance_m=None, max_pair_distance_m=None,
                  coincident_pair_count=0, spacing_definition="mean nearest-neighbour WGS84 surface distance",
                  points=points, nearest_neighbor_m={})
    if n < 2:
        return result
    lat = np.array([p[2] for p in points]); lon = np.array([p[3] for p in points])
    nearest = np.full(n, np.inf); total = 0.; count = 0; low = np.inf; high = 0.
    for i in range(n - 1):
        _, _, d = GEOD.inv(np.full(n-i-1, lon[i]), np.full(n-i-1, lat[i]), lon[i+1:], lat[i+1:])
        d = np.asarray(d)
        nearest[i] = min(nearest[i], float(d.min()))
        nearest[i+1:] = np.minimum(nearest[i+1:], d)
        total += float(d.sum()); count += len(d)
        low = min(low, float(d.min())); high = max(high, float(d.max()))
        result["coincident_pair_count"] += int(np.sum(d < .001))
    result.update(mean_nearest_neighbor_m=float(nearest.mean()), mean_pair_distance_m=total/count,
                  min_pair_distance_m=low, max_pair_distance_m=high,
                  nearest_neighbor_m={f"{p[0]}.{p[1]}": float(d) for p, d in zip(points, nearest)})
    return result


def fetch_stationxml(client, row, start, end):
    """Preserve the exact response-level StationXML returned by the server."""
    buf = BytesIO()
    inv = client.get_stations(network=row.network, station=row.station,
                              location=row.location or "--", channel=row.channel,
                              starttime=start, endtime=end, level="response",
                              format="xml", filename=buf)
    payload = buf.getvalue()
    if not payload and inv is not None:  # Custom/test clients returning Inventory.
        inv.write(buf, format="STATIONXML")
        payload = buf.getvalue()
    if not payload:
        raise ValueError("Provider returned empty response-level StationXML")
    return payload, read_inventory(BytesIO(payload), format="STATIONXML")


def matching_channel(inv, tr):
    candidates = [ch for net in inv if net.code == tr.stats.network
                  for sta in net if sta.code == tr.stats.station
                  for ch in sta if ch.code == tr.stats.channel and ch.location_code == tr.stats.location
                  and (ch.start_date is None or ch.start_date.ns <= tr.stats.starttime.ns)
                  and (ch.end_date is None or ch.end_date.ns >= tr.stats.endtime.ns)]
    if len(candidates) != 1:
        raise ValueError(f"{tr.id}: unique metadata epoch must cover every actual sample; split at response changes")
    ch = candidates[0]
    if ch.response is None or not ch.response.response_stages or ch.response.instrument_sensitivity is None:
        raise ValueError(f"{tr.id}: complete instrument response is missing")
    if not np.isclose(ch.sample_rate, tr.stats.sampling_rate, atol=1e-6, rtol=0):
        raise ValueError(f"{tr.id}: metadata/native sampling-rate mismatch")
    if ch.latitude is None or ch.longitude is None or not np.isfinite([ch.latitude, ch.longitude]).all():
        raise ValueError(f"{tr.id}: missing or invalid channel coordinates")
    if not (-90 <= ch.latitude <= 90 and -180 <= ch.longitude <= 180):
        raise ValueError(f"{tr.id}: invalid channel coordinates")
    # Raw FDSN MiniSEED must be digital counts. Do not invent a unit declaration.
    out_units = str(ch.response.instrument_sensitivity.output_units or "").upper()
    if out_units not in ("COUNTS", "COUNT", "DIGITAL COUNTS"):
        raise ValueError(f"{tr.id}: response output units {out_units!r} are not digital counts")
    return ch


def export_h5(mseed_filename, h5_filename, *, stationxml=None, response_xml="",
              provider="", request=None, array_geometry=None, strict=False, verbose=False):
    """Atomic conversion. strict requires metadata; legacy conversion is archival.

    A legacy file without response remains clearly incomplete and cannot pass
    Studio's import audit. Values and original timing/encoding are preserved.
    """
    source = Path(mseed_filename); target = Path(h5_filename)
    raw = source.read_bytes(); st = read(BytesIO(raw), format="MSEED")
    if not st:
        raise ValueError("MiniSEED has no traces")
    xml = stationxml.encode("utf-8") if isinstance(stationxml, str) else stationxml
    inv = read_inventory(BytesIO(xml), format="STATIONXML") if xml else None
    if strict and inv is None:
        raise ValueError("Complete StationXML is required for dispersion-ready H5")
    channels = [matching_channel(inv, tr) if inv else None for tr in st]
    temp = target.with_name(target.name + ".part")
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        with h5py.File(temp, "w") as h:
            h.attrs.update(h5data_schema=SCHEMA, producer="h5data 2.2", created_utc=utc9(UTCDateTime()),
                           provider=str(provider), raw_units="counts", metadata_complete=inv is not None,
                           preprocessing="none", source_mseed_sha256=hashlib.sha256(raw).hexdigest())
            meta = h.create_group("metadata"); meta.attrs["role"] = "metadata"
            archive = meta.create_dataset("original_mseed", data=np.frombuffer(raw, dtype=np.uint8), compression="gzip")
            archive.attrs["role"] = "binary_archive"
            text_type = h5py.string_dtype("utf-8")
            if xml:
                meta.create_dataset("stationxml", data=xml.decode("utf-8"), dtype=text_type)
                h.attrs["stationxml_sha256"] = hashlib.sha256(xml).hexdigest()
            for name, value in (("download_request", request or {}), ("array_geometry", array_geometry or {})):
                meta.create_dataset(name, data=json.dumps(json_value(value), ensure_ascii=False, allow_nan=False), dtype=text_type)
            for i, (tr, ch) in enumerate(zip(st, channels)):
                name = tr.id if tr.id not in h else f"{tr.id}_{i}"
                d = h.create_dataset(name, data=tr.data, compression="gzip", shuffle=True)
                d.attrs.update(role="waveform", network=tr.stats.network, station=tr.stats.station,
                               station_label=f"{tr.stats.network}_{tr.stats.station}",
                               location=tr.stats.location, channel=tr.stats.channel, seed_id=tr.id,
                               sampling_rate=tr.stats.sampling_rate, sampling_rate_hz=tr.stats.sampling_rate,
                               starttime=utc9(tr.stats.starttime), start_utc=utc9(tr.stats.starttime),
                               start_utc_ns=tr.stats.starttime.ns, endtime=utc9(tr.stats.endtime),
                               end_utc_ns=tr.stats.endtime.ns, npts=tr.stats.npts, delta=tr.stats.delta,
                               time_format="iso_utc", units="counts", precorrected=False,
                               polarity=1, polarity_verified=False, clock_offset_s=0.,
                               clock_offset_verified=False, response_xml=str(response_xml),
                               metadata_complete=ch is not None, usable_for_dispersion=tr.stats.channel.endswith("Z") and ch is not None,
                               trace_stats_json=json.dumps(json_value(dict(tr.stats)), ensure_ascii=False))
                if xml:
                    d.attrs["response_embedded_dataset"] = XML_DATASET
                if ch is not None:
                    optional = dict(elevation_m=ch.elevation, depth_m=ch.depth,
                                    azimuth_deg=ch.azimuth, dip_deg=ch.dip)
                    d.attrs.update({k: v for k, v in optional.items() if v is not None})
                    d.attrs.update(latitude=ch.latitude, longitude=ch.longitude,
                                   response_input_units=str(ch.response.instrument_sensitivity.input_units),
                                   response_output_units=str(ch.response.instrument_sensitivity.output_units),
                                   instrument_sensitivity=ch.response.instrument_sensitivity.value,
                                   sensor_description=str(getattr(ch.sensor, "description", "") or ""),
                                   epoch_start_utc=utc9(ch.start_date) if ch.start_date else "",
                                   epoch_end_utc=utc9(ch.end_date) if ch.end_date else "")
            h.flush()
        os.replace(temp, target)
    except Exception:
        temp.unlink(missing_ok=True)
        raise
    if verbose:
        print(f"Saved {target}: {len(st)} trace(s), metadata_complete={inv is not None}")
    return len(st)


def h5_is_complete(path):
    try:
        with h5py.File(path, "r") as h:
            if h.attrs.get("h5data_schema") != SCHEMA or not bool(h.attrs.get("metadata_complete")) or XML_DATASET not in h:
                return False
            xml = h[XML_DATASET][()]
            if isinstance(xml, str):
                xml = xml.encode("utf-8")
            return hashlib.sha256(xml).hexdigest() == h.attrs.get("stationxml_sha256")
    except (OSError, ValueError):
        return False


def save_metadata(client, row, start, end, root, cache):
    """One full response request per channel/epoch; content-addressed sidecar."""
    key = (row.seed_id, start.ns, end.ns)
    if key not in cache:
        payload, _ = fetch_stationxml(client, row, start, end)
        digest = hashlib.sha256(payload).hexdigest()
        path = Path(root) / "metadata" / (digest + ".xml")
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists():
            # Different channel requests can return identical StationXML.
            # Each writer needs its own temporary file before atomic publish.
            import tempfile
            with tempfile.NamedTemporaryFile(dir=path.parent, prefix=digest+'.', suffix='.part', delete=False) as handle:
                handle.write(payload); tmp = Path(handle.name)
            try:
                os.replace(tmp, path)
            finally:
                tmp.unlink(missing_ok=True)
        cache[key] = (payload, path.relative_to(root).as_posix())
    return cache[key]


def write_run_manifest(root, jobs, results, provider="", **options):
    rows = [j.row for j in jobs]
    paths = []
    for r in results:
        paths.extend(str(Path(p).resolve().relative_to(Path(root).resolve())) for p in r.h5_files)
    manifest = dict(schema=SCHEMA, provider=str(provider), options=options,
                    geometry=geometry(rows), time_convention="UTC; actual samples preserved; no gap interpolation",
                    jobs=[dict(seed_id=j.row.seed_id, start_utc=utc9(j.start), end_utc=utc9(j.end), row=json_value(vars(j.row))) for j in jobs],
                    results=json_value([vars(r) for r in results]), written_h5_files=paths)
    path = Path(root)/"download_manifest.json"
    temp = path.with_suffix(".json.part")
    temp.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    os.replace(temp, path)
    return path
