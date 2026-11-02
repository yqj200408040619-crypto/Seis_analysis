"""Probe: run select_array + aperture/frequency gate for a chosen candidate.

Live FDSN station-metadata query only (no bulk waveform download). Writes a
UTF-8 JSON summary so the result can be read directly regardless of console
codepage. Reuses the program's own candidate definitions and gate.
"""
import sys, json, argparse, traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ndstudio.array_search_e2e import CANDIDATES, resolvable_band, _norm_utc
from ndstudio.acquisition import select_array


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--name', default='Parkfield_NC')
    ap.add_argument('--start', default='2023-11-01')
    ap.add_argument('--end', default='2023-11-04')
    ap.add_argument('--out', default=str(HERE / 'runs' / 'array_search'))
    a = ap.parse_args()

    out_root = Path(a.out).resolve()
    cand = next(c for c in CANDIDATES if c.name == a.name)
    root = out_root / cand.name
    root.mkdir(parents=True, exist_ok=True)
    selected = root / 'array_selection.json'
    summary = {'candidate': cand.name, 'provider': cand.provider, 'network': cand.network,
               'start': a.start, 'end': a.end}
    try:
        res = select_array(
            cand.provider, cand.network, a.start, a.end, selected,
            station=cand.station, channel=cand.channel, size=cand.size,
            min_spacing_m=cand.min_spacing_m, max_spacing_m=cand.max_spacing_m,
            max_aperture_m=cand.max_aperture_m,
            min_stations=cand.min_stations, max_stations=cand.max_stations,
            minlatitude=cand.minlatitude, maxlatitude=cand.maxlatitude,
            minlongitude=cand.minlongitude, maxlongitude=cand.maxlongitude)
        geom = res['geometry']
        band = resolvable_band(geom, cand.velocity_min_m_s, cand.velocity_max_m_s)
        summary.update(
            ok=True,
            eligible_stations=res['eligible_stations'],
            metadata_rows_queried=res['metadata_rows_queried'],
            selected_station_count=len(res['stations']),
            selected_stations=[s['seed_id'] for s in res['stations']],
            anchor=res['anchor_seed_id'],
            aperture_m=geom['max_pair_distance_m'],
            min_nn_m=min(geom['nearest_neighbor_m'].values()),
            max_nn_m=max(geom['nearest_neighbor_m'].values()),
            mean_nn_m=geom['mean_nearest_neighbor_m'],
            band=band,
            velocity_range_m_s=[cand.velocity_min_m_s, cand.velocity_max_m_s],
        )
    except Exception as exc:
        summary.update(ok=False, error=f'{type(exc).__name__}: {exc}',
                       traceback=traceback.format_exc())

    out_json = root / 'select_probe_summary.json'
    out_json.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    print('WROTE', out_json)


if __name__ == '__main__':
    main()
