"""Step 1: high-frequency dispersion via NDS 'traditional route' (scan).

Re-preprocess the already-downloaded LA dense-core waveforms at 20 Hz (so 5 Hz < Nyquist)
and run the scan (the program's traditional FTAN/multi-route group+phase route) over a
HIGH-frequency band (user: ~2 s to 5 Hz). Then read conditional_dispersion.csv to find a
stable continuous band to anchor the later extrapolation.

Reuses the LA-core 14-day download under runs/core_la_14d/download (100 Hz native MiniSEED
embedded in the H5), so NO new network download is needed; only ~2 days are analysed.
"""
import sys, json, os, argparse, traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import ndstudio
import h5data as download
from dispersion.common import write_json
from ndstudio.acquisition import audit_acquisition, screen_array
from ndstudio.project import new_project
from ndstudio.real_compare import preprocess
from obspy import UTCDateTime
import dataclasses

PROVIDER, NETWORK = 'SCEDC', 'CI'
CORE = ['BHP', 'GR2', 'HLL', 'LCG', 'LGB', 'PDR', 'SMF2', 'USC']
# ~2 days
START = os.environ.get('HF_START', '2023-11-01')
END = os.environ.get('HF_END', '2023-11-03')


def pick_rows(stations, start, end):
    cli = download.make_client(PROVIDER)
    rows = download.list_stations(PROVIDER, NETWORK, station=','.join(stations), channel='HHZ,BHZ',
                                  starttime=UTCDateTime(start), endtime=UTCDateTime(end), client=cli)
    by = {}
    for r in rows:
        if (r.channel.endswith('Z') and r.latitude is not None and r.sample_rate and r.sample_rate >= 20
                and not (r.start and r.start > UTCDateTime(start)) and not (r.end and r.end < UTCDateTime(end))):
            by.setdefault(r.station, []).append(r)
    chosen = []
    for s in stations:
        cs = by.get(s, [])
        if cs:
            cs.sort(key=lambda r: (0 if r.channel == 'HHZ' else 1, -(r.sample_rate or 0)))
            chosen.append(cs[0])
    return chosen


def build_selection(out_path, start, end):
    rows = pick_rows(CORE, start, end); rows.sort(key=lambda r: r.seed_id)
    g = download.array_geometry(rows)
    ser = lambda r: {k: (str(v) if isinstance(v, UTCDateTime) else v) for k, v in dataclasses.asdict(r).items()}
    result = {'schema': 'ndstudio.array_selection/1', 'provider': PROVIDER, 'network': NETWORK,
              'start_utc': str(UTCDateTime(start)), 'end_utc': str(UTCDateTime(end)),
              'anchor_seed_id': rows[0].seed_id,
              'constraints': {'station_count': len(rows), 'minimum_station_count': 6,
                              'maximum_station_count': 50, 'minimum_nearest_neighbor_m': 0.,
                              'maximum_nearest_neighbor_m': 1e9, 'maximum_pair_distance_m': 1e9},
              'metadata_rows_queried': len(rows), 'eligible_stations': len(rows),
              'station_pattern': ','.join(CORE), 'station_filters': {},
              'matching_station_catalog': download.station_catalog(rows),
              'channel_filter': 'HHZ,BHZ', 'rejected_metadata': [], 'geometry': g,
              'stations': [ser(r) for r in rows], 'selection_uses_waveforms_or_previous_curves': False}
    write_json(out_path, result); return result, g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workdir', default=str(HERE/'runs'/'hf_core'))
    ap.add_argument('--source-download', default=str(HERE/'runs'/'core_la_14d'/'download'),
                    help='existing download dir to reuse (copied in by caller)')
    a = ap.parse_args()
    root = Path(a.workdir).resolve(); root.mkdir(parents=True, exist_ok=True)
    data = root/'download'; project = root/'project'; selected = root/'array_selection.json'

    result, g = build_selection(selected, START, END)
    print('HF selection: n=%d aperture=%.0f m, band target 0.2-5 Hz' % (len(result['stations']), g['max_pair_distance_m']), flush=True)

    try:
        from ndstudio.acquisition import fetch_array
        fetch_array(selected, data, workers=4, chunk_seconds=21600., allow_missing_data=True)
        audit_acquisition(data)
        screen_array(selected, data, project=None)
    except Exception as exc:
        (root/'error.json').write_text(json.dumps({'stage':'audit/screen','error':f'{type(exc).__name__}: {exc}','tb':traceback.format_exc()},ensure_ascii=False,indent=2),encoding='utf-8')
        print('ERROR audit/screen (error.json)', flush=True); return

    # 20 Hz so that scan_high 5 Hz < Nyquist 10 Hz. Bandpass encloses [0.2,5].
    fs = 20.0
    if not (project/'project.json').exists():
        new_project(project)
        spec = json.loads((project/'project.json').read_text(encoding='utf-8-sig'))
        spec.update(input_folder=str(data), automatic_geometry=True, exclude_increasing_velocity=True)
        write_json(project/'project.json', spec)
        c = json.loads((project/'config.json').read_text(encoding='utf-8-sig'))
        c.update(sample_rate_hz=fs, bandpass_hz=[0.1, 8.0], window_s=600.0, max_lag_s=120.0,
                 min_clean_segment_s=150.0)
        c['response'].update(pre_filt_hz=[0.03, 0.08, 8.5, 9.5], table_step_hz=0.005)
        # traditional-route HF scan band (2 s = 0.5 Hz ... up toward 5 Hz); start at 0.2 Hz
        c['scan'].update(min_hz=0.2, max_hz=5.0, step_hz=0.02, min_band_hz=0.2, min_points=8,
                         group_velocity_bounds_m_s=[200.0, 5000.0],
                         arrival_min_s=0.2, arrival_max_s=60.0, noise_lag_s=[70.0, 110.0])
        c['phase'].update(velocity_bounds_m_s=[150.0, 5000.0])
        c['cohorts'].update(feature_band_hz=[0.5, 5.0])
        write_json(project/'config.json', c)

    try:
        preprocess(root)   # runs acquisition-audit, import, prepare, scan, controls, families
    except Exception as exc:
        (root/'error.json').write_text(json.dumps({'stage':'preprocess','error':f'{type(exc).__name__}: {exc}','tb':traceback.format_exc()},ensure_ascii=False,indent=2),encoding='utf-8')
        print('ERROR preprocess (error.json)', flush=True); return
    print('HF scan done; output under', project/'output', flush=True)


if __name__ == '__main__':
    main()
