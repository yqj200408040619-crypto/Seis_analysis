"""Extended-aperture direction experiment (user design):
 - dense 20-30 km core (LA basin CI) for velocity refinement
 - a few SPARSE extenders (different stations, same CI net via SCEDC) to widen the
   overall aperture and sharpen direction resolution near 0.2 Hz (5 s).

Builds a selection JSON in the program's own schema using the program's geometry
helper (no hand-rolled geometry), fetches via SCEDC, and runs the MUSIC/front
direction stage in DIAGNOSTIC mode so per-frequency direction estimates and array
response are produced even if the strict filter/joint gates reject points.
"""
import sys, json, argparse
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import ndstudio
import h5data as download
from dispersion.common import write_json
from ndstudio.acquisition import fetch_array, audit_acquisition, screen_array
from ndstudio.project import new_project
from ndstudio.real_compare import preprocess
from ndstudio.pipeline import execute
from obspy import UTCDateTime
import dataclasses

import os
PROVIDER = 'SCEDC'
NETWORK = 'CI'
START = os.environ.get('EXT_START', '2023-11-01')
END = os.environ.get('EXT_END', '2023-11-04')

# dense core (8) + sparse azimuthal extenders (4). Prefer HHZ; fall back to BHZ.
CORE = ['BHP', 'GR2', 'HLL', 'LCG', 'LGB', 'PDR', 'SMF2', 'USC']
EXTENDERS = ['MWC', 'RPV', 'RIO']   # NE, S, E  (~30-35 km out; WSS dropped: mixed native grid)
ALL_STATIONS = CORE if os.environ.get('EXT_CORE_ONLY') == '1' else CORE + EXTENDERS


def pick_rows(stations):
    """One vertical channel per station, prefer HHZ then BHZ, full-coverage only."""
    cli = download.make_client(PROVIDER)
    rows = download.list_stations(PROVIDER, NETWORK, station=','.join(stations),
                                  channel='HHZ,BHZ', starttime=UTCDateTime(START),
                                  endtime=UTCDateTime(END), client=cli)
    by_sta = {}
    for r in rows:
        if not r.channel.endswith('Z'):
            continue
        if r.latitude is None or r.longitude is None:
            continue
        if r.sample_rate is None or r.sample_rate < 20:
            continue
        if (r.start and r.start > UTCDateTime(START)) or (r.end and r.end < UTCDateTime(END)):
            continue
        by_sta.setdefault(r.station, []).append(r)
    chosen = []
    for s in stations:
        cands = by_sta.get(s, [])
        if not cands:
            continue
        cands.sort(key=lambda r: (0 if r.channel == 'HHZ' else 1, -(r.sample_rate or 0)))
        chosen.append(cands[0])
    return chosen


def build_selection(out_path):
    rows = pick_rows(ALL_STATIONS)
    rows.sort(key=lambda r: r.seed_id)
    g = download.array_geometry(rows)
    def ser(r):
        return {k: (str(v) if isinstance(v, UTCDateTime) else v) for k, v in dataclasses.asdict(r).items()}
    result = {
        'schema': 'ndstudio.array_selection/1', 'provider': PROVIDER, 'network': NETWORK,
        'start_utc': str(UTCDateTime(START)), 'end_utc': str(UTCDateTime(END)),
        'anchor_seed_id': rows[0].seed_id,
        'constraints': {'station_count': len(rows), 'minimum_station_count': 7,
                        'maximum_station_count': 50,
                        'minimum_nearest_neighbor_m': 0., 'maximum_nearest_neighbor_m': 1e9,
                        'maximum_pair_distance_m': 1e9},
        'metadata_rows_queried': len(rows), 'eligible_stations': len(rows),
        'station_pattern': ','.join(ALL_STATIONS), 'station_filters': {},
        'matching_station_catalog': download.station_catalog(rows),
        'channel_filter': 'HHZ,BHZ', 'rejected_metadata': [], 'geometry': g,
        'stations': [ser(r) for r in rows],
        'selection_uses_waveforms_or_previous_curves': False,
    }
    write_json(out_path, result)
    return result, g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workdir', default=str(HERE/'runs'/'extended_la'))
    ap.add_argument('--reuse-download', action='store_true')
    a = ap.parse_args()
    root = Path(a.workdir).resolve(); root.mkdir(parents=True, exist_ok=True)
    data = root/'download'; project = root/'project'
    selected = root/'array_selection.json'

    result, g = build_selection(selected)
    summary = {'stations': [s['seed_ids'] if 'seed_ids' in s else s for s in result['matching_station_catalog']],
               'n': len(result['stations']),
               'aperture_m': g['max_pair_distance_m'],
               'min_pair_m': g['min_pair_distance_m'],
               'mean_nn_m': g['mean_nearest_neighbor_m']}
    write_json(root/'extended_geometry.json', summary)
    print('SELECTION built: n=%d aperture=%.0f m' % (summary['n'], summary['aperture_m']), flush=True)

    if not a.reuse_download:
        fetch_array(selected, data, workers=4, chunk_seconds=21600., allow_missing_data=True)
    import traceback as _tb
    try:
        audit_acquisition(data)
        screen_array(selected, data, project=None)
    except Exception as exc:
        (root/'error.json').write_text(json.dumps(
            {'stage':'audit/screen','error':f'{type(exc).__name__}: {exc}','tb':_tb.format_exc()},
            ensure_ascii=False, indent=2), encoding='utf-8')
        print('ERROR at audit/screen (see error.json)', flush=True); return

    # project: band around 0.2 Hz (5 s), broad enough for direction diagnostics
    f_min, f_max = 0.08, 0.40
    if not (project/'project.json').exists():
        new_project(project)
        spec = json.loads((project/'project.json').read_text(encoding='utf-8-sig'))
        spec.update(input_folder=str(data), automatic_geometry=True, exclude_increasing_velocity=True)
        write_json(project/'project.json', spec)
        c = json.loads((project/'config.json').read_text(encoding='utf-8-sig'))
        c.update(sample_rate_hz=10., bandpass_hz=[0.04, 0.8], max_lag_s=400.)
        c['response'].update(pre_filt_hz=[0.01, 0.02, 3.0, 4.0], table_step_hz=0.002)
        c['scan'].update(min_hz=f_min, max_hz=f_max, step_hz=0.005, min_band_hz=0.03,
                         group_velocity_bounds_m_s=[1500., 4500.])
        c['phase'].update(velocity_bounds_m_s=[1000., 4500.])
        write_json(project/'config.json', c)

    # Authorized relaxation of the SCAN validation gate (route SNR, block wave-packet,
    # arrival). Applied to the existing config even on reuse. These are NON-DEFAULT,
    # NON-CALIBRATED thresholds set deliberately with the user's permission.
    if os.environ.get('EXT_RELAX') == '1':
        c = json.loads((project/'config.json').read_text(encoding='utf-8-sig'))
        v = c.setdefault('validation', {})
        v.update(
            snr_min=float(os.environ.get('EXT_SNR_MIN', '3.0')),
            block_snr_min=float(os.environ.get('EXT_BLOCK_SNR_MIN', '3.0')),
            block_fraction_min=float(os.environ.get('EXT_BLOCK_FRAC_MIN', '0.5')),
            block_phase_max_rad=float(os.environ.get('EXT_BLOCK_PHASE_MAX', '0.7')),
            # arrival kept reasonable: only modest loosening
            arrival_fraction_max=float(os.environ.get('EXT_ARR_FRAC_MAX', '0.2')),
            leaveout_arrival_fraction_max=float(os.environ.get('EXT_LO_ARR_FRAC_MAX', '0.2')),
            block_arrival_fraction_max=float(os.environ.get('EXT_BLK_ARR_FRAC_MAX', '0.25')),
        )
        write_json(project/'config.json', c)
        write_json(root/'relaxed_validation_used.json', v)

    # DIAGNOSTIC direction mode: outputs per-frequency direction + array response,
    # without the strict filter/joint accept gate masking everything.
    _cycles = float(os.environ.get('EXT_FRONT_CYCLES', '200'))
    _minwin = int(os.environ.get('EXT_MIN_WINDOWS', '3'))
    direction = dict(mode='diagnostic', max_pair_angle_deg=45., direction_step_deg=2.,
                     velocity_points=128, velocity_min_m_s=1000., velocity_max_m_s=4500.,
                     frequency_half_width_hz=0., frequency_bins=1,
                     music_front_cycles=_cycles,
                     minimum_independent_windows=_minwin,
                     minimum_accepted_windows=_minwin)
    write_json(project/'phase_settings.json', {'direction': direction})

    try:
        preprocess(root)
        execute(project, 'music')
    except Exception as exc:
        import traceback as _tb2
        (root/'error.json').write_text(json.dumps(
            {'stage':'preprocess/music','error':f'{type(exc).__name__}: {exc}','tb':_tb2.format_exc()},
            ensure_ascii=False, indent=2), encoding='utf-8')
        print('ERROR at preprocess/music (see error.json)', flush=True); return
    print('MUSIC diagnostic done; see', project/'output'/'music_front', flush=True)


if __name__ == '__main__':
    main()
