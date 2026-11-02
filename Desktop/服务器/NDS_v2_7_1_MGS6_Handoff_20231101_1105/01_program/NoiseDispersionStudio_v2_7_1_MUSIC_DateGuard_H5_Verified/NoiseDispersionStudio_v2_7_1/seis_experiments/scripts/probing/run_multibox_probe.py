"""Scan several candidate dense-array boxes against LIVE FDSN metadata and report,
for each, how many broadband vertical channels fully cover the interval and the
compact-subset geometry (min/max nearest-neighbor, aperture) the program would pick.

Goal: find a region with >6 broadband vertical stations inside a ~20 km aperture,
every nearest-neighbor within [min_spacing, max_spacing], aperture <= max_aperture.

Writes UTF-8 JSON.
"""
import sys, json, argparse, traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import ndstudio  # sets core path
import h5data as download
from ndstudio.acquisition import select_array
from ndstudio.array_search_e2e import resolvable_band
from obspy import UTCDateTime

# (name, provider, network, station, channel, box lat0,lat1,lon0,lon1)
# Boxes chosen around classic dense broadband arrays with ambient-noise literature.
BOXES = [
    # Anza, S. California (CI network) - dense broadband, ~tens of km
    ('Anza_CI',        'IRIS', 'CI', '*', 'HHZ', 33.3, 33.65, -116.80, -116.35),
    # Long Valley caldera (NC) broadband
    ('LongValley_NC',  'IRIS', 'NC', '*', 'HHZ', 37.55, 37.80, -119.10, -118.75),
    # Yellowstone (WY/US)
    ('Yellowstone_WY', 'IRIS', 'WY', '*', 'HHZ', 44.35, 44.70, -110.95, -110.40),
    # Imperial Valley / Salton (CI)
    ('Salton_CI',      'IRIS', 'CI', '*', 'HHZ', 33.00, 33.35, -115.80, -115.45),
    # Swiss Valais dense (CH)
    ('Valais_CH',      'ETH',  'CH', '*', 'HHZ', 46.10, 46.45, 7.30, 7.75),
    # Northern Italy / INGV dense around Po plain or Alps
    ('CentralItaly_IV','INGV', 'IV', '*', 'HHZ', 42.55, 42.90, 13.00, 13.45),
    # Iceland (GEOFON/various) - skip, provider mix
    # Utah (UU) Wasatch
    ('Wasatch_UU',     'IRIS', 'UU', '*', 'HHZ', 40.55, 40.90, -111.95, -111.55),
    # Southern California LB (CI) metro dense
    ('LA_CI',          'IRIS', 'CI', '*', 'HHZ', 34.00, 34.25, -118.40, -118.05),
]

GATE = dict(size=7, min_stations=7, max_stations=20,
            min_spacing_m=1500., max_spacing_m=12000., max_aperture_m=25000.)
VMIN, VMAX = 1500., 4500.


def count_bbz(provider, net, sta, chan, box, t0, t1):
    client = download.make_client(provider)
    rows = download.list_stations(provider, net, station=sta, channel=chan,
                                  starttime=UTCDateTime(t0), endtime=UTCDateTime(t1), client=client,
                                  minlatitude=box[0], maxlatitude=box[1],
                                  minlongitude=box[2], maxlongitude=box[3])
    # count stations whose vertical channel fully covers interval and rate>=20
    good = {}
    for r in rows:
        if not r.channel.endswith('Z'):
            continue
        if r.latitude is None or r.longitude is None:
            continue
        if r.sample_rate is None or r.sample_rate < 20:
            continue
        if (r.start is not None and r.start > UTCDateTime(t0)) or (r.end is not None and r.end < UTCDateTime(t1)):
            continue
        good.setdefault((r.network, r.station), []).append(r.seed_id)
    unambiguous = {k: v for k, v in good.items() if len(v) == 1}
    return len(rows), len(good), len(unambiguous), sorted('.'.join(k) for k in unambiguous)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--start', default='2023-11-01')
    ap.add_argument('--end', default='2023-11-04')
    ap.add_argument('--out', default=str(HERE / 'runs' / 'array_search' / 'multibox_probe.json'))
    a = ap.parse_args()
    results = []
    for name, provider, net, sta, chan, *bb in BOXES:
        box = tuple(bb)
        entry = dict(name=name, provider=provider, network=net, channel=chan, box=box)
        try:
            nrows, nsta, nun, names = count_bbz(provider, net, sta, chan, box, a.start, a.end)
            entry.update(metadata_rows=nrows, stations_with_vertical=nsta,
                         unambiguous_full_coverage=nun, station_names=names)
            # Try the actual compact-subset selection if there are enough
            if nun >= GATE['size']:
                sel_path = Path(a.out).parent / f'_sel_{name}.json'
                try:
                    res = select_array(provider, net, a.start, a.end, sel_path,
                                       station=sta, channel=chan,
                                       size=GATE['size'], min_stations=GATE['min_stations'],
                                       max_stations=GATE['max_stations'],
                                       min_spacing_m=GATE['min_spacing_m'],
                                       max_spacing_m=GATE['max_spacing_m'],
                                       max_aperture_m=GATE['max_aperture_m'],
                                       minlatitude=box[0], maxlatitude=box[1],
                                       minlongitude=box[2], maxlongitude=box[3])
                    g = res['geometry']
                    band = resolvable_band(g, VMIN, VMAX)
                    entry['selection'] = dict(
                        ok=True, selected=[s['seed_id'] for s in res['stations']],
                        aperture_m=g['max_pair_distance_m'],
                        min_nn_m=min(g['nearest_neighbor_m'].values()),
                        max_nn_m=max(g['nearest_neighbor_m'].values()),
                        band=band)
                except Exception as exc:
                    entry['selection'] = dict(ok=False, error=f'{type(exc).__name__}: {exc}')
        except Exception as exc:
            entry.update(error=f'{type(exc).__name__}: {exc}',
                         traceback=traceback.format_exc().splitlines()[-3:])
        results.append(entry)
        # incremental write so partial progress survives
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        Path(a.out).write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding='utf-8')
        print('done', name, entry.get('unambiguous_full_coverage'), flush=True)
    print('WROTE', a.out)


if __name__ == '__main__':
    main()
