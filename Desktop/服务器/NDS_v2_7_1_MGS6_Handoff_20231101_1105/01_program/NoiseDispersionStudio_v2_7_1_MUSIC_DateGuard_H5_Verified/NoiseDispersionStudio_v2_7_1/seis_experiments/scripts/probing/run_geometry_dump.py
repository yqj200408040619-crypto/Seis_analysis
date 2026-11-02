"""Dump coordinates + full pairwise geometry for the two promising boxes, and try a
grid of (min_spacing, max_spacing, aperture, size) settings to find a compact subset
of >6 broadband vertical stations near a ~20 km aperture. Writes UTF-8 JSON.
"""
import sys, json, itertools
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import ndstudio
import h5data as download
from ndstudio.acquisition import select_array
from ndstudio.array_search_e2e import resolvable_band
from obspy import UTCDateTime
from pyproj import Geod

START, END = '2023-11-01', '2023-11-04'
GEOD = Geod(ellps='WGS84')

TARGETS = [
    ('LongValley_NC', 'IRIS', 'NC', 37.55, 37.80, -119.10, -118.75),
    ('LA_CI',         'IRIS', 'CI', 34.00, 34.25, -118.40, -118.05),
]


def eligible(provider, net, box):
    client = download.make_client(provider)
    rows = download.list_stations(provider, net, station='*', channel='HHZ',
                                  starttime=UTCDateTime(START), endtime=UTCDateTime(END), client=client,
                                  minlatitude=box[0], maxlatitude=box[1],
                                  minlongitude=box[2], maxlongitude=box[3])
    good = {}
    for r in rows:
        if (r.channel.endswith('Z') and r.latitude is not None and r.longitude is not None
                and r.sample_rate and r.sample_rate >= 20
                and not (r.start and r.start > UTCDateTime(START))
                and not (r.end and r.end < UTCDateTime(END))):
            good.setdefault((r.network, r.station), []).append(r)
    unamb = [v[0] for v in good.values() if len(v) == 1]
    return unamb


def pairwise(stations):
    d = []
    for a, b in itertools.combinations(stations, 2):
        _, _, dist = GEOD.inv(a.longitude, a.latitude, b.longitude, b.latitude)
        d.append((a.station, b.station, round(dist)))
    return d


def nn_stats(stations):
    import numpy as np
    n = len(stations)
    mins = []
    for i, a in enumerate(stations):
        dmin = None
        for j, b in enumerate(stations):
            if i == j:
                continue
            _, _, dist = GEOD.inv(a.longitude, a.latitude, b.longitude, b.latitude)
            dmin = dist if dmin is None else min(dmin, dist)
        mins.append(dmin)
    return mins


def main():
    out = {}
    for name, provider, net, *bb in TARGETS:
        box = tuple(bb)
        stns = eligible(provider, net, box)
        coords = [dict(seed=s.seed_id, sta=s.station, lat=s.latitude, lon=s.longitude) for s in stns]
        pw = pairwise(stns)
        aperture_full = max(p[2] for p in pw) if pw else None
        nn = nn_stats(stns)
        entry = dict(count=len(stns), coords=coords,
                     full_aperture_m=aperture_full,
                     nn_min_m=round(min(nn)) if nn else None,
                     nn_max_m=round(max(nn)) if nn else None,
                     pairwise_m=pw)
        # grid of selection attempts
        attempts = []
        for size in (7, 8):
            for maxnn in (10000., 15000., 20000.):
                for ap in (20000., 25000., 30000.):
                    for minnn in (500., 1000., 1500.):
                        if not (minnn <= maxnn <= ap):
                            continue
                        try:
                            res = select_array(provider, net, START, END,
                                               Path(HERE) / 'runs' / 'array_search' / f'_grid_{name}.json',
                                               station='*', channel='HHZ', size=size,
                                               min_stations=7, max_stations=20,
                                               min_spacing_m=minnn, max_spacing_m=maxnn,
                                               max_aperture_m=ap,
                                               minlatitude=box[0], maxlatitude=box[1],
                                               minlongitude=box[2], maxlongitude=box[3])
                            g = res['geometry']
                            band = resolvable_band(g, 1500., 4500.)
                            attempts.append(dict(size=size, minnn=minnn, maxnn=maxnn, aperture_cap=ap,
                                                 ok=True, selected=[s['seed_id'] for s in res['stations']],
                                                 aperture_m=round(g['max_pair_distance_m']),
                                                 min_nn_m=round(min(g['nearest_neighbor_m'].values())),
                                                 max_nn_m=round(max(g['nearest_neighbor_m'].values())),
                                                 band=band))
                        except Exception as exc:
                            attempts.append(dict(size=size, minnn=minnn, maxnn=maxnn, aperture_cap=ap,
                                                 ok=False, error=str(exc)[:80]))
        entry['attempts_ok'] = [a for a in attempts if a.get('ok')]
        out[name] = entry
    p = Path(HERE) / 'runs' / 'array_search' / 'geometry_dump.json'
    p.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
    print('WROTE', p)


if __name__ == '__main__':
    main()
