"""Diagnose the Long Valley compact-subset selection: reproduce select_array for a
range of bounds and print the geometry dict and full error for the best one.
"""
import sys, json, traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import ndstudio
import h5data as download
from obspy import UTCDateTime
from pyproj import Geod

START, END = '2023-11-01', '2023-11-04'
BOX = (37.55, 37.80, -119.10, -118.75)
GEOD = Geod(ellps='WGS84')


def eligible():
    client = download.make_client('IRIS')
    rows = download.list_stations('IRIS', 'NC', station='*', channel='HHZ',
                                  starttime=UTCDateTime(START), endtime=UTCDateTime(END), client=client,
                                  minlatitude=BOX[0], maxlatitude=BOX[1],
                                  minlongitude=BOX[2], maxlongitude=BOX[3])
    good = {}
    for r in rows:
        if (r.channel.endswith('Z') and r.latitude is not None and r.longitude is not None
                and r.sample_rate and r.sample_rate >= 20
                and not (r.start and r.start > UTCDateTime(START))
                and not (r.end and r.end < UTCDateTime(END))):
            good.setdefault((r.network, r.station), []).append(r)
    return [v[0] for v in good.values() if len(v) == 1]


def geom_of(subset):
    g = download.array_geometry(subset)
    return g


def main():
    stns = eligible()
    out = {'n': len(stns)}
    # For each center, list the 7 and 8 nearest and their geometry
    reports = []
    for center in stns:
        ranked = sorted(stns, key=lambda r: (
            abs(GEOD.inv(center.longitude, center.latitude, r.longitude, r.latitude)[2]), r.seed_id))
        for size in (7, 8):
            subset = ranked[:size]
            g = geom_of(subset)
            nn = list(g['nearest_neighbor_m'].values())
            reports.append(dict(center=center.station, size=size,
                                members=[s.station for s in subset],
                                min_nn=round(min(nn)), max_nn=round(max(nn)),
                                aperture=round(g['max_pair_distance_m']),
                                coord_conflicts=g.get('coordinate_conflicts')))
    out['subsets'] = reports
    p = HERE / 'runs' / 'array_search' / 'longvalley_diag.json'
    p.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
    print('WROTE', p)


if __name__ == '__main__':
    main()
