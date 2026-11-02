"""Find sparse broadband vertical stations in a WIDE box around the dense LA-basin
subarray, across whatever networks SCEDC/IRIS serve, to EXTEND the overall aperture
for better direction resolution near 0.2 Hz (5 s). We then (separately) keep the dense
20-30 km subarray for velocity refinement.

Strategy:
 - dense core (verified): CI BHP,GR2,HLL,LCG,LGB,PDR,SMF2,USC (~29.6 km aperture).
 - look OUTSIDE the core box, 30-90 km from the core centroid, in several azimuths,
   for a few broadband Z stations (any net) that cover the interval and download.
Writes UTF-8 JSON (metadata only; cheap).
"""
import sys, json, itertools
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import ndstudio
import h5data as download
from obspy import UTCDateTime
from pyproj import Geod

START, END = '2023-11-01', '2023-11-04'
GEOD = Geod(ellps='WGS84')

# Wide box around LA basin core (core ~34.0-34.3N, -118.45..-117.95W).
WIDE_BOX = (33.4, 34.8, -119.1, -117.1)

# Networks that may host broadband Z around S. California, served by SCEDC and/or IRIS.
QUERIES = [
    ('SCEDC', 'CI'),
    ('IRIS',  'CI'),
    ('IRIS',  'II'),   # global GSN
    ('IRIS',  'IU'),   # global GSN
    ('IRIS',  'US'),   # US backbone
    ('IRIS',  'BK'),   # Berkeley (N. Cal but some S)
    ('IRIS',  'NP'),   # NSMP
    ('IRIS',  'CE'),   # CGS strong motion (may be BB)
    ('IRIS',  'AZ'),   # Anza (SE)
    ('IRIS',  'PB'),   # PBO/NOTA borehole
]


def core_centroid():
    cli = download.make_client('SCEDC')
    rows = download.list_stations('SCEDC', 'CI', station='*', channel='HHZ',
                                  starttime=UTCDateTime(START), endtime=UTCDateTime(END), client=cli,
                                  minlatitude=33.95, maxlatitude=34.30,
                                  minlongitude=-118.45, maxlongitude=-117.95)
    pts = {}
    for r in rows:
        if r.channel.endswith('Z') and r.latitude and r.longitude:
            pts.setdefault(r.station, (r.latitude, r.longitude))
    core = {k: pts[k] for k in ('BHP','GR2','HLL','LCG','LGB','PDR','SMF2','USC') if k in pts}
    lat = sum(v[0] for v in core.values())/len(core)
    lon = sum(v[1] for v in core.values())/len(core)
    return (lat, lon), core


def main():
    (clat, clon), core = core_centroid()
    out = {'core_centroid': [clat, clon], 'core_stations': core, 'networks': {}}
    for provider, net in QUERIES:
        key = f'{provider}:{net}'
        try:
            cli = download.make_client(provider)
            rows = download.list_stations(provider, net, station='*', channel='?HZ',
                                          starttime=UTCDateTime(START), endtime=UTCDateTime(END), client=cli,
                                          minlatitude=WIDE_BOX[0], maxlatitude=WIDE_BOX[1],
                                          minlongitude=WIDE_BOX[2], maxlongitude=WIDE_BOX[3])
            cands = []
            seen = set()
            for r in rows:
                if not r.channel.endswith('Z'):
                    continue
                if r.latitude is None or r.longitude is None:
                    continue
                if r.sample_rate is None or r.sample_rate < 20:
                    continue
                if (r.start and r.start > UTCDateTime(START)) or (r.end and r.end < UTCDateTime(END)):
                    continue
                k = (r.network, r.station)
                if k in seen:
                    continue
                seen.add(k)
                _, _, dist = GEOD.inv(clon, clat, r.longitude, r.latitude)
                # want SPARSE extenders 25-90 km from core centroid
                if 25000 <= dist <= 90000:
                    az, _, _ = GEOD.inv(clon, clat, r.longitude, r.latitude)
                    cands.append(dict(seed=r.seed_id, sta=r.station, chan=r.channel,
                                      lat=r.latitude, lon=r.longitude, rate=r.sample_rate,
                                      dist_km=round(dist/1000, 1), az_deg=round(az, 1)))
            cands.sort(key=lambda c: c['dist_km'])
            out['networks'][key] = dict(count=len(cands), stations=cands)
        except Exception as exc:
            out['networks'][key] = dict(error=f'{type(exc).__name__}: {exc}'[:160])
    Path(HERE/'runs'/'extend_probe.json').write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
    print('WROTE extend_probe.json')


if __name__ == '__main__':
    main()
