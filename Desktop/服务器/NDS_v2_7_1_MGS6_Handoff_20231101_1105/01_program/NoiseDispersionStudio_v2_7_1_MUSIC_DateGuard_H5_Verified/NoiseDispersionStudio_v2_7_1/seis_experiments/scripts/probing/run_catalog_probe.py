"""Catalog probe: list stations in a geographic box for a network/channel over an
interval, WITHOUT the compact-subset selection, so we can see what is actually
available before gating. Writes UTF-8 JSON.
"""
import sys, json, argparse
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import ndstudio  # sets up core path
import h5data as download
from obspy import UTCDateTime


def probe(provider, network, station, channel, box, t0, t1):
    client = download.make_client(provider)
    rows = download.list_stations(provider, network, station=station, channel=channel,
                                  starttime=UTCDateTime(t0), endtime=UTCDateTime(t1), client=client,
                                  minlatitude=box[0], maxlatitude=box[1],
                                  minlongitude=box[2], maxlongitude=box[3])
    out = []
    for r in rows:
        out.append(dict(seed_id=r.seed_id, net=r.network, sta=r.station, chan=r.channel,
                        lat=r.latitude, lon=r.longitude, rate=r.sample_rate,
                        start=str(r.start) if r.start else None,
                        end=str(r.end) if r.end else None))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--start', default='2023-11-01')
    ap.add_argument('--end', default='2023-11-04')
    ap.add_argument('--out', default=str(HERE / 'runs' / 'array_search' / 'catalog_probe.json'))
    a = ap.parse_args()
    box = (35.78, 36.00, -120.63, -120.38)
    results = {}
    queries = [
        ('NCEDC', 'NC', '*', 'HHZ'),
        ('NCEDC', 'NC', '*', 'BHZ'),
        ('NCEDC', 'BK', '*', 'HHZ'),
        ('NCEDC', 'BK', '*', 'BHZ'),
        ('NCEDC', 'BP', '*', 'HHZ'),   # Parkfield HRSN borehole array
        ('NCEDC', 'BP', '*', 'DPZ'),
        ('NCEDC', '*', '*', 'HHZ'),
    ]
    for provider, net, sta, chan in queries:
        key = f'{provider}:{net}:{chan}'
        try:
            rows = probe(provider, net, sta, chan, box, a.start, a.end)
            results[key] = dict(count=len(rows), stations=rows)
        except Exception as exc:
            results[key] = dict(error=f'{type(exc).__name__}: {exc}')
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding='utf-8')
    print('WROTE', a.out)


if __name__ == '__main__':
    main()
