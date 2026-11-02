"""Scan NEW candidate dense-array regions against LIVE FDSN metadata, targeting
>=8 broadband vertical stations within a 20-30 km aperture (< 30 km). For each box
we count unambiguous full-coverage vertical channels and, if >=8, try the program's
compact-subset select_array at size 8/9/10 with aperture cap 30 km. Writes UTF-8 JSON.
"""
import sys, json, traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import ndstudio
import h5data as download
from ndstudio.acquisition import select_array
from ndstudio.array_search_e2e import resolvable_band
from obspy import UTCDateTime

START, END = '2023-11-01', '2023-11-04'

# (name, provider, network, channel, lat0,lat1,lon0,lon1)
# NOTE: provider = where WAVEFORMS live. NC->NCEDC, CI/AZ/others->IRIS/EARTHSCOPE,
# CH->ETH, IV->INGV, GE/various->GFZ. We test metadata via the same provider.
BOXES = [
    # San Jacinto / Anza wider (CI + AZ dense) - classic dense broadband
    ('SanJacinto_CI', 'IRIS', 'CI', 'HHZ', 33.30, 33.75, -116.90, -116.30),
    ('Anza_AZ',       'IRIS', 'AZ', 'HHZ', 33.20, 33.75, -116.90, -116.30),
    # Coso geothermal (CI) eastern CA
    ('Coso_CI',       'IRIS', 'CI', 'HHZ', 35.90, 36.25, -117.95, -117.55),
    # Salton Sea / Imperial wider (CI)
    ('SaltonWide_CI', 'IRIS', 'CI', 'HHZ', 32.90, 33.45, -115.95, -115.35),
    # Mammoth / Long Valley WIDER box (NCEDC waveforms) -> include MBM, 8 stns ~23km
    ('LongValleyWide_NC','NCEDC','NC','HHZ', 37.50, 37.90, -119.15, -118.70),
    # Yellowstone wider: WY + US + IW federated on EarthScope
    ('YellowstoneWide_WY','IRIS','WY','HHZ', 44.30, 44.85, -111.10, -110.30),
    # Utah Wasatch wider (UU)
    ('WasatchWide_UU','IRIS','UU','HHZ', 40.40, 40.95, -112.10, -111.45),
    # Swiss plateau / Wallis wider (CH, ETH)
    ('SwissWide_CH',  'ETH',  'CH', 'HHZ', 46.00, 46.55, 7.10, 7.95),
    # Central Apennines wider (IV, INGV)
    ('ApenninesWide_IV','INGV','IV','HHZ', 42.40, 42.95, 12.90, 13.60),
    # Southern California Los Angeles wider box (CI) aiming 20-30km subset
    ('LAbasin_CI',    'IRIS', 'CI', 'HHZ', 33.95, 34.30, -118.45, -117.95),
    # Hawaii Big Island (HV, IRIS) volcanic dense
    ('Kilauea_HV',    'IRIS', 'HV', 'HHZ', 19.25, 19.60, -155.45, -155.05),
    # Mount St Helens / Cascades (UW/CC)
    ('Cascades_UW',   'IRIS', 'UW', 'HHZ', 46.00, 46.40, -122.40, -121.90),
]

VMIN, VMAX = 1500., 4500.


def eligible(provider, net, chan, box):
    client = download.make_client(provider)
    rows = download.list_stations(provider, net, station='*', channel=chan,
                                  starttime=UTCDateTime(START), endtime=UTCDateTime(END), client=client,
                                  minlatitude=box[0], maxlatitude=box[1],
                                  minlongitude=box[2], maxlongitude=box[3])
    good = {}
    for r in rows:
        if (r.channel.endswith('Z') and r.latitude is not None and r.longitude is not None
                and r.sample_rate and r.sample_rate >= 20
                and not (r.start and r.start > UTCDateTime(START))
                and not (r.end and r.end < UTCDateTime(END))):
            good.setdefault((r.network, r.station), []).append(r.seed_id)
    unamb = {k: v for k, v in good.items() if len(v) == 1}
    return len(rows), len(good), sorted('.'.join(k) for k in unamb), len(unamb)


def main():
    out = []
    for name, provider, net, chan, *bb in BOXES:
        box = tuple(bb)
        entry = dict(name=name, provider=provider, network=net, box=box)
        try:
            nrows, nsta, names, nun = eligible(provider, net, chan, box)
            entry.update(metadata_rows=nrows, unambiguous_full_coverage=nun, stations=names)
            if nun >= 8:
                best = None
                for size in (8, 9, 10):
                    if size > nun:
                        break
                    for ap in (25000., 30000.):
                        for maxnn in (12000., 18000.):
                            for minnn in (500., 1000.):
                                try:
                                    res = select_array(provider, net, START, END,
                                        HERE/'runs'/'array_search'/f'_newsel_{name}.json',
                                        station='*', channel=chan, size=size,
                                        min_stations=7, max_stations=25,
                                        min_spacing_m=minnn, max_spacing_m=maxnn,
                                        max_aperture_m=ap,
                                        minlatitude=box[0], maxlatitude=box[1],
                                        minlongitude=box[2], maxlongitude=box[3])
                                    g = res['geometry']
                                    band = resolvable_band(g, VMIN, VMAX)
                                    cand = dict(size=size, minnn=minnn, maxnn=maxnn, ap_cap=ap,
                                                selected=[s['station'] for s in res['stations']],
                                                aperture_m=round(g['max_pair_distance_m']),
                                                min_nn_m=round(min(g['nearest_neighbor_m'].values())),
                                                max_nn_m=round(max(g['nearest_neighbor_m'].values())),
                                                mean_nn_m=round(g['mean_nearest_neighbor_m']),
                                                band=band)
                                    # prefer largest size, then aperture closest to ~25km but <30km
                                    key = (size, -abs(cand['aperture_m']-25000))
                                    if best is None or key > best[0]:
                                        best = (key, cand)
                                except Exception:
                                    pass
                if best:
                    entry['best_selection'] = best[1]
                else:
                    entry['best_selection'] = None
        except Exception as exc:
            entry.update(error=f'{type(exc).__name__}: {exc}')
        out.append(entry)
        Path(HERE/'runs'/'array_search'/'newregion_probe.json').write_text(
            json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
        print('done', name, entry.get('unambiguous_full_coverage'),
              (entry.get('best_selection') or {}).get('aperture_m') if entry.get('best_selection') else None,
              flush=True)
    print('WROTE newregion_probe.json')


if __name__ == '__main__':
    main()
