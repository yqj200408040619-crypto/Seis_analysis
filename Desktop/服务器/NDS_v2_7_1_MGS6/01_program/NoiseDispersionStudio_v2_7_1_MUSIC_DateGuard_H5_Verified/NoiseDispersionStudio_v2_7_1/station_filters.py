"""Station-name and geographic filters shared by browse, preview and acquisition.

Bounds include their endpoints. West > east means a box across the dateline.
No selection limit is imposed here: every matching channel epoch is retained.
"""
import csv
import fnmatch
import math
import re
from collections import defaultdict


def normalize_station_filters(*, station_prefix='', minlatitude=None, maxlatitude=None,
                              minlongitude=None, maxlongitude=None):
    prefixes=[]
    for raw in str(station_prefix or '').split(','):
        value=raw.strip().upper()
        if not value: continue
        if value.endswith('*'): value=value[:-1]
        if not value or not re.fullmatch(r'[A-Z0-9_-]+',value):
            raise ValueError('Station prefix: use initial letters/digits, e.g. QS or QS,QT (no internal wildcards).')
        if value not in prefixes: prefixes.append(value)
    result={'station_prefix':','.join(prefixes)}
    for key,value,limit in [('minlatitude',minlatitude,90),('maxlatitude',maxlatitude,90),
                            ('minlongitude',minlongitude,180),('maxlongitude',maxlongitude,180)]:
        if value is None or str(value).strip()=='': result[key]=None; continue
        try: number=float(value)
        except (TypeError,ValueError): raise ValueError(f'{key} must be a number') from None
        if not math.isfinite(number) or not -limit<=number<=limit:
            raise ValueError(f'{key} must lie in {-limit}..{limit}')
        result[key]=number
    if result['minlatitude'] is not None and result['maxlatitude'] is not None and result['minlatitude']>result['maxlatitude']:
        raise ValueError('Minimum latitude must not exceed maximum latitude')
    return result


def geographic_queries(filters):
    bounds={key:value for key,value in filters.items() if key!='station_prefix' and value is not None}
    west,east=filters['minlongitude'],filters['maxlongitude']
    if west is not None and east is not None and west>east:
        return [dict(bounds,maxlongitude=180.),dict(bounds,minlongitude=-180.)]
    return [bounds]


def station_query_pattern(station,filters):
    """Narrow a broad server query, then use the local intersection as authority."""
    pattern=str(station or '*').strip().upper() or '*'
    prefixes=filters['station_prefix'].split(',') if filters['station_prefix'] else []
    return ','.join(p+'*' for p in prefixes) if pattern=='*' and prefixes else pattern


def filter_station_rows(rows, *, station='*', **options):
    filters=normalize_station_filters(**options)
    prefixes=tuple(filters['station_prefix'].split(',')) if filters['station_prefix'] else ()
    patterns=[x.strip().upper() for x in str(station or '*').split(',') if x.strip()] or ['*']
    geographic=any(filters[key] is not None for key in filters if key!='station_prefix')
    result=[]
    for row in rows:
        name=row.station.upper()
        if not any(fnmatch.fnmatchcase(name,p) for p in patterns): continue
        if prefixes and not name.startswith(prefixes): continue
        if geographic:
            lat,lon=row.latitude,row.longitude
            if lat is None or lon is None: continue
            if not math.isfinite(lat) or not math.isfinite(lon) or not -90<=lat<=90 or not -180<=lon<=180: continue
            if filters['minlatitude'] is not None and lat<filters['minlatitude']: continue
            if filters['maxlatitude'] is not None and lat>filters['maxlatitude']: continue
            west,east=filters['minlongitude'],filters['maxlongitude']
            if west is not None and east is not None and west>east:
                if not (lon>=west or lon<=east): continue
            else:
                if west is not None and lon<west: continue
                if east is not None and lon>east: continue
        result.append(row)
    return result


def unique_channel_rows(rows):
    result=[]; seen=set()
    for row in rows:
        key=repr(row)
        if key not in seen: seen.add(key); result.append(row)
    return result


def station_catalog(rows):
    """One record per physical NET.STA; retain channels/epochs and coordinate conflicts."""
    grouped=defaultdict(list)
    for row in rows: grouped[(row.network,row.station)].append(row)
    result=[]
    for (network,station),group in sorted(grouped.items()):
        coordinates=sorted({(r.latitude,r.longitude) for r in group if r.latitude is not None and r.longitude is not None
                            and math.isfinite(r.latitude) and math.isfinite(r.longitude)})
        lat,lon=coordinates[0] if coordinates else (None,None)
        starts=[r.start for r in group if r.start is not None]
        ends=[r.end for r in group if r.end is not None]
        result.append(dict(network=network,station=station,latitude=lat,longitude=lon,
            channels=','.join(sorted({r.channel for r in group})),
            seed_ids=','.join(sorted({r.seed_id for r in group})),
            sample_rates_hz=','.join(str(x) for x in sorted({r.sample_rate for r in group if r.sample_rate is not None})),
            channel_epoch_count=len(group),coordinate_conflict=len(coordinates)>1,
            start_utc=str(min(starts)) if starts else '',
            end_utc='' if any(r.end is None for r in group) or not ends else str(max(ends))))
    return result


def export_station_catalog(rows,filename):
    records=station_catalog(rows)
    columns=['network','station','latitude','longitude','channels','seed_ids','sample_rates_hz',
             'channel_epoch_count','coordinate_conflict','start_utc','end_utc']
    with open(filename,'w',newline='',encoding='utf-8-sig') as handle:
        writer=csv.DictWriter(handle,fieldnames=columns);writer.writeheader();writer.writerows(records)
    return records
