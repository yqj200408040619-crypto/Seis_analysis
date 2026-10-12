"""Program-owned FDSN array selection and lossless, resumable H5 acquisition."""
from pathlib import Path
from dataclasses import asdict
import json
from contextlib import contextmanager
import numpy as np
from obspy import UTCDateTime
from pyproj import Geod
import h5data as download
from h5data_handoff import json_value, h5_is_complete, write_run_manifest
from dispersion.common import write_json


def select_array(provider, network, start, end, output, *, size=10,
                 min_spacing_m=200., max_spacing_m=2000., max_aperture_m=6000.,
                 anchor=None, channel='HHZ', client=None,min_stations=6,max_stations=49,
                 station='*',station_prefix='',minlatitude=None,maxlatitude=None,minlongitude=None,maxlongitude=None):
    """Select one native vertical channel per station, before looking at signals.

    Station count is strictly 5<n<50. Bounds apply to EVERY nearest-neighbor
    distance, rather than accepting a deceptively small network-wide average.
    The maximum pair distance bounds the aperture. All channel epochs must
    cover the complete requested interval; larger networks are reduced to a
    deterministic compact subset, without using old dispersion observations.
    """
    if not 5<=min_stations<=max_stations<=100:raise ValueError('Station-count limits must lie in 5–100')
    if not min_stations <= size <= max_stations:
        raise ValueError('Array size must satisfy 5 < n < 50' if (min_stations,max_stations)==(6,49)
                         else f'Array size must lie in {min_stations}–{max_stations}')
    if not 0 < min_spacing_m <= max_spacing_m <= max_aperture_m:
        raise ValueError('Require 0 < minimum spacing <= maximum spacing <= aperture')
    t0, t1 = UTCDateTime(start), UTCDateTime(end)
    if t1 <= t0:
        raise ValueError('End UTC must follow start UTC')
    client = client or download.make_client(provider)
    rows = download.list_stations(provider, network, station=station,channel=channel,
                                  starttime=t0, endtime=t1, client=client,
                                  station_prefix=station_prefix,minlatitude=minlatitude,maxlatitude=maxlatitude,
                                  minlongitude=minlongitude,maxlongitude=maxlongitude)
    station_filters=download.normalize_station_filters(station_prefix=station_prefix,minlatitude=minlatitude,
                    maxlatitude=maxlatitude,minlongitude=minlongitude,maxlongitude=maxlongitude)
    rows=download.filter_station_rows(rows,station=station,**station_filters)
    candidates, rejected = {}, []
    for row in rows:
        reason = None
        if not row.channel.endswith('Z'):
            reason = 'not_vertical'
        elif row.latitude is None or row.longitude is None:
            reason = 'missing_coordinates'
        elif not np.isfinite([row.latitude, row.longitude]).all():
            reason = 'invalid_coordinates'
        elif row.sample_rate is None or row.sample_rate < 20:
            reason = 'native_rate_below_analysis_rate'
        elif (row.start is not None and row.start > t0) or (row.end is not None and row.end < t1):
            reason = 'epoch_does_not_cover_full_interval'
        if reason:
            rejected.append({'seed_id': row.seed_id, 'reason': reason})
        else:
            candidates.setdefault((row.network, row.station), []).append(row)
    usable = []
    for key, channels in sorted(candidates.items()):
        if len(channels) != 1:
            rejected.append({'station': '.'.join(key), 'reason': 'ambiguous_channel_or_epoch',
                             'seed_ids': [r.seed_id for r in channels]})
        else:
            usable.append(channels[0])
    if len(usable) < size:
        raise ValueError(f'Only {len(usable)} unambiguous eligible stations; require {size}')
    geod = Geod(ellps='WGS84')
    choices = []
    for center in usable:
        if anchor and center.station != anchor and center.seed_id != anchor:
            continue
        ranked = sorted(usable, key=lambda r: (
            abs(geod.inv(center.longitude, center.latitude, r.longitude, r.latitude)[2]), r.seed_id))
        subset = ranked[:size]
        g = download.array_geometry(subset)
        nearest = list(g['nearest_neighbor_m'].values())
        valid = (len(nearest) == size and not g['coordinate_conflicts'] and
                 min(nearest) >= min_spacing_m and max(nearest) <= max_spacing_m and
                 g['max_pair_distance_m'] <= max_aperture_m)
        if valid:
            choices.append((g['max_pair_distance_m'], g['mean_nearest_neighbor_m'],
                            center.seed_id, subset, g))
    if not choices:
        raise ValueError('No compact subset satisfies every spacing/aperture constraint; change network or explicit bounds')
    _, _, center, selected, g = min(choices, key=lambda x: x[:3])
    selected.sort(key=lambda r: r.seed_id)
    result = {'schema': 'ndstudio.array_selection/1', 'provider': provider, 'network': network,
              'start_utc': str(t0), 'end_utc': str(t1), 'anchor_seed_id': center,
              'constraints': {'station_count': size,'minimum_station_count':min_stations,
                              'maximum_station_count':max_stations,'minimum_nearest_neighbor_m': min_spacing_m,
                              'maximum_nearest_neighbor_m': max_spacing_m,
                              'maximum_pair_distance_m': max_aperture_m},
              'metadata_rows_queried': len(rows), 'eligible_stations': len(usable),
              'station_pattern':station,'station_filters':station_filters,'matching_station_catalog':download.station_catalog(rows),
              'channel_filter': channel, 'rejected_metadata': rejected, 'geometry': g,
              'stations': [{k: str(v) if isinstance(v, UTCDateTime) else v
                            for k, v in asdict(r).items()} for r in selected],
              'selection_uses_waveforms_or_previous_curves': False}
    write_json(output, result)
    preview_array(result, Path(output).with_suffix('.png'))
    return result


def list_matching_stations(provider,network,output,*,station='*',channel='*',start=None,end=None,
                           client=None,**station_filters):
    """Complete station/epoch catalogue, without automatic compact-subset selection."""
    import csv
    filters=download.normalize_station_filters(**station_filters)
    t0,t1=UTCDateTime(start) if start else None,UTCDateTime(end) if end else None
    if t0 is not None and t1 is not None and t1<=t0:raise ValueError('End UTC must follow start UTC')
    rows=download.list_stations(provider,network,station,channel=channel,
                               starttime=t0,endtime=t1,client=client,**filters)
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    catalog=download.export_station_catalog(rows,output/'stations.csv')
    serialized=[{k:str(v) if isinstance(v,UTCDateTime) else v for k,v in asdict(r).items()} for r in rows]
    with (output/'channels.csv').open('w',newline='',encoding='utf-8-sig') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(download.ChannelRow.__dataclass_fields__))
        writer.writeheader();writer.writerows(serialized)
    result={'schema':'ndstudio.station_catalog/1','provider':provider,'network':network,
            'station_pattern':station,'channel_filter':channel,'start_utc':str(t0) if t0 else None,
            'end_utc':str(t1) if t1 else None,'station_filters':filters,
            'matching_station_count':len(catalog),'matching_channel_epochs':len(rows),
            'automatic_subset_selection_applied':False,'stations':catalog,'channel_epochs':serialized}
    write_json(output/'station_query.json',result)
    return result


def preview_array(selection, filename):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from pyproj import CRS, Transformer
    points = selection['geometry']['points']
    lon0 = np.mean([p[3] for p in points]); lat0 = np.mean([p[2] for p in points])
    crs = CRS.from_proj4(f'+proj=aeqd +lat_0={lat0} +lon_0={lon0} +datum=WGS84 +units=m')
    transform = Transformer.from_crs('EPSG:4326', crs, always_xy=True)
    fig, ax = plt.subplots(figsize=(7, 6))
    for net, sta, lat, lon, *_ in points:
        x, y = transform.transform(lon, lat)
        ax.scatter(x, y, marker='^', color='#146c94', s=75)
        ax.annotate(f'{net}.{sta}', (x, y), xytext=(5, 5), textcoords='offset points', fontsize=9)
    g = selection['geometry']
    ax.set_title(f"{selection['network']} | {len(points)} stations\nMean nearest-neighbor spacing: {g['mean_nearest_neighbor_m']:.1f} m")
    ax.set(xlabel='East (m)', ylabel='North (m)', aspect='equal')
    ax.grid(alpha=.25); fig.tight_layout(); fig.savefig(filename, dpi=160); plt.close(fig)


@contextmanager
def acquisition_lock(output):
    """OS-owned lock survives no crashed process and excludes another writer."""
    Path(output).mkdir(parents=True, exist_ok=True)
    handle = (Path(output)/'.acquisition.lock').open('a+b')
    handle.seek(0)
    if not handle.read(1):
        handle.write(b'0'); handle.flush()
    handle.seek(0)
    try:
        try:
            if __import__('os').name == 'nt':
                import msvcrt
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise RuntimeError('Another acquisition is writing this directory; wait for it to finish') from exc
        yield
    finally:
        handle.close()


def fetch_array(selection_file, output, *, chunk_seconds=21600., resume=True, client=None, workers=4, allow_missing_data=False):
    with acquisition_lock(output):
        return _fetch_array(selection_file, output, chunk_seconds=chunk_seconds,
                            resume=resume, client=client, workers=workers, allow_missing_data=allow_missing_data)


def _fetch_array(selection_file, output, *, chunk_seconds, resume, client, workers, allow_missing_data):
    selection = json.loads(Path(selection_file).read_text(encoding='utf-8-sig'))
    if selection.get('schema') != 'ndstudio.array_selection/1':
        raise ValueError('Expected an array-select output file')
    jobs = []
    for item in selection['stations']:
        row = dict(item)
        for key in ('start', 'end'):
            row[key] = UTCDateTime(row[key]) if row[key] else None
        jobs.append(download.DownloadJob(download.ChannelRow(**row),
                    UTCDateTime(selection['start_utc']), UTCDateTime(selection['end_utc'])))
    client = client or download.make_client(selection['provider'])
    Path(output).mkdir(parents=True, exist_ok=True)
    write_json(Path(output)/'array_selection.json', selection)
    if not isinstance(workers, int) or not 1 <= workers <= 4:
        raise ValueError('Download workers must be an integer from 1 to 4')
    from concurrent.futures import ThreadPoolExecutor
    from threading import Lock
    progress_lock = Lock()
    def progress(done, total, message):
        with progress_lock:
            print(f'fetch: {done}/{total} {message}', flush=True)
    def acquire(job):
        return download.run_downloads([job], output, client, chunk_seconds=chunk_seconds,
                 convert_to_h5=True, layout='array', resume=resume,
                 max_consecutive_errors=3, provider=selection['provider'], write_manifest=False,
                 array_geometry_override=selection['geometry'],
                 progress=progress)[0]
    # One job per physical station; no two workers write the same waveform or
    # response file. Publish the combined manifest only after all jobs finish.
    with ThreadPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(acquire, jobs))
    write_run_manifest(output, jobs, results, provider=selection['provider'],
                       layout='array', chunk_seconds=chunk_seconds, resume=resume)
    all_h5 = (len(results) == len(jobs) and all(not r.errors and not r.cancelled and
                not r.aborted and all(h5_is_complete(
                    download.job_out_dir(output, job.row, 'array') /
                    Path(download.chunk_filename(job.row, a, b)).with_suffix('.h5'))
                    for a, b in download.split_into_chunks(job.start, job.end, chunk_seconds))
                for job, r in zip(jobs, results)))
    complete=all_h5 or (allow_missing_data and len(results)==len(jobs) and
                        all(not r.errors and not r.cancelled and not r.aborted and
                            r.chunks_done==r.chunks_total for r in results))
    write_json(Path(output)/'acquisition_status.json', {'complete': complete,
               'all_requested_chunks_have_h5':all_h5,'allow_missing_data':allow_missing_data,
               'results': json_value([asdict(r) for r in results])})
    if not complete:
        raise RuntimeError('Acquisition incomplete; inspect acquisition_status.json and rerun fetch to resume')
    return results


def audit_acquisition(output):
    """Verify actual archives against embedded original FDSN MiniSEED bytes."""
    import hashlib
    import h5py
    from io import BytesIO
    from obspy import read, read_inventory
    from h5data_handoff import matching_channel
    root=Path(output);records=[]
    for path in sorted(root.rglob('*.h5')):
        if not h5_is_complete(path):
            raise ValueError('Incomplete H5 archive: '+str(path))
        with h5py.File(path) as h:
            raw=h['metadata/original_mseed'][()].tobytes()
            if hashlib.sha256(raw).hexdigest()!=h.attrs['source_mseed_sha256']:
                raise ValueError('Embedded MiniSEED checksum mismatch: '+str(path))
            original=path.with_suffix('.mseed')
            if original.exists() and original.read_bytes()!=raw:
                raise ValueError('MiniSEED sidecar differs from embedded bytes: '+str(path))
            xml=h['metadata/stationxml'][()]
            inventory=read_inventory(BytesIO(xml))
            stream=read(BytesIO(raw),format='MSEED')
            datasets=[obj for obj in h.values() if isinstance(obj,h5py.Dataset) and obj.attrs.get('role')=='waveform']
            if len(datasets)!=len(stream):
                raise ValueError('Trace count mismatch: '+str(path))
            remaining=list(datasets)
            for tr in stream:
                # HDF keys are sorted, so identify each actual trace by NSLC
                # and nanosecond start rather than assuming key order.
                matches=[q for q in remaining if q.attrs['seed_id']==tr.id and
                         int(q.attrs['start_utc_ns'])==tr.stats.starttime.ns and
                         len(q)==tr.stats.npts and float(q.attrs['sampling_rate_hz'])==tr.stats.sampling_rate]
                ds=next((q for q in matches if np.array_equal(q[()],tr.data,equal_nan=True)),None)
                if ds is None:raise ValueError('No identical archived trace: '+tr.id)
                remaining=[q for q in remaining if q.name!=ds.name]
                matching_channel(inventory,tr)
                if float(ds.attrs['sampling_rate_hz'])!=tr.stats.sampling_rate or int(ds.attrs['end_utc_ns'])!=tr.stats.endtime.ns:
                    raise ValueError('Native UTC/rate changed: '+str(path))
            records.append({'file':str(path.relative_to(root)), 'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                            'traces':len(stream), 'samples':sum(t.stats.npts for t in stream),
                            'native_values_utc_rate_equal_to_original_mseed':True,
                            'stationxml_sha256':str(h.attrs['stationxml_sha256'])})
        print('acquisition audit:',path.name,'passed',flush=True)
    if not records:raise ValueError('No H5 files to audit')
    result={'checks_passed':True,'h5_files':len(records),'native_samples':sum(r['samples'] for r in records),
            'records':records,'original_mseed_and_full_stationxml_embedded':True}
    write_json(root/'acquisition_audit.json',result)
    return result


def screen_array(selection_file, folder, project=None):
    """Conservative waveform-input screening, before phase/dispersion scanning.

    A station is excluded for extreme counts DC/robust fluctuation (>10,000)
    in every downloaded trace, or for no finite variation. This is an explicit
    screening policy, not a diagnosis of hardware saturation. Keep all H5
    archives and report exclusions. Then remove geographically isolated
    stations to keep the original nearest-neighbor constraints satisfied.
    """
    import h5py
    selection=json.loads(Path(selection_file).read_text(encoding='utf-8-sig'))
    root=Path(folder);accepted=[];exclusions=[];quality=[]
    for item in selection['stations']:
        metrics=[]
        for path in sorted((root/item['network']/item['station']).glob('*.h5')):
            with h5py.File(path) as h:
                for ds in h.values():
                    if not isinstance(ds,h5py.Dataset) or ds.attrs.get('role')!='waveform':continue
                    x=ds[()];x=x[np.isfinite(x)]
                    median=float(np.median(x)) if len(x) else 0.
                    mad=float(1.4826*np.median(abs(x.astype(float)-median))) if len(x) else 0.
                    ratio=abs(median)/mad if mad else None
                    suspect=not len(x) or not mad or ratio>10000
                    metrics.append({'file':path.name,'dataset':ds.name,'samples':len(x),
                                    'median_counts':median,'robust_scale_counts':mad,
                                    'minimum_counts':float(x.min()) if len(x) else None,
                                    'maximum_counts':float(x.max()) if len(x) else None,
                                    'absolute_dc_to_robust_scale':ratio,'screening_flag':bool(suspect)})
        bad=not metrics or all(m['screening_flag'] for m in metrics)
        quality.append({'seed_id':'.'.join([item['network'],item['station'],item['location'],item['channel']]),
                        'traces':metrics,'excluded_by_waveform_screen':bad})
        if bad:exclusions.append({'station':item['network']+'_'+item['station'],
                                  'reason':'extreme_DC_relative_to_variation_in_all_records_or_no_variation'})
        else:accepted.append(item)
    def geometry(items):
        return download.array_geometry([download.ChannelRow(**(r|{'start':None,'end':None})) for r in items])
    bounds=selection['constraints']
    minimum_count=selection['constraints'].get('minimum_station_count',6)
    maximum_count=selection['constraints'].get('maximum_station_count',49)
    while len(accepted)>=minimum_count:
        g=geometry(accepted)
        isolated=[key for key,value in g['nearest_neighbor_m'].items()
                  if value>bounds['maximum_nearest_neighbor_m']]
        if not isolated:break
        accepted=[r for r in accepted if r['network']+'.'+r['station'] not in isolated]
        exclusions.extend({'station':key.replace('.','_'),
                           'reason':'nearest_neighbor_exceeds_selected_bound_after_waveform_exclusion'} for key in isolated)
    if not minimum_count<=len(accepted)<=maximum_count:raise ValueError('Too few usable stations after waveform/geometric screening')
    g=geometry(accepted)
    if min(g['nearest_neighbor_m'].values())<bounds['minimum_nearest_neighbor_m']:
        raise ValueError('Screened stations violate the minimum spacing constraint')
    analysis=selection|{'stations':accepted,'geometry':g,
                        'constraints':bounds|{'station_count':len(accepted)},'requested_constraints':bounds,
                        'original_station_count':len(selection['stations']),
                        'screened_before_dispersion':True,'exclusions':exclusions}
    write_json(root/'analysis_array_selection.json',analysis)
    preview_array(analysis,root/'analysis_array_selection.png')
    report={'screened_stations':len(accepted),'geometry':g,'exclusions':exclusions,
            'quality':quality,'all_original_h5_retained':True,'dc_ratio_screening_threshold':10000,
            'screening_is_not_a_verified_hardware_fault_diagnosis':True}
    write_json(root/'array_screening.json',report)
    if project:
        specfile=Path(project)/'project.json'
        spec=json.loads(specfile.read_text(encoding='utf-8-sig'))
        spec['excluded_stations']=[r['station'] for r in exclusions]
        write_json(specfile,spec)
    return report
