"""Native single-station download entry using Studio's existing FDSN backend."""
from pathlib import Path
from dataclasses import asdict
import json
from obspy import UTCDateTime
import h5data as download
from h5data_handoff import json_value
from .hvsr import discover,write_json
from .acquisition import audit_acquisition


def fetch_single_station(provider,network,station,start,end,output,*,channel='HH?',
                         instrument_id=None,catalog_file=None,chunk_seconds=1800.,resume=True,client=None):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    (output/'hvsr_download_complete.json').unlink(missing_ok=True)
    start,end=download.validate_time_window(start,end)
    client=client or download.make_client(provider)
    if catalog_file:
        catalog=json.loads(Path(catalog_file).read_text(encoding='utf-8-sig'))
        if catalog.get('provider')!=provider:raise ValueError('目录与下载服务商不一致')
        rows=[]
        for item in catalog['channel_epochs']:
            value=dict(item)
            for k in ('start','end'):value[k]=UTCDateTime(value[k]) if value.get(k) else None
            rows.append(download.ChannelRow(**value))
    else:
        rows=download.list_stations(provider,network,station,channel=channel,starttime=start,endtime=end,client=client)
    rows=download.filter_channels(download.filter_station_rows(rows,station=station),channel)
    rows=[r for r in rows if r.network==network]
    groups={}
    for row in rows:groups.setdefault('.'.join([row.network,row.station,row.location,row.channel[:-1]]),[]).append(row)
    valid={}
    for key,items in groups.items():
        suffix={r.channel[-1] for r in items}
        if len(items)==3 and (suffix==set('NEZ') or suffix==set('12Z')):valid[key]=items
    candidates=[key for key in valid if not instrument_id or key==instrument_id]
    if len(candidates)!=1:raise ValueError('请选择唯一三分量仪器 --instrument-id；可用：'+','.join(sorted(valid)))
    key=candidates[0];selected=valid[key]
    if len({r.sample_rate for r in selected})!=1:raise ValueError('三分量采样率不同，不能自动重采样')
    if any((r.start and r.start>start) or (r.end and r.end<end) for r in selected):raise ValueError('请求跨越通道历元，须缩短或分段下载分析')
    serial=[{k:str(v) if isinstance(v,UTCDateTime) else v for k,v in asdict(r).items()} for r in selected]
    selection=dict(provider=provider,instrument_id=key,start_utc=str(start),end_utc=str(end),channels=serial,
                   transport='Studio h5data.run_downloads; original MiniSEED + complete StationXML + H5',catalog_file=str(catalog_file) if catalog_file else None)
    write_json(output/'hvsr_download_selection.json',selection)
    jobs,skipped=download.build_jobs(selected,'window',start,end)
    if skipped or len(jobs)!=3:raise ValueError('三分量下载任务未全部覆盖请求范围')
    def progress(done,total,message):print(f'Download {done}/{total}: {message}',flush=True)
    results=download.run_downloads(jobs,output,client,chunk_seconds=chunk_seconds,convert_to_h5=True,
                                  layout='array',resume=resume,provider=provider,progress=progress,max_consecutive_errors=3)
    state=download.summarize_results(results,jobs,chunk_seconds)
    write_json(output/'hvsr_download_status.json',state)
    failures=[dict(seed_id=r.seed_id,errors=r.errors,skipped=r.skipped,aborted=r.aborted,cancelled=r.cancelled) for r in results if r.errors or r.aborted or r.cancelled]
    if len(results)!=3 or failures:raise ValueError('三分量下载未完成：'+json.dumps(json_value(failures),ensure_ascii=False))
    audit=audit_acquisition(output)
    inventory=discover(output);ready=[r for r in inventory['instruments'] if r['ready'] and r['instrument_id']==key]
    if len(ready)!=1:raise ValueError('下载后未发现完整可分析三分量，查看inventory.json')
    write_json(output/'hvsr_inventory.json',inventory)
    result=dict(selection=selection,download_status=state,audit=audit,source='live FDSN requests via application backend')
    write_json(output/'hvsr_download_complete.json',result)
    print('Single station download complete: '+key,flush=True)
    return result
