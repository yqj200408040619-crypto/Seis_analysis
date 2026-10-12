"""FDSN byte transport → complete H5 → movable Studio input, offline reproducible."""
from pathlib import Path
from io import BytesIO
import json
import shutil
import copy
from dataclasses import replace
import numpy as np
import h5py
import pytest
from obspy import UTCDateTime, Trace, Stream, read
from obspy.core.inventory import Inventory, Network, Station, Channel, Site, Equipment
from obspy.core.inventory.response import Response
import h5data as download
import h5data_handoff as bridge
from ndstudio.hdf import inspect_folder, read_waveform, MappingError
from ndstudio.project import new_project, import_hdf
from dispersion.common import load_config
from dispersion.prepare import read_station, response_table, correct_response

DATE = UTCDateTime('2026-01-01T00:00:00.123456Z')


def inventory(stations=('A','B','C'), channels=('HHZ',), fs=20):
    response = Response.from_paz([], [-2*np.pi*3], stage_gain=1000.,
                                input_units='M/S', output_units='COUNTS', normalization_factor=1.)
    rows = []
    for i, sta in enumerate(stations):
        chs = [Channel(code=c, location_code='', latitude=48., longitude=11.+.0015*i,
                       elevation=500.+i, depth=1., azimuth=0., dip=-90. if c.endswith('Z') else 0.,
                       sample_rate=fs, start_date=DATE-86400, end_date=DATE+86400*40,
                       response=copy.deepcopy(response), sensor=Equipment(type='Broadband', description='TEST sensor')) for c in channels]
        rows.append(Station(code=sta, latitude=48., longitude=11.+.0015*i, elevation=500.+i,
                            site=Site(name='Test field; metadata must be retained'), channels=chs))
    return Inventory([Network(code='XX', stations=rows, description='Offline transport fixture')], source='Fixture')


class FakeFDSN:
    def __init__(self, inv=None, fail_metadata=False, gaps=False):
        self.inv=inv or inventory();self.fail_metadata=fail_metadata;self.gaps=gaps
        self.waveform_calls=0;self.station_calls=0;self.payloads=[]
    def get_stations(self, filename=None, **kwargs):
        self.station_calls+=1
        assert kwargs['level']=='response'
        if self.fail_metadata:raise RuntimeError('response service unavailable')
        buf=BytesIO();self.inv.write(buf,format='STATIONXML');payload=buf.getvalue()
        if filename is not None:filename.write(payload);return None
        return self.inv.copy()
    def get_waveforms(self, net, sta, loc, cha, start, end, filename=None):
        self.waveform_calls+=1;fs=20.
        lo=int(np.ceil((start.ns-DATE.ns)/1e9*fs-1e-5));hi=int(np.floor((end.ns-DATE.ns)/1e9*fs+1e-5))+1
        # Deterministic integer counts preserve exact overlap at chunk boundaries.
        k=np.arange(lo,hi);delay=(ord(sta)-ord('A'))*.2
        x=np.round(10000*(np.sin(2*np.pi*2.3*(k/fs-delay))+.3*np.cos(2*np.pi*3.1*(k/fs-delay)))).astype(np.int32)
        header=dict(network=net,station=sta,location=loc,channel=cha,sampling_rate=fs,starttime=DATE+lo/fs)
        tr=Trace(x,header);st=Stream([tr])
        if self.gaps and len(x)>50:st=Stream([tr.slice(tr.stats.starttime,tr.stats.starttime+19/fs),tr.slice(tr.stats.starttime+40/fs,tr.stats.endtime)])
        buf=BytesIO();st.write(buf,format='MSEED');payload=buf.getvalue();self.payloads.append(payload)
        if filename is not None:filename.write(payload);return None
        return st


def jobs(inv=None,duration=120):
    return [download.DownloadJob(r,DATE,DATE+duration) for r in download.inventory_to_rows(inv or inventory())]


def studio(root, folder):
    p=new_project(root);spec=json.loads((p/'project.json').read_text());spec.update(input_folder=str(folder),automatic_geometry=False)
    (p/'project.json').write_text(json.dumps(spec));return p


@pytest.mark.parametrize('changed_response',[False,True])
def test_resumed_xml_headers_are_preserved_and_real_response_changes_rejected(tmp_path,changed_response):
    raw=tmp_path/'raw';client=FakeFDSN()
    download.run_downloads(jobs(duration=60),raw,client,60,True,layout='array')
    newer=client.inv.copy();newer.source='Service regenerated XML header'
    if changed_response:newer[0][0][0].response.instrument_sensitivity.value+=1
    later=[download.DownloadJob(r,DATE+60,DATE+120) for r in download.inventory_to_rows(newer)]
    download.run_downloads(later,raw,FakeFDSN(newer),60,True,layout='array')
    p=studio(tmp_path/'project',raw)
    if changed_response:
        with pytest.raises(MappingError,match='响应或单位不同'):import_hdf(p)
    else:
        import_hdf(p)
        assert len(list((p/'metadata').glob('*.xml')))==2
        audit=json.loads((p/'hdf_import_audit.json').read_text())['records']
        assert len({r['embedded_stationxml_sha256'] for r in audit})==2
        assert len({r['canonical_response_xml'] for r in audit})==1


def test_nested_identical_record_overlaps_preserve_every_unique_sample(tmp_path):
    inv=inventory(stations=('A','B'));buf=BytesIO();inv.write(buf,format='STATIONXML')
    raw=tmp_path/'raw';raw.mkdir();expected=np.arange(250,dtype=np.int32)
    for sta in ('A','B'):
        for number,(start,lengths) in enumerate([(0,(120,150)),(100,(120,150))]):
            st=Stream([Trace(expected[start:start+n].copy(),dict(network='XX',station=sta,location='',channel='HHZ',
                       sampling_rate=20,starttime=DATE+start/20)) for n in lengths])
            path=raw/f'{sta}_{number}.mseed';st.write(str(path),format='MSEED')
            bridge.export_h5(path,path.with_suffix('.h5'),stationxml=buf.getvalue(),strict=True)
    p=studio(tmp_path/'project',raw);import_hdf(p)
    c=load_config(p/'config.json')
    import pandas as pd
    for row in pd.read_csv(p/'stations.csv',keep_default_na=False).to_dict('records'):
        st,_=read_station(row,c)
        assert len(st)==1 and np.array_equal(st[0].data,expected)


@pytest.mark.parametrize('microseconds,conflict,accepted',[(2,False,True),(-2,False,True),(3,False,False),(-3,False,False),(1,True,False)])
def test_native_chunk_header_rounding_requires_identical_samples_and_preserves_utc(tmp_path,microseconds,conflict,accepted):
    inv=inventory(stations=('A','B'),fs=200);buf=BytesIO();inv.write(buf,format='STATIONXML')
    raw=tmp_path/'raw';raw.mkdir();expected=np.arange(300,dtype=np.int32)
    for sta in ('A','B'):
        for number,start in enumerate((0,100)):
            data=expected[start:start+200].copy()
            if number and conflict:data[3]+=1
            stamp=DATE+start/200+(microseconds*1e-6 if number else 0)
            tr=Trace(data,dict(network='XX',station=sta,location='',channel='HHZ',sampling_rate=200,starttime=stamp))
            path=raw/f'{sta}_{number}.mseed';Stream([tr]).write(str(path),format='MSEED')
            bridge.export_h5(path,path.with_suffix('.h5'),stationxml=buf.getvalue(),strict=True)
    p=studio(tmp_path/'project',raw)
    if not accepted:
        with pytest.raises(MappingError,match='重叠'):import_hdf(p)
        return
    import_hdf(p)
    for sta in ('A','B'):
        with np.load(p/f'canonical/XX_{sta}_00001.npz') as z:
            np.testing.assert_array_equal(z['data'],expected[200:])
            assert UTCDateTime(z['start_utc'].item()).ns==(DATE+1+microseconds*1e-6).ns
    audit=json.loads((p/'hdf_import_audit.json').read_text())['records']
    assert max(r['maximum_identical_overlap_time_residual_ns'] for r in audit)==2000
    assert all(not r['original_native_utc_shifted'] for r in audit)
    import_hdf(p)
    assert all(r['canonical_cache_reused'] for r in json.loads((p/'hdf_import_audit.json').read_text())['records'])


def test_geometry_is_metres_unique_stations_and_not_sorted_adjacent():
    from pyproj import Geod
    g=Geod(ellps='WGS84');origin=(11.,48.)
    rows=[]
    for name,x in [('Z',0),('A',100),('M',200)]:
        lon,lat,_=g.fwd(*origin,90,x)
        rows.append(download.ChannelRow('XX',name,'','HHZ',20,lat,lon,DATE,None))
    rows.append(replace(rows[0],channel='HHN'))
    info=bridge.geometry(rows)
    assert info['stations']==3
    assert info['mean_nearest_neighbor_m']==pytest.approx(100,abs=1e-6)
    assert info['mean_pair_distance_m']==pytest.approx(400/3,abs=1e-6)
    assert info['min_pair_distance_m']==pytest.approx(100,abs=1e-6)


def test_geometry_missing_single_colocated_and_dateline():
    r=download.ChannelRow('XX','A','','HHZ',20,0.,179.999,DATE,None)
    assert bridge.geometry([r])['mean_nearest_neighbor_m'] is None
    info=bridge.geometry([r,replace(r,station='B',longitude=-179.999)])
    assert 222<info['mean_nearest_neighbor_m']<223
    info=bridge.geometry([r,replace(r,station='B'),replace(r,station='C',latitude=None)])
    assert info['missing_coordinates']==1 and info['coincident_pair_count']==1


def test_geometry_reports_coordinate_epoch_conflict():
    r=download.inventory_to_rows(inventory())[0]
    info=bridge.geometry([r,replace(r,latitude=r.latitude+1),replace(r,station='B')])
    assert info['coordinate_conflicts']==['XX.A'] and info['stations_with_coordinates']==1


def test_full_download_is_byte_lossless_and_importable_after_move(tmp_path):
    client=FakeFDSN();raw=tmp_path/'download';res=download.run_downloads(jobs(),raw,client,60,True,layout='array',provider='OFFLINE')
    assert not any(r.errors for r in res);assert client.station_calls==3
    scan=inspect_folder(raw);assert scan['ready'];assert len(scan['items'])==6
    for item in scan['items']:
        with h5py.File(item['file'],'r') as h:
            source=Path(item['file']).with_suffix('.mseed').read_bytes()
            assert h['metadata/original_mseed'][()].tobytes()==source
            assert source in client.payloads
            assert b'Test field' in h['metadata/stationxml'][()]
            d=h[item['dataset']];original=read(BytesIO(source))[0]
            np.testing.assert_array_equal(d[:],original.data)
            assert d.attrs['start_utc_ns']==original.stats.starttime.ns
            assert d.attrs['units']=='counts' and not d.attrs['precorrected']
            assert not d.attrs['polarity_verified'] and not d.attrs['clock_offset_verified']
            assert d.attrs['dip_deg']==-90.
            assert json.loads(d.attrs['trace_stats_json'])['mseed']
    manifest=json.loads((raw/'download_manifest.json').read_text());assert manifest['geometry']['stations']==3
    moved=tmp_path/'moved';moved.mkdir()
    for item in scan['items']:shutil.copy2(item['file'],moved/Path(item['file']).name)
    # Only H5 files: no XML sidecars or original MiniSEED files carried over.
    project=studio(tmp_path/'project',moved);import_hdf(project)
    table=__import__('pandas').read_csv(project/'stations.csv',keep_default_na=False)
    c=load_config(project/'config.json')
    for row in table.to_dict('records'):
        assert (project/row['response_xml']).is_file()
        st,_=read_station(row,c);_,audit=response_table(row,st,c)
        assert audit['response_status']=='nominal_response_corrected'
    assert all(x['retained_utc_and_values_preserved'] for x in json.loads((project/'hdf_import_audit.json').read_text())['records'])


def test_horizontal_data_is_retained_but_not_used_as_vertical(tmp_path):
    inv=inventory(channels=('HHZ','HHN','HHE'));client=FakeFDSN(inv)
    download.run_downloads(jobs(inv,60),tmp_path/'download',client,60,True)
    scan=inspect_folder(tmp_path/'download');assert scan['ready']
    assert len(scan['items'])==9 and sum(i['enabled'] for i in scan['items'])==3
    p=studio(tmp_path/'project',tmp_path/'download');import_hdf(p)


def test_resume_does_not_request_completed_chunks_or_response(tmp_path):
    client=FakeFDSN();download.run_downloads(jobs(),tmp_path,client,60,True)
    calls=(client.waveform_calls,client.station_calls)
    results=download.run_downloads(jobs(),tmp_path,client,60,True,resume=True)
    assert calls==(client.waveform_calls,client.station_calls)
    assert sum(len(r.existing) for r in results)==sum(download.count_chunks(j.start,j.end,60) for j in jobs())


def test_resume_upgrades_legacy_h5_and_missing_h5_without_redownload(tmp_path):
    client=FakeFDSN();js=jobs(duration=60);res=download.run_downloads(js,tmp_path,client,60,False)
    paths=[r.files[0] for r in res]
    download.convert_mseed_to_h5_strict(paths[0],paths[0].with_suffix('.h5'))
    calls=client.waveform_calls
    results=download.run_downloads(js,tmp_path,client,60,True,resume=True)
    assert not any(r.errors for r in results);assert client.waveform_calls==calls
    assert all(bridge.h5_is_complete(p.with_suffix('.h5')) for p in paths)


def test_metadata_failure_keeps_raw_files_and_does_not_invent_h5(tmp_path):
    client=FakeFDSN(fail_metadata=True)
    res=download.run_downloads(jobs(duration=60),tmp_path,client,60,True)
    assert all(r.errors for r in res);assert len(list(tmp_path.glob('*.mseed')))==sum(download.count_chunks(j.start,j.end,60) for j in jobs(duration=60))
    assert not list(tmp_path.glob('*.h5'))
    assert not list(tmp_path.rglob('*.part'))


@pytest.mark.parametrize('defect',['missing_response','epoch','rate','units'])
def test_invalid_metadata_is_rejected_without_overwriting_existing_h5(tmp_path,defect):
    client=FakeFDSN();r=jobs()[0].row;mseed=tmp_path/'a.mseed'
    client.get_waveforms('XX','A','','HHZ',DATE,DATE+10,filename=(buf:=BytesIO()))
    mseed.write_bytes(buf.getvalue());inv=inventory();ch=inv[0][0][0]
    if defect=='missing_response':ch.response=None
    elif defect=='epoch':ch.end_date=DATE+1
    elif defect=='rate':ch.sample_rate=40.
    else:ch.response.instrument_sensitivity.output_units='V'
    xml=BytesIO();inv.write(xml,format='STATIONXML')
    target=tmp_path/'a.h5';target.write_bytes(b'previous complete file')
    with pytest.raises(ValueError):bridge.export_h5(mseed,target,stationxml=xml.getvalue(),strict=True)
    assert target.read_bytes()==b'previous complete file'


def test_gaps_are_preserved_as_distinct_native_traces(tmp_path):
    client=FakeFDSN(gaps=True);download.run_downloads(jobs(duration=60),tmp_path,client,60,True)
    scan=inspect_folder(tmp_path);assert scan['ready'] and len(scan['items'])==6
    p=studio(tmp_path/'project',tmp_path);import_hdf(p)
    table=__import__('pandas').read_csv(p/'stations.csv',keep_default_na=False);c=load_config(p/'config.json')
    st,_=read_station(table.iloc[0].to_dict(),c);assert len(st)==2
    assert st[1].stats.starttime.ns-st[0].stats.endtime.ns==1050000000


def test_embedded_response_missing_or_wrong_epoch_blocks_import(tmp_path):
    client=FakeFDSN();download.run_downloads(jobs(duration=60),tmp_path/'download',client,60,True)
    path=next((tmp_path/'download').glob('*.h5'))
    with h5py.File(path,'a') as h:del h['metadata/stationxml']
    p=studio(tmp_path/'project',tmp_path/'download')
    with pytest.raises(MappingError,match='内嵌'):import_hdf(p)


def test_xml_transport_retains_exact_server_payload(tmp_path):
    client=FakeFDSN();row=jobs()[0].row
    xml,_=bridge.fetch_stationxml(client,row,DATE,DATE+60)
    source=BytesIO();client.inv.write(source,format='STATIONXML')
    assert xml==source.getvalue()


def test_subsecond_request_filenames_do_not_collide():
    row=jobs()[0].row
    assert download.chunk_filename(row,DATE,DATE+1)!=download.chunk_filename(row,DATE+.1,DATE+1.1)


def test_identical_station_codes_in_different_networks_remain_distinct(tmp_path):
    inv=inventory(stations=('A',));other=copy.deepcopy(inv[0]);other.code='YY';inv.networks.append(other)
    client=FakeFDSN(inv)
    download.run_downloads(jobs(inv,60),tmp_path/'download',client,60,True)
    project=studio(tmp_path/'project',tmp_path/'download');import_hdf(project)
    table=__import__('pandas').read_csv(project/'stations.csv')
    assert set(table.station)=={'XX_A','YY_A'} and set(table.seed_id)=={'XX.A..HHZ','YY.A..HHZ'}


def test_embedded_xml_checksum_rejects_changed_payload(tmp_path):
    client=FakeFDSN();download.run_downloads(jobs(duration=60),tmp_path/'download',client,60,True)
    path=next((tmp_path/'download').glob('*.h5'))
    with h5py.File(path,'a') as h:
        xml=h['metadata/stationxml'][()].replace(b'Offline transport fixture',b'Changed transport fixture')
        del h['metadata/stationxml'];h['metadata'].create_dataset('stationxml',data=xml.decode(),dtype=h5py.string_dtype('utf-8'))
    project=studio(tmp_path/'project',tmp_path/'download')
    with pytest.raises(MappingError,match='校验失败'):import_hdf(project)
