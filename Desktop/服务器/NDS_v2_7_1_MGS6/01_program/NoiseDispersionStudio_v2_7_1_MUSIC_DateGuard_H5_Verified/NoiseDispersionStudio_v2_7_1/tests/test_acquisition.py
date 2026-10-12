from dataclasses import replace
import json
import pytest
from pyproj import Geod
import h5data as d
from ndstudio.acquisition import select_array, fetch_array, acquisition_lock
from test_download_handoff import FakeFDSN, inventory, DATE


def rows(n=7):
    geod=Geod(ellps='WGS84');result=[]
    for i in range(n):
        lon,lat,_=geod.fwd(11,48,90,500*i)
        result.append(d.ChannelRow('XX',chr(65+i),'','HHZ',20,lat,lon,DATE-10,DATE+1000))
    return result


def select(tmp_path, monkeypatch, source, **kwargs):
    monkeypatch.setattr(d,'list_stations',lambda *a,**k: source)
    return select_array('FAKE','XX',DATE,DATE+120,tmp_path/'selection.json',
                        client=object(),size=6,**kwargs)


def test_selection_bounds_each_station_and_full_epoch(tmp_path,monkeypatch):
    source=rows();source.append(replace(source[0],station='BAD',start=DATE+10))
    r=select(tmp_path,monkeypatch,source)
    assert len(r['stations'])==6 and r['metadata_rows_queried']==8
    assert r['geometry']['mean_nearest_neighbor_m']==pytest.approx(500,abs=.01)
    assert not r['selection_uses_waveforms_or_previous_curves']
    assert r['rejected_metadata'][0]['reason']=='epoch_does_not_cover_full_interval'
    assert (tmp_path/'selection.png').is_file()
    with pytest.raises(ValueError,match='No compact subset'):
        select(tmp_path,monkeypatch,source,max_spacing_m=400)


def test_ambiguous_channel_and_count_are_not_silently_accepted(tmp_path,monkeypatch):
    source=rows(6);source.append(replace(source[0],location='01'))
    with pytest.raises(ValueError,match='Only 5'):
        select(tmp_path,monkeypatch,source)
    with pytest.raises(ValueError,match='5 < n < 50'):
        select_array('FAKE','XX',DATE,DATE+120,tmp_path/'a.json',size=50,client=object())

def test_explicit_station_names_prevent_nearby_station_substitution(tmp_path,monkeypatch):
    source=rows(7)
    monkeypatch.setattr(d,'list_stations',lambda *a,**k:source)
    result=select_array('FAKE','XX',DATE,DATE+120,tmp_path/'selection.json',size=5,
                         min_stations=5,station='A,B,C,D,E',client=object())
    assert {r['station'] for r in result['stations']}==set('ABCDE')
    assert result['eligible_stations']==5 and result['station_pattern']=='A,B,C,D,E'


def test_fetch_resume_and_parallel_manifest_preserve_full_geometry(tmp_path,monkeypatch):
    from h5data_handoff import geometry
    inv=inventory(stations=tuple('ABCDEF'));source=d.inventory_to_rows(inv)
    monkeypatch.setattr(d,'list_stations',lambda *a,**k: source)
    select_array('FAKE','XX',DATE,DATE+120,tmp_path/'selection.json',size=6,
                 min_spacing_m=100,max_spacing_m=200,client=object())
    client=FakeFDSN(inv);out=tmp_path/'download'
    fetch_array(tmp_path/'selection.json',out,chunk_seconds=60,client=client,workers=4)
    calls=client.waveform_calls
    fetch_array(tmp_path/'selection.json',out,chunk_seconds=60,client=client,workers=4)
    assert client.waveform_calls==calls==18
    assert json.loads((out/'acquisition_status.json').read_text())['complete']
    manifest=json.loads((out/'download_manifest.json').read_text())
    assert len(manifest['jobs'])==6
    import h5py
    for path in out.rglob('*.h5'):
        with h5py.File(path) as h:
            assert json.loads(h['metadata/array_geometry'][()])['stations']==6


def test_acquisition_directory_excludes_second_writer_and_releases(tmp_path):
    with acquisition_lock(tmp_path):
        with pytest.raises(RuntimeError,match='Another acquisition'):
            with acquisition_lock(tmp_path):pass
    with acquisition_lock(tmp_path):pass


@pytest.mark.parametrize('count',[5,100])
def test_explicit_station_limits_support_second_test_range(tmp_path,monkeypatch,count):
    import numpy as np
    geod=Geod(ellps='WGS84');source=[]
    for i in range(count):
        lon,lat,_=geod.fwd(11,48,i*360/count,4000)
        source.append(d.ChannelRow('XX',f'S{i:03d}','','HHZ',100,lat,lon,DATE-10,DATE+1000))
    monkeypatch.setattr(d,'list_stations',lambda *a,**k:source)
    r=select_array('FAKE','XX',DATE,DATE+120,tmp_path/'selection.json',client=object(),
                   size=count,min_stations=5,max_stations=100,max_spacing_m=5000,max_aperture_m=8001)
    assert r['geometry']['stations']==count


def test_temporal_station_screen_uses_availability_and_preserves_spacing(tmp_path,monkeypatch):
    import numpy as np
    import pandas as pd
    from ndstudio.project import new_project
    from ndstudio.diagnostics import screen_prepared_array
    source=rows(7);monkeypatch.setattr(d,'list_stations',lambda *a,**k:source)
    data=tmp_path/'data';data.mkdir()
    selection=select_array('FAKE','XX',DATE,DATE+120,data/'analysis_array_selection.json',size=7,client=object())
    project=new_project(tmp_path/'project');cache=project/'output/cache';cache.mkdir(parents=True)
    pd.DataFrame({'station':['XX_'+r.station for r in source]}).to_csv(project/'stations.csv',index=False)
    valid=np.ones((7,24000),bool);valid[0,:20000]=False
    np.save(cache/'aligned_valid.npy',valid)
    pd.DataFrame({'accepted':[False]}).to_csv(project/'output/window_quality.csv',index=False)
    result=screen_prepared_array(project,data)
    assert len(result['stations'])==6
    assert [r['station'] for r in result['stations']]==list('BCDEFG')
    spec=json.loads((project/'project.json').read_text())
    assert spec['excluded_stations']==['XX_A']
    assert not result['temporal_screening_uses_phase_or_dispersion']
