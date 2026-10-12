"""Station query intersection, complete catalogues and GUI download scope."""
from dataclasses import replace
import csv
import json
import sys
import numpy as np
import pytest
import h5data as d
from obspy.clients.fdsn.header import FDSNNoDataException
from test_download_handoff import inventory,DATE
from station_filters import normalize_station_filters,filter_station_rows


def hundred_station_inventory():
    return inventory(stations=tuple([f'QS{i:02d}' for i in range(20)]+[f'AA{i:02d}' for i in range(80)]),
                     channels=('HHZ','HHN','HHE'))


class MetadataClient:
    """Intentionally ignores station/geographic filters, to test local enforcement."""
    def __init__(self,inv=None):self.inv=inv or hundred_station_inventory();self.calls=[]
    def get_stations(self,**kwargs):
        self.calls.append(kwargs)
        return self.inv.select(channel=kwargs.get('channel','*'))


def test_100_station_array_returns_all_20_qs_stations():
    client=MetadataClient()
    rows=d.list_array_channels('TEST','XX',client=client,station_prefix='qs')
    assert client.calls[0]['station']=='QS*'
    assert len(rows)==60
    assert [x['station'] for x in d.station_catalog(rows)]==[f'QS{i:02d}' for i in range(20)]


def test_prefix_and_geography_and_existing_pattern_intersect():
    client=MetadataClient()
    rows=d.list_stations('TEST','XX','QS1?',client=client,station_prefix='QS',
                        minlatitude=48,maxlatitude=48,minlongitude=11.0225,maxlongitude=11.0285)
    assert {r.station for r in rows}=={f'QS{i:02d}' for i in range(15,20)}
    assert client.calls[0]['minlatitude']==48
    assert client.calls[0]['maxlongitude']==11.0285


def test_prefix_matches_start_not_middle_and_supports_multiple():
    r=d.inventory_to_rows(inventory(stations=('QS1','XQS1','QT1','AA1')))
    assert {x.station for x in filter_station_rows(r,station_prefix='qs*,qt')}=={'QS1','QT1'}


def test_unknown_coordinates_retained_without_box_excluded_with_box():
    r=d.inventory_to_rows(inventory(stations=('QS1',)))[0]
    rows=[replace(r,latitude=None),replace(r,latitude=np.nan),replace(r,latitude=48)]
    assert len(filter_station_rows(rows,station_prefix='QS'))==3
    assert len(filter_station_rows(rows,station_prefix='QS',minlatitude=47))==1


@pytest.mark.parametrize('opts',[{'minlatitude':91},{'minlongitude':181},{'maxlongitude':'nan'},
    {'minlatitude':49,'maxlatitude':48},{'station_prefix':'Q?S'},{'maxlatitude':'north'}])
def test_invalid_filters_fail_before_query(opts):
    client=MetadataClient()
    with pytest.raises(ValueError):d.list_stations('TEST','XX',client=client,**opts)
    assert not client.calls


def test_dateline_two_queries_deduplicate_epochs():
    inv=inventory(stations=('QS1','QS2','QS3'))
    for station,lon in zip(inv[0].stations,[179.,-179.,0.]):
        station.longitude=lon
        for channel in station:channel.longitude=lon
    client=MetadataClient(inv)
    rows=d.list_stations('TEST','XX',client=client,minlongitude=170,maxlongitude=-170)
    assert len(client.calls)==2
    assert (client.calls[0]['minlongitude'],client.calls[0]['maxlongitude'])==(170,180)
    assert (client.calls[1]['minlongitude'],client.calls[1]['maxlongitude'])==(-180,-170)
    assert {r.station for r in rows}=={'QS1','QS2'} and len(rows)==2


def test_no_data_in_one_dateline_half_does_not_hide_other_half():
    inv=inventory(stations=('QS1',));inv[0][0].longitude=-179.;inv[0][0][0].longitude=-179.
    class Client(MetadataClient):
        def get_stations(self,**kwargs):
            if kwargs['minlongitude']>=0:raise FDSNNoDataException('empty western half')
            return super().get_stations(**kwargs)
    rows=d.list_stations('TEST','XX',client=Client(inv),minlongitude=170,maxlongitude=-170)
    assert len(rows)==1 and rows[0].longitude==-179.


def test_network_discovery_forwards_prefix_and_bounds():
    client=MetadataClient()
    found=d.list_networks('TEST','*',client=client,station_prefix='QS',minlatitude=47,maxlatitude=49)
    assert found[0].code=='XX'
    assert client.calls[0]['level']=='network' and client.calls[0]['station']=='QS*'
    assert client.calls[0]['minlatitude']==47


def test_cli_catalogue_preserves_every_matching_station_and_epoch(tmp_path,monkeypatch):
    import importlib.util
    from pathlib import Path
    spec=importlib.util.spec_from_file_location('studio_cli',Path(__file__).resolve().parents[1]/'run.py')
    run=importlib.util.module_from_spec(spec);spec.loader.exec_module(run)
    client=MetadataClient();monkeypatch.setattr(d,'make_client',lambda *_:client)
    monkeypatch.setattr(sys,'argv',['run.py','station-list','--provider','TEST','--network','XX',
        '--station-prefix','QS','--min-latitude','47','--max-latitude','49',
        '--channel','*','--output',str(tmp_path)])
    assert run.main()==0
    with (tmp_path/'stations.csv').open(encoding='utf-8-sig') as handle:rows=list(csv.DictReader(handle))
    assert len(rows)==20
    spec=json.loads((tmp_path/'station_query.json').read_text())
    assert len(spec['channel_epochs'])==60 and not spec['automatic_subset_selection_applied']


def test_array_selection_operates_only_on_filtered_stations(tmp_path,monkeypatch):
    from ndstudio.acquisition import select_array
    client=MetadataClient()
    result=select_array('TEST','XX',DATE,DATE+120,tmp_path/'selection.json',client=client,
                       size=20,station_prefix='QS',min_spacing_m=50,max_spacing_m=500,max_aperture_m=4000)
    assert len(result['stations'])==20
    assert all(r['station'].startswith('QS') for r in result['stations'])
    assert len(result['matching_station_catalog'])==20
    assert result['station_filters']['station_prefix']=='QS'
