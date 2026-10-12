"""Single-instrument native CLI acquisition contract and resume audit."""
import json
import pytest
from test_download_handoff import FakeFDSN,inventory,DATE
from ndstudio.hvsr_acquisition import fetch_single_station


class BrowseFakeFDSN(FakeFDSN):
    def get_stations(self,filename=None,**kwargs):
        if kwargs.get('level')=='channel':
            self.station_calls+=1
            return self.inv.copy()
        return super().get_stations(filename=filename,**kwargs)


def client():
    inv=inventory(stations=('A',),channels=('HHN','HHE','HHZ'),fs=20)
    inv[0][0][1].azimuth=90.
    return BrowseFakeFDSN(inv)


def test_native_single_station_fetch_and_resume_preserve_complete_archive(tmp_path):
    source=client()
    result=fetch_single_station('MOCK','XX','A',DATE,DATE+120,tmp_path,chunk_seconds=60,client=source)
    assert result['selection']['instrument_id']=='XX.A..HH'
    assert result['audit']['checks_passed'] and result['audit']['h5_files']==9
    assert result['download_status']['errors']==0
    before=source.waveform_calls
    resumed=fetch_single_station('MOCK','XX','A',DATE,DATE+120,tmp_path,chunk_seconds=60,client=source)
    assert source.waveform_calls==before
    assert resumed['download_status']['existing']==9
    assert (tmp_path/'hvsr_download_complete.json').exists()


def test_missing_component_is_rejected_before_any_waveform_request(tmp_path):
    source=BrowseFakeFDSN(inventory(stations=('A',),channels=('HHZ',),fs=20))
    with pytest.raises(ValueError,match='三分量'):fetch_single_station('MOCK','XX','A',DATE,DATE+120,tmp_path,client=source)
    assert source.waveform_calls==0


def test_request_across_metadata_epoch_is_not_silently_clipped(tmp_path):
    source=client();source.inv[0][0][0].end_date=DATE+90
    with pytest.raises(ValueError,match='历元'):fetch_single_station('MOCK','XX','A',DATE,DATE+120,tmp_path,client=source)
    assert source.waveform_calls==0
