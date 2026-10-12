"""CLI configuration precedence and single-station channel defaults."""
import json
import sys
import importlib.util
from pathlib import Path
spec=importlib.util.spec_from_file_location("studio_cli",Path(__file__).resolve().parents[1]/"run.py")
run=importlib.util.module_from_spec(spec);spec.loader.exec_module(run)
from ndstudio import hvsr,hvsr_acquisition


def test_cli_date_overrides_settings_and_resolves_relative_source(tmp_path,monkeypatch):
    p=tmp_path/'hvsr_settings.json'
    p.write_text(json.dumps(dict(folder='../download',instrument_id='BW.MGS01..HH',
                                settings=dict(start='2020-01-01',end='2020-01-02'))))
    seen={}
    def fake(folder,output,instrument_id,config):
        seen.update(folder=folder,output=output,instrument_id=instrument_id,config=config)
        return dict(instrument_id=instrument_id,selected_peak=None,structure={})
    monkeypatch.setattr(hvsr,'run_hvsr',fake)
    monkeypatch.setattr(sys,'argv',['run.py','hvsr','--project',str(tmp_path),'--start','2020-11-01','--end','2020-11-02'])
    assert run.main()==0
    assert seen['config']['start']=='2020-11-01' and seen['config']['end']=='2020-11-02'
    assert seen['folder']==tmp_path/'../download'


def test_single_station_download_defaults_to_three_components(monkeypatch,tmp_path):
    seen={}
    def fake(*args,**kwargs):seen.update(kwargs);return {}
    monkeypatch.setattr(hvsr_acquisition,'fetch_single_station',fake)
    monkeypatch.setattr(sys,'argv',['run.py','hvsr-fetch','--provider','LMU','--network','BW','--station','MGS01',
                                   '--start','2020-11-01','--end','2020-11-02','--output',str(tmp_path)])
    assert run.main()==0 and seen['channel']=='HH?'
