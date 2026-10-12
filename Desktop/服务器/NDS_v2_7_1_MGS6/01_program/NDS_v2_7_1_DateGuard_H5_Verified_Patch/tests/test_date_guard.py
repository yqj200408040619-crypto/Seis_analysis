import json
import h5py
import pandas as pd
import pytest
from ndstudio.date_guard import bind_workdir, validate_download

START='2023-11-01T00:00:00Z'
END='2023-11-05T00:00:00Z'


def archive(folder, *, request_start=START, actual_start=START, seed='BW.MGS01..HHZ'):
    folder.mkdir(parents=True, exist_ok=True)
    path=folder/'trace.h5'
    with h5py.File(path,'w') as h:
        req=dict(requested_start_utc=request_start, requested_end_utc=END,
                 chunk_start_utc=START,chunk_end_utc='2023-11-01T06:00:00Z',seed_id=seed)
        h.create_dataset('metadata/download_request',data=json.dumps(req))
        ds=h.create_dataset('trace',data=[1,2])
        ns=pd.Timestamp(actual_start).value
        ds.attrs.update(role='waveform',start_utc_ns=ns,end_utc_ns=ns+3*10**9,seed_id=seed)
    return path


def test_orphan_old_h5_cannot_enter_fresh_workdir(tmp_path):
    archive(tmp_path/'download',request_start='2022-11-01T00:00:00Z')
    with pytest.raises(ValueError,match='时段防混用'):
        bind_workdir(tmp_path,START,END)
    assert not (tmp_path/'run_period.json').exists()


def test_forged_request_cannot_hide_old_actual_waveform(tmp_path):
    archive(tmp_path,actual_start='2022-11-01T00:00:00Z')
    with pytest.raises(ValueError,match='实际波形'):
        validate_download(tmp_path,START,END)


def test_boundary_padding_and_equivalent_timezone_are_valid(tmp_path):
    archive(tmp_path/'download',actual_start='2023-10-31T23:59:58Z')
    bind_workdir(tmp_path,'2023-11-01T01:00:00+01:00',END)
    assert validate_download(tmp_path/'download',START,END)['h5_files']==1
    with pytest.raises(ValueError,match='时段防混用'):
        bind_workdir(tmp_path,'2022-11-01', '2022-11-05')


def test_manifest_and_project_wrong_dates_are_rejected(tmp_path):
    d=tmp_path/'download';d.mkdir()
    (d/'download_manifest.json').write_text(json.dumps(dict(jobs=[dict(start_utc='2022-11-01',end_utc='2022-11-05')])))
    with pytest.raises(ValueError,match='时段防混用'):
        bind_workdir(tmp_path,START,END)
    (d/'download_manifest.json').unlink()
    p=tmp_path/'project';p.mkdir()
    (p/'project.json').write_text(json.dumps(dict(requested_utc=['2022-11-01','2022-11-05'])))
    with pytest.raises(ValueError,match='时段防混用'):
        bind_workdir(tmp_path,START,END)


def test_old_aligned_cache_is_rejected_even_with_current_raw_data(tmp_path):
    archive(tmp_path/'download')
    p=tmp_path/'project/output/music_front';p.mkdir(parents=True)
    with h5py.File(p/'aligned_waveforms.h5','w') as h:
        h.attrs.update(start_utc='2022-11-01',end_utc_exclusive='2022-11-05')
    with pytest.raises(ValueError,match='校准缓存'):
        bind_workdir(tmp_path,START,END)


@pytest.mark.parametrize('value',[None,'','NaT'])
def test_missing_time_never_becomes_today(tmp_path,value):
    with pytest.raises(ValueError):
        bind_workdir(tmp_path,value,END)
