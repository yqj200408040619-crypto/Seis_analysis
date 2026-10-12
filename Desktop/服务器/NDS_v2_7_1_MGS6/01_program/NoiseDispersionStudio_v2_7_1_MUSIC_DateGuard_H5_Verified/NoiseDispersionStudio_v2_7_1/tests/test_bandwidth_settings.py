"""Band gates, preparation compatibility and historical-result provenance."""
from copy import deepcopy
from pathlib import Path
import json
import numpy as np
import pandas as pd
import pytest
from dispersion.common import DEFAULTS,bands,ensure_ready,same_preparation_config,sha,write_json
from dispersion.scan import track_peaks
from ndstudio.scan_diagnostics import write_scan_diagnostics


def test_bandwidth_change_alters_tracks_and_bands_including_exact_boundary():
    c=deepcopy(DEFAULTS);c['scan'].update(min_band_hz=.1,min_points=8)
    f=np.round(.27+.005*np.arange(15),9)
    nodes=[[dict(arrival=2.,snr=100.,edge=False)] for _ in f]
    assert not track_peaks(nodes,f,c) and not bands(np.ones(15,bool),f,8,.1)
    c['scan']['min_band_hz']=.05
    assert track_peaks(nodes,f,c) and bands(np.ones(15,bool),f,8,.05)
    c['scan']['min_band_hz']=.07
    assert track_peaks(nodes,f,c) and bands(np.ones(15,bool),f,8,.07)
    assert not bands(np.ones(15,bool),f,16,.05)


def test_only_bandwidth_and_point_count_can_reuse_preparation():
    c=deepcopy(DEFAULTS);d=deepcopy(c);d['scan'].update(min_band_hz=.05,min_points=12)
    assert same_preparation_config(c,d)
    for section,key,val in [(None,'window_s',900),('scan','step_hz',.005),('scan','arrival_max_s',28)]:
        q=deepcopy(d)
        if section:q[section][key]=val
        else:q[key]=val
        assert not same_preparation_config(c,q)


def test_new_scan_snapshot_overrides_old_preparation_for_diagnostics(tmp_path):
    c=deepcopy(DEFAULTS);c['output']=str(tmp_path);c['scan'].update(step_hz=.005,min_band_hz=.1)
    (tmp_path/'cache').mkdir();write_json(tmp_path/'cache/provenance.json',dict(config=c))
    c['scan']['min_band_hz']=.05
    pd.DataFrame(dict(gate_id=['P1']*15,frequency_hz=.27+.005*np.arange(15),robust_point=[True]*15)).to_csv(tmp_path/'point_diagnostics.csv',index=False)
    pd.DataFrame([dict(gate_id='P1',pair='A--B',cohort='all',branch='positive')]).to_csv(tmp_path/'candidate_gates.csv',index=False)
    write_json(tmp_path/'scan_complete.json',dict(config=c,point_diagnostics_sha256=sha(tmp_path/'point_diagnostics.csv'),candidate_gates_sha256=sha(tmp_path/'candidate_gates.csv')))
    result=write_scan_diagnostics(c)
    assert result['minimum_bandwidth_hz']==.05 and result['accepted_continuous_runs']==1
    assert result['threshold_source']=='scan_completion' and result['thresholds_match_current_config']
    c['scan']['min_band_hz']=.01;result=write_scan_diagnostics(c)
    assert result['minimum_bandwidth_hz']==.05 and not result['thresholds_match_current_config']


def test_downstream_rejects_old_scan_after_bandwidth_change(tmp_path):
    import dispersion.common as common
    c=deepcopy(DEFAULTS);c['output']=str(tmp_path)
    cache=tmp_path/'cache';cache.mkdir()
    write_json(cache/'provenance.json',dict(config=c,code_sha256={p.name:sha(p) for p in sorted(Path(common.__file__).parent.glob('*.py'))},sources=[]))
    write_json(tmp_path/'prepare_complete.json',{})
    write_json(tmp_path/'scan_complete.json',dict(config=c))
    newer=deepcopy(c);newer['scan']['min_band_hz']=.05
    ensure_ready(newer,'prepare_complete.json')
    with pytest.raises(ValueError,match='重新扫描'):ensure_ready(newer,'scan_complete.json')
    write_json(tmp_path/'scan_complete.json',dict(config=newer))
    ensure_ready(newer,'scan_complete.json')
