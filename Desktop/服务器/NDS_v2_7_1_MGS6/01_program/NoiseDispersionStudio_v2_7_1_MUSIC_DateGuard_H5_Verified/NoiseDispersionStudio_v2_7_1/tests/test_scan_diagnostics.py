import copy
import json
import numpy as np
import pandas as pd
import pytest
from dispersion.common import DEFAULTS
from ndstudio.scan_diagnostics import write_scan_diagnostics

def fixture(tmp_path,frequencies,flags):
    c=copy.deepcopy(DEFAULTS);c['output']=str(tmp_path)
    c['scan'].update(step_hz=.005,min_band_hz=.1,min_points=8)
    pd.DataFrame({'gate_id':['P1']*len(frequencies),'frequency_hz':frequencies,'robust_point':flags}).to_csv(tmp_path/'point_diagnostics.csv',index=False)
    pd.DataFrame([{'gate_id':'P1','cohort':'all','pair':'A--B','branch':'positive'}]).to_csv(tmp_path/'candidate_gates.csv',index=False)
    return c

def test_short_band_with_all_points_passing_is_reported(tmp_path,capsys):
    c=fixture(tmp_path,.27+.005*np.arange(15),[True]*15)
    result=write_scan_diagnostics(c)
    q=pd.read_csv(tmp_path/'scan_band_diagnostics.csv')
    assert result['robust_contiguous_runs']==1 and result['accepted_continuous_runs']==0
    assert q.points.iloc[0]==15 and q.bandwidth_hz.iloc[0]==pytest.approx(.07)
    assert q.all_point_checks_pass.iloc[0] and not q.passes_minimum_bandwidth.iloc[0]
    assert '带宽不足(0.070<0.100 Hz)' in capsys.readouterr().out

def test_frequency_gap_and_false_string_split_intervals(tmp_path):
    c=fixture(tmp_path,[.27,.275,.285,.29,.295,.3],['true','true','true','true','False','true'])
    write_scan_diagnostics(c)
    q=pd.read_csv(tmp_path/'scan_band_diagnostics.csv')
    assert sorted(q.points.tolist())==[1,2,2]
    assert q.bandwidth_hz.max()==pytest.approx(.005)

def test_empty_scan_still_exports_readable_diagnostic_headers(tmp_path):
    c=fixture(tmp_path,[],[])
    result=write_scan_diagnostics(c)
    assert result['robust_contiguous_runs']==result['candidate_gates']==0
    assert 'reason' in pd.read_csv(tmp_path/'scan_band_diagnostics.csv')

def test_edited_config_does_not_relabel_historical_scan_as_accepted(tmp_path):
    c=fixture(tmp_path,.27+.005*np.arange(15),[True]*15)
    (tmp_path/'cache').mkdir();(tmp_path/'cache/provenance.json').write_text(json.dumps({'config':c}))
    c['scan']['min_band_hz']=.05
    result=write_scan_diagnostics(c)
    assert not result['thresholds_match_current_config'] and result['accepted_continuous_runs']==0

def test_candidate_bottlenecks_report_nan_and_phase_failures(tmp_path):
    c=fixture(tmp_path,[.27,.275,.28],[True,False,False])
    q=pd.read_csv(tmp_path/'point_diagnostics.csv')
    q['gate_change_rad']=[.1,.9,np.nan];q.to_csv(tmp_path/'point_diagnostics.csv',index=False)
    write_scan_diagnostics(c)
    gates=pd.read_csv(tmp_path/'scan_candidate_diagnostics.csv')
    assert gates.failed_gate_change_rad.iloc[0]==2
    assert '门窗扰动相位:2/3点' in gates.main_bottlenecks.iloc[0]
