"""Regression: arbitrary curve counts and project metadata in phase reports."""
import json
import re
import numpy as np
import pandas as pd
import pytest
from ndstudio.phase_velocity import figures,run_phase_velocity


def sample_curves(count=1):
    summary=[];targets=[]
    for cycle in range(1,count+1):
        for f in (.27,.30,.34):
            velocity=360*f*1000/(40+360*cycle)
            summary.append(dict(curve_id=f'test_N{cycle}',frequency_hz=f,segment_id='test',
                pair='A--B',branch='positive',after_count=2,after_raw_phase_deg=20.,
                after_unit_phase_deg=20.,predicted_phase_deg=20.,phase_offset_low_deg=-5.,
                phase_offset_high_deg=5.,accepted_windows=2,valid_windows=5))
            targets.append(dict(curve_id=f'test_N{cycle}',frequency_hz=f,constant_cycle=cycle,
                source_phase_rad=np.pi/3,distance_m=1000.,predicted_unwrapped_phase_rad=np.radians(20.),
                phase_velocity_m_s=velocity,velocity_min_m_s=velocity*.99,velocity_max_m_s=velocity*1.01))
    return pd.DataFrame(summary),pd.DataFrame(targets)


def test_nine_curves_all_visible(tmp_path,monkeypatch):
    from ndstudio.phase_velocity import calculate
    import matplotlib.pyplot as plt
    s,t=sample_curves(9);details,_=calculate(s,t)
    captured=[]
    def close(fig):
        if hasattr(fig,'axes'):captured.append(fig)
    original_close=plt.close
    monkeypatch.setattr(plt,'close',close)
    figures(details,tmp_path)
    shown=[ax.get_title().split('\n')[0] for ax in captured[0].axes if ax.axison]
    assert shown==[f'test_N{n}' for n in range(1,10)]
    assert len(captured[0].axes)>=10
    for fig in captured:original_close(fig)


def test_report_uses_actual_source_and_period(tmp_path):
    folder=tmp_path/'stack';folder.mkdir();output=tmp_path/'comparison'
    s,t=sample_curves();s.to_csv(folder/'phase_stack_summary.csv',index=False)
    t.to_csv(folder/'selected_phase_targets.csv',index=False)
    (folder/'phase_stack_complete.json').write_text(json.dumps(dict(project_name='BW.MGS test',
        start_utc='2020-11-01',end_utc_exclusive='2020-11-05',range_source='measured',
        tolerance_deg=0,repeat_cycles=20,slide_cycles=1)))
    audit=run_phase_velocity(folder,output,verify_windows=False)
    report=(output/'report.html').read_text()
    report=re.sub(r'<img[^>]*>','',report)
    assert 'BW.MGS test' in report and '2020-11-01' in report
    assert 'DARE' not in report and '301 个曲线频点' not in report
    assert audit['source_phase_deg']==pytest.approx(60.)
    assert len(pd.read_csv(output/'curve_comparison.csv'))==3
    s.loc[1,'after_count']=0;s.to_csv(folder/'phase_stack_summary.csv',index=False)
    with pytest.raises(ValueError,match='Empty stack'):run_phase_velocity(folder,output)
    assert not (output/'report.html').exists()
    assert not (output/'fv_curves.png').exists()


def test_incomplete_window_write_preserves_previous_archive(tmp_path,monkeypatch):
    import zipfile
    from ndstudio.phase_stack import save_windows
    path=tmp_path/'windows.npz'
    np.savez(path,cross=np.array([1+2j]))
    before=path.read_bytes()
    def partial(handle,**arrays):handle.write(b'PK incomplete data')
    monkeypatch.setattr(np,'savez_compressed',partial)
    with pytest.raises(zipfile.BadZipFile):save_windows(path,cross=np.array([3+4j]))
    assert path.read_bytes()==before
    assert not path.with_suffix('.npz.part').exists()
