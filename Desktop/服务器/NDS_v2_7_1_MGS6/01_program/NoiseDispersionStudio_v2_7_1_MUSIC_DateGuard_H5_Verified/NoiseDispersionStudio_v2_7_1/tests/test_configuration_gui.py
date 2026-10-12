"""Actual Qt saves/step buttons and a real native scan child process."""
from pathlib import Path
import json,os,time
import numpy as np
import pandas as pd
import pytest
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
from PySide6.QtWidgets import QApplication,QMessageBox
from PySide6.QtCore import QProcess
from ndstudio.gui import MainWindow
from ndstudio.project import new_project
from dispersion.common import load_config,write_json,sha


@pytest.fixture
def view(tmp_path,monkeypatch):
    app=QApplication.instance() or QApplication([]);errors=[]
    monkeypatch.setattr(QMessageBox,'warning',lambda *a:errors.append(a[-1]))
    root=new_project(tmp_path/'project');window=MainWindow(root);window.show()
    yield window,root,app,errors
    window.close()


def edit_json(window,**changes):
    c=json.loads(window.configeditor.toPlainText());c['scan'].update(changes)
    window.configeditor.setPlainText(json.dumps(c));window.save_advanced()


def test_advanced_json_saves_syncs_and_survives_startup_save(view):
    w,root,app,errors=view
    edit_json(w,min_band_hz=.05,min_points=12,step_hz=.00001)
    assert w.configfields['scan.min_band_hz'].value()==pytest.approx(.05)
    assert w.configfields['scan.min_points'].value()==12
    w.persist_common();c=load_config(root/'config.json')
    assert c['scan']['min_band_hz']==.05 and c['scan']['step_hz']==.00001
    w.load(root);w.persist_common()
    assert load_config(root/'config.json')['scan']['step_hz']==.00001
    assert not errors


def test_step_buttons_use_frequency_fraction_and_duration_units(view):
    w,root,app,errors=view
    for key,delta in [('scan.min_band_hz',.005),('scan.step_hz',.001),('scan.min_hz',.005),
                      ('min_clean_fraction',.01),('spikes.max_repair_s',.001)]:
        field=w.configfields[key];before=field.value();field.stepUp()
        assert field.value()-before==pytest.approx(delta)
    assert w.configfields['window_s'].singleStep()==1
    assert w.hvsr.fields['overlap'].singleStep()==.01
    assert w.phasetol.singleStep()==.1
    assert not errors


def test_saving_common_keeps_unsaved_json_edits_and_worker_geometry_updates(view):
    w,root,app,errors=view
    # An import worker may adjust geometry on disk while the editor remains open.
    c=json.loads((root/'config.json').read_text());c['scan']['arrival_max_s']=28
    write_json(root/'config.json',c)
    edited=json.loads(w.configeditor.toPlainText());edited['scan']['min_band_hz']=.05
    w.configeditor.setPlainText(json.dumps(edited))
    w.configfields['scan.step_hz'].setValue(.005);w.persist_common()
    c=load_config(root/'config.json')
    assert c['scan']['arrival_max_s']==28 and c['scan']['min_band_hz']==.05 and c['scan']['step_hz']==.005
    assert not errors


def prepared_scan_fixture(root):
    c=json.loads((root/'config.json').read_text())
    c.update(bandpass_hz=[.05,7],discard_first_s=0,window_s=60,min_clean_segment_s=15,max_lag_s=10)
    c['response']['pre_filt_hz']=[.005,.02,7.5,9]
    c['scan'].update(min_hz=.27,max_hz=.34,step_hz=.005,min_band_hz=.1,min_points=8,arrival_max_s=4,noise_lag_s=[6,8])
    c['cohorts'].update(clock_groups=False,states=0)
    c['validation'].update(block_s=180,min_windows_each_split=2,min_windows_per_block=2,min_eligible_blocks=2,bootstrap_draws=0)
    c['routes']=c['routes'][:2];write_json(root/'config.json',c);c=load_config(root/'config.json')
    out=root/'output';cache=out/'cache';cache.mkdir(parents=True,exist_ok=True)
    import dispersion.common as common
    write_json(cache/'provenance.json',dict(config=c,code_sha256={p.name:sha(p) for p in sorted(Path(common.__file__).parent.glob('*.py'))},sources=[],response_corrected=True))
    write_json(out/'prepare_complete.json',dict(test_fixture=True))
    # Controlled no-signal CCF cache: tests configuration propagation, not recovery.
    x=np.zeros((8,1,401));np.savez(cache/'window_routes.npz',ram5=x,ram20=x,
        starts=np.arange(8)*1200,timestamps=np.arange(8)*60,training=np.arange(8)<4)
    pd.DataFrame([dict(pair_index=0,pair='A--B',distance_m=100,azimuth_deg=90)]).to_csv(out/'pairs.csv',index=False)
    return c


def test_actual_qprocess_scan_uses_new_bandwidth_without_rewriting_preparation(view):
    w,root,app,errors=view;prepared_scan_fixture(root);w.load(root)
    edit_json(w,min_band_hz=.05)
    w.start('scan');assert w.current_job=='scan'
    deadline=time.monotonic()+20
    while w.process.state()!=QProcess.ProcessState.NotRunning and time.monotonic()<deadline:
        app.processEvents();time.sleep(.01)
    app.processEvents()
    assert w.process.state()==QProcess.ProcessState.NotRunning and w.process.exitCode()==0,w.buffer
    marker=json.loads((root/'output/scan_complete.json').read_text())
    diagnostic=json.loads((root/'output/scan_diagnostics_summary.json').read_text())
    assert marker['config']['scan']['min_band_hz']==.05
    assert diagnostic['minimum_bandwidth_hz']==.05 and diagnostic['threshold_source']=='scan_completion'
    assert diagnostic['thresholds_match_current_config']
    assert json.loads((root/'output/cache/provenance.json').read_text())['config']['scan']['min_band_hz']==.1
    assert '最小连续带宽 0.05 Hz' in w.buffer and not errors



def test_advanced_save_preserves_worker_changes_not_edited_in_json(view):
    w,root,app,errors=view
    c=json.loads((root/'config.json').read_text());c['scan']['arrival_max_s']=28
    write_json(root/'config.json',c)
    edit_json(w,min_band_hz=.05)
    c=load_config(root/'config.json')
    assert c['scan']['arrival_max_s']==28 and c['scan']['min_band_hz']==.05
    assert w.config_snapshot['scan']['arrival_max_s']==28 and not errors
