"""Actual Qt actions and download handoff without a live FDSN service."""
import json
import os
import time
from types import SimpleNamespace
from pathlib import Path
import pytest
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
from PySide6.QtWidgets import QApplication
from PySide6.QtCore import QProcess
from ndstudio.gui import MainWindow
from ndstudio.project import new_project
from test_hvsr import simple


def test_qt_single_station_child_process_and_primary_button_do_not_run_pair_pipeline(tmp_path):
    app=QApplication.instance() or QApplication([])
    raw=tmp_path/'raw';simple(raw);root=new_project(tmp_path/'project')
    window=MainWindow(root);window.show();window.tabs.setCurrentWidget(window.hvsr)
    window.hvsr.folder.setText(str(raw));window.hvsr.scan()
    assert window.hvsr.instrument.currentData()=='XX.A..HH'
    window.hvsr.fields['fmin_hz'].setValue(.3);window.hvsr.fields['fmax_hz'].setValue(7)
    before=(root/'config.json').read_bytes()
    assert window.runbutton.text()=='一键 H/V 分析' and not window.autoinvert.isEnabled()
    window.primary_analysis()
    assert window.current_job=='hvsr'
    deadline=time.monotonic()+25
    while window.process.state()!=QProcess.ProcessState.NotRunning and time.monotonic()<deadline:
        app.processEvents();time.sleep(.01)
    app.processEvents()
    assert window.process.exitCode()==0
    assert (root/'output/hvsr/hvsr_complete.json').exists()
    assert (root/'config.json').read_bytes()==before
    assert not (root/'stations.csv').exists() and not (root/'canonical').exists()
    assert not window.hvsr.image.pixmap().isNull()
    saved=json.loads((root/'hvsr_settings.json').read_text());assert saved['instrument_id']=='XX.A..HH'
    window.close()


def test_download_button_opens_single_station_project_and_preserves_choice(tmp_path,monkeypatch):
    import h5data_gui
    raw=tmp_path/'download';simple(raw)
    calls=[];errors=[]
    monkeypatch.setattr(h5data_gui.subprocess,'Popen',lambda *args,**kwargs:calls.append((args,kwargs)))
    view=SimpleNamespace(_busy=False,outdir_var=SimpleNamespace(get=lambda:str(raw)),_log=lambda message:None,_error=lambda *e:errors.append(e))
    h5data_gui.H5DataApp.open_hvsr(view)
    assert not errors and len(calls)==1
    command=calls[0][0][0];assert command[-1]=='--hvsr-tab'
    root=Path(command[command.index('--project')+1]);cfg=root/'hvsr_settings.json'
    saved=json.loads(cfg.read_text());saved['settings']['window_s']=120;cfg.write_text(json.dumps(saved))
    h5data_gui.H5DataApp.open_hvsr(view)
    assert json.loads(cfg.read_text())['settings']['window_s']==120
    assert json.loads((root/'project.json').read_text())['input_folder']==str(raw.resolve())


def test_download_hv_handoff_gives_actionable_missing_component_error(tmp_path,monkeypatch):
    import h5data_gui
    raw=tmp_path/'vertical_only';simple(raw,channels=('HHZ',),locations=('',))
    errors=[];calls=[];monkeypatch.setattr(h5data_gui.subprocess,'Popen',lambda *a,**k:calls.append(a))
    view=SimpleNamespace(_busy=False,outdir_var=SimpleNamespace(get=lambda:str(raw)),_log=lambda message:None,_error=lambda *e:errors.append(e))
    h5data_gui.H5DataApp.open_hvsr(view)
    assert not calls and errors and 'three' not in errors[0][0].lower()
    assert 'HH?' in errors[0][1] and '3-component' in errors[0][1]
