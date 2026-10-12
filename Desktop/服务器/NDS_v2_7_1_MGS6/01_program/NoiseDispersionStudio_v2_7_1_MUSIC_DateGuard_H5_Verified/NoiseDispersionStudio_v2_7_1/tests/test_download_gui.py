"""Real Tk widgets, map backends and threaded download/handoff actions."""
import json
import os
import time
import tkinter as tk
from pathlib import Path
import pytest
import h5data
import h5data_gui
from test_download_handoff import FakeFDSN, inventory, DATE


@pytest.fixture(scope='module')
def tkroot():
    try:
        root=tk.Tk()
    except tk.TclError as exc:
        pytest.skip('Tk display unavailable: '+str(exc))
    root.geometry('1550x1080')
    yield root
    root.destroy()


def pump(root, condition, timeout=15):
    deadline=time.monotonic()+timeout
    while not condition() and time.monotonic()<deadline:
        root.update();time.sleep(.01)
    root.update()
    assert condition(), 'GUI action did not finish'


@pytest.fixture
def app(tkroot,tmp_path):
    client=FakeFDSN(inventory(channels=('HHZ','HHN','HHE')))
    frame=h5data_gui.H5DataApp(tkroot,client_factory=lambda *_a,**_k:client,
                             ask_confirm=lambda *_a,**_k:True,use_matplotlib=True,
                             show_error=lambda *_a,**_k:None)
    frame.pack(fill='both',expand=True)
    frame.start_var.set(str(DATE));frame.end_var.set(str(DATE+120));frame.chunk_var.set('60')
    frame.outdir_var.set(str(tmp_path/'download'))
    frame._show_preview('XX',h5data.inventory_to_rows(client.inv))
    pump(tkroot,lambda:'calculating' not in frame.spacing_var.get())
    frame._requested_array_code = 'XX'
    yield frame
    frame.destroy()


def test_filtered_array_preview_matches_download_stations(tkroot,app):
    assert app.convert_var.get()
    assert app.array_chan_var.get()=='HHZ'
    assert '111.' in app.spacing_var.get() and ' m ' in app.spacing_var.get()
    assert len(app.station_map.points)==3
    assert app.station_map.ax.get_xlabel()=='East (m)'
    app.array_chan_var.set('NONE')
    pump(tkroot,lambda:'calculating' not in app.spacing_var.get())
    assert not app.station_map.points
    assert '--' in app.spacing_var.get()


def test_channel_tab_has_selection_map_and_meter_spacing(tkroot,app):
    app._show_rows(h5data.inventory_to_rows(inventory()))
    app.tree.selection_set('0','1')
    app._update_channel_preview()
    pump(tkroot,lambda:'calculating' not in app.channel_spacing_var.get())
    assert len(app.channel_map.points)==2
    assert ' m ' in app.channel_spacing_var.get()


def test_real_download_button_worker_h5_and_open_project(tkroot,app,monkeypatch):
    app.download_array()
    pump(tkroot,lambda:not app.busy,30)
    folder=Path(app.outdir_var.get());assert len(list(folder.rglob('*.h5')))==9
    assert 'Done' in app.status_var.get()
    calls=[];monkeypatch.setattr(h5data_gui.subprocess,'Popen',lambda *a,**kw:calls.append((a,kw)))
    app.open_dispersion();assert calls
    root=Path(calls[0][0][0][-1]);spec=json.loads((root/'project.json').read_text())
    assert Path(spec['input_folder'])==folder.resolve()
    assert calls[0][0][0][-3:] == ['gui','--project',str(root)]


def test_canvas_fallback_and_dateline_points(tkroot):
    frame=h5data_gui.StationMap(tkroot,use_matplotlib=False);frame.pack()
    frame.plot([('XX','A',0,179.999,0),('XX','B',0,-179.999,0)],title='Dateline')
    tkroot.update()
    assert len(frame.canvas.find_withtag('station'))==2
    frame.destroy()


def test_provider_change_clears_old_selection(tkroot,app):
    app.provider_var.set('OTHER')
    assert not app.array_rows and not app.rows and not app.networks
    pump(tkroot,lambda:'calculating' not in app.spacing_var.get())
    assert not app.station_map.points
