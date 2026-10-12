"""Actual Tk widgets and asynchronous station queries, without external services."""
from dataclasses import replace
import csv
import pytest
import h5data as d
import h5data_gui
from test_download_gui import tkroot,pump
from test_station_filters import MetadataClient
from test_download_handoff import DATE


@pytest.fixture
def filter_app(tkroot,tmp_path):
    client=MetadataClient();errors=[]
    app=h5data_gui.H5DataApp(tkroot,client_factory=lambda *_:client,
        ask_confirm=lambda *_a,**_kw:True,show_error=lambda *a,**kw:errors.append(a),use_matplotlib=True)
    app.pack(fill='both',expand=True)
    app.start_var.set(str(DATE));app.end_var.set(str(DATE+120));app.chunk_var.set('60')
    app.outdir_var.set(str(tmp_path/'download'));app._requested_array_code='XX'
    app._show_preview('XX',d.inventory_to_rows(client.inv))
    yield app,client,errors
    app.destroy()


def test_array_filter_lists_20_and_downloads_only_20(tkroot,filter_app,monkeypatch,tmp_path):
    app,client,errors=filter_app
    app.station_prefix_var.set('QS')
    app.geo_vars['minlatitude'].set('47');app.geo_vars['maxlatitude'].set('49')
    assert not app.array_rows and not app.preview_tree.get_children()
    app.apply_station_filters()
    pump(tkroot,lambda:len(app.array_rows)==60 and 'calculating' not in app.spacing_var.get())
    assert len(app.preview_tree.get_children())==60
    assert '20 matching stations' in app.array_summary_var.get()
    assert len(app.station_map.points)==20
    captured=[]
    monkeypatch.setattr(app,'_start_download',lambda jobs,*args:captured.extend(jobs))
    app.download_array()
    assert len(captured)==20
    assert {j.row.station for j in captured}=={f'QS{i:02d}' for i in range(20)}
    assert all(j.row.channel=='HHZ' for j in captured)
    filename=tmp_path/'matching.csv'
    monkeypatch.setattr(h5data_gui.filedialog,'asksaveasfilename',lambda **kw:str(filename))
    app.export_stations('array')
    with filename.open(encoding='utf-8-sig') as handle:catalog=list(csv.DictReader(handle))
    assert len(catalog)==20 and all(r['station'].startswith('QS') for r in catalog)
    assert not errors
    old_token=app._preview_token
    app.station_prefix_var.set('ZZ')
    app._handle(('preview',old_token,'XX',d.inventory_to_rows(client.inv)))
    assert not app.array_rows and not app.preview_tree.get_children()


def test_channels_filter_and_stale_result_is_rejected(tkroot,filter_app,monkeypatch):
    app,client,errors=filter_app
    app.notebook.select(app.channels_tab);app.network_var.set('XX')
    app.station_prefix_var.set('QS');app.apply_station_filters()
    pump(tkroot,lambda:not app.busy and len(app.rows)==20)
    old_token=app._browse_token
    app.tree.selection_set(*app.tree.get_children())
    captured=[];monkeypatch.setattr(app,'_start_download',lambda jobs,*args:captured.extend(jobs))
    app.download()
    assert len(captured)==20 and all(j.row.station.startswith('QS') for j in captured)
    app.station_prefix_var.set('ZZ')
    app._handle(('rows',old_token,d.inventory_to_rows(client.inv)))
    assert not app.rows and not app.tree.get_children()
    app.apply_station_filters();pump(tkroot,lambda:not app.busy)
    assert not app.rows and '0 matching stations' in app.count_var.get()
    assert not errors


def test_every_matching_table_row_is_shown_not_capped(filter_app,monkeypatch):
    app,client,_=filter_app
    monkeypatch.setattr(app,'_request_geometry',lambda *args:None)
    row=d.inventory_to_rows(client.inv)[0]
    rows=[replace(row,station=f'S{i:05d}') for i in range(h5data_gui.MAX_TREE_ROWS+3)]
    app._show_rows(rows)
    assert len(app.tree.get_children())==len(rows)
    app._show_preview('XX',rows)
    assert len(app.preview_tree.get_children())==len(rows)
