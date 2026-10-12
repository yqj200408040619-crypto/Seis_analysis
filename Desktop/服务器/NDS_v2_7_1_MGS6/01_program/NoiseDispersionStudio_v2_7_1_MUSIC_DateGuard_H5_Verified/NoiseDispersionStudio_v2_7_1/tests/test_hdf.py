import json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import pytest
from ndstudio.hdf import inspect_folder,read_waveform,utc,truth,MappingError
from ndstudio.project import new_project,import_hdf

def record(path,station='A',x=None,start='2026-01-01T00:00:00.123456789Z',**attrs):
    x=np.arange(20.,dtype=float) if x is None else x
    with h5py.File(path,'w') as h:
        h.attrs.update(station=station,channel='BHZ',fs=20.,start_utc=start,units='m/s',x_m=0 if station=='A' else 1000,y_m=0)
        h.attrs.update(attrs);h.create_dataset('waveform',data=x)
    return path

def project(tmp_path):
    root=new_project(tmp_path/'project');s=json.loads((root/'project.json').read_text(encoding='utf-8-sig'));s.update(input_folder='../h5',automatic_geometry=False)
    (root/'project.json').write_text(json.dumps(s));raw=tmp_path/'h5';raw.mkdir();return root,raw

def test_utc_precision_and_mask_are_preserved(tmp_path):
    path=record(tmp_path/'A.h5')
    with h5py.File(path,'a') as h:
        mask=np.ones(20,bool);mask[7:9]=False;h.create_dataset('valid',data=mask);h['waveform'].attrs['valid_dataset']='/valid'
    data=inspect_folder(tmp_path);assert data['ready'];i=data['items'][0];m=i['resolved']
    assert m['start_utc']=='2026-01-01T00:00:00.123456789Z'
    x=read_waveform(i,m);assert np.isnan(x[7:9]).all();assert x[6]==6 and x[9]==9

@pytest.mark.parametrize('axis',[0,1])
def test_two_array_orders_extract_declared_vertical(tmp_path,axis):
    x=np.column_stack([np.arange(20)+10,np.arange(20)+20,np.arange(20)+30])
    record(tmp_path/'A.h5',x=x if axis==0 else x.T,components='BHN,BHE,BHZ')
    s=inspect_folder(tmp_path);assert s['ready'];i=s['items'][0];assert i['resolved']['sample_axis']==axis
    np.testing.assert_array_equal(read_waveform(i,i['resolved']),np.arange(20)+30)

def test_undeclared_order_and_missing_utc_do_not_use_file_time(tmp_path):
    record(tmp_path/'A.h5',x=np.ones((20,3)))
    assert not inspect_folder(tmp_path)['ready']
    with h5py.File(tmp_path/'A.h5','a') as h:h['waveform'].attrs['components']='BHN BHE BHZ';del h.attrs['start_utc']
    data=inspect_folder(tmp_path);assert not data['ready'];assert 'UTC' in data['items'][0]['issues'][0]

@pytest.mark.parametrize('fmt,value',[('unix_s',1767225600),('unix_ms',1767225600000),('unix_us',1767225600000000),('unix_ns',1767225600123456789)])
def test_numeric_time_units_must_be_explicit(fmt,value):
    with pytest.raises(MappingError):utc(value)
    assert utc(value,fmt).startswith('2026-01-01T00:00:00.')
    if fmt=='unix_ns':assert utc(value,fmt).endswith('123456789Z')

@pytest.mark.parametrize('value,expected',[('False',False),('true',True),(0,False),(1,True)])
def test_boolean_declarations(value,expected):assert truth(value) is expected

def test_mapping_rules_use_raw_attributes(tmp_path):
    record(tmp_path/'A.h5')
    with h5py.File(tmp_path/'A.h5','a') as h:del h.attrs['fs'];h['waveform'].attrs['my_rate']=20.
    assert not inspect_folder(tmp_path)['ready']
    profile={'rules':[{'file_pattern':'*.h5','dataset_pattern':'/waveform','values':{'attribute_map':{'sampling_rate_hz':'my_rate'}}}]}
    assert inspect_folder(tmp_path,profile)['ready']

def test_identical_overlap_is_trimmed_and_native_time_kept(tmp_path):
    root,raw=project(tmp_path);record(raw/'a1.h5',x=np.arange(20.),start='2026-01-01T00:00:00Z')
    record(raw/'a2.h5',x=np.arange(10.,30.),start='2026-01-01T00:00:00.5Z');record(raw/'b.h5',station='B')
    import_hdf(root)
    with np.load(root/'canonical/A_00001.npz') as z:
        np.testing.assert_array_equal(z['data'],np.arange(20.,30.));assert z['start_utc'].item()=='2026-01-01T00:00:01.000000000Z'
    audit=json.loads((root/'hdf_import_audit.json').read_text(encoding='utf-8-sig'));assert audit['records'][1]['identical_overlap_removed_samples']==10

def test_conflicting_overlap_is_rejected(tmp_path):
    root,raw=project(tmp_path);record(raw/'a1.h5',start='2026-01-01T00:00:00Z')
    record(raw/'a2.h5',x=np.arange(10.,30.)+1,start='2026-01-01T00:00:00.5Z');record(raw/'b.h5',station='B')
    with pytest.raises(MappingError,match='重叠'):import_hdf(root)

def test_generic_hdf_does_not_inherit_mseed_header_tolerance(tmp_path):
    root,raw=project(tmp_path)
    record(raw/'a1.h5',x=np.arange(20.),start='2026-01-01T00:00:00Z')
    record(raw/'a2.h5',x=np.arange(10.,30.),start='2026-01-01T00:00:00.500002Z')
    record(raw/'b.h5',station='B')
    with pytest.raises(MappingError,match='重叠'):import_hdf(root)

def test_counts_without_response_are_not_used(tmp_path):
    root,raw=project(tmp_path);record(raw/'a.h5',units='counts');record(raw/'b.h5',station='B')
    with pytest.raises(MappingError,match='StationXML'):import_hdf(root)

def test_wrong_vertical_and_unknown_units_are_rejected(tmp_path):
    record(tmp_path/'a.h5',channel='BHN',units='volts')
    assert not inspect_folder(tmp_path)['ready']

def test_rate_change_is_not_silently_joined(tmp_path):
    root,raw=project(tmp_path);record(raw/'a1.h5',start='2026-01-01T00:00:00Z')
    record(raw/'a2.h5',start='2026-01-01T00:00:02Z',fs=40.);record(raw/'b.h5',station='B')
    with pytest.raises(MappingError,match='采样率'):import_hdf(root)

def test_fractional_component_index_is_not_truncated(tmp_path):
    record(tmp_path/'A.h5',x=np.ones((20,3)))
    scan=inspect_folder(tmp_path,{'defaults':{'component_index':1.5}})
    assert not scan['ready'] and '整数' in scan['items'][0]['issues'][0]

def test_unused_nan_attributes_do_not_break_inventory(tmp_path):
    record(tmp_path/'A.h5',unused_calibration=np.nan,latitude=np.nan,longitude=np.nan)
    scan=inspect_folder(tmp_path);assert scan['ready'];json.dumps(scan,allow_nan=False)

def test_explicit_null_mapping_can_remove_wrong_identity(tmp_path):
    record(tmp_path/'A.h5',seed_id='XX.BAD..BHZ')
    scan=inspect_folder(tmp_path,{'defaults':{'seed_id':None}})
    assert scan['ready'] and scan['items'][0]['resolved']['station']=='A'
