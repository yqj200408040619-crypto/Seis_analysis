"""Reject mixed acquisition periods before downloads or cached processing."""
from pathlib import Path
import json
import h5py
import pandas as pd
from dispersion.common import write_json


def utc(value):
    if value is None or str(value).strip() == '':
        raise ValueError('缺少数据时段，不能确认目录归属')
    t = pd.Timestamp(value)
    if pd.isna(t):
        raise ValueError('无效数据时段')
    return t.tz_localize('UTC') if t.tzinfo is None else t.tz_convert('UTC')


def window(start, end):
    a, b = utc(start), utc(end)
    if b <= a:
        raise ValueError('结束时间必须晚于开始时间')
    return [a.isoformat(), b.isoformat()]


def check_window(values, expected, source):
    if len(values) != 2 or window(*values) != expected:
        raise ValueError(f'时段防混用：{source} 与请求时段不一致；请使用新的工作目录')


def validate_download(root, start, end):
    """Check request provenance AND actual trace epochs; preserve FDSN padding."""
    root = Path(root)
    expected = window(start, end)
    for path in (root/'array_selection.json',):
        if path.exists():
            r = json.loads(path.read_text(encoding='utf-8-sig'))
            check_window([r.get('start_utc'), r.get('end_utc')], expected, path)
    manifest = root/'download_manifest.json'
    if manifest.exists():
        for job in json.loads(manifest.read_text(encoding='utf-8-sig'))['jobs']:
            check_window([job.get('start_utc'), job.get('end_utc')], expected, manifest)
    records = []
    for path in sorted(root.rglob('*.h5')):
        with h5py.File(path, 'r') as h:
            if 'metadata/download_request' not in h:
                raise ValueError(f'时段防混用：H5缺少下载请求 {path}')
            req = json.loads(h['metadata/download_request'][()])
            check_window([req.get('requested_start_utc'), req.get('requested_end_utc')], expected, path)
            a, b = utc(req.get('chunk_start_utc')), utc(req.get('chunk_end_utc'))
            if not utc(start) <= a < b <= utc(end):
                raise ValueError(f'时段防混用：分块超出请求范围 {path}')
            traces = 0
            for ds in h.values():
                if not isinstance(ds, h5py.Dataset) or ds.attrs.get('role') != 'waveform':
                    continue
                lo = pd.Timestamp(int(ds.attrs['start_utc_ns']), tz='UTC')
                hi = pd.Timestamp(int(ds.attrs['end_utc_ns']), tz='UTC')
                # FDSN returns complete MiniSEED records around boundaries.
                padding = pd.Timedelta(seconds=60)
                if hi < lo or hi < a or lo > b or lo < a-padding or hi > b+padding:
                    raise ValueError(f'时段防混用：实际波形时间不属于下载分块 {path}')
                if ds.attrs.get('seed_id') != req.get('seed_id'):
                    raise ValueError(f'时段防混用：波形台站与请求不一致 {path}')
                traces += 1
            if not traces:
                raise ValueError(f'时段防混用：没有波形 {path}')
            records.append(dict(file=str(path.relative_to(root)), traces=traces))
    return dict(requested_utc=expected, h5_files=len(records), records=records,
                requests_and_actual_epochs_checked=True, fdsn_boundary_padding_limit_s=60)


def bind_workdir(root, start, end):
    root = Path(root)
    expected = window(start, end)
    marker = root/'run_period.json'
    if marker.exists():
        check_window(json.loads(marker.read_text(encoding='utf-8-sig'))['requested_utc'], expected, marker)
    for path in (root/'array_selection.json', root/'project/project.json'):
        if path.exists():
            r = json.loads(path.read_text(encoding='utf-8-sig'))
            values = r.get('requested_utc') if path.name == 'project.json' else [r.get('start_utc'), r.get('end_utc')]
            if values is not None:
                check_window(values, expected, path)
    validate_download(root/'download', start, end)
    for path in (root/'project').rglob('aligned_waveforms.h5'):
        with h5py.File(path, 'r') as h:
            a, b = utc(h.attrs['start_utc']), utc(h.attrs['end_utc_exclusive'])
            padding = pd.Timedelta(seconds=60)
            if a < utc(start)-padding or b > utc(end)+padding or b <= utc(start) or a >= utc(end):
                raise ValueError(f'时段防混用：校准缓存不属于本次请求 {path}')
    write_json(marker, dict(schema='ndstudio.run_period/1', requested_utc=expected))
