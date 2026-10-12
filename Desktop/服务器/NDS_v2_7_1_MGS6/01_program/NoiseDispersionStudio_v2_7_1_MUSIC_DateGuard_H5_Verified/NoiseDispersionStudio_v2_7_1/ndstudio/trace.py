"""Readable stage logs and machine-readable provenance, including failures."""
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
import json, sys, time, uuid, math, re, shutil, shutil, shutil, shutil, shutil
from dispersion.common import write_json, sha


class Tee:
    def __init__(self, console, handle): self.console, self.handle = console, handle
    def write(self, text):
        self.handle.write(text); self.handle.flush()
        return self.console.write(text)
    def flush(self): self.handle.flush(); self.console.flush()
    def __getattr__(self, name): return getattr(self.console, name)


def clean(value):
    if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    if isinstance(value,float) and not math.isfinite(value):return None
    if isinstance(value,Path):return str(value)
    if hasattr(value,'tolist'):return clean(value.tolist())
    return value


def _repair_json_text(text):
    """Repair legacy JSON that contains raw Windows backslashes inside strings.

    Old/generated project markers from some Windows runs can contain paths such
    as C:\\Users\\... written with single backslashes.  json.loads then raises
    ``Invalid \\escape`` (and valid-looking sequences such as ``\\t`` may be
    misread as control characters).  This function only runs after normal JSON
    parsing has failed.  It repairs path-like string tokens conservatively and
    leaves already-valid escaped backslashes unchanged.
    """
    token_re = re.compile(r'"(?:\\.|[^"\\])*"')

    def fix_token(match):
        token = match.group(0)
        body = token[1:-1]
        is_drive_path = bool(re.match(r'^[A-Za-z]:\\', body))
        is_unc_path = body.startswith('\\\\')
        if is_drive_path or is_unc_path:
            if is_unc_path and not body.startswith(r'\\\\'):
                body=r'\\'+body
            out=[]; i=0
            while i < len(body):
                ch=body[i]
                if ch != '\\':
                    out.append(ch); i += 1; continue
                if i+1 < len(body) and body[i+1] == '\\':
                    # Already JSON-escaped backslash; keep the pair.
                    out.append('\\\\'); i += 2; continue
                if i+1 < len(body) and body[i+1] == '"':
                    # Preserve an escaped quote if one ever occurs.
                    out.append('\\\\"'); i += 2; continue
                # Raw Windows separator -> JSON escaped separator.
                out.append('\\\\'); i += 1
            return '"'+''.join(out)+'"'
        # Non-path syntax errors must not be guessed or silently repaired.
        return token

    return token_re.sub(fix_token, text)


def safe_json_load(path, *, repair_generated=False):
    """Load JSON without letting a malformed generated summary abort science stages.

    Normal JSON is returned unchanged.  If a legacy Windows path produced an
    invalid escape, one conservative repair is attempted.  Generated files can
    then be rewritten atomically in canonical JSON form, with the original kept
    beside it as ``*.invalid_escape_backup``.
    """
    path=Path(path)
    raw=path.read_text(encoding='utf-8-sig')
    try:
        return json.loads(raw), None
    except json.JSONDecodeError as first:
        repaired=_repair_json_text(raw)
        try:
            data=json.loads(repaired)
        except json.JSONDecodeError:
            raise first
        warning=f'已修复旧JSON中的Windows路径转义: {path.name}'
        if repair_generated:
            backup=path.with_name(path.name+'.invalid_escape_backup')
            if not backup.exists():
                shutil.copy2(path,backup)
            write_json(path,clean(data))
        return data, warning


class Trace:
    def __init__(self, directory, settings, kind='workflow'):
        self.directory = Path(directory); self.directory.mkdir(parents=True, exist_ok=True)
        self.kind = kind; self.settings = clean(settings)
        self.run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')+'_'+uuid.uuid4().hex[:8]
        self.rows = []; self.started = time.time(); self.active = None

    def __enter__(self):
        self.handle = (self.directory/(self.kind+'_processing.log')).open('w', encoding='utf-8')
        self.console = sys.stdout; sys.stdout = Tee(self.console, self.handle)
        write_json(self.directory/(self.kind+'_settings_used.json'), dict(run_id=self.run_id, settings=self.settings))
        print(f'【流程开始】{self.kind} · run={self.run_id}', flush=True)
        print('【实际参数】'+json.dumps(self.settings, ensure_ascii=False, default=str), flush=True)
        self.save(False)
        return self

    def save(self, complete):
        write_json(self.directory/(self.kind+'_trace.json'), clean(dict(run_id=self.run_id, complete=complete,
            elapsed_s=time.time()-self.started, stages=self.rows, settings_file=self.kind+'_settings_used.json',
            log_file=self.kind+'_processing.log')))

    def note(self, stage, message, **details):
        print(f'【{stage}】{message}', flush=True)
        row=dict(stage=stage, status='info', message=message, details=details, utc=datetime.now(timezone.utc).isoformat())
        self.rows.append(row); self.save(False)

    @contextmanager
    def stage(self, name, parameters=None):
        started=time.time(); row=dict(stage=name,status='running',parameters=parameters or {})
        self.rows.append(row); self.active=name; self.save(False)
        print(f'【阶段开始】{name} · '+json.dumps(parameters or {}, ensure_ascii=False, default=str), flush=True)
        print('@@EVENT '+json.dumps(dict(stage=name,status='running'),ensure_ascii=False),flush=True)
        result={}
        try:
            yield result
        except Exception as exc:
            row.update(status='failed',error=str(exc),elapsed_s=time.time()-started)
            print(f'【阶段失败】{name} · {exc}',flush=True); self.save(False); raise
        else:
            row.update(status=result.pop('status','complete'),results=result,elapsed_s=time.time()-started)
            print(f'【阶段结果】{name} · 耗时 {row["elapsed_s"]:.2f} s · '+json.dumps(result,ensure_ascii=False,default=str),flush=True)
            print('@@EVENT '+json.dumps(dict(stage=name,status=row['status']),ensure_ascii=False),flush=True)
            self.save(False)
        finally: self.active=None

    def __exit__(self, cls, exc, tb):
        if exc is not None:
            self.note('流程中断',str(exc), error_type=cls.__name__)
        self.save(exc is None)
        print(f'【流程结束】{"完成" if exc is None else "失败"} · 追溯文件：{self.directory/(self.kind+"_trace.json")}',flush=True)
        sys.stdout=self.console; self.handle.close()


def file_summary(root, stage):
    """A bounded summary; CSV row counts are parsed without loading waveforms."""
    import pandas as pd
    root=Path(root); out=root/'output'
    markers={'import':root/'hdf_import_audit.json','prepare':out/'prepare_complete.json',
        'scan':out/'scan_complete.json','controls':out/'controls_complete.json',
        'families':out/'families_complete.json','verify':out/'verification.json',
        'invert':out/'inversion_complete.json'}
    files={'import':['input_station_qc.csv'],'prepare':['pairs.csv','window_quality.csv'],'scan':['measured_phase_segments.csv','scan_candidate_intervals.csv'],
        'controls':['timeshift_summary.csv'],'families':['conditional_dispersion.csv','curve_catalog.csv'],'review':['candidate_review.csv'],
        'invert':['structure_layers.csv','inversion_attempts.csv']}
    result={}
    marker=markers.get(stage)
    if marker and marker.exists():
        try:
            data,warning=safe_json_load(marker,repair_generated=True); result['完成标记']=str(marker)
            result['摘要']={k:v for k,v in data.items() if isinstance(v,(bool,int,float,str))}
            if warning: result['追溯警告']=warning
        except (json.JSONDecodeError,UnicodeDecodeError) as exc:
            # A summary/provenance parse error must not invalidate a completed
            # waveform stage.  Preserve the marker for diagnosis and continue.
            result['追溯警告']=f'无法读取完成标记 {marker.name}: {exc}'
    for name in files.get(stage,[]):
        path=(root if stage=='import' else out)/name
        if path.exists():
            try:
                q=pd.read_csv(path); summary={'行数':len(q),'文件':str(path)}
                for c in ['station','pair','curve_id','frequency_hz']:
                    if c in q: summary[c+'_数量']=int(q[c].nunique())
                if 'frequency_hz' in q and len(q): summary['频率范围_Hz']=[float(q.frequency_hz.min()),float(q.frequency_hz.max())]
                result[name]=summary
            except pd.errors.EmptyDataError: result[name]={'行数':0,'文件':str(path)}
    if stage=='autoselect':
        q,_=safe_json_load(root/'selection.json',repair_generated=True); result={'假设组数':len(q.get('cases',[])),'文件':str(root/'selection.json')}
    if stage=='report':
        result['报告文件']=[str(out/name) for name in ('studio_report.html','report.html','studio_summary.json') if (out/name).exists()]
    return result


def stage_parameters(root,stage,config):
    root=Path(root)
    if stage=='import':
        spec,_=safe_json_load(root/'project.json');return {k:spec.get(k) for k in ('input_folder','excluded_stations','automatic_geometry')}
    if stage=='prepare':return {k:config[k] for k in ('sample_rate_hz','bandpass_hz','filter_guard_s','window_s','discard_first_s','minimum_windows','response')}
    if stage in ('scan','scan-diagnostics'):return dict(scan=config['scan'],validation=config['validation'])
    if stage=='controls':return config['controls']
    if stage=='families':return config['phase']
    if stage=='phase-stack':
        path=root/'phase_settings.json';return safe_json_load(path)[0] if path.exists() else {'direction_mode':'off (legacy CLI default)','time_mode':'full'}
    if stage=='invert':
        path=root/'selection.json';data=safe_json_load(path,repair_generated=True)[0] if path.exists() else {}
        return dict(settings=data.get('settings',{}),cases=len(data.get('cases',[])),selection_file=str(path))
    return dict(config_sha256=sha(root/'config.json'))
