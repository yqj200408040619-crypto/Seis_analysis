"""Header-only UTC extent discovery for full-record phase screening."""
from pathlib import Path
import zipfile
import numpy as np
from dispersion.common import load_config,parse_utc
from dispersion.prepare import station_table,source_files


def record_headers(row,c):
    """Do not inflate a multi-day waveform merely to discover its time axis."""
    from obspy import read
    records=[]
    for path in source_files(row,c):
        if Path(path).suffix.lower()=='.npz':
            with zipfile.ZipFile(path) as archive,archive.open('data.npy') as handle:
                version=np.lib.format.read_magic(handle)
                reader={(1,0):np.lib.format.read_array_header_1_0,(2,0):np.lib.format.read_array_header_2_0}.get(version)
                if reader is None:raise ValueError('Unsupported waveform NPY header version')
                shape,_,_=reader(handle)
            if len(shape)!=1:raise ValueError('NPZ data must be one-dimensional')
            with np.load(path,allow_pickle=False) as z:
                fs=float(z['fs_hz']);start=parse_utc(str(z['start_utc'].item()))
            traces=[(start,fs,shape[0])]
        else:
            traces=[(tr.stats.starttime,float(tr.stats.sampling_rate),tr.stats.npts)
                    for tr in read(path,headonly=True).select(id=row.seed_id)]
        for start,fs,count in traces:
            if not np.isfinite(fs) or fs<=0 or count<1:raise ValueError('Invalid waveform time axis: '+path)
            start-=float(row.clock_offset_s)
            records.append(dict(path=path,start=start,end=start+count/fs,fs_hz=fs,samples=count))
    if not records:raise ValueError('No matching waveform records: '+row.station)
    records.sort(key=lambda r:r['start'].ns)
    if len({round(r['fs_hz'],9) for r in records})!=1:raise ValueError('Within-station sample-rate changes require separate input epochs')
    if any(b['start']<=a['end']-1/a['fs_hz'] for a,b in zip(records,records[1:])):
        raise ValueError('Overlapping waveform samples: resolve duplicates or mask the overlap explicitly')
    return records


def available_phase_extent(root,station_names):
    """Union timeline: each pair retains ALL jointly valid windows, including tails.

    A global intersection of every selected station would needlessly drop data
    from other pairs. Missing samples remain NaN; no gap is concatenated away.
    CCF scan start/end/cohort cuts do not shorten this independent time traversal.
    """
    c=load_config(Path(root)/'config.json');table,_,_=station_table(c)
    table=table[table.station.isin(station_names)]
    if set(table.station)!=set(station_names):raise ValueError('Missing requested station')
    bounds=[]
    for _,row in table.iterrows():
        records=record_headers(row,c)
        bounds.append(dict(station=row.station,start_utc=str(records[0]['start']),
                           end_utc_exclusive=str(max(r['end'] for r in records)),files=len({r['path'] for r in records})))
    return min(parse_utc(r['start_utc']) for r in bounds),max(parse_utc(r['end_utc_exclusive']) for r in bounds),bounds


def phase_window_geometry(frequency,fs,repeat_cycles,slide_cycles):
    values=np.array([frequency,fs,repeat_cycles,slide_cycles],float)
    if not np.isfinite(values).all() or (values<=0).any() or frequency>=fs/2:
        raise ValueError('Frequency must be positive and below Nyquist; cycles/sample rate must be positive')
    n=int(np.floor(repeat_cycles*fs/frequency+.5));step=max(1,int(np.floor(slide_cycles*fs/frequency+.5)))
    if n<2:raise ValueError('Frequency-adaptive window must contain at least two samples')
    return dict(window_samples=n,slide_samples=step,window_duration_s=n/fs,slide_duration_s=step/fs,
                requested_window_s=repeat_cycles/frequency,requested_slide_s=slide_cycles/frequency)


def full_window_starts(sample_count,n,step,include_last=False):
    """All complete windows; optional end anchor, uniform hops match MATLAB."""
    if n>sample_count:raise ValueError('Record is shorter than the frequency-adaptive window')
    starts=np.arange(0,sample_count-n+1,step,dtype=np.int64)
    last=sample_count-n
    return np.r_[starts,last] if include_last and starts[-1]!=last else starts


class NativePhaseReader:
    """Inflate a large canonical NPZ once onto disk; read only the needed slice."""
    def __init__(self,row,c,records,workdir):
        self.row=row;self.records=records;self.workdir=Path(workdir);self.cache={}

    def read(self,start,end):
        from obspy import Stream,Trace,read
        import shutil,hashlib
        result=Stream();paths=sorted({r['path'] for r in self.records if r['start']<=end and r['end']>start})
        for path in paths:
            if Path(path).suffix.lower()=='.npz':
                r=next(q for q in self.records if q['path']==path)
                if path not in self.cache:
                    key=hashlib.sha256(path.encode()).hexdigest()[:24];arrays={}
                    with zipfile.ZipFile(path) as archive:
                        for name in ('data','valid'):
                            if name+'.npy' not in archive.namelist():continue
                            target=self.workdir/(key+'_'+name+'.npy')
                            with archive.open(name+'.npy') as source,target.open('wb') as dest:shutil.copyfileobj(source,dest,1024*1024)
                            arrays[name]=np.load(target,mmap_mode='r',allow_pickle=False)
                    if arrays['data'].ndim!=1 or ('valid' in arrays and arrays['valid'].shape!=arrays['data'].shape):raise ValueError('Invalid waveform data/mask')
                    self.cache[path]=arrays
                arrays=self.cache[path];fs=r['fs_hz']
                first=max(0,int(np.ceil(float(start-r['start'])*fs-1e-7)))
                last=min(r['samples'],int(np.floor(float(end-r['start'])*fs+1e-7))+1)
                if last<=first:continue
                data=np.array(arrays['data'][first:last],dtype=float)
                if 'valid' in arrays:data[~np.asarray(arrays['valid'][first:last],bool)]=np.nan
                ids=self.row.seed_id.split('.')
                result+=Trace(data,dict(network=ids[0],station=ids[1],location=ids[2],channel=ids[3],sampling_rate=fs,starttime=r['start']+first/fs))
            else:
                offset=float(self.row.clock_offset_s)
                st=read(path,starttime=start+offset,endtime=end+offset).select(id=self.row.seed_id)
                for tr in st:
                    tr.stats.starttime-=offset;tr.data=np.asarray(np.ma.filled(tr.data,np.nan),float)
                result+=st
        result.sort(keys=['starttime']);joined=Stream();current=[]
        for tr in result:
            if current:
                prev=current[-1];residual=float(tr.stats.starttime-prev.stats.endtime)*tr.stats.sampling_rate-1
                if tr.stats.starttime.ns<=prev.stats.endtime.ns:raise ValueError('Overlapping waveform samples')
                if abs(residual)>1e-5:
                    q=current[0].copy();q.data=np.concatenate([t.data for t in current]);joined+=q;current=[]
            current.append(tr)
        if current:
            q=current[0].copy();q.data=np.concatenate([t.data for t in current]);joined+=q
        return joined

    def close(self):
        for arrays in self.cache.values():
            for data in arrays.values():data._mmap.close()
        self.cache.clear()
