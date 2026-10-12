"""MUSIC direction correction of CCF/FTAN candidates before review and phase gates."""
from pathlib import Path
import json,tempfile
import numpy as np
import pandas as pd
from scipy.signal.windows import hann
from dispersion.common import load_config,csv,write_json,sha


def direction_settings(root):
    from .direction import settings
    path=Path(root)/'phase_settings.json'
    phase=json.loads(path.read_text(encoding='utf-8-sig')) if path.exists() else {}
    return settings(phase.get('direction') or {'mode':'filter'})


def signature(root):
    root=Path(root)
    return dict(config=sha(root/'config.json'),stations=sha(root/'stations.csv'),
        curves=sha(root/'output/conditional_dispersion.csv'),direction=direction_settings(root),
        code={name:sha(Path(__file__).with_name(name)) for name in ('music_front.py','music.py','direction.py','phase_stack.py','phase_records.py')},
        canonical={str(p.relative_to(root)):sha(p) for p in sorted((root/'canonical').glob('*')) if p.is_file()})


def current(root):
    p=Path(root)/'output/music_front/music_complete.json'
    if not p.exists():return False
    marker=json.loads(p.read_text(encoding='utf-8-sig'))
    if marker.get('signature')!=signature(root):return False
    return bool(marker.get('output_sha256')) and all((p.parent/name).exists() and sha(p.parent/name)==digest for name,digest in marker['output_sha256'].items())


def frequency_reference(result,starts,n,delta,c,branch,angle_gate=False):
    """Choose a measured direction at THIS frequency, without a curve velocity prior."""
    from .direction import pair_gate,independent_indices
    projection=delta[0]*np.sin(np.radians(result['azimuth_deg']))+delta[1]*np.cos(np.radians(result['azimuth_deg']))
    sign=1 if branch=='positive' else -1;distance=np.linalg.norm(delta)
    good=result['quality']&(sign*projection>distance*.1)
    if angle_gate:good &= pair_gate(result,*delta,c)['direction_valid']
    if not good.any():return None
    groups,counts=np.unique(result['direction_group'][good],return_counts=True);group=int(groups[np.argmax(counts)])
    selected=good&(result['direction_group']==group);ids=independent_indices(starts,n,selected)
    if len(ids)<c['minimum_independent_windows']:return None
    az=float(np.degrees(np.angle(np.mean(np.exp(1j*np.radians(result['azimuth_deg'][ids])))))%360)
    reference=float(delta@np.array([np.sin(np.radians(az)),np.cos(np.radians(az))]))
    uncertainty=float(np.max(result['direction_uncertainty_deg'][ids]))
    # A near-transverse projection crossing zero has no stable speed correction.
    angle=abs((az-np.degrees(np.arctan2(delta[0],delta[1]))+90)%180-90)
    if angle+uncertainty>=85:return None
    return dict(azimuth_deg=az,selected_direction_group=group,projected_distance_m=reference,distance_m=float(distance),
        independent_windows=len(ids),direction_uncertainty_deg=uncertainty,music_velocity_m_s=float(np.median(result['velocity_m_s'][ids])),
        indices=ids,eligible=selected)


def correct_candidates(curves,pairs,references,scan):
    from .phase_bands import assign_bands
    from .reconstruction import derived_group
    rows=[]
    for r in curves.to_dict('records'):
        ref=references.get((r['pair'],r['branch'],round(r['frequency_hz'],9)))
        row=dict(r);row.update(accepted=bool(ref),direction_method='music',projection_applied=bool(ref),
            velocity_type='wave_vector_phase_velocity',apparent_pair_velocity_m_s=r['phase_velocity_m_s'],
            apparent_pair_group_velocity_m_s=r.get('group_velocity_m_s',np.nan),parent_curve_id=r['curve_id'])
        if ref:
            row.update({k:v for k,v in ref.items() if k not in ('indices','eligible')})
            factor=abs(ref['projected_distance_m'])/ref['distance_m']
            row.update(phase_velocity_m_s=r['phase_velocity_m_s']*factor,propagation_azimuth_deg=ref['azimuth_deg'],confirmed_azimuth_deg=ref['azimuth_deg'],
                wavenumber_rad_m=2*np.pi*r['frequency_hz']/(r['phase_velocity_m_s']*factor),
                distance_over_wavelength=ref['distance_m']*r['frequency_hz']/r['phase_velocity_m_s'])
        else:
            row.update(phase_velocity_m_s=np.nan,propagation_azimuth_deg=np.nan,projected_distance_m=np.nan,
                distance_m=float(pairs.loc[r['pair_index'],'distance_m']),music_velocity_m_s=np.nan)
        rows.append(row)
    audit=pd.DataFrame(rows)
    if not len(audit):
        audit=curves.iloc[:0].copy()
        for key in ('accepted','projection_applied','propagation_azimuth_deg','projected_distance_m','distance_m','parent_curve_id','velocity_type'):audit[key]=pd.Series(dtype=object)
        return audit,audit,pd.DataFrame(columns=['curve_id','parent_curve_id','segment_id','cohort','pair_index','pair','branch','constant_cycle','source_phase_rad','fmin_hz','fmax_hz','points','background_discriminated','compatible_other_pairs','measurement_method'])
    # Direction is estimated at each frequency. Bin-edge changes alone cannot break a band.
    band_input=audit.drop(columns=['selected_direction_group'],errors='ignore')
    checked=assign_bands(band_input,np.sort(curves.frequency_hz.unique()),scan)
    audit['band_id']=checked.band_id;audit['band_reason']=checked.band_reason;audit['accepted']=checked.accepted
    final=audit[audit.accepted].copy();catalog=[];parts=[]
    for band,g in final.groupby('band_id',sort=False):
        g=g.sort_values('frequency_hz').copy();g['curve_id']=band
        # Frequency-varying projection: differentiate corrected f/c, never multiply u by one mean cosine.
        g['group_velocity_m_s']=derived_group(g.frequency_hz,g.phase_velocity_m_s);g['group_velocity_is_derived']=True
        if not np.isfinite(g.group_velocity_m_s).all() or (g.group_velocity_m_s<=0).any():
            audit.loc[g.index,'accepted']=False;audit.loc[g.index,'band_reason']='修正后派生群速度非正';continue
        parts.append(g);catalog.append(dict(curve_id=band,parent_curve_id=g.parent_curve_id.iloc[0],segment_id=g.segment_id.iloc[0],
            cohort=g.cohort.iloc[0],pair_index=int(g.pair_index.iloc[0]),pair=g.pair.iloc[0],branch=g.branch.iloc[0],
            constant_cycle=int(g.constant_cycle.iloc[0]),source_phase_rad=float(g.source_phase_rad.iloc[0]),
            fmin_hz=float(g.frequency_hz.min()),fmax_hz=float(g.frequency_hz.max()),points=len(g),
            background_discriminated=bool(g.background_discriminated.all()),compatible_other_pairs=0,measurement_method='music_projected_ccf'))
    final=pd.concat(parts,ignore_index=True) if parts else final.iloc[:0]
    cat=pd.DataFrame(catalog,columns=['curve_id','parent_curve_id','segment_id','cohort','pair_index','pair','branch','constant_cycle','source_phase_rad','fmin_hz','fmax_hz','points','background_discriminated','compatible_other_pairs','measurement_method'])
    return final,audit,cat


def run_music_front(project,aligned_h5=None):
    from .direction import DirectionSession
    from .phase_stack import align_project
    from .phase_records import available_phase_extent,phase_window_geometry
    from .trace import Trace
    from dispersion.prepare import station_table
    from dispersion.common import parse_utc
    root=Path(project).resolve();dest=root/'output/music_front';dest.mkdir(parents=True,exist_ok=True)
    wanted=signature(root)
    if aligned_h5 is None and current(root):
        print('【MUSIC前段】复用已校验的逐频点方向/速度',flush=True);return json.loads((dest/'music_complete.json').read_text())
    curves=pd.read_csv(root/'output/conditional_dispersion.csv');pairs=pd.read_csv(root/'output/pairs.csv').set_index('pair_index')
    c=load_config(root/'config.json');dc=direction_settings(root);table,_,_=station_table(c);xy=table.set_index('station');names=set(table.station)
    resources=[];references={};rows=[]
    # Reuse aligned bytes only when their canonical input and calibration hashes still match.
    existing=dest/'music_complete.json'
    if aligned_h5 is None and existing.exists() and (dest/'aligned_waveforms.h5').exists():
        previous=json.loads(existing.read_text(encoding='utf-8-sig'))
        old=previous.get('signature',{})
        if old.get('canonical')==wanted['canonical'] and old.get('config')==wanted['config'] and old.get('stations')==wanted['stations'] and previous.get('output_sha256',{}).get('aligned_waveforms.h5')==sha(dest/'aligned_waveforms.h5'):
            aligned_h5=dest/'aligned_waveforms.h5'
    (dest/'music_complete.json').unlink(missing_ok=True)
    with tempfile.TemporaryDirectory(prefix='.music_work_',dir=dest) as work,Trace(dest,dict(method='music',direction=dc),'music') as trace:
        try:
            if aligned_h5:
                import h5py,shutil
                path=Path(aligned_h5)
                with h5py.File(path) as h:
                    if not h.attrs.get('response_corrected'):raise ValueError('MUSIC需要响应校正波形')
                    fs=float(h.attrs['fs_hz']);origin=parse_utc(h.attrs['start_utc'])
                    if fs!=c['sample_rate_hz']:raise ValueError('MUSIC波形采样率与项目不符')
                    arrays={s:np.where(h['stations/'+s+'/valid'][:],h['stations/'+s+'/data'][:],np.nan) for s in names}
                    limits={s:(int(h['stations/'+s].attrs.get('record_start_sample',0)),int(h['stations/'+s].attrs.get('record_end_sample_exclusive',len(arrays[s])))) for s in names}
                if path.resolve()!=(dest/'aligned_waveforms.h5').resolve():shutil.copyfile(path,dest/'aligned_waveforms.h5')
            else:
                first,last,bounds=available_phase_extent(root,names)
                arrays,fs,origin,_=align_project(root,names,str(first),str(last),dest/'aligned_waveforms.h5',work,resources)
                limits={r['station']:(round(float(parse_utc(r['start_utc'])-origin)*fs),round(float(parse_utc(r['end_utc_exclusive'])-origin)*fs)) for r in bounds}
            session=DirectionSession(arrays,xy,fs,origin,limits,dest,work,dc,0.,trace,True)
            for f,points in curves.groupby('frequency_hz',sort=True):
                geometry=phase_window_geometry(float(f),fs,dc['music_front_cycles'],dc['music_front_cycles']);n=geometry['window_samples'];step=geometry['slide_samples']
                result=session.get(float(f),float(f),n,step,hann(n,sym=True))
                for (pair,branch),_ in points.groupby(['pair','branch']):
                    a,b=pair.split('--');delta=xy.loc[b,['x_m','y_m']].to_numpy(float)-xy.loc[a,['x_m','y_m']].to_numpy(float)
                    ref=frequency_reference(result,result['start_sample'],n,delta,dc,branch,angle_gate=False)
                    references[(pair,branch,round(f,9))]=ref
                    rows.append(dict(pair=pair,branch=branch,frequency_hz=f,accepted=bool(ref),**({k:v for k,v in ref.items() if k not in ('indices','eligible')} if ref else {})))
                trace.note('MUSIC前段逐频修正',f'{f:g} Hz · 阵列合格{int(result["quality"].sum())}/{len(result["quality"])}窗；前段不使用45°台对门槛')
            final,audit,cat=correct_candidates(curves,pairs,references,c['scan'])
            links=compatible_links(final,c);csv(dest/'cross_pair_compatible_curves.csv',links)
            if len(cat):
                indices=cat.set_index('curve_id').pair_index.to_dict();neighbors={cid:set() for cid in cat.curve_id}
                for r in links.itertuples():neighbors[r.curve_a].add(indices[r.curve_b]);neighbors[r.curve_b].add(indices[r.curve_a])
                cat['compatible_other_pairs']=cat.curve_id.map(lambda cid:len(neighbors[cid]))
            csv(dest/'dispersion.csv',final);csv(dest/'point_audit.csv',audit);csv(dest/'catalog.csv',cat)
            csv(dest/'frequency_directions.csv',rows,columns=None if rows else ['pair','branch','frequency_hz','accepted'])
            session.finish(np.sort(curves.frequency_hz.unique()),c['scan'])
            marker=dict(complete=True,method='music',signature=wanted,source_velocity_type='apparent_pair_velocity',
                velocity_type='wave_vector_phase_velocity',front_pair_angle_limit_applied=False,
                corrected_points=len(final),corrected_curves=len(cat),no_qualified_direction_no_true_velocity=True)
            marker['output_sha256']={name:sha(dest/name) for name in ('dispersion.csv','catalog.csv','point_audit.csv','frequency_directions.csv','aligned_waveforms.h5')}
            plot_front(final,dest)
            (dest/'report.html').write_text('<meta charset="utf-8"><h1>MUSIC逐频点方向与速度修正</h1><p>前段不使用45°台对门槛。c = c_app |ℓ(f)| / d；方向、投影距离和校正系数逐频点计算。近直线台阵、无信号、多波源及混叠窗不输出真实速度。</p><p><a href="direction/report.html">MUSIC诊断</a> · <a href="frequency_directions.csv">逐频点方向</a> · <a href="dispersion.csv">修正频散</a></p><img src="dispersion_direction.png" style="max-width:100%">'+pd.DataFrame(rows).to_html(index=False),encoding='utf-8')
            write_json(dest/'music_complete.json',marker)
            print(f'【MUSIC前段完成】{len(final)}点 / {len(cat)}条连续候选；{dest / "report.html"}',flush=True)
            return marker
        finally:
            for item in reversed(resources):
                if isinstance(item,np.memmap):item._mmap.close()
                else:item.close()


def compatible_links(curves,c):
    """Cross-pair links use measured MUSIC propagation and corrected velocities."""
    import itertools
    columns=['curve_a','curve_b','segment_a','segment_b','cohort','branch','overlap_min_hz','overlap_max_hz','common_points','phase_velocity_median_fraction','phase_velocity_max_fraction','group_velocity_median_fraction','same_rayleigh_mode_verified']
    rows=[];groups={cid:g for cid,g in curves.groupby('curve_id')}
    for ca,cb in itertools.combinations(groups,2):
        a,b=groups[ca],groups[cb]
        if a.pair_index.iloc[0]==b.pair_index.iloc[0] or a.cohort.iloc[0]!=b.cohort.iloc[0]:continue
        q=a.merge(b,on='frequency_hz',suffixes=('_a','_b'))
        if len(q)<c['scan']['min_points'] or q.frequency_hz.max()-q.frequency_hz.min()<c['scan']['min_band_hz']-1e-9:continue
        angle=abs((q.propagation_azimuth_deg_a-q.propagation_azimuth_deg_b+180)%360-180)
        if (angle>c['phase']['max_direction_difference_deg']).any():continue
        dv=abs(q.phase_velocity_m_s_a-q.phase_velocity_m_s_b)/((q.phase_velocity_m_s_a+q.phase_velocity_m_s_b)/2)
        du=abs(q.group_velocity_m_s_a-q.group_velocity_m_s_b)/((q.group_velocity_m_s_a+q.group_velocity_m_s_b)/2)
        if dv.median()<=c['phase']['link_median_fraction'] and dv.max()<=c['phase']['link_max_fraction'] and du.median()<=c['phase']['link_group_median_fraction']:
            rows.append(dict(curve_a=ca,curve_b=cb,segment_a=a.segment_id.iloc[0],segment_b=b.segment_id.iloc[0],cohort=a.cohort.iloc[0],branch=a.branch.iloc[0],overlap_min_hz=q.frequency_hz.min(),overlap_max_hz=q.frequency_hz.max(),common_points=len(q),phase_velocity_median_fraction=dv.median(),phase_velocity_max_fraction=dv.max(),group_velocity_median_fraction=du.median(),same_rayleigh_mode_verified=False))
    return pd.DataFrame(rows,columns=columns)


def plot_front(curves,dest):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    fig,axes=plt.subplots(2,1,figsize=(10,7),constrained_layout=True)
    for i,(cid,q) in enumerate(curves.groupby('curve_id')):
        if i>=24:break
        q=q.sort_values('frequency_hz');axes[0].plot(q.frequency_hz,q.apparent_pair_velocity_m_s,'--',alpha=.25,color='#64748b')
        axes[0].plot(q.frequency_hz,q.phase_velocity_m_s,'-',alpha=.65,color='#087f8c')
        axes[1].scatter(q.frequency_hz,q.propagation_azimuth_deg,s=12,alpha=.6,color='#087f8c')
    axes[0].set(xlabel='Frequency / Hz',ylabel='Phase velocity / m/s',title='MUSIC projection: corrected (solid), original apparent (dashed)')
    axes[1].set(xlabel='Frequency / Hz',ylabel='Propagation azimuth / deg',ylim=(0,360))
    if curves.empty:axes[0].text(.5,.5,'No qualified continuous MUSIC band',ha='center',transform=axes[0].transAxes)
    for ax in axes:ax.grid(alpha=.2)
    fig.savefig(Path(dest)/'dispersion_direction.png',dpi=150);plt.close(fig)
