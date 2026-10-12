"""Program-owned diagnostics for failed real-data preprocessing."""
from pathlib import Path
import numpy as np
import pandas as pd
from dispersion.common import load_config, intervals, csv, write_json


def diagnose_preparation(project):
    root=Path(project);c=load_config(root/'config.json');out=root/'output'
    valid=np.load(out/'cache/aligned_valid.npy',mmap_mode='r')
    stations=pd.read_csv(root/'stations.csv');rows=[]
    for i,row in stations.iterrows():
        lengths=np.asarray([b-a for a,b in intervals(valid[i])])/c['sample_rate_hz']
        rows.append({'station':row.station,'valid_fraction':float(valid[i].mean()),
                     'valid_runs':len(lengths),'maximum_run_s':float(lengths.max()) if len(lengths) else 0.,
                     'runs_at_least_minimum_segment':int(sum(lengths>=c['min_clean_segment_s']))})
    csv(out/'preparation_station_diagnostics.csv',rows)
    lengths=np.asarray([b-a for a,b in intervals(np.all(valid,axis=0))])/c['sample_rate_hz']
    summary={'stations':rows,'joint_valid_fraction':float(np.all(valid,axis=0).mean()),
             'joint_maximum_run_s':float(lengths.max()) if len(lengths) else 0.,
             'joint_long_runs':int(sum(lengths>=c['min_clean_segment_s'])),
             'minimum_clean_segment_s':c['min_clean_segment_s'],
             'window_quality':pd.read_csv(out/'window_quality.csv').accepted.value_counts().to_dict()}
    # JSON object keys require built-in bool strings.
    summary['window_quality']={str(k):int(v) for k,v in summary['window_quality'].items()}
    write_json(out/'preparation_diagnostics.json',summary)
    print(pd.DataFrame(rows).to_string(index=False),flush=True)
    print('joint:',summary['joint_valid_fraction'],summary['joint_maximum_run_s'],summary['joint_long_runs'],flush=True)
    return summary


def screen_prepared_array(project, folder, minimum_usable_fraction=.5,target_stations=None):
    """Remove chronically unavailable stations using only time-domain masks."""
    import json
    from h5data import ChannelRow, array_geometry
    from .acquisition import preview_array
    root=Path(project);data=Path(folder);report=diagnose_preparation(root)
    excluded=[r['station'] for r in report['stations'] if r['valid_fraction']<minimum_usable_fraction]
    if target_stations is not None:
        remaining=sorted([r for r in report['stations'] if r['station'] not in excluded],
                          key=lambda r:(-r['valid_fraction'],r['station']))
        if not isinstance(target_stations,int) or target_stations<1:raise ValueError('Target station count must be positive')
        excluded.extend(r['station'] for r in remaining[target_stations:])
    if not excluded:raise ValueError('No station meets the fixed low-availability exclusion rule; inspect data/configuration')
    path=data/'analysis_array_selection.json';selection=json.loads(path.read_text(encoding='utf-8-sig'))
    kept=[r for r in selection['stations'] if r['network']+'_'+r['station'] not in excluded]
    if not selection['constraints'].get('minimum_station_count',6)<=len(kept)<=selection['constraints'].get('maximum_station_count',49):raise ValueError('Too few stations after temporal-availability screening')
    g=array_geometry([ChannelRow(**(r|{'start':None,'end':None})) for r in kept])
    bounds=selection['constraints'];nearest=list(g['nearest_neighbor_m'].values())
    if min(nearest)<bounds['minimum_nearest_neighbor_m'] or max(nearest)>bounds['maximum_nearest_neighbor_m']:
        raise ValueError('Low-availability exclusion would violate the specified station spacing')
    for name in excluded:selection.setdefault('exclusions',[]).append({'station':name,'reason':'preprocessed_time_availability_screen',
                                                       'minimum_fraction':minimum_usable_fraction,'target_station_count':target_stations})
    selection.update(stations=kept,geometry=g,constraints=bounds|{'station_count':len(kept)},
                     temporal_screening_uses_phase_or_dispersion=False)
    write_json(path,selection);preview_array(selection,path.with_suffix('.png'))
    specpath=root/'project.json';spec=json.loads(specpath.read_text(encoding='utf-8-sig'))
    spec['excluded_stations']=sorted(set(spec.get('excluded_stations',[])+excluded));write_json(specpath,spec)
    write_json(data/'temporal_screening.json',{'excluded':excluded,'minimum_usable_fraction':minimum_usable_fraction,
               'target_station_count':target_stations,
               'station_diagnostics':report['stations'],'remaining_stations':len(kept),'geometry':g,
               'phase_or_dispersion_used_in_selection':False,'quality_thresholds_relaxed':False})
    print('temporal screening: excluded',excluded,'remaining',len(kept),flush=True)
    return selection
