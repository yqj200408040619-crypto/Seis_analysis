"""Synthetic HDF5 waveform and known-layered-model inversion examples."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import h5py
from dispersion.common import write_json,merge,DEFAULTS
from .project import new_project
from .candidates import INVERSION_DEFAULTS


def waveform_example(directory,field='plane'):
 from dispersion.demo import generate
 root=Path(directory).resolve();raw=root/'h5';raw.mkdir(parents=True,exist_ok=True);native=root/'_synthetic_native'
 generate(native,field=field);station=pd.read_csv(native/'stations.csv')
 for row in station.to_dict('records'):
  with np.load(native/row['file']) as z:
   path=raw/(row['station']+'.h5')
   with h5py.File(path,'w') as h:
    # Alternate array order deliberately to test both N×3 and 3×N mappings.
    x=z['data'];values=np.column_stack((x*.05,x*.03,x))
    if row['station'] in ('D02','D04'):values=values.T
    ds=h.create_dataset('waveform',data=values,compression='gzip',shuffle=True)
    ds.attrs['components']=np.asarray(['HHN','HHE','HHZ'],dtype='S3')
    for k in ['station','seed_id','x_m','y_m','units','precorrected','polarity','polarity_verified','clock_offset_s']:ds.attrs[k]=row[k]
    ds.attrs['start_utc']=z['start_utc'].item();ds.attrs['sampling_rate_hz']=float(z['fs_hz'])
    h.create_dataset('valid',data=z['valid']);ds.attrs['valid_dataset']='/valid'
 # The valid mask is an explicit mapping, not guessed from an arbitrary array.
 project=new_project(root/'project');spec=json.loads((project/'project.json').read_text(encoding='utf-8-sig'));spec.update(input_folder=str(raw),automatic_geometry=False);write_json(project/'project.json',spec)
 profile=json.loads((project/'mapping.json').read_text(encoding='utf-8-sig'));profile['defaults']['valid_dataset']='/valid';write_json(project/'mapping.json',profile)
 c=json.loads((native/'config.json').read_text(encoding='utf-8-sig'));c.update(data_root='canonical',stations_csv='stations.csv',output='output');write_json(project/'config.json',merge(DEFAULTS,c))
 write_json(root/'truth.json',json.loads((native/'truth.json').read_text(encoding='utf-8-sig')))
 print('HDF_EXAMPLE',project,flush=True);return project


def inversion_example(directory):
 from dispersion.inversion import forward
 root=Path(directory).resolve();root.mkdir(parents=True,exist_ok=True)
 # 35 m slow finite layer above an 800 m/s halfspace; fixed Vp/Vs=2 and rho=2.
 model=np.array([[.035,.7,.35,2.],[0.,1.6,.8,2.]])
 f=np.linspace(1.5,5.,45);c=forward(model,f,0);u=forward(model,f,0,True);rows=[]
 for freq,phase,group in zip(f,c,u):
  rows.append(dict(curve_id='KNOWN_N0',segment_id='SYNTHETIC',cohort='synthetic',pair_index=0,pair='A--B',branch='positive',
    frequency_hz=freq,phase_velocity_m_s=phase,group_velocity_m_s=group,constant_cycle=0,source_phase_rad=0.,distance_over_wavelength=2000*freq/phase,
    instrument_response_corrected=True,background_discriminated=True,source_phase_verified=True,is_model_prediction=False,
    synthetic_observation=True,absolute_cycle_resolved=True,rayleigh_mode_verified=True))
 pd.DataFrame(rows).to_csv(root/'curves.csv',index=False)
 settings=dict(INVERSION_DEFAULTS);settings.update(layers=[1],starts=5,max_evaluations=180,thickness_bounds_m=[10.,150.],vs_bounds_m_s=[200.,1200.],require_independent_solver=False)
 write_json(root/'selection.json',{'settings':settings,'cases':[{'name':'KNOWN_FUNDAMENTAL','bands':[{'curve_id':'KNOWN_N0','mode':0}],
   'assumptions':{'synthetic_truth_case':True,'common_1d_structure':True}}]})
 write_json(root/'truth.json',{'finite_layer_thickness_m':[35.],'vs_m_s':[350.,800.],'vp_vs_ratio':2.,'density_g_cm3':2.,
   'is_real_earth_observation':False,'generated_by':'disba; independent solver validates forward values when available'})
 return root
