from pathlib import Path
import math
import numpy as np
import pandas as pd
import pytest

from ndstudio.geothermal import save_settings, DEFAULT_SETTINGS
from ndstudio.tomography import (
    Grid2D, _segment_lengths, _prepare_observations, run_phase_velocity_tomography,
    TomographyError,
)


def test_segment_lengths_conserve_distance():
    grid=Grid2D(np.linspace(0,4000,5),np.linspace(0,4000,5),np.array([500,1500,2500,3500]),np.array([500,1500,2500,3500]),4,4)
    row=_segment_lengths(0,0,4000,4000,grid)
    assert np.isclose(row.sum(),math.sqrt(2)*4000,rtol=1e-12)
    assert np.count_nonzero(row) == 4


def test_projected_wavevector_converted_to_pair_equivalent_and_oblique_rejected():
    base=dict(DEFAULT_SETTINGS,tomo_quality_mode='all_valid',tomo_max_direction_angle_deg=15,tomo_duplicate_policy='error')
    df=pd.DataFrame([
        dict(pair='A--B',frequency_hz=1.0,phase_velocity_m_s=2000,distance_m=1000,
             velocity_type='wave_vector_phase_velocity',measurement_method='direction_projected_goertzel',
             pair_azimuth_deg=0,confirmed_azimuth_deg=10,projected_distance_m=900),
        dict(pair='A--C',frequency_hz=1.0,phase_velocity_m_s=2000,distance_m=1000,
             velocity_type='wave_vector_phase_velocity',measurement_method='direction_projected_goertzel',
             pair_azimuth_deg=0,confirmed_azimuth_deg=40,projected_distance_m=760),
    ])
    q,w=_prepare_observations(df,base)
    assert len(q)==1
    assert q.iloc[0].pair=='A--B'
    assert np.isclose(q.iloc[0].tomography_velocity_m_s,1000*2000/900)
    assert any('剔除' in x for x in w)


def test_duplicate_pair_frequency_is_error_by_default():
    df=pd.DataFrame([
        dict(pair='A--B',frequency_hz=1.,phase_velocity_m_s=1000,distance_m=1000,default_quality_pass=True),
        dict(pair='A--B',frequency_hz=1.,phase_velocity_m_s=1200,distance_m=1000,default_quality_pass=True),
    ])
    with pytest.raises(TomographyError,match='多个候选'):
        _prepare_observations(df,dict(DEFAULT_SETTINGS))


def _synthetic_paths():
    pts={
        'A':(0.,0.),'B':(2000.,0.),'C':(4000.,0.),'D':(4000.,2000.),
        'E':(4000.,4000.),'F':(2000.,4000.),'G':(0.,4000.),'H':(0.,2000.),
    }
    pairs=[]
    names=list(pts)
    # A dense set of chords through the square.
    for i,a in enumerate(names):
        for b in names[i+1:]:
            x0,y0=pts[a];x1,y1=pts[b]
            d=math.hypot(x1-x0,y1-y0)
            if d<2500: continue
            pairs.append((a,b,x0,y0,x1,y1,d))
    # 4x4 true slowness checker, same geometry as configured inversion grid.
    grid=Grid2D(np.linspace(0,4000,5),np.linspace(0,4000,5),np.array([500,1500,2500,3500]),np.array([500,1500,2500,3500]),4,4)
    pattern=np.array([1 if ((ix//2+iy//2)%2==0) else -1 for iy in range(4) for ix in range(4)],float)
    rows=[]
    for freq,base_v in [(1.0,2000.),(2.0,1700.)]:
        s0=1000/base_v
        s=s0*(1+0.05*pattern)
        for k,(a,b,x0,y0,x1,y1,d) in enumerate(pairs):
            G=_segment_lengths(x0,y0,x1,y1,grid)/1000
            t=G@s
            v=d/t
            rows.append(dict(curve_id=f'C{freq}_{k}',pair=f'{a}--{b}',frequency_hz=freq,
                phase_velocity_m_s=v,distance_m=d,x_a_m=x0,y_a_m=y0,x_b_m=x1,y_b_m=y1,
                default_quality_pass=True,velocity_type='conditional_pair_velocity',measurement_method='synthetic'))
    return pd.DataFrame(rows)


def test_end_to_end_straight_ray_tomography_and_local_dispersion(tmp_path: Path):
    root=tmp_path/'project';root.mkdir()
    settings=dict(DEFAULT_SETTINGS,tomo_nx=4,tomo_ny=4,tomo_margin_fraction=0.0,tomo_damping=0.01,tomo_smoothing=0.05,
                  tomo_min_paths_per_frequency=6,tomo_min_cell_rays=1,tomo_checkerboard=True,tomo_checkerboard_block_cells=2,
                  tomo_max_figures=2,tomo_quality_mode='default_pass')
    save_settings(root,settings)
    source=root/'paths.csv';_synthetic_paths().to_csv(source,index=False)
    result=run_phase_velocity_tomography(root,paths_file=source)
    assert result['complete'] is True
    assert result['frequencies_inverted']==2
    out=root/'output/geothermal/tomography2d'
    models=pd.read_csv(out/'phase_velocity_tomography.csv')
    local=pd.read_csv(out/'local_dispersion.csv')
    cat=pd.read_csv(out/'frequency_catalog.csv')
    assert len(models)==32
    assert set(local.frequency_hz)=={1.0,2.0}
    assert (models.ray_count>0).sum()>8
    assert (models.resolution.between(0,1)).all()
    assert (cat.status=='accepted').all()
    assert (out/'tomography_2d_report.html').is_file()
    assert len(list((out/'figures').glob('*.png')))==2
    # Synthetic data are exact; regularized fit should still be much better than 0.1 s RMS.
    assert float(cat.travel_time_rms_s.max()) < 0.1


def test_selected_curve_ids_supports_selection_bands(tmp_path: Path):
    import json
    from ndstudio.geothermal import _selected_curve_ids
    root=tmp_path/'p';root.mkdir()
    (root/'selection.json').write_text(json.dumps({'cases':[{'name':'CASE_A','bands':[{'curve_id':'C1','mode':0},{'curve_id':'C2','mode':0}]},{'name':'CASE_B','curve_ids':['C3']}]}),encoding='utf-8')
    assert _selected_curve_ids(root,'CASE_A')=={'C1','C2'}
    assert _selected_curve_ids(root,'CASE_B')=={'C3'}
    assert _selected_curve_ids(root)=={'C1','C2','C3'}
