"""Recheck stored phase gates and sample transforms against aligned waveform inputs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import h5py
from scipy import signal
from .phase_screening import screen_cross
from .phase_stack import phase_interval_mask,goertzel_windows,wrap
from dispersion.common import sha,write_json

def audit(folder):
    dest=Path(folder);m=json.loads((dest/'phase_stack_complete.json').read_text(encoding='utf-8-sig'));gates=pd.read_csv(dest/'phase_gate_summary.csv')
    if m['module_sha256']!=sha(Path(__file__).with_name('phase_stack.py')):raise ValueError('Phase computation source changed')
    if m.get('phase_records_sha256') and m['phase_records_sha256']!=sha(Path(__file__).with_name('phase_records.py')):raise ValueError('Phase record/window geometry source changed')
    if m.get('phase_screening_sha256') and m['phase_screening_sha256']!=sha(Path(__file__).with_name('phase_screening.py')):raise ValueError('Phase screening source changed')
    if m.get('direction_sha256') and m['direction_sha256']!=sha(Path(__file__).with_name('direction.py')):raise ValueError('Direction computation source changed')
    if m.get('reconstruction_module_sha256') and m['reconstruction_module_sha256']!=sha(Path(__file__).with_name('reconstruction.py')):raise ValueError('Projection reconstruction source changed')
    if m.get('music_module_sha256') and m['music_module_sha256']!=sha(Path(__file__).with_name('music.py')):raise ValueError('MUSIC source changed')
    if m['aligned_h5_sha256']!=sha(dest/'aligned_waveforms.h5'):raise ValueError('Aligned waveform fingerprint mismatch')
    geometry=None
    if m.get('station_geometry_sha256'):
        if m['station_geometry_sha256']!=sha(dest/'direction/station_geometry.csv'):raise ValueError('Projection geometry changed')
        geometry=pd.read_csv(dest/'direction/station_geometry.csv').set_index('station')
    used=json.loads((dest/'phase_settings_used.json').read_text(encoding='utf-8-sig'))['settings']
    checked=0;maximum_error=0.;accepted=0;valid_count=0
    with h5py.File(dest/'aligned_waveforms.h5') as h:
        arrays={s:np.where(np.asarray(g['valid'],bool),np.asarray(g['data'],float),np.nan) for s,g in h['stations'].items()}
        fs=float(h.attrs['fs_hz'])
        for r in gates.itertuples():
            with np.load(dest/r.windows_file,allow_pickle=False) as z:
                cross=z['cross'];valid=z['valid'];mask=z['accepted'];starts=z['start_sample'];n=int(z['window_samples'])
                expected=np.isfinite(cross.real)&np.isfinite(cross.imag)&(abs(cross)>0)
                if not np.array_equal(valid,expected):raise ValueError('Stored base mask mismatch')
                mode=str(z['screening_mode']) if 'screening_mode' in z else 'curve'
                split=bool(z['split_branches']) if 'split_branches' in z else False
                amp=bool(z['amplitude_rejection']) if 'amplitude_rejection' in z else False
                # MATLAB mode's branch choice is the original curve center, not the preliminary center.
                curve_center=np.radians(r.predicted_phase_deg)
                eligible=z['direction_eligible'] if 'direction_applied' in z and bool(z['direction_applied']) else None
                checked_masks=screen_cross(cross,curve_center,float(z['phase_offset_low_rad']),float(z['phase_offset_high_rad']),mode,amp,split,eligible)
                if eligible is not None:
                    c=m['direction'];angle=z['pair_angle_deg'];upper=z['pair_angle_upper_deg']
                    if np.any(eligible&(angle>c['max_pair_angle_deg']+1e-8)):raise ValueError('Accepted pair angle exceeds threshold')
                    if c['conservative_angle'] and np.any(eligible&(upper>c['max_pair_angle_deg']+1e-8)):raise ValueError('Accepted angle envelope exceeds threshold')
                    if np.any(eligible&(z['angle_velocity_error_fraction']>c['max_angle_velocity_error_fraction']+1e-8)):raise ValueError('Accepted angular velocity error exceeds budget')
                expected=checked_masks['accepted']
                if not np.array_equal(mask,expected):raise ValueError('Stored phase mask mismatch')
                for sign in ('positive','negative'):
                    if 'accepted_'+sign in z and not np.array_equal(z['accepted_'+sign],checked_masks['accepted_'+sign]):raise ValueError('Stored branch mask mismatch')
                if int(mask.sum())!=r.accepted_windows or int(valid.sum())!=r.valid_windows:raise ValueError('Window count mismatch')
                if float(z['frequency_hz'])!=r.frequency_hz:raise ValueError('Frequency mismatch')
                a,b=str(z['pair']).split('--');pick=np.unique(np.linspace(0,len(starts)-1,min(3,len(starts)),dtype=int));ix=starts[pick]
                taper=np.ones(n) if m['taper']=='none' else signal.windows.hann(n,sym=True)
                actual_f=float(z['actual_frequency_hz']) if 'actual_frequency_hz' in z else r.frequency_hz
                ga=goertzel_windows(arrays[a],ix,n,actual_f,fs,taper);gb=goertzel_windows(arrays[b],ix,n,actual_f,fs,taper)
                if m.get('direction',{}).get('mode','off')!='off':
                    stable=int(np.floor(float(m.get('stabilization_s',used['stabilization_s']))*fs+.5))
                    for name,coeff in ((a,ga),(b,gb)):
                        g=h['stations/'+name];first=int(g.attrs.get('record_start_sample',0));last=int(g.attrs.get('record_end_sample_exclusive',len(arrays[name])))
                        good=(ix>=max(0,first+stable))&(ix+n<=min(len(arrays[name]),last));coeff[~good]=np.nan+1j*np.nan
                rec=gb*np.conj(ga)
                if 'raw_cross' in z:
                    from .reconstruction import align_cross
                    np.testing.assert_allclose(rec,z['raw_cross'][pick],rtol=1e-11,atol=1e-30,equal_nan=True)
                    reference=float(z['reference_projected_distance_m']);v=float(z['prior_true_velocity_m_s'])
                    if np.isfinite(reference) and np.isfinite(v):
                        delta=geometry.loc[b,['x_m','y_m']].to_numpy(float)-geometry.loc[a,['x_m','y_m']].to_numpy(float)
                        az=np.radians(z['propagation_azimuth_deg']);projected=delta[0]*np.sin(az)+delta[1]*np.cos(az)
                        np.testing.assert_allclose(projected,z['window_projected_distance_m'],rtol=1e-12,atol=1e-10,equal_nan=True)
                        azref=np.radians(float(z['confirmed_azimuth_deg']))
                        np.testing.assert_allclose(reference,delta[0]*np.sin(azref)+delta[1]*np.cos(azref),rtol=1e-12,atol=1e-10)
                        rebuilt,rotation=align_cross(z['raw_cross'],actual_f,z['window_projected_distance_m'],reference,v)
                        np.testing.assert_allclose(rotation,z['alignment_rotation_rad'],rtol=1e-12,atol=1e-12)
                        np.testing.assert_allclose(rebuilt,cross,rtol=1e-11,atol=1e-30,equal_nan=True)
                        if np.any(eligible & (z['window_projected_distance_m']*reference<=0)):raise ValueError('Projection sign disagrees with propagation')
                        from .direction import angular_distance
                        if np.any(eligible & (angular_distance(z['propagation_azimuth_deg'],float(z['confirmed_azimuth_deg']))>m['direction']['direction_group_width_deg']+1e-8)):raise ValueError('Window direction differs from confirmed direction')
                    rec=rec*np.exp(1j*z['alignment_rotation_rad'][pick])
                np.testing.assert_allclose(rec,cross[pick],rtol=1e-11,atol=1e-30,equal_nan=True)
                finite=np.isfinite(rec)&(abs(rec)>0)
                if finite.any():maximum_error=max(maximum_error,float(abs(wrap(np.angle(rec[finite])-np.angle(cross[pick][finite]))).max()))
                checked+=len(pick);accepted+=int(mask.sum());valid_count+=int(valid.sum())
    if accepted!=m['accepted_window_evaluations'] or valid_count!=m['valid_window_evaluations']:raise ValueError('Completion count mismatch')
    result=dict(pass_audit=True,gates_checked=len(gates),waveform_cross_recomputed=checked,accepted_window_evaluations=accepted,valid_window_evaluations=valid_count,
                max_recomputed_phase_error_rad=maximum_error,module_sha256=m['module_sha256'],aligned_h5_sha256=m['aligned_h5_sha256'])
    write_json(dest/'phase_audit.json',result);return result
