"""Reject incompatible phase resumes and stale paired inversion inputs."""
import json
import pandas as pd
import pytest
from ndstudio import real_compare
from ndstudio.direction import settings
from dispersion.common import sha


def test_phase_resume_rejects_changed_window_parameters(tmp_path):
    options=dict(repeat_cycles=200.,slide_cycles=200.,taper='hann',frequency_method='exact',
        monotonic_only=True,direction={'mode':'filter'},prior={'enabled':True,'extension_hz':.1})
    waveform=tmp_path/'aligned_waveforms.h5';waveform.write_bytes(b'fingerprint-fixture')
    recorded=dict(options,direction=settings(options['direction']),
        module_sha256=sha(real_compare.Path(real_compare.__file__).with_name('phase_stack.py')),
        direction_sha256=sha(real_compare.Path(real_compare.__file__).with_name('direction.py')),
        music_module_sha256=sha(real_compare.Path(real_compare.__file__).with_name('music.py')),
        aligned_h5_sha256=sha(waveform))
    (tmp_path/'phase_stack_complete.json').write_text(json.dumps(recorded))
    assert real_compare.completed_phase(tmp_path,options)
    with pytest.raises(ValueError,match='参数'):
        real_compare.completed_phase(tmp_path,dict(options,repeat_cycles=20.))
    waveform.write_bytes(b'changed')
    with pytest.raises(ValueError,match='波形'):
        real_compare.completed_phase(tmp_path,options)


def test_pair_inversion_recomputes_after_curve_change(tmp_path,monkeypatch):
    from ndstudio import inversion
    calls=[]
    def fake_invert(table,selection,dest,minimum_points,minimum_band_hz):
        calls.append(table.phase_velocity_m_s.tolist())
        (dest/'inversion_complete.json').write_text(json.dumps({'compatible_models':1}))
    monkeypatch.setattr(inversion,'invert_tables',fake_invert)
    q=pd.DataFrame({'frequency_hz':[1.,1.1],'phase_velocity_m_s':[700.,650.]})
    scan=dict(min_points=6,min_band_hz=.25);selection={'cases':[]}
    real_compare.paired_inversion(q,selection,tmp_path,scan)
    real_compare.paired_inversion(q,selection,tmp_path,scan)
    assert len(calls)==1
    changed=q.copy();changed.loc[0,'phase_velocity_m_s']=710.
    real_compare.paired_inversion(changed,selection,tmp_path,scan)
    assert len(calls)==2
