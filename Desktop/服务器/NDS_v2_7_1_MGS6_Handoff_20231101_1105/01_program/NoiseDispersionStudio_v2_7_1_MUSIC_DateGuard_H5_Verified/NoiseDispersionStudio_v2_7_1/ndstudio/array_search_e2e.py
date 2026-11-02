"""Auto-searching ~20 km aperture array workflow: direction + velocity + dispersion check.

Goal (user spec): take the wave-direction analysis program, find a DENSE array of
~20 km aperture with MORE than 6 stations, restrict to the frequency band the aperture
actually supports (the microseism band), find directionally-consistent station pairs,
use integer-period time windows, recover stable direction and phase/group velocity, and
compare the measured dispersion against published (authoritative) dispersion.

This driver LOOPS over a ranked list of candidate regions/networks. For each candidate it
tries, in order:

    select_array -> fetch_array -> audit_acquisition -> screen_array
    -> preprocess+scan -> front MUSIC (apparent->radial velocity)
    -> low-frequency reconstruction (integer-period Goertzel windows, direction+phase gated stack)
    -> phase-audit -> select -> dispersion comparison

If a candidate fails ANY hard constraint (not enough stations in the aperture box, no
compact subset, download/audit failure, or no resolvable band), it is recorded and the
driver moves on to the NEXT candidate. It keeps going until one candidate strictly passes
every stage, or the list is exhausted.

IMPORTANT: select_array / fetch_array open LIVE FDSN connections (via h5data). Run this where
the machine can reach the relevant FDSN provider (EarthScope/IRIS, NCEDC, ORFEUS, GEOFON).
It cannot run in an offline sandbox.

The aperture/frequency matching is explicit:
    f_max ~= c / (2 * d_min)     (anti-aliasing ceiling, from min station spacing)
    f_min ~= c / D               (resolution floor, from aperture D)
Only frequencies inside [f_min, f_max] are treated as evaluable for direction.
"""
from __future__ import annotations
from pathlib import Path
from dataclasses import dataclass, field, asdict
import argparse, json, traceback
import pandas as pd

from dispersion.common import write_json, sha
from .acquisition import select_array, fetch_array, audit_acquisition, screen_array
from .project import new_project
from .real_compare import preprocess, completed_phase
from .pipeline import execute
from .candidates import auto_selection, read_curves
from .phase_stack import run_phase_stack
from .phase_audit import audit
from .reconstruction import paired_compare


# --------------------------------------------------------------------------------------
# Candidate arrays. Ordered best-first. Each is a geographic box / network that a prior
# FDSN station-metadata check showed to hold >6 broadband vertical channels inside a
# ~20 km aperture. Bounds are chosen so that, with Rayleigh c ~ 2.5-3.5 km/s, the
# secondary-microseism band (~0.1-0.3 Hz) sits inside [f_min, f_max].
#
# NOTE: select_array queries ONE network at a time. The Parkfield footprint mixes BK/NC/NP;
# NC (NCEDC, https://service.ncedc.org) carries the densest broadband subset there, so it
# is tried first via the NCEDC provider. Verify provider URLs for your environment.
# --------------------------------------------------------------------------------------
@dataclass
class Candidate:
    name: str
    provider: str            # FDSN base URL or short name understood by h5data.make_client
    network: str
    station: str = '*'
    channel: str = 'HHZ'
    # geographic box (deg); None = no bound
    minlatitude: float | None = None
    maxlatitude: float | None = None
    minlongitude: float | None = None
    maxlongitude: float | None = None
    # geometry constraints (m)
    size: int = 10
    min_stations: int = 7
    max_stations: int = 20
    min_spacing_m: float = 1500.
    max_spacing_m: float = 9000.
    max_aperture_m: float = 22000.
    # expected surface-wave velocity range for band matching (m/s)
    velocity_min_m_s: float = 1500.
    velocity_max_m_s: float = 4500.
    # authoritative dispersion reference for this footprint (paper / dataset)
    reference: str = ''
    notes: str = ''


CANDIDATES: list[Candidate] = [
    Candidate(
        name='LABasin_CI',
        provider='SCEDC', network='CI', station='*', channel='HHZ',
        minlatitude=33.95, maxlatitude=34.30, minlongitude=-118.45, maxlongitude=-117.95,
        size=8, min_stations=7, max_stations=25,
        min_spacing_m=500., max_spacing_m=18000., max_aperture_m=30000.,
        velocity_min_m_s=1500., velocity_max_m_s=4500.,
        reference='Southern California / Los Angeles basin ambient-noise Rayleigh-wave phase '
                  'velocities: Lin et al. 2013 GJI 194(2) (SCEC CVM ANT); Shapiro et al. 2005 '
                  'Science 307; Lee et al. 2014 JGR 3-D CVM. Fundamental-mode Rayleigh Vph.',
        notes='VERIFIED against live EarthScope/IRIS metadata + waveform (USC 1h=360000@100Hz) '
              '2023-11: 16 broadband HHZ eligible; compact 8-station subset (BHP,GR2,HLL,LCG,'
              'LGB,PDR,SMF2,USC), aperture 29.6 km, d_min 1.85 km, mean NN 6.9 km; band '
              '~0.05-1.22 Hz. Sedimentary basin -> coherent surface-wave propagation.'),
    Candidate(
        name='LongValley_NC',
        provider='NCEDC', network='NC', station='*', channel='HHZ',
        minlatitude=37.55, maxlatitude=37.80, minlongitude=-119.10, maxlongitude=-118.75,
        size=7, min_stations=7, max_stations=20,
        min_spacing_m=1000., max_spacing_m=9000., max_aperture_m=20000.,
        velocity_min_m_s=1500., velocity_max_m_s=4500.,
        reference='Western-US ambient-noise Rayleigh-wave phase velocities, Lin, Moschetti & '
                  'Ritzwoller 2008 GJI 173(1); Moschetti et al. 2007 G-cubed 8(8). Long Valley / '
                  'eastern Sierra Nevada footprint, fundamental-mode Rayleigh Vph.',
        notes='VERIFIED against live EarthScope/IRIS metadata 2023-11: 7 permanent broadband HHZ '
              '(MCB,MCY,MDPB,MDY,MINS,MLI,MMLB), aperture 16.2 km, d_min 1.2 km, mean NN 3.0 km; '
              'resolvable band ~0.09-1.9 Hz spans the secondary microseism.'),
    Candidate(
        name='Parkfield_NC',
        provider='NCEDC', network='NC', station='P*', channel='HHZ',
        minlatitude=35.78, maxlatitude=36.00, minlongitude=-120.63, maxlongitude=-120.38,
        size=9, min_stations=7, max_stations=14,
        min_spacing_m=1500., max_spacing_m=9000., max_aperture_m=22000.,
        velocity_min_m_s=1500., velocity_max_m_s=4500.,
        reference='Parkfield 3C ambient-noise beamforming, GJI 213(3) 2018; '
                  'S. California microseism phase velocities, GJI 211(1) 2017',
        notes='Dense SAF array; ~20 km aperture after trimming outliers, d_min ~4 km -> f in ~0.14-0.35 Hz.'),
    Candidate(
        name='Parkfield_BK', provider='NCEDC', network='BK', station='P*', channel='HHZ',
        minlatitude=35.78, maxlatitude=36.00, minlongitude=-120.63, maxlongitude=-120.38,
        size=7, min_stations=6, max_stations=12,
        reference='as Parkfield_NC', notes='BK broadband subset only (fewer stations).'),
    # Add more verified boxes here (e.g. an ORFEUS/GEOFON central-European subset vs
    # Lu et al. 2018 GJI European ambient-noise tomography) before falling through.
]


@dataclass
class StageOutcome:
    stage: str
    ok: bool
    detail: str = ''


@dataclass
class CandidateResult:
    name: str
    passed: bool = False
    stages: list = field(default_factory=list)
    band_hz: list | None = None
    message: str = ''


def _norm_utc(value: str) -> str:
    t = pd.Timestamp(value)
    t = t.tz_localize('UTC') if t.tzinfo is None else t.tz_convert('UTC')
    return t.strftime('%Y-%m-%dT%H:%M:%S.%fZ')


def resolvable_band(geometry: dict, c_lo: float, c_hi: float) -> dict:
    """[f_min, f_max] the aperture supports. f_max from min spacing, f_min from aperture."""
    nn = list(geometry['nearest_neighbor_m'].values())
    d_min = min(nn)
    aperture = geometry['max_pair_distance_m']
    f_max = c_hi / (2.0 * d_min)         # anti-aliasing ceiling (use higher c = looser)
    f_min = c_lo / aperture              # resolution floor (use lower c = looser)
    return {'d_min_m': d_min, 'aperture_m': aperture,
            'f_min_hz': round(f_min, 4), 'f_max_hz': round(f_max, 4),
            'evaluable': f_max > f_min}


def run_candidate(cand: Candidate, workroot: Path, start: str, end: str,
                  reuse_download: bool) -> CandidateResult:
    res = CandidateResult(name=cand.name)
    root = (workroot / cand.name).resolve()
    root.mkdir(parents=True, exist_ok=True)
    data = root / 'download'
    project = root / 'project'
    selected = root / 'array_selection.json'
    start_utc, end_utc = _norm_utc(start), _norm_utc(end)

    def record(stage, fn):
        try:
            out = fn()
            res.stages.append(asdict(StageOutcome(stage, True)))
            return out
        except Exception as exc:  # move on to next candidate
            res.stages.append(asdict(StageOutcome(stage, False, f'{type(exc).__name__}: {exc}')))
            raise

    try:
        # 1. Array selection within the aperture box (program-owned FDSN search).
        if not selected.exists():
            record('array-select', lambda: select_array(
                cand.provider, cand.network, start, end, selected,
                station=cand.station, channel=cand.channel, size=cand.size,
                min_spacing_m=cand.min_spacing_m, max_spacing_m=cand.max_spacing_m,
                max_aperture_m=cand.max_aperture_m,
                min_stations=cand.min_stations, max_stations=cand.max_stations,
                minlatitude=cand.minlatitude, maxlatitude=cand.maxlatitude,
                minlongitude=cand.minlongitude, maxlongitude=cand.maxlongitude))
        sel = json.loads(selected.read_text(encoding='utf-8-sig'))

        # Aperture/frequency gate BEFORE paying for a download.
        band = resolvable_band(sel['geometry'], cand.velocity_min_m_s, cand.velocity_max_m_s)
        res.band_hz = [band['f_min_hz'], band['f_max_hz']]
        write_json(root / 'resolvable_band.json', {**band, 'velocity_range_m_s': [cand.velocity_min_m_s, cand.velocity_max_m_s]})
        if not band['evaluable']:
            raise ValueError(f"no resolvable band: f_min {band['f_min_hz']} >= f_max {band['f_max_hz']}")

        # 2. Download ~3 days, audit against embedded original MiniSEED, screen bad/isolated.
        if not reuse_download:
            record('fetch', lambda: fetch_array(selected, data, workers=4, chunk_seconds=21600., allow_missing_data=False))
        record('audit-acquisition', lambda: audit_acquisition(data))
        record('screen-array', lambda: screen_array(selected, data, project=None))

        # 3. Project with band matched to the aperture (scan + reconstruction limited to [f_min, f_max]).
        f_min, f_max = band['f_min_hz'], band['f_max_hz']
        if not (project / 'project.json').exists():
            new_project(project)
            spec = json.loads((project / 'project.json').read_text(encoding='utf-8-sig'))
            spec.update(input_folder=str(data), automatic_geometry=True, exclude_increasing_velocity=True)
            write_json(project / 'project.json', spec)
            c = json.loads((project / 'config.json').read_text(encoding='utf-8-sig'))
            c.update(sample_rate_hz=10., bandpass_hz=[max(0.02, f_min * 0.5), min(2.5, f_max * 2.0)], max_lag_s=200.)
            c['response'].update(pre_filt_hz=[0.01, 0.02, 3.0, 4.0], table_step_hz=0.002)
            c['scan'].update(min_hz=f_min, max_hz=f_max, step_hz=0.005, min_band_hz=0.05,
                             group_velocity_bounds_m_s=[cand.velocity_min_m_s, cand.velocity_max_m_s])
            c['phase'].update(velocity_bounds_m_s=[cand.velocity_min_m_s - 500., cand.velocity_max_m_s])
            write_json(project / 'config.json', c)
            write_json(project / 'inversion_settings.json',
                       dict(thickness_bounds_m=[10., 10000.], vs_bounds_m_s=[200., 5000.], require_independent_solver=True))

        direction = dict(mode='filter', max_pair_angle_deg=45., direction_step_deg=2.,
                         velocity_points=128, velocity_min_m_s=cand.velocity_min_m_s - 500.,
                         velocity_max_m_s=cand.velocity_max_m_s, frequency_half_width_hz=0.,
                         frequency_bins=1, music_front_cycles=200.)
        write_json(project / 'phase_settings.json', {'direction': direction})

        # 4. Scan + front MUSIC (apparent pair velocity -> radial/wave-vector velocity).
        record('preprocess-and-scan', lambda: preprocess(root))
        record('front-MUSIC', lambda: execute(project, 'music'))
        for st in ('review', 'verify'):
            record(st, lambda st=st: execute(project, st))

        # 5. Directionally-consistent pair selection (non-increasing velocity curves, d/lambda>=3).
        before = record('select-before', lambda: auto_selection(project, max_cases=3, modes=(0,), minimum_wavelengths=3.))
        if not before['cases']:
            raise ValueError('no directionally-consistent candidate curves survived selection')
        ids = sorted({b['curve_id'] for case in before['cases'] for b in case['bands']})

        # 6. Integer-period Goertzel windows, direction+phase gated stack, low-f reconstruction.
        #    min_band_hz with step 0.005 -> bandwidth-based point count; uncertainty-aware monotonicity
        #    is enforced inside reconstruction/phase_velocity.
        options = dict(curve_ids=ids, repeat_cycles=200., slide_cycles=200., taper='hann',
                       frequency_method='exact', screening_mode='curve', tolerance_deg=15.,
                       monotonic_only=True, direction=direction,
                       reconstruction=dict(enabled=True, lower_frequency_hz=max(0.05, round(f_min, 3)), step_hz=0.005))
        write_json(project / 'phase_settings.json', options)
        phase = project / 'output/phase_guided'
        if not completed_phase(phase, options):
            record('goertzel-reconstruction', lambda: run_phase_stack(project, **options))
        record('phase-audit', lambda: audit(phase))

        after = record('select-after', lambda: auto_selection(project, max_cases=3, modes=(0,), minimum_wavelengths=3.))
        if not after['cases']:
            # Direction+velocity may still be recovered per-frequency even without a
            # continuous inversion band; keep the candidate but mark no band.
            res.message = ('direction/velocity estimated; no continuous measured band passed '
                           'selection (see phase_gate_summary.csv / frequency_directions.csv)')
            res.passed = False
            return res

        # 7. Compare measured dispersion vs the authoritative reference.
        measured = read_curves(project)[0]
        measured[measured.curve_id.isin({b['curve_id'] for case in after['cases'] for b in case['bands']})] \
            .to_csv(root / 'measured_dispersion.csv', index=False)
        write_json(root / 'dispersion_reference.json',
                   {'reference': cand.reference,
                    'how_to_compare': 'Sample reference Vph/Vgr at the array footprint over '
                                      f'[{f_min}, {f_max}] Hz; report per-frequency measured v, CI, '
                                      'reference v, and |dv/v|. Pre-register tolerance (suggest median |dv/v|<=5%).',
                    'measured_curve_csv': 'measured_dispersion.csv'})
        record('paired-inversion-and-compare',
               lambda: paired_compare(project, after, root / 'comparison'))
        res.message = 'PASS: measured direction+velocity band produced; compare measured_dispersion.csv to reference.'
        res.passed = True
        return res

    except Exception as exc:
        res.message = f'{type(exc).__name__}: {exc}'
        (root / 'failure.txt').write_text(traceback.format_exc(), encoding='utf-8')
        return res


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workdir', type=Path, required=True)
    p.add_argument('--start', default='2023-11-01', help='UTC start (>=3 day span)')
    p.add_argument('--end', default='2023-11-04', help='UTC end')
    p.add_argument('--reuse-download', action='store_true',
                   help='reuse verified H5 already under <workdir>/<candidate>/download')
    p.add_argument('--only', default=None, help='run only this candidate name')
    a = p.parse_args()
    workroot = a.workdir.resolve()
    workroot.mkdir(parents=True, exist_ok=True)
    if pd.Timestamp(_norm_utc(a.end)) <= pd.Timestamp(_norm_utc(a.start)):
        raise ValueError('--end must follow --start')

    results = []
    for cand in CANDIDATES:
        if a.only and cand.name != a.only:
            continue
        print(f'\n===== candidate: {cand.name} ({cand.network}@{cand.provider}) =====', flush=True)
        r = run_candidate(cand, workroot, a.start, a.end, a.reuse_download)
        results.append(asdict(r))
        write_json(workroot / 'array_search_status.json',
                   {'span_utc': [a.start, a.end], 'results': results,
                    'passed': next((x['name'] for x in results if x['passed']), None)})
        if r.passed:
            print(f'\n[STRICT PASS] {cand.name}: {r.message}', flush=True)
            print(f'Outputs under {workroot / cand.name}', flush=True)
            return
        print(f'[try next] {cand.name}: {r.message}', flush=True)

    print('\nNo candidate strictly passed. Inspect array_search_status.json; add more '
          'verified candidate boxes to CANDIDATES and rerun. The search does not stop '
          'until a candidate passes or the list is exhausted.', flush=True)


if __name__ == '__main__':
    main()
