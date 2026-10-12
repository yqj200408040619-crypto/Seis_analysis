"""The uploaded MATLAB branch/MAD framework with an optional curve phase gate."""
import numpy as np


def wrap(a):return np.angle(np.exp(1j*np.asarray(a)))


def amplitude_flags(amplitude):
    good=np.isfinite(amplitude);vals=amplitude[good]
    if not len(vals):return np.zeros(len(amplitude),bool),np.zeros(len(amplitude),bool),np.nan,np.nan
    med=np.median(vals);scale=1.4826*np.median(abs(vals-med));threshold=med+5*scale if scale>0 else 5*med
    delta=np.r_[0.,abs(np.diff(amplitude))]
    # MATLAB excludes the artificial first-window zero from delta statistics.
    dv=delta[1:];dv=dv[np.isfinite(dv)]
    dm=np.median(dv) if len(dv) else 0.;ds=1.4826*np.median(abs(dv-dm)) if len(dv) else 0.
    dt=dm+5*ds if ds>0 else 5*dm
    return amplitude>threshold,delta>dt,float(threshold),float(dt)


def screen_cross(cross,center,low,high,mode='curve',amplitude_rejection=True,split_branches=True,eligible=None):
    from .phase_stack import phase_interval_mask
    if mode not in ('curve','matlab'):raise ValueError('screening_mode must be curve or matlab')
    phase=np.angle(cross);amplitude=abs(cross)
    valid=np.isfinite(cross.real)&np.isfinite(cross.imag)&(amplitude>0)
    large,jump,ath,dth=amplitude_flags(amplitude)
    amplitude_valid=valid&~(large|jump) if amplitude_rejection else valid.copy()
    eligible=np.ones(len(cross),bool) if eligible is None else np.asarray(eligible,bool)
    if eligible.shape!=valid.shape:raise ValueError('Direction eligibility mask shape mismatch')
    positive=amplitude_valid&eligible&(phase>=0);negative=amplitude_valid&eligible&(phase<0)
    preliminary=[]
    for mask in (positive,negative):
        total=cross[mask].sum();preliminary.append(float(np.angle(total)) if abs(total)>0 else np.nan)
    if mode=='matlab':
        accepted_positive=positive&(abs(wrap(phase-preliminary[0]))<=np.pi/2)
        accepted_negative=negative&(abs(wrap(phase-preliminary[1]))<=np.pi/2)
    else:
        gate=phase_interval_mask(phase,center,low,high)
        accepted_positive=positive&gate;accepted_negative=negative&gate
    branch='nonnegative' if float(wrap(center))>=0 else 'negative'
    accepted=(accepted_positive if branch=='nonnegative' else accepted_negative) if split_branches or mode=='matlab' else accepted_positive|accepted_negative
    return dict(valid=valid,amplitude_valid=amplitude_valid,accepted=accepted,accepted_positive=accepted_positive,
                accepted_negative=accepted_negative,amplitude_large=large,amplitude_jump=jump,
                preliminary_positive_phase_rad=preliminary[0],preliminary_negative_phase_rad=preliminary[1],
                amplitude_threshold=ath,delta_threshold=dth,phase_sign_branch=branch if split_branches or mode=='matlab' else 'combined')


def geometric_fraction(center,low,high,branch='combined'):
    """Uniform circular-phase measure of a gate, including the sign restriction."""
    width=high-low
    if branch=='combined':return min(1.,width/(2*np.pi))
    if width>=2*np.pi:return .5
    a,b=center+low,center+high
    bounds=(0.,np.pi) if branch=='nonnegative' else (-np.pi,0.)
    first=int(np.floor((a-bounds[1])/(2*np.pi)))-1;last=int(np.ceil((b-bounds[0])/(2*np.pi)))+1
    length=sum(max(0.,min(b,bounds[1]+k*2*np.pi)-max(a,bounds[0]+k*2*np.pi)) for k in range(first,last+1))
    return length/(2*np.pi)


def accumulation_series(cross,mask,limit=20000):
    """Display a bounded set of points from the TRUE all-window cumulative sum."""
    count=int(mask.sum())
    if not count:return dict(count=np.empty(0,int),raw=np.empty(0,complex),unit=np.empty(0,complex))
    pick=np.unique(np.linspace(0,count-1,min(count,limit),dtype=np.int64))
    raw_out=np.empty(len(pick),complex);unit_out=np.empty(len(pick),complex);offset=0;raw=0j;unit=0j
    for lo in range(0,len(cross),65536):
        z=cross[lo:lo+65536][mask[lo:lo+65536]]
        if not len(z):continue
        sums=np.cumsum(z)+raw;unit_sums=np.cumsum(z/abs(z))+unit
        left=np.searchsorted(pick,offset);right=np.searchsorted(pick,offset+len(z))
        ix=pick[left:right]-offset;raw_out[left:right]=sums[ix];unit_out[left:right]=unit_sums[ix]
        raw=sums[-1];unit=unit_sums[-1];offset+=len(z)
    return dict(count=pick+1,raw=raw_out,unit=unit_out)
