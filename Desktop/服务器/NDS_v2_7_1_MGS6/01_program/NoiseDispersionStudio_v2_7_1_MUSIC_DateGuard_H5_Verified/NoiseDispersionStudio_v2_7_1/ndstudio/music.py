"""Narrowband MUSIC on signed horizontal slowness, shared by both stages.

DFT convention: exp(-i wt); propagation steering: exp(-i 2 pi f x.s).
Production spectra have independent DPSS tapered snapshots at the SAME target
frequency. A single-vector input is supported for noiseless numerical checks.
Grid envelopes and temporal quantiles are resolution diagnostics, not CIs.
"""
import numpy as np
from scipy.optimize import minimize


def source_order(eigenvalues, snapshots, c):
    maximum=min(c['music_max_sources'],len(eigenvalues)-1)
    if c['music_sources']:
        return min(c['music_sources'],len(eigenvalues)-1)
    if snapshots<2:
        return 1
    e=np.maximum(eigenvalues[::-1],max(eigenvalues[-1],1e-30)*1e-8)
    scores=[]
    for k in range(maximum+1):
        tail=e[k:];ratio=np.exp(np.log(tail).mean())/tail.mean()
        scores.append(-snapshots*len(tail)*np.log(max(ratio,1e-30))+.5*k*(2*len(e)-k)*np.log(snapshots))
    order=int(np.argmin(scores))
    # Finite DPSS bandwidth gives a weak second eigenvalue even for one delayed broadband wave.
    # MDL alone treats that deterministic leakage as a second physical wavefront.
    strong=int(np.sum(e/np.sum(e)>=c['music_min_source_fraction']))
    return min(order,max(1,strong)) if order else 0


def estimate(coefficients, xy, frequencies, c):
    from .direction import grid,geometry_ratio,peak_diagnostics
    z=np.asarray(coefficients,complex);xy=np.asarray(xy,float);freq=np.asarray(frequencies,float)
    if z.ndim==3:z=z[...,None]
    if z.ndim!=4 or z.shape[:2]!=(len(freq),len(xy)) or not np.isfinite(xy).all():
        raise ValueError('MUSIC频谱形状须为(频点,台站,窗,快照)，坐标须为有限数值')
    if len(freq)<1 or not np.isfinite(freq).all() or (freq<=0).any():raise ValueError('MUSIC频率必须为正有限值')
    nw=z.shape[2];sv,az,slow=grid(c);na=len(az);nv=len(slow);central=len(freq)//2
    keys=['azimuth_deg','velocity_m_s','slow_e_s_m','slow_n_s_m','coherence','secondary_ratio',
          'direction_uncertainty_deg','velocity_envelope_min_m_s','velocity_envelope_max_m_s','phase_rms_deg',
          'grid_half_cell_phase_rad','music_peak','music_eigenvalue_ratio']
    out={k:np.full(nw,np.nan) for k in keys}
    out.update(reason=np.full(nw,'few_stations',dtype='U32'),station_count=np.zeros(nw,int),geometry_ratio=np.zeros(nw),
               music_source_count=np.zeros(nw,int),music_snapshot_count=np.full(nw,z.shape[3],int))
    mask=np.all(np.isfinite(z.real)&np.isfinite(z.imag),axis=(0,3))&np.all(np.sum(abs(z)**2,axis=3)>0,axis=0)
    patterns,inverse=np.unique(mask.T,axis=0,return_inverse=True);examples=[]
    for pattern,valid in enumerate(patterns):
        ids=np.flatnonzero(inverse==pattern);m=int(valid.sum());out['station_count'][ids]=m
        if m<c['min_stations']:continue
        coords=xy[valid];ratio=geometry_ratio(coords);out['geometry_ratio'][ids]=ratio
        if ratio<c['min_geometry_ratio']:out['reason'][ids]='geometry';continue
        coords=coords-coords.mean(axis=0);aperture=np.max(np.linalg.norm(coords[:,None]-coords[None,:],axis=-1))
        # Work one bounded window batch at a time, never build window x grid x station cubes.
        steering=[np.exp(-2j*np.pi*f*(sv@coords.T))/np.sqrt(m) for f in freq]
        for start in range(0,len(ids),64):
            ix=ids[start:start+64];spectra=z[:,:,ix,:][:,valid]
            covariances=[];eigvecs=[];orders=[];eigvals=[]
            for fi in range(len(freq)):
                x=spectra[fi].transpose(1,0,2)
                energy=np.mean(abs(x)**2,axis=-1)
                x=x/np.sqrt(np.maximum(energy,1e-300))[...,None]
                R=np.einsum('wsk,wrk->wsr',x,x.conj())/z.shape[3]
                R=(R+R.conj().transpose(0,2,1))/2
                # Scalar loading stabilises MDL noise estimates without changing eigenvectors.
                R+=np.eye(m)[None]*c['music_diagonal_loading']
                values,U=np.linalg.eigh(R)
                order=np.array([source_order(e,z.shape[3],c) for e in values])
                covariances.append(R);eigvecs.append(U);eigvals.append(values);orders.append(order)
            power=np.zeros((len(ix),len(sv)))
            for fi,a in enumerate(steering):
                # Signal-space complement equals the noise-space norm. Group model
                # orders to scan many windows in one BLAS call, with bounded memory.
                for q in np.unique(orders[fi]):
                    group=np.flatnonzero(orders[fi]==q)
                    if q==0:
                        power[group]+=1.;continue
                    U=eigvecs[fi][group,:,-q:].transpose(1,0,2).reshape(m,-1)
                    projection=(a.conj()@U).reshape(len(sv),len(group),q)
                    noise=np.maximum(1.-np.sum(abs(projection)**2,axis=2),1e-12)
                    power[group]+=(1/noise).T
            power/=len(freq)
            peak=power.argmax(axis=1);best=power[np.arange(len(ix)),peak]
            refined=sv[peak].copy();coherence=np.empty(len(ix));rms=np.empty(len(ix));refpeak=best.copy()
            # Refine the MUSIC noise-subspace objective, never fall back to phase-line fitting.
            scale=c['velocity_min_m_s'];Q=[]
            for w in range(len(ix)):
                matrices=[]
                for fi in range(len(freq)):
                    q=orders[fi][w];En=eigvecs[fi][w,:,:m-q];matrices.append(En@En.conj().T)
                def objective(u):
                    value=0.;derivative=np.zeros(2)
                    for fi,P in enumerate(matrices):
                        a=np.exp(-2j*np.pi*freq[fi]*(coords@(u/scale)))/np.sqrt(m)
                        Pa=P@a;value+=float(np.vdot(a,Pa).real)
                        da=(-2j*np.pi*freq[fi]/scale)*coords*a[:,None]
                        derivative+=2*np.real(da.conj().T@Pa)
                    return value/len(freq),derivative/len(freq)
                fit=minimize(objective,refined[w]*scale,jac=True,method='BFGS',options=dict(gtol=1e-10,maxiter=60))
                candidate=fit.x/scale;v=1/max(np.linalg.norm(candidate),1e-30)
                if c['velocity_min_m_s']<v<c['velocity_max_m_s'] and objective(fit.x)[0]<=objective(refined[w]*scale)[0]+1e-12:
                    refined[w]=candidate
                a=np.exp(-2j*np.pi*freq[central]*(coords@refined[w]))/np.sqrt(m)
                R=covariances[central][w]-np.eye(m)*c['music_diagonal_loading']
                coherence[w]=np.clip(np.vdot(a,R@a).real/m,0,1)
                u=eigvecs[central][w,:,-1];res=np.angle(u*np.conj(a));res=np.angle(np.exp(1j*(res-np.angle(np.sum(np.exp(1j*res))))))
                rms[w]=np.degrees(np.sqrt(np.mean(res**2))) if orders[central][w]<=1 else 0.
                refpeak[w]=1/max(objective(refined[w]*scale)[0],1e-12)
            velocity=1/np.linalg.norm(refined,axis=1);theta=np.degrees(np.arctan2(refined[:,0],refined[:,1]))%360
            uncertainty,vmin,vmax,secondary=peak_diagnostics(np.ascontiguousarray(power),az,slow,theta,velocity,best,
                c['envelope_fraction'],max(15.,2*c['direction_step_deg']))
            uncertainty=np.maximum(uncertainty,c['direction_step_deg']/2)
            vmin=np.minimum(vmin,velocity);vmax=np.maximum(vmax,velocity)
            halfcell=2*np.pi*max(freq)*aperture*np.hypot((slow[1]-slow[0])/2,slow[peak//na]*np.radians(c['direction_step_deg'])/2)
            reason=np.full(len(ix),'accepted',dtype='U32')
            reason[uncertainty>c['max_direction_uncertainty_deg']]='uncertain'
            reason[secondary>=c['max_secondary_peak_ratio']]='ambiguous'
            reason[(velocity<=c['velocity_min_m_s']*1.001)|(velocity>=c['velocity_max_m_s']*.999)]='boundary'
            reason[halfcell>np.pi/2]='grid_resolution'
            reason[(coherence<c['min_coherence'])|(rms>c['phase_rms_max_deg'])]='low_coherence'
            reason[orders[central]==0]='no_signal'
            # A single phase correction must not silently combine equally strong wavefronts.
            reason[orders[central]>1]='multiple_sources'
            values=eigvals[central];eigenratio=values[:,-1]/np.maximum(values[:,:-1].mean(axis=1),1e-30)
            for key,value in dict(azimuth_deg=theta,velocity_m_s=velocity,slow_e_s_m=refined[:,0],slow_n_s_m=refined[:,1],
                coherence=coherence,secondary_ratio=secondary,direction_uncertainty_deg=uncertainty,
                velocity_envelope_min_m_s=vmin,velocity_envelope_max_m_s=vmax,phase_rms_deg=rms,
                grid_half_cell_phase_rad=halfcell,music_peak=refpeak,music_eigenvalue_ratio=eigenratio,
                music_source_count=orders[central]).items():out[key][ix]=value
            out['reason'][ix]=reason
            if len(examples)<2:
                examples.append(dict(window_index=int(ix[0]),power=(power[0]/best[0]).reshape(nv,na),azimuth=az,velocity=1/slow,
                    reference_slowness=refined[0],coordinates=coords,frequencies=freq))
    out['quality']=out['reason']=='accepted'
    width=c['direction_group_width_deg'];out['direction_group']=np.floor((np.nan_to_num(out['azimuth_deg'])+width/2)/width).astype(int)%round(360/width)
    return out,examples
