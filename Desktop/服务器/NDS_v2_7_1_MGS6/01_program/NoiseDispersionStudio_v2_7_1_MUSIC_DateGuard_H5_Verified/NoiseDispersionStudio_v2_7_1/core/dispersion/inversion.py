"""Optional fresh multistart conditional Rayleigh inversion; no previous models."""
from pathlib import Path
import json,warnings
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from .common import ensure_ready,csv,write_json

def model_from_log(x,nl,ratio,density):
    values=np.exp(x);h=np.r_[values[:nl],0];vs=values[nl:]
    return np.column_stack((h/1000,vs*ratio/1000,vs/1000,np.full(nl+1,density)))

def forward(model,f,mode,group=False):
    from disba import PhaseDispersion,GroupDispersion,DispersionError
    periods=1/np.asarray(f);order=np.argsort(periods)
    calc=(GroupDispersion if group else PhaseDispersion)(*model.T,algorithm="dunkin",dc=.001)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore");result=calc(periods[order],mode=mode,wave="rayleigh")
    except DispersionError as exc:
        raise ValueError("Missing modal root: "+str(exc)) from exc
    if len(result.velocity)!=len(f) or not np.isfinite(result.velocity).all(): raise ValueError("Missing modal root")
    v=np.empty(len(f));v[order]=result.velocity*1000
    if (v<=0).any() or (not group and (v>=model[-1,2]*1000).any()): raise ValueError("Invalid/untrapped modal root")
    return v

def independent_phase(model,f,mode):
    from pysurf96 import surf96
    periods=1/np.asarray(f);order=np.argsort(periods)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # pysurf96 1.0.1 accepts at most 60 periods per Fortran call. Evaluate
        # every requested root in sorted batches, never interpolate missing roots.
        ordered=periods[order]
        v=np.concatenate([surf96(*model.T,ordered[k:k+60],wave="rayleigh",mode=mode+1,
            velocity="phase",flat_earth=True) for k in range(0,len(ordered),60)])*1000
    if not np.isfinite(v).all() or (v<=0).any(): raise ValueError("Missing independent modal root")
    result=np.empty(len(f));result[order]=v;return result

def invert(c,selection_path):
    out,cache,meta=ensure_ready(c,"families_complete.json")
    selection=json.loads(Path(selection_path).read_text(encoding="utf-8-sig"))
    settings={"layers":[1,2],"starts":12,"max_evaluations":250,"vp_vs_ratio":2.0,"density_g_cm3":2.0,
        "thickness_bounds_m":[2.,1500.],"vs_bounds_m_s":[120.,3500.],"phase_rms_max_fraction":.03,
        "phase_point_max_fraction":.05,"group_median_max_fraction":.15,"require_independent_solver":True,
        "independent_max_difference_m_s":.5,"allow_background_unverified":False,"seed":9075,
        "finite_difference_log_step":.003}
    extra=set(selection)-{"settings","cases"}
    if extra: raise ValueError("Unknown selection keys: "+str(extra))
    for k,v in selection.get("settings",{}).items():
        if k not in settings: raise ValueError("Unknown inversion setting: "+k)
        settings[k]=v
    if not selection.get("cases"): raise ValueError("Specify at least one explicit cycle/mode case in selection.json; no automatic unique-cycle assumption")
    if any(not isinstance(n,int) or not 1<=n<=4 for n in settings["layers"]): raise ValueError("layers: finite layer counts 1..4")
    if settings["vp_vs_ratio"]<=np.sqrt(4/3) or settings["density_g_cm3"]<=0 or settings["starts"]<1: raise ValueError("Invalid elastic parameters")
    try:
        import pysurf96
        independent_available=True
    except ImportError: independent_available=False
    if settings["require_independent_solver"] and not independent_available:
        raise ImportError("Independent solver pysurf96 is required. Install requirements-inversion.txt or explicitly request exploratory unchecked fits")
    curves=pd.read_csv(out/"conditional_dispersion.csv");rng=np.random.default_rng(settings["seed"]);models=[];observed=[];predictions=[];caseids=set();attempts=[]
    for case in selection["cases"]:
        caseid=case["name"]
        if caseid in caseids: raise ValueError("Duplicate case name")
        caseids.add(caseid);ds=[];modes=[]
        if set(case)-{"name","bands"}: raise ValueError("Unknown case key")
        for band in case["bands"]:
            if set(band)-{"curve_id","mode","fmin_hz","fmax_hz"}: raise ValueError("Unknown band key")
            q=curves[curves.curve_id==band["curve_id"]].sort_values("frequency_hz")
            if "fmin_hz" in band: q=q[q.frequency_hz>=band["fmin_hz"]]
            if "fmax_hz" in band: q=q[q.frequency_hz<=band["fmax_hz"]]
            if len(q)<c["scan"]["min_points"] or q.frequency_hz.max()-q.frequency_hz.min()<c["scan"]["min_band_hz"]-1e-8: raise ValueError("Missing/too short selected curve")
            if not q.instrument_response_corrected.all(): raise ValueError("Structural inversion requires response-corrected inputs")
            if not settings["allow_background_unverified"] and not q.background_discriminated.all(): raise ValueError("Selected curve lacks >=20 valid wrong-time trials / background separation")
            mode=band["mode"]
            if not isinstance(mode,int) or mode<0: raise ValueError("Explicit nonnegative Rayleigh mode required (0=fundamental)")
            for old in ds:
                if old.pair_index.iloc[0]==q.pair_index.iloc[0] and set(np.round(old.frequency_hz,8))&set(np.round(q.frequency_hz,8)):
                    raise ValueError("Overlapping same-pair observations cannot be counted twice; select one fixed cohort over the band")
            ds.append(q);modes.append(mode)
        if not ds: raise ValueError("Case has no measured bands")
        y=np.concatenate([q.phase_velocity_m_s.to_numpy() for q in ds]);weights=np.concatenate([np.full(len(q),1/np.sqrt(len(ds)*len(q))) for q in ds])
        for nl in settings["layers"]:
            hlo,hhi=settings["thickness_bounds_m"];vlo,vhi=settings["vs_bounds_m_s"]
            bounds=(np.log([hlo]*nl+[vlo]*(nl+1)),np.log([hhi]*nl+[vhi]*(nl+1)))
            for attempt in range(settings["starts"]):
                vbase=float(np.median(y));vs=vbase*rng.uniform(.55,1.6,nl+1)
                vs[-1]=max(float(y.max())*rng.uniform(1.05,1.5),vs[-1])
                thick=np.exp(rng.uniform(np.log(max(hlo,5)),np.log(min(hhi,800)),nl))
                x0=np.clip(np.log(np.r_[thick,vs]),bounds[0]+1e-5,bounds[1]-1e-5)
                def residual(x):
                    try:
                        model=model_from_log(x,nl,settings["vp_vs_ratio"],settings["density_g_cm3"])
                        pred=np.concatenate([forward(model,q.frequency_hz.to_numpy(),mode) for q,mode in zip(ds,modes)])
                        return (pred/y-1)*weights
                    except (ValueError,ZeroDivisionError,RuntimeError): return np.full(len(y),10.)*weights
                fit=least_squares(residual,x0,bounds=bounds,max_nfev=settings["max_evaluations"],
                    diff_step=settings["finite_difference_log_step"],ftol=2e-7,xtol=2e-7,gtol=2e-7)
                model=model_from_log(fit.x,nl,settings["vp_vs_ratio"],settings["density_g_cm3"])
                error=residual(fit.x)/weights;rms=float(np.sqrt(np.sum((error*weights)**2)));maxerr=float(np.max(abs(error)))
                diagnostic=dict(case=caseid,finite_layers=nl,attempt=attempt,optimizer_status=int(fit.status),evaluations=int(fit.nfev),
                    phase_weighted_rms_fraction=rms,phase_max_fraction=maxerr,accepted=False,reason="phase_criterion_or_missing_root")
                attempts.append(diagnostic)
                if rms>settings["phase_rms_max_fraction"] or maxerr>settings["phase_point_max_fraction"]: continue
                pred=[];groups=[];inderror=np.nan;gerror=[];indok=False
                try:
                    for q,mode in zip(ds,modes):
                        pp=forward(model,q.frequency_hz.to_numpy(),mode);gg=forward(model,q.frequency_hz.to_numpy(),mode,True)
                        pred.append(pp);groups.append(gg);gerror.append(float(np.median(abs(gg/q.group_velocity_m_s.to_numpy()-1))))
                    if independent_available:
                        inderror=max(float(np.max(abs(independent_phase(model,q.frequency_hz.to_numpy(),mode)-pp))) for q,mode,pp in zip(ds,modes,pred))
                        indok=inderror<=settings["independent_max_difference_m_s"]
                except (ValueError,ZeroDivisionError,RuntimeError) as exc:
                    diagnostic.update(reason="group_or_independent_root_missing",numerical_failure=type(exc).__name__+": "+str(exc));continue
                groupok=max(gerror)<=settings["group_median_max_fraction"]
                diagnostic.update(maximum_group_median_fraction=max(gerror),independent_max_difference_m_s=inderror,
                    reason="group_criterion" if not groupok else "independent_solver_criterion")
                if not groupok or (settings["require_independent_solver"] and not indok): continue
                # Keep distinct structures; small residuals are not posterior probabilities.
                if any(r["case"]==caseid and r["finite_layers"]==nl and np.max(abs(np.log(np.asarray(r["parameters"])/np.exp(fit.x))))<.025 for r in models):
                    diagnostic["reason"]="duplicate_structure";continue
                mid=f"M{len(models)+1:05d}"
                diagnostic.update(accepted=True,reason="conditional_compatible",model_id=mid)
                models.append(dict(model_id=mid,case=caseid,finite_layers=nl,parameters=np.exp(fit.x).tolist(),
                    thickness_m=(model[:-1,0]*1000).tolist(),vs_m_s=(model[:,2]*1000).tolist(),vp_vs_ratio=settings["vp_vs_ratio"],density_g_cm3=settings["density_g_cm3"],
                    phase_weighted_rms_fraction=rms,phase_max_fraction=maxerr,maximum_group_median_fraction=max(gerror),
                    independent_max_difference_m_s=inderror if np.isfinite(inderror) else None,independent_solver_pass=indok,
                    verified_conditional_example=bool(groupok and indok),source_phase_verified=c["phase"]["source_phase_verified"],
                    absolute_cycles_or_modes_verified=False,posterior_probability_calibrated=False))
                for bandno,(q,mode,pp,gg) in enumerate(zip(ds,modes,pred,groups)):
                    for record in q.to_dict("records"):
                        observed.append(dict(model_id=mid,case=caseid,band_index=bandno,assumed_mode=mode,**record))
                fgrid=np.linspace(c["scan"]["min_hz"],c["scan"]["max_hz"],141)
                for mode in sorted(set(modes)):
                    # Evaluate each frequency independently to keep missing modal roots missing.
                    for freq in fgrid:
                        value=np.nan
                        try: value=float(forward(model,np.asarray([freq]),mode)[0])
                        except (ValueError,ZeroDivisionError,RuntimeError): pass
                        covered=any(m==mode and q.frequency_hz.min()<=freq<=q.frequency_hz.max() for q,m in zip(ds,modes))
                        predictions.append(dict(model_id=mid,case=caseid,frequency_hz=freq,assumed_mode=mode,phase_velocity_m_s=value,
                            is_observation=False,within_selected_observed_band=covered,root_available=bool(np.isfinite(value))))
            print("invert:",caseid,"layers",nl,"compatible examples",len(models),flush=True)
    csv(out/"inversion_attempts.csv",attempts)
    write_json(out/"conditional_models.json",dict(settings=settings,models=models,source="fresh_multistart_no_previous_models",probabilities_calibrated=False))
    csv(out/"selected_observed_dispersion.csv",observed,columns=None if observed else ["model_id","curve_id","frequency_hz","phase_velocity_m_s"])
    csv(out/"model_predictions.csv",predictions,columns=["model_id","case","frequency_hz","assumed_mode","phase_velocity_m_s","is_observation","within_selected_observed_band","root_available"])
    write_json(out/"inversion_complete.json",dict(cases=len(selection["cases"]),compatible_examples=len(models),independent_solver_available=independent_available,
               observed_and_predicted_points_in_separate_files=True,unique_structure_identified=False))
