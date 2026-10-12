"""Explicit source-phase/cycle families, observed points only; no gap filling."""
import itertools
import numpy as np
import pandas as pd
from .common import csv,write_json,ensure_ready

CURVE_COLUMNS=["curve_id","segment_id","cohort","pair_index","pair","branch","frequency_hz","phase_velocity_m_s",
    "group_velocity_m_s","wavenumber_rad_m","distance_over_wavelength","constant_cycle","source_phase_rad",
    "source_phase_verified","absolute_cycle_resolved","rayleigh_mode_verified","instrument_response_corrected",
    "background_discriminated","conditional_interpretation","is_model_prediction"]

def families(c):
    out,cache,meta=ensure_ready(c,"scan_complete.json")
    measured=pd.read_csv(out/"measured_phase_segments.csv");cat=pd.read_csv(out/"measured_segment_catalog.csv")
    pairs=pd.read_csv(out/"pairs.csv").set_index("pair_index")
    control=pd.read_csv(out/"timeshift_summary.csv").set_index("segment_id") if (out/"controls_complete.json").exists() else None
    rows=[];catalog=[];prepared={};p=c["phase"];low,high=p["velocity_bounds_m_s"]
    for sid,q in measured.groupby("segment_id",sort=False):
        q=q.sort_values("frequency_hz");f=q.frequency_hz.to_numpy();theta=q.measured_unwrapped_phase_rad.to_numpy();r=float(pairs.loc[int(q.pair_index.iloc[0]),"distance_m"])
        background=bool(control.loc[sid,"background_threshold_exceeded"]) if control is not None and sid in control.index else False
        for iphi,phi in enumerate(p["source_phase_rad"]):
            for cycle in range(int(p["cycle_min"]),int(p["cycle_max"])+1):
                kr=phi-theta+2*np.pi*cycle
                if (kr<=0).any(): continue
                velocity=2*np.pi*f*r/kr
                if not np.isfinite(velocity).all() or (velocity<low).any() or (velocity>high).any(): continue
                cid=f"{sid}_S{iphi}_N{cycle}";prepared[cid]=(f,velocity,q.group_velocity_m_s.to_numpy(),q)
                catalog.append(dict(curve_id=cid,segment_id=sid,cohort=q.cohort.iloc[0],pair_index=int(q.pair_index.iloc[0]),branch=q.branch.iloc[0],
                    source_phase_rad=phi,constant_cycle=cycle,fmin_hz=float(f[0]),fmax_hz=float(f[-1]),points=len(f),
                    background_discriminated=background,instrument_response_corrected=meta["response_corrected"],
                    source_phase_verified=p["source_phase_verified"],absolute_cycle_resolved=False,rayleigh_mode_verified=False,
                    median_distance_over_wavelength=float(np.median(kr/(2*np.pi))),compatible_other_pairs=0))
                for k,v,phasekr,u in zip(f,velocity,kr,q.group_velocity_m_s):
                    rows.append(dict(curve_id=cid,segment_id=sid,cohort=q.cohort.iloc[0],pair_index=int(q.pair_index.iloc[0]),pair=q.pair.iloc[0],branch=q.branch.iloc[0],
                        frequency_hz=k,phase_velocity_m_s=v,group_velocity_m_s=u,wavenumber_rad_m=phasekr/r,distance_over_wavelength=phasekr/(2*np.pi),
                        constant_cycle=cycle,source_phase_rad=phi,source_phase_verified=p["source_phase_verified"],absolute_cycle_resolved=False,
                        rayleigh_mode_verified=False,instrument_response_corrected=meta["response_corrected"],background_discriminated=background,
                        conditional_interpretation=True,is_model_prediction=False))
    frame=pd.DataFrame(catalog);links=[];connections=[]
    # Compare entire overlapping bands; no per-frequency choice of N or windows.
    if len(frame):
        groups={sid:q.curve_id.tolist() for sid,q in frame.groupby("segment_id")}
        for a,b in itertools.combinations(cat.to_dict("records"),2):
            if a["cohort"]!=b["cohort"]: continue
            if a["pair_index"]==b["pair_index"]:
                if a["branch"]!=b["branch"]:continue
            lo=max(a["fmin_hz"],b["fmin_hz"]);hi=min(a["fmax_hz"],b["fmax_hz"])
            qa=measured[measured.segment_id==a["segment_id"]].sort_values("frequency_hz")
            qb=measured[measured.segment_id==b["segment_id"]].sort_values("frequency_hz")
            common=qa.merge(qb,on="frequency_hz",suffixes=("_a","_b"))
            overlap=hi-lo>=c["scan"]["min_band_hz"]-1e-8 and len(common)>=c["scan"]["min_points"]
            if a["pair_index"]==b["pair_index"]:
                shift=0;error=np.nan;passed=False
                if overlap:
                    shift=int(np.rint(np.median((common.measured_unwrapped_phase_rad_a-common.measured_unwrapped_phase_rad_b)/(2*np.pi))))
                    error=float(np.max(abs(common.measured_unwrapped_phase_rad_b+2*np.pi*shift-common.measured_unwrapped_phase_rad_a)))
                    arrival=float(np.median(abs(common.arrival_s_a-common.arrival_s_b)/common.arrival_s_a))
                    passed=error<=c["validation"]["phase_route_max_rad"] and arrival<=c["validation"]["arrival_fraction_max"]
                connections.append(dict(segment_a=a["segment_id"],segment_b=b["segment_id"],same_cohort=True,overlap_points=len(common),
                    phase_alignment_add_to_b_rad=2*np.pi*shift,cycle_b_equals_cycle_a_minus=shift,maximum_aligned_phase_error_rad=error,
                    measured_overlap_connection_pass=passed,gap_interpolated=False))
                continue
            if not overlap: continue
            for ca,cb in itertools.product(groups.get(a["segment_id"],[]),groups.get(b["segment_id"],[])):
                fa,va,ua,_=prepared[ca];fb,vb,ub,_=prepared[cb]
                use=(fa>=lo-1e-8)&(fa<=hi+1e-8);fc=fa[use]
                diff=abs(va[use]-np.interp(fc,fb,vb))/((va[use]+np.interp(fc,fb,vb))/2)
                gd=abs(ua[use]-np.interp(fc,fb,ub))/((ua[use]+np.interp(fc,fb,ub))/2)
                md=float(np.median(diff));mx=float(diff.max());gm=float(np.median(gd))
                if md<=p["link_median_fraction"] and mx<=p["link_max_fraction"] and gm<=p["link_group_median_fraction"]:
                    links.append(dict(curve_a=ca,curve_b=cb,segment_a=a["segment_id"],segment_b=b["segment_id"],cohort=a["cohort"],branch=a["branch"],
                        overlap_min_hz=lo,overlap_max_hz=hi,common_points=len(fc),phase_velocity_median_fraction=md,phase_velocity_max_fraction=mx,
                        group_velocity_median_fraction=gm,same_rayleigh_mode_verified=False))
        if links:
            counts={cid:set() for cid in prepared};indices=frame.set_index("curve_id").pair_index.to_dict()
            for l in links:
                counts[l["curve_a"]].add(indices[l["curve_b"]]);counts[l["curve_b"]].add(indices[l["curve_a"]])
            frame["compatible_other_pairs"]=frame.curve_id.map(lambda cid:len(counts[cid]))
    csv(out/"conditional_dispersion.csv",rows,CURVE_COLUMNS)
    csv(out/"curve_catalog.csv",frame if len(frame) else pd.DataFrame(columns=["curve_id","segment_id","cohort","pair_index","branch","source_phase_rad","constant_cycle","fmin_hz","fmax_hz","points","background_discriminated","compatible_other_pairs"]))
    csv(out/"cross_pair_compatible_curves.csv",links,columns=["curve_a","curve_b","segment_a","segment_b","cohort","branch","overlap_min_hz","overlap_max_hz","common_points","phase_velocity_median_fraction","phase_velocity_max_fraction","group_velocity_median_fraction","same_rayleigh_mode_verified"])
    csv(out/"segment_overlap_connections.csv",connections,columns=["segment_a","segment_b","same_cohort","overlap_points","phase_alignment_add_to_b_rad","cycle_b_equals_cycle_a_minus","maximum_aligned_phase_error_rad","measured_overlap_connection_pass","gap_interpolated"])
    write_json(out/"families_complete.json",dict(conditional_curves=len(catalog),observed_curve_points=len(rows),compatible_cross_pair_links=len(links),
        selected_unique_physical_curve=False,cycle_search_range=[p["cycle_min"],p["cycle_max"]],source_phase_hypotheses_rad=p["source_phase_rad"],
        response_corrected=meta["response_corrected"],source_phase_verified=p["source_phase_verified"],absolute_cycles_or_modes_verified=False,
        gap_predictions_mixed_into_observations=False))
