#!/usr/bin/env python3
"""Usage: python run.py all --config config.json; python run.py demo --directory demo."""
import os
for variable in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"):
    os.environ.setdefault(variable,"1")
from pathlib import Path
import argparse,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from dispersion.common import load_config,DEFAULTS,write_json

def all_steps(c):
    from dispersion.prepare import prepare
    from dispersion.scan import scan
    from dispersion.controls import controls
    from dispersion.families import families
    from dispersion.products import report,verify
    for function in (prepare,scan,controls,families,report,verify): function(c)

def main():
    ap=argparse.ArgumentParser(description="可复用的传统环境噪声频散提取")
    ap.add_argument("stage",choices=["init","all","prepare","scan","controls","families","report","verify","invert","demo","check-demo"])
    ap.add_argument("--config",type=Path);ap.add_argument("--directory",type=Path,default=Path("demo"));ap.add_argument("--selection",type=Path)
    ap.add_argument("--field",choices=["plane","noise"],default="plane",help="Synthetic recovery or independent-noise rejection")
    a=ap.parse_args()
    if a.stage=="init":
        a.directory.mkdir(parents=True,exist_ok=True);target=a.directory/"config.json"
        if target.exists() or (a.directory/"stations.csv").exists(): raise FileExistsError("Input/config already exists; choose a new directory")
        write_json(target,DEFAULTS)
        (a.directory/"stations.csv").write_text("station,file,seed_id,x_m,y_m,response_xml,precorrected,units,polarity,polarity_verified,clock_offset_s\nA,A.mseed,XX.A..HHZ,0,0,response.xml,false,counts,1,false,0\nB,B.mseed,XX.B..HHZ,1000,0,response.xml,false,counts,1,false,0\n",encoding="utf-8")
        print("Template created:",target);return
    if a.stage in ("demo","check-demo"):
        from dispersion.demo import generate,check
        if a.stage=="demo": all_steps(load_config(generate(a.directory,field=a.field)))
        print(check(a.directory).to_string(index=False));return
    if a.config is None: ap.error("--config is required")
    c=load_config(a.config)
    if a.stage=="all": all_steps(c)
    elif a.stage=="prepare":
        from dispersion.prepare import prepare;prepare(c)
    elif a.stage=="scan":
        from dispersion.scan import scan;scan(c)
    elif a.stage=="controls":
        from dispersion.controls import controls;controls(c)
    elif a.stage=="families":
        from dispersion.families import families;families(c)
    elif a.stage in ("report","verify"):
        from dispersion.products import report,verify
        (report if a.stage=="report" else verify)(c)
    elif a.stage=="invert":
        if a.selection is None: ap.error("--selection is required for explicit cycle/mode cases")
        from dispersion.inversion import invert
        invert(c,a.selection)

if __name__=="__main__": main()
