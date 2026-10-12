"""Reproduce the real BW.MGS01 download/HVSR test using only Studio modules.

From the program directory:
  python examples/hvsr_live_test.py --directory live_hvsr
  python examples/hvsr_live_test.py --directory /path/to/unpacked/test --offline
"""
from pathlib import Path
import argparse,sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ndstudio.hvsr import DEFAULTS,run_hvsr,write_json
from ndstudio.hvsr_acquisition import fetch_single_station
from ndstudio.hvsr_validation import CASES,real_validation_report
from ndstudio.project import new_project


def main():
    p=argparse.ArgumentParser();p.add_argument('--directory',type=Path,default=Path('live_hvsr'))
    p.add_argument('--offline',action='store_true',help='Use downloaded MiniSEED/H5/StationXML; no network requests')
    a=p.parse_args();root=a.directory.resolve();download=root/'download'
    if not a.offline:
        fetch_single_station('LMU','BW','MGS01','2020-11-01T00:00:00Z','2020-11-01T02:00:00Z',download,chunk_seconds=1800)
    elif not (download/'acquisition_audit.json').exists():p.error('--offline requires the unpacked real-data test download directory')
    project=new_project(root/'project')
    spec=json.loads((project/'project.json').read_text());spec.update(input_folder='../download',name='BW.MGS01 real HVSR test')
    write_json(project/'project.json',spec)
    for case,sub,file in CASES:
        config=dict(DEFAULTS,start='2020-11-01T00:00:00Z',end='2020-11-01T02:00:00Z',window_s=120)
        if case=='same_window':config.update(fmin_hz=.1)
        elif case=='low_frequency':config.update(fmin_hz=.05,window_s=300)
        write_json(project/file,dict(folder='../download',instrument_id='BW.MGS01..HH',settings=config))
        run_hvsr(download,project/'output'/sub,'BW.MGS01..HH',config)
    result=real_validation_report(root)
    print(json.dumps(result,ensure_ascii=False,indent=2))
    print('Report:',root/'validation/real_report.html')


if __name__=='__main__':main()
