"""Actual Qt/QProcess smoke test; run with QT_QPA_PLATFORM=offscreen in CI."""
import argparse,json,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import ndstudio
from ndstudio.gui import MainWindow,MappingDialog,APP_ROOT,readjson
from ndstudio.project import new_project
from ndstudio.candidates import import_curve_csv
from ndstudio.examples import inversion_example
from dispersion.common import write_json
from PySide6.QtWidgets import QApplication
from PySide6.QtGui import QFontDatabase,QFont
from PySide6.QtCore import QTimer,QProcess

p=argparse.ArgumentParser();p.add_argument('--directory',type=Path,required=True);p.add_argument('--folder',type=Path,required=True);a=p.parse_args()
root=new_project(a.directory/'project');known=inversion_example(a.directory/'known')
spec=readjson(root/'project.json');spec.update(input_folder=str(a.folder.resolve()));write_json(root/'project.json',spec)
app=QApplication([]);fid=QFontDatabase.addApplicationFont(str(APP_ROOT/'assets/NotoSansSC.otf'));families=QFontDatabase.applicationFontFamilies(fid)
if families:app.setFont(QFont(families[0],10))
w=MainWindow(root);w.show();beats=[];timer=QTimer();timer.timeout.connect(lambda:beats.append(time.monotonic()));timer.start(20)
def wait_job():
    deadline=time.monotonic()+120
    while w.process and w.process.state()!=QProcess.ProcessState.NotRunning and time.monotonic()<deadline:app.processEvents();time.sleep(.01)
    app.processEvents();assert w.process.exitCode()==0,w.buffer
    assert w.tabs.isEnabled() and not w.stopbutton.isEnabled()
w.configfields['discard_first_s'].setValue(601)
w.start('inspect');assert not w.openbutton.isEnabled();wait_job();assert w.hdftable.rowCount()==4 and w.inventory['ready'] and w.openbutton.isEnabled()
assert readjson(root/'config.json')['discard_first_s']==601
images=a.directory/'screenshots';images.mkdir(exist_ok=True)
w.grab().save(str(images/'01_hdf_mapping.png'))
dialog=MappingDialog(w.inventory['items'][0],readjson(root/'mapping.json'),w);dialog.show();app.processEvents();dialog.grab().save(str(images/'05_mapping_dialog.png'));dialog.close()
import_curve_csv(root,known/'curves.csv');s=readjson(known/'selection.json');write_json(root/'selection.json',s);write_json(root/'inversion_settings.json',s['settings']);w.load(root)
w.tabs.setCurrentIndex(2);w.candidatetable.selectRow(0);app.processEvents();w.grab().save(str(images/'02_candidates.png'))
w.tabs.setCurrentIndex(3);start=len(beats);w.start('invert');wait_job();assert len(beats)-start>=5,'GUI event loop blocked during inversion'
marker=readjson(root/'output/inversion_complete.json');assert marker['compatible_models']>=1
w.start('report');wait_job();w.tabs.setCurrentIndex(3);app.processEvents();w.grab().save(str(images/'03_inversion.png'))
w.tabs.setCurrentIndex(4);app.processEvents();w.grab().save(str(images/'04_structure.png'))
result={'hdf_records_displayed':4,'implicit_parameter_save_passed':True,'qt_child_process_success':True,
        'ui_timer_ticks_during_inversion':len(beats)-start,'conditional_models':marker['compatible_models'],'actual_ui_screenshots':5}
write_json(a.directory/'gui_smoke_result.json',result);print(json.dumps(result));w.close()
