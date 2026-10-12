"""Chinese Qt desktop workflow; all expensive pipeline stages run in a child process."""
from pathlib import Path
import json,sys,codecs,datetime,os
import numpy as np
import pandas as pd
from PySide6.QtCore import Qt,QProcess,QProcessEnvironment,QTimer
from PySide6.QtGui import QPixmap,QDesktopServices,QFontDatabase,QFont
from PySide6.QtCore import QUrl
from PySide6.QtWidgets import (QApplication,QMainWindow,QWidget,QVBoxLayout,QHBoxLayout,QFormLayout,QTabWidget,
 QLineEdit,QPushButton,QLabel,QFileDialog,QMessageBox,QTableWidget,QTableWidgetItem,QHeaderView,QSplitter,
 QPlainTextEdit,QComboBox,QSpinBox,QDoubleSpinBox,QCheckBox,QProgressBar,QDialog,QDialogButtonBox,QListWidget,QAbstractItemView,QSizePolicy,QScrollArea)
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg,NavigationToolbar2QT
from .project import new_project
from .hdf import override_for,truth
from .candidates import INVERSION_DEFAULTS,manual_selection,auto_selection,import_curve_csv
from dispersion.common import write_json
APP_ROOT=Path(__file__).resolve().parent.parent


def readjson(path,default=None):return json.loads(Path(path).read_text(encoding='utf-8-sig')) if Path(path).exists() else default

def button(text,callback):
 w=QPushButton(text);w.clicked.connect(callback);return w

def spin(value,lo,hi,decimals=3,step=None):
 w=QDoubleSpinBox();w.setDecimals(decimals);w.setRange(lo,hi);w.setSingleStep(step if step is not None else (1. if decimals==0 else .01));w.setValue(value);return w

def fill_table(table,frame,columns=None):
 columns=columns or list(frame.columns);table.setRowCount(len(frame));table.setColumnCount(len(columns));table.setHorizontalHeaderLabels(columns)
 for i,row in enumerate(frame[columns].itertuples(index=False,name=None)):
  for j,value in enumerate(row):
   item=QTableWidgetItem(str(value));item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable);table.setItem(i,j,item)
 table.resizeColumnsToContents();table.horizontalHeader().setStretchLastSection(True)


class JsonDialog(QDialog):
 def __init__(self,title,data,parent=None):
  super().__init__(parent);self.setWindowTitle(title);self.resize(850,680);layout=QVBoxLayout(self)
  self.editor=QPlainTextEdit(json.dumps(data,ensure_ascii=False,indent=2));layout.addWidget(self.editor)
  box=QDialogButtonBox(QDialogButtonBox.StandardButton.Save|QDialogButtonBox.StandardButton.Cancel);box.accepted.connect(self.check);box.rejected.connect(self.reject);layout.addWidget(box)
 def check(self):
  try:self.value=json.loads(self.editor.toPlainText());self.accept()
  except Exception as exc:QMessageBox.warning(self,'JSON格式错误',str(exc))


class MappingDialog(QDialog):
 def __init__(self,item,profile,parent=None):
  super().__init__(parent);self.setWindowTitle('波形字段映射');self.resize(820,780);self.fields={};self.item=item
  m=dict(item['metadata']);m.update(override_for(item,profile));m.update(item.get('resolved',{}));self.old=override_for(item,profile)
  layout=QVBoxLayout(self);layout.addWidget(QLabel(Path(item['file']).name+'  '+item['dataset']+'  '+str(item['shape'])))
  form=QFormLayout();layout.addLayout(form)
  names=[('station','原始台站代码'),('station_label','项目台站标签（跨台网区分）'),('seed_id','NSLC（NET.STA.LOC.CHA）'),('channel','垂直通道名'),('sampling_rate_hz','原始采样率 Hz'),('start_utc','记录起始UTC'),
    ('units','单位：counts / m/s / m/s^2 / m'),('latitude','纬度（或下方本地坐标）'),('longitude','经度'),('x_m','本地x / m'),('y_m','本地y / m'),
    ('sample_axis','二维采样轴 0/1；留空自动'),('component_index','二维垂直分量索引，从0起'),('valid_dataset','有效掩码数据集（可空）'),
    ('response_xml','StationXML路径（本版内嵌响应自动读取）'),('clock_offset_s','记录UTC−真实UTC / s'),('polarity','极性 +1/-1')]
  for key,label in names:
   value=m.get(key,'');field=QLineEdit('' if value is None else str(value));self.fields[key]=field;form.addRow(label,field)
  self.timefmt=QComboBox();self.timefmt.addItems(['iso_utc','unix_s','unix_ms','unix_us','unix_ns']);self.timefmt.setCurrentText(m.get('time_format','iso_utc'));form.addRow('时间编码（数值必须明确单位）',self.timefmt)
  self.vertical=QCheckBox('确认非Z结尾通道确为垂直分量');self.vertical.setChecked(truth(m.get('vertical_confirmed',False)));form.addRow(self.vertical)
  self.polarity=QCheckBox('现场极性已独立核实');self.polarity.setChecked(truth(m.get('polarity_verified',False)));form.addRow(self.polarity)
  layout.addWidget(QLabel('缺少UTC或采样率时不能计算。经纬度与本地坐标选一套；已校正物理单位代表你的输入声明。'))
  box=QDialogButtonBox(QDialogButtonBox.StandardButton.Save|QDialogButtonBox.StandardButton.Cancel);box.accepted.connect(self.save);box.rejected.connect(self.reject);layout.addWidget(box)
 def save(self):
  try:
   value=dict(self.old)
   for key,w in self.fields.items():
    text=w.text().strip()
    if text:
     value[key]=float(text) if key in ('sampling_rate_hz','latitude','longitude','x_m','y_m','clock_offset_s') else int(text) if key in ('sample_axis','component_index','polarity') else text
    else:value[key]=None
   value.update(time_format=self.timefmt.currentText(),vertical_confirmed=self.vertical.isChecked(),polarity_verified=self.polarity.isChecked())
   if 'units' in value:value['precorrected']=value['units']!='counts'
   from .hdf import resolve
   resolve(self.item,value);self.value=value;self.accept()
  except Exception as exc:QMessageBox.warning(self,'字段尚未完整',str(exc))


class MainWindow(QMainWindow):
 def __init__(self,project=None):
  super().__init__();from . import __version__;self.setWindowTitle(f'Noise Dispersion Studio {__version__} · 单台 H/V、频散与浅层结构');self.resize(1400,940)
  self.config_snapshot={};self.config_dirty_fields=set()
  self.root=None;self.inventory=None;self.curves=None;self.catalog=pd.DataFrame();self.process=None;self.current_job='';self.buffer='';self.pending='';self.worker_log=None;self.cancelled=False
  widget=QWidget();self.setCentralWidget(widget);layout=QVBoxLayout(widget)
  heading=QLabel('噪声频散与条件反演');heading.setStyleSheet('font-size:24px;font-weight:600;color:#174963;');layout.addWidget(heading)
  bar=QHBoxLayout();self.newbutton=button('新建项目',self.new);self.openbutton=button('打开项目',self.open);bar.addWidget(self.newbutton);bar.addWidget(self.openbutton);bar.addWidget(button('下载数据 / 阵列预览',self.download_gui));self.pathlabel=QLabel('选择项目后导入H5文件夹');bar.addWidget(self.pathlabel,1);bar.addWidget(button('打开输出目录',self.open_output));layout.addLayout(bar)
  self.tabs=QTabWidget();layout.addWidget(self.tabs,1);self.input_tab();self.config_tab();self.candidate_tab();self.inversion_tab();self.figures_tab();self.phase_tab();self.hvsr_tab();self.geothermal_tab();self.tabs.currentChanged.connect(self.refresh_image_size)
  bottom=QHBoxLayout();self.autoinvert=QCheckBox('分析完成后自动组合并反演（条件假设）');self.autoinvert.setChecked(True);bottom.addWidget(self.autoinvert)
  self.runbutton=button('一键分析',self.primary_analysis);bottom.addWidget(self.runbutton)
  self.stopbutton=button('停止',self.stop);self.stopbutton.setEnabled(False);bottom.addWidget(self.stopbutton);self.status=QLabel('就绪');bottom.addWidget(self.status,1);self.progress=QProgressBar();self.progress.setMaximumWidth(220);self.progress.setValue(0);bottom.addWidget(self.progress);layout.addLayout(bottom)
  self.log=QPlainTextEdit();self.log.setReadOnly(True);self.log.setMaximumBlockCount(10000);self.log.setMaximumHeight(160);layout.addWidget(self.log)
  self.setStyleSheet('QPushButton{padding:7px 13px;} QLineEdit,QComboBox,QSpinBox,QDoubleSpinBox{min-height:25px;} QTabWidget::pane{border:1px solid #d8e1e9;}')
  self.tabs.currentChanged.connect(self.primary_tab_changed)
  if project:self.load(project)
 def primary_tab_changed(self,*_):
  hv=self.tabs.currentWidget() is self.hvsr;phase=self.tabs.currentWidget() is self.phasepage;geo=self.tabs.currentWidget() is self.geothermal
  self.runbutton.setText('一键 H/V 分析' if hv else '一键相位筛选与叠加' if phase else '生成地热成果' if geo else '一键分析');self.autoinvert.setEnabled(not (hv or phase or geo));self.autoinvert.setVisible(not (phase or geo))
 def primary_analysis(self):
  hv=self.tabs.currentWidget() is self.hvsr;phase=self.tabs.currentWidget() is self.phasepage;geo=self.tabs.currentWidget() is self.geothermal
  self.start('hvsr' if hv else 'phase-stack' if phase else 'geothermal' if geo else 'all',auto=False if hv or phase or geo else self.autoinvert.isChecked())
 def require_root(self):
  if self.root is None:raise ValueError('请先新建或打开项目')
 def hvsr_tab(self):
  from .hvsr_gui import HvsrPanel
  self.hvsr=HvsrPanel(self);self.tabs.addTab(self.hvsr,'7 单台 H/V 与浅层')
 def geothermal_tab(self):
  from .geothermal_gui import GeothermalPanel
  self.geothermal=GeothermalPanel(self);self.tabs.addTab(self.geothermal,'8 地热与空间层析')
 def phase_tab(self):
  from .phase_gui import build_phase_page
  build_phase_page(self,spin,button)
 def phase_settings(self):
  from .phase_gui import collect_direction_settings
  custom=self.phasetimemode.currentData()=='custom'
  saved=readjson(self.root/'phase_settings.json',{}) if self.root else {}
  extras={k:saved[k] for k in ('monotonic_only','prior') if k in saved}
  reconstruction=dict(saved.get('reconstruction',{}),enabled=self.phasereconstruct.isChecked(),lower_frequency_hz=self.phaselow.value(),step_hz=self.phaselowstep.value())
  return dict(curve_ids=self.phaseids.text().strip() or None,time_mode='custom' if custom else 'full',
   start=self.phasestart.text().strip() if custom else None,end=self.phaseend.text().strip() if custom else None,
   repeat_cycles=self.phaserepeat.value(),slide_cycles=self.phaseslide.value(),tolerance_deg=self.phasetol.value(),
   range_source=self.phaserange.currentData(),taper=self.phasetaper.currentData(),stabilization_s=self.phasestabilization.value(),
   frequency_method=self.phasefrequency.currentData(),amplitude_rejection=self.phaseamp.isChecked(),split_branches=True,screening_mode=self.phasescreening.currentData(),direction=collect_direction_settings(self),reconstruction=reconstruction,**extras)
 def download_gui(self):
  import subprocess
  subprocess.Popen([sys.executable,str(APP_ROOT/'h5data_gui.py')],cwd=str(APP_ROOT))
 def guard(self,fn):
  try:return fn()
  except Exception as exc:QMessageBox.warning(self,'操作未完成',str(exc));return None
 def new(self):
  folder=QFileDialog.getExistingDirectory(self,'选择新项目目录')
  if folder:self.guard(lambda:self.load(new_project(folder)))
 def open(self):
  folder=QFileDialog.getExistingDirectory(self,'打开包含project.json的项目目录')
  if folder:self.guard(lambda:self.load(folder))
 def load(self,folder):
  root=Path(folder).resolve()
  if self.process and self.process.state()!=QProcess.ProcessState.NotRunning:raise ValueError('计算时不能切换项目，请先停止任务')
  if not (root/'project.json').exists():raise ValueError('目录中没有project.json；请选择项目或先新建')
  self.root=root;self.pathlabel.setText(str(root));spec=readjson(root/'project.json');input_folder=spec.get('input_folder','');self.folder.setText(str((root/input_folder).resolve()) if input_folder else '');self.excluded.setText(','.join(spec.get('excluded_stations',[])));self.autogeo.setChecked(spec.get('automatic_geometry',True))
  phase=readjson(root/'phase_settings.json',{})
  from .phase_gui import load_direction_settings
  load_direction_settings(self,phase.get('direction',{'mode':'filter'}))
  ids=phase.get('curve_ids') or '';self.phaseids.setText(','.join(ids) if isinstance(ids,list) else ids)
  self.phasestart.setText(phase.get('start') or '');self.phaseend.setText(phase.get('end') or '')
  self.phasetimemode.setCurrentIndex(max(0,self.phasetimemode.findData(phase.get('time_mode','full'))));self.phaseextras.setChecked(phase.get('time_mode')=='custom')
  self.phaserepeat.setValue(phase.get('repeat_cycles',20));self.phaseslide.setValue(phase.get('slide_cycles',1));self.phasetol.setValue(phase.get('tolerance_deg',15))
  self.phaserange.setCurrentIndex(max(0,self.phaserange.findData(phase.get('range_source','fixed'))));self.phasetaper.setCurrentIndex(max(0,self.phasetaper.findData(phase.get('taper','none'))))
  self.phasestabilization.setValue(phase.get('stabilization_s',600));self.phaseamp.setChecked(phase.get('amplitude_rejection',True))
  for field,key,default in ((self.phasefrequency,'frequency_method','matlab_bin'),(self.phasescreening,'screening_mode','curve')):field.setCurrentIndex(max(0,field.findData(phase.get(key,default))))
  reconstruction=phase.get('reconstruction',{});blocked=self.phasereconstruct.blockSignals(True)
  self.phasereconstruct.setChecked(reconstruction.get('enabled',False));self.phasereconstruct.blockSignals(blocked)
  self.phaselow.setValue(reconstruction.get('lower_frequency_hz',.2));self.phaselowstep.setValue(reconstruction.get('step_hz',.05))
  from dispersion.common import merge,DEFAULTS
  self.sync_config(merge(DEFAULTS,readjson(root/'config.json')))
  settings=dict(INVERSION_DEFAULTS)|readjson(root/'inversion_settings.json',{});self.layerfield.setText(','.join(map(str,settings['layers'])))
  for key,w in self.invfields.items():w.setValue(settings[key])
  self.independent.setChecked(settings['require_independent_solver']);self.invjson.setPlainText(json.dumps(readjson(root/'selection.json',{'settings':settings,'cases':[]}),ensure_ascii=False,indent=2))
  self.hvsr.load(root,self.folder.text());self.geothermal.load(root);self.inventory=None;self.hdftable.setRowCount(0);self.refresh()
  from .phase_gui import update_phase_preview
  update_phase_preview(self)
 def input_tab(self):
  page=QWidget();v=QVBoxLayout(page);row=QHBoxLayout();self.folder=QLineEdit();row.addWidget(QLabel('HDF5文件夹'));row.addWidget(self.folder,1);row.addWidget(button('选择文件夹',self.choose_folder));row.addWidget(button('扫描与识别',self.scan));v.addLayout(row)
  row=QHBoxLayout();row.addWidget(button('编辑所选映射',self.edit_mapping));row.addWidget(button('批量映射模板 JSON',self.edit_profile));row.addWidget(button('导入映射模板',self.load_profile));row.addWidget(button('导出映射模板',self.export_profile));row.addWidget(button('统一设置StationXML',self.response));row.addStretch();v.addLayout(row)
  self.hdftable=QTableWidget();self.hdftable.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows);self.hdftable.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection);self.hdftable.cellDoubleClicked.connect(lambda *_:self.edit_mapping());v.addWidget(self.hdftable,1)
  self.hdfstatus=QLabel('支持 .h5 / .hdf5 / .hdf，递归扫描。识别失败的字段可映射；文件以只读方式打开。');self.hdfstatus.setWordWrap(True);v.addWidget(self.hdfstatus)
  row=QHBoxLayout();self.excluded=QLineEdit();row.addWidget(QLabel('排除台站（逗号分隔）'));row.addWidget(self.excluded,1);row.addWidget(button('仅导入与核查',lambda:self.start('import')));row.addWidget(button('查看输入质量表',self.input_quality));row.addWidget(button('从已有候选CSV开始',self.import_csv));v.addLayout(row);self.tabs.addTab(page,'1 输入与映射')
 def choose_folder(self):
  path=QFileDialog.getExistingDirectory(self,'选择H5文件夹')
  if path:self.folder.setText(path)
 def input_quality(self):
  def act():
   self.require_root();path=self.root/'input_station_qc.csv'
   if not path.exists():raise ValueError('请先导入与核查波形')
   dialog=QDialog(self);dialog.setWindowTitle('输入台站质量：异常提示需结合仪器记录核实');dialog.resize(1100,560);layout=QVBoxLayout(dialog)
   layout.addWidget(QLabel('大直流/近常数提示可能来自ADC轨或偏置；不凭此自动排站。可在输入页指定排除台站后重跑。'))
   table=QTableWidget();fill_table(table,pd.read_csv(path));layout.addWidget(table);dialog.exec()
  self.guard(act)
 def scan(self):self.start('inspect')
 def populate_inventory(self,data):
  self.inventory=data;rows=[]
  for item in data['items']:
   m=item.get('resolved',item['metadata']);rows.append({'使用':item.get('enabled',True),'文件':Path(item['file']).name,'数据集':item['dataset'],'形状':str(item['shape']),
    '台站':m.get('station',''),'通道':m.get('channel',''),'采样率':m.get('sampling_rate_hz',''),'UTC':m.get('start_utc',''),'单位':m.get('units',''),'状态':'; '.join(item['issues']) or '已识别'})
  fill_table(self.hdftable,pd.DataFrame(rows,columns=['使用','文件','数据集','形状','台站','通道','采样率','UTC','单位','状态']))
  for i,item in enumerate(data['items']):
   cell=self.hdftable.item(i,0);cell.setFlags(cell.flags()|Qt.ItemFlag.ItemIsUserCheckable);cell.setCheckState(Qt.CheckState.Checked if item.get('enabled',True) else Qt.CheckState.Unchecked)
  self.hdfstatus.setText(f"{data['files']}个HDF5文件，{len(rows)}个波形数据集；"+('必需字段完整' if data['ready'] else '存在待补字段或文件错误')+'。'+str(data.get('file_errors',[])))
 def save_dataset_checks(self):
  if not self.inventory:return
  profile=readjson(self.root/'mapping.json');over=profile.setdefault('overrides',{})
  for i,item in enumerate(self.inventory['items']):
   key=item['file']+'::'+item['dataset'];over.setdefault(key,{})['enabled']=self.hdftable.item(i,0).checkState()==Qt.CheckState.Checked
  write_json(self.root/'mapping.json',profile)
 def edit_mapping(self):
  def act():
   self.require_root();row=self.hdftable.currentRow()
   if row<0 or not self.inventory:raise ValueError('请先扫描并选择一个数据集')
   item=self.inventory['items'][row];profile=readjson(self.root/'mapping.json');dialog=MappingDialog(item,profile,self)
   if dialog.exec():
    profile.setdefault('overrides',{})[item['file']+'::'+item['dataset']]=dialog.value;write_json(self.root/'mapping.json',profile);self.scan()
  self.guard(act)
 def edit_profile(self):
  def act():
   self.require_root();dialog=JsonDialog('可复用HDF5映射：defaults / rules / overrides',readjson(self.root/'mapping.json'),self)
   if dialog.exec():write_json(self.root/'mapping.json',dialog.value);self.scan()
  self.guard(act)
 def load_profile(self):
  def act():
   self.require_root();path,_=QFileDialog.getOpenFileName(self,'导入映射模板','','JSON (*.json)')
   if path:write_json(self.root/'mapping.json',readjson(path));self.scan()
  self.guard(act)
 def export_profile(self):
  def act():
   self.require_root();path,_=QFileDialog.getSaveFileName(self,'导出映射模板','hdf_mapping.json','JSON (*.json)')
   if path:write_json(Path(path),readjson(self.root/'mapping.json'))
  self.guard(act)
 def response(self):
  def act():
   self.require_root();path,_=QFileDialog.getOpenFileName(self,'选择完整StationXML','','StationXML (*.xml);;All (*)')
   if path:
    profile=readjson(self.root/'mapping.json');profile['defaults']['response_xml']=path;write_json(self.root/'mapping.json',profile);self.scan()
  self.guard(act)
 def import_csv(self):
  def act():
   self.require_root();path,_=QFileDialog.getOpenFileName(self,'导入程序生成的候选频散表','','CSV (*.csv)')
   if path:import_curve_csv(self.root,path);self.refresh();self.tabs.setCurrentIndex(2)
  self.guard(act)
 def config_tab(self):
  page=QWidget();layout=QHBoxLayout(page);left=QWidget();form=QFormLayout(left);scroll=QScrollArea();scroll.setWidgetResizable(True);scroll.setWidget(left);scroll.setMinimumWidth(500);layout.addWidget(scroll);self.configfields={}
  entries=[('scan.min_hz','频率下限 Hz',.000001,100,.005),('scan.max_hz','频率上限 Hz',.000001,100,.005),
   ('scan.step_hz','扫描步长 Hz',.000001,1,.001),('scan.min_band_hz','最小连续带宽 Hz',.000001,100,.005),
   ('sample_rate_hz','分析采样率 Hz',1,1000,1),('discard_first_s','共同去前 / s',0,100000,1),('window_s','分析窗 / s',1,100000,1),
   ('min_clean_segment_s','连续干净片段 / s',.001,100000,1),('min_clean_fraction','共同干净比例',.001,1,.01),
   ('validation.snr_min','开发/验证SNR门槛',.001,1000,.5),('spikes.threshold_mad','尖峰MAD阈值',.001,100,.5),
   ('spikes.max_repair_s','最大尖峰修复时长 / s',.000001,1,.001)]
  from dispersion.common import DEFAULTS
  for path,label,lo,hi,step in entries:
   val=DEFAULTS
   for key in path.split('.'):val=val[key]
   w=spin(val,lo,hi,6,step);self.configfields[path]=w;form.addRow(label,w)
   w.valueChanged.connect(lambda _,key=path:self.config_dirty_fields.add(key))
  w=QSpinBox();w.setRange(3,1000000);w.setSingleStep(1);w.setValue(DEFAULTS['scan']['min_points']);self.configfields['scan.min_points']=w
  w.valueChanged.connect(lambda _:self.config_dirty_fields.add('scan.min_points'));form.addRow('最少连续频点数',w)
  self.confignote=QLabel('带宽=末频率−首频率；须同时满足最少连续频点数。修改后点击“仅扫描”重新筛选。');self.confignote.setWordWrap(True);self.confignote.setMaximumWidth(420);form.addRow(self.confignote)
  self.autogeo=QCheckBox('导入后自动核对几何并调整到时/噪声区');self.autogeo.setChecked(True);form.addRow(self.autogeo)
  form.addRow(button('保存常用参数',self.save_common));form.addRow(button('仅预处理',lambda:self.start('prepare')));form.addRow(button('仅扫描',lambda:self.start('scan')));form.addRow(button('错时对照',lambda:self.start('controls')));form.addRow(button('候选枚举',lambda:self.start('families')))
  right=QVBoxLayout();layout.addLayout(right,1);right.addWidget(QLabel('高级配置 JSON：带通、响应、UTC分界、验证块、频散周期等'))
  self.configeditor=QPlainTextEdit();right.addWidget(self.configeditor,1);right.addWidget(button('保存高级JSON并同步常用参数',self.save_advanced));self.tabs.addTab(page,'2 分析设置')
 def sync_config(self,c):
  from copy import deepcopy
  from decimal import Decimal
  self.config_snapshot=deepcopy(c)
  for path,w in self.configfields.items():
   value=c
   for key in path.split('.'):value=value[key]
   blocked=w.blockSignals(True)
   try:
    if isinstance(w,QDoubleSpinBox):w.setDecimals(max(w.decimals(),min(12,max(0,-Decimal(str(value)).as_tuple().exponent))))
    w.setRange(min(w.minimum(),value),max(w.maximum(),value));w.setValue(value)
   finally:w.blockSignals(blocked)
  self.config_dirty_fields.clear();self.configeditor.setPlainText(json.dumps(c,ensure_ascii=False,indent=2))
  scan=c['scan'];self.confignote.setText(f"带宽=末频率−首频率；当前门限 {scan['min_band_hz']:g} Hz、{scan['min_points']}点。修改后点击“仅扫描”重新筛选。")
 def persist_common(self):
  self.require_root()
  for w in self.configfields.values():w.interpretText()
  from dispersion.common import merge,DEFAULTS,validate_config
  from .configuration import apply_editor_changes
  current=merge(DEFAULTS,readjson(self.root/'config.json'))
  edited=merge(DEFAULTS,json.loads(self.configeditor.toPlainText()))
  c=apply_editor_changes(self.config_snapshot,current,edited)
  for path in self.config_dirty_fields:
   obj=c;keys=path.split('.')
   for key in keys[:-1]:obj=obj[key]
   obj[keys[-1]]=self.configfields[path].value()
  validate_config(c);write_json(self.root/'config.json',c);self.sync_config(c)
 def settings_saved(self):
  c=self.config_snapshot['scan'];message=f"参数已保存：最小连续带宽 {c['min_band_hz']:g} Hz，最少 {c['min_points']} 点；旧结果保留，重新扫描后应用新门限。"
  self.status.setText(message);self.log.appendPlainText(message)
  self.confignote.setText(f"带宽=末频率−首频率；当前门限 {c['min_band_hz']:g} Hz、{c['min_points']}点。旧结果保留，修改后点击“仅扫描”重新筛选。")
 def save_common(self):
  def act():self.persist_common();self.settings_saved()
  self.guard(act)
 def save_advanced(self):
  def act():
   self.require_root();c=json.loads(self.configeditor.toPlainText());from dispersion.common import merge,DEFAULTS,validate_config
   from .configuration import apply_editor_changes
   full=apply_editor_changes(self.config_snapshot,merge(DEFAULTS,readjson(self.root/'config.json')),merge(DEFAULTS,c))
   validate_config(full);write_json(self.root/'config.json',full)
   self.sync_config(full);self.settings_saved()
  self.guard(act)
 def candidate_tab(self):
  page=QWidget();v=QVBoxLayout(page);row=QHBoxLayout();self.search=QLineEdit();self.search.setPlaceholderText('搜索curve_id / 台对 / 时间组');self.search.textChanged.connect(self.filter_candidates);row.addWidget(self.search,1)
  self.displaylimit=QSpinBox();self.displaylimit.setRange(20,20000);self.displaylimit.setValue(300);self.displaylimit.valueChanged.connect(self.filter_candidates);row.addWidget(QLabel('最多显示'));row.addWidget(self.displaylimit);row.addWidget(button('刷新结果',self.refresh));v.addLayout(row)
  split=QSplitter();self.candidatetable=QTableWidget();self.candidatetable.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows);self.candidatetable.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection);self.candidatetable.itemSelectionChanged.connect(self.plot_candidates);split.addWidget(self.candidatetable)
  self.figure=Figure(figsize=(6,5));self.canvas=FigureCanvasQTAgg(self.figure);right=QWidget();rl=QVBoxLayout(right);rl.addWidget(NavigationToolbar2QT(self.canvas,self));rl.addWidget(self.canvas,1);split.addWidget(right);split.setSizes([750,550]);v.addWidget(split,1)
  self.candidateinfo=QLabel('候选为固定整数周期和源相位下的条件相速度；点击或多选绘图。');self.candidateinfo.setWordWrap(True);v.addWidget(self.candidateinfo)
  row=QHBoxLayout();self.manualmode=QSpinBox();self.manualmode.setRange(0,5);row.addWidget(QLabel('假设Rayleigh模态（0=基阶）'));row.addWidget(self.manualmode);row.addWidget(button('所选各自建组',lambda:self.add_cases(False)));row.addWidget(button('所选联合建组',lambda:self.add_cases(True)));row.addStretch();v.addLayout(row);self.tabs.addTab(page,'3 候选曲线')
 def selected_ids(self):return [self.candidatetable.item(i.row(),0).text() for i in self.candidatetable.selectionModel().selectedRows()]
 def filter_candidates(self,*_):
  if self.catalog.empty:fill_table(self.candidatetable,pd.DataFrame(columns=['curve_id']));return
  q=self.catalog;term=self.search.text().lower().strip()
  if term:q=q[q.apply(lambda row:term in ' '.join(map(str,row.values)).lower(),axis=1)]
  cols=[k for k in ['curve_id','pair','cohort','branch','constant_cycle','fmin_hz','fmax_hz','minimum_phase_velocity_m_s','maximum_phase_velocity_m_s','compatible_other_pairs','review_priority'] if k in q]
  self.candidatetable.blockSignals(True);fill_table(self.candidatetable,q.head(self.displaylimit.value()),cols);self.candidatetable.blockSignals(False)
  self.candidateinfo.setText(f'共有{len(self.catalog)}条条件曲线；当前匹配{len(q)}条，显示前{min(len(q),self.displaylimit.value())}条。排序不是概率，完整结果仍保留。')
 def plot_candidates(self):
  ids=self.selected_ids()
  if not ids or self.curves is None:return
  self.figure.clear();ax=self.figure.add_subplot(111)
  for cid in ids[:30]:
   q=self.curves[self.curves.curve_id==cid].sort_values('frequency_hz');ax.plot(q.frequency_hz,q.phase_velocity_m_s,label=cid)
   if len(ids)==1:ax.plot(q.frequency_hz,q.group_velocity_m_s,'k:',label='Measured group speed')
  ax.set(xlabel='Frequency (Hz)',ylabel='Velocity (m/s)',title='Conditional phase-speed candidates');ax.grid(alpha=.2)
  if len(ids)<=8:ax.legend(fontsize=7)
  self.figure.tight_layout();self.canvas.draw_idle()
 def add_cases(self,joint):
  def act():
   self.require_root();manual_selection(self.root,self.selected_ids(),self.manualmode.value(),joint);self.invjson.setPlainText(json.dumps(readjson(self.root/'selection.json'),ensure_ascii=False,indent=2));self.tabs.setCurrentIndex(3)
  self.guard(act)
 def inversion_tab(self):
  page=QWidget();v=QVBoxLayout(page);row=QHBoxLayout();self.maxcases=QSpinBox();self.maxcases.setRange(1,1000);self.maxcases.setValue(12);self.automodes=QLineEdit('0');self.automodes.setMaximumWidth(90);row.addWidget(QLabel('自动组数上限'));row.addWidget(self.maxcases);row.addWidget(QLabel('自动假设模态，如0,1'));row.addWidget(self.automodes);row.addWidget(button('自动组合',self.auto_cases));row.addWidget(button('运行选定组反演',lambda:self.start('invert')));row.addStretch();v.addLayout(row)
  split=QSplitter();left=QWidget();form=QFormLayout(left);self.layerfield=QLineEdit('1,2');form.addRow('有限层数（另加半空间）',self.layerfield);self.invfields={}
  for key,label,lo,hi,dec in [('starts','每组/层数起点数',1,1000,0),('max_evaluations','每起点最大评价次数',10,10000,0),('vp_vs_ratio','固定Vp/Vs',1.16,5,3),('density_g_cm3','固定密度 g/cm³',.5,6,3),('phase_rms_max_fraction','相速度RMS上限比例',.001,.5,3),('phase_point_max_fraction','相速度最大误差比例',.001,1,3),('group_median_max_fraction','群速度中位误差比例',.001,1,3)]:
   w=spin(INVERSION_DEFAULTS[key],lo,hi,dec);self.invfields[key]=w;form.addRow(label,w)
  self.independent=QCheckBox('要求pysurf96独立求解器通过');form.addRow(self.independent);form.addRow(button('保存反演设置与组JSON',self.save_selection));note=QLabel('每组分别拟合多种结构。基阶/高阶为假设；拟合良好不能独立证明周期或模态。参数上下界、群速度权重可在组JSON中编辑。');note.setWordWrap(True);form.addRow(note);split.addWidget(left)
  self.invjson=QPlainTextEdit();split.addWidget(self.invjson);split.setSizes([450,800]);v.addWidget(split,1)
  self.modeltable=QTableWidget();v.addWidget(self.modeltable,1);self.tabs.addTab(page,'4 条件反演')
 def inv_settings(self,base=None):
  data=dict(INVERSION_DEFAULTS)|(base if base is not None else readjson(self.root/'inversion_settings.json',{}) if self.root else {})
  for key,w in self.invfields.items():data[key]=int(w.value()) if key in ('starts','max_evaluations') else w.value()
  data['layers']=[int(x.strip()) for x in self.layerfield.text().split(',') if x.strip()];data['require_independent_solver']=self.independent.isChecked();return data
 def persist_selection(self):
  self.require_root();s=json.loads(self.invjson.toPlainText());s['settings']=self.inv_settings(s.get('settings',{}));write_json(self.root/'selection.json',s);write_json(self.root/'inversion_settings.json',s['settings']);self.invjson.setPlainText(json.dumps(s,ensure_ascii=False,indent=2))
 def save_selection(self):self.guard(self.persist_selection)
 def auto_cases(self):
  def act():
   self.require_root();write_json(self.root/'inversion_settings.json',self.inv_settings());selection=auto_selection(self.root,max_cases=self.maxcases.value(),modes=tuple(int(x) for x in self.automodes.text().split(',')));self.invjson.setPlainText(json.dumps(selection,ensure_ascii=False,indent=2))
  self.guard(act)
 def figures_tab(self):
  page=QWidget();v=QVBoxLayout(page);row=QHBoxLayout();row.addWidget(button('生成/刷新报告',lambda:self.start('report')));row.addWidget(button('打开离线报告',self.open_report));row.addWidget(button('计算自洽检查',lambda:self.start('verify')));row.addStretch();v.addLayout(row)
  split=QSplitter();self.figurelist=QListWidget();self.figurelist.currentTextChanged.connect(self.show_figure);split.addWidget(self.figurelist);self.image=QLabel('完成分析后显示图件');self.image.setAlignment(Qt.AlignmentFlag.AlignCenter);self.image.setSizePolicy(QSizePolicy.Policy.Ignored,QSizePolicy.Policy.Ignored);self.image.setMinimumSize(100,100);split.addWidget(self.image);split.setSizes([300,1000]);v.addWidget(split,1);self.tabs.addTab(page,'5 图件与报告')
 def refresh_image_size(self,*_):
  if hasattr(self,'image') and self.tabs.currentIndex()==4:
   QTimer.singleShot(0,lambda:self.show_figure(self.figurelist.currentItem().text() if self.figurelist.currentItem() else ''))
 def resizeEvent(self,event):
  super().resizeEvent(event);self.refresh_image_size()
 def show_figure(self,name):
  if self.root and name:
   p=QPixmap(str(self.root/'output/figures'/name));self.image.setPixmap(p.scaled(max(400,self.image.width()-20),max(300,self.image.height()-20),Qt.AspectRatioMode.KeepAspectRatio,Qt.TransformationMode.SmoothTransformation))
 def refresh(self):
  if not self.root:return
  path=self.root/'output/candidate_review.csv';self.catalog=pd.read_csv(path) if path.exists() else pd.DataFrame()
  spec=readjson(self.root/'project.json');curvepath=self.root/(spec.get('external_curves_file','')) if spec.get('analysis_mode')=='curves_only' else self.root/'output/conditional_dispersion.csv'
  self.curves=pd.read_csv(curvepath) if curvepath.is_file() else None;self.filter_candidates()
  if (self.root/'output/families_complete.json').exists() and ((self.root/'output/phase_guided/phase_stack_complete.json').exists() or (self.root/'output/music_front/music_complete.json').exists()):
   try:
    from .candidates import read_curves,review
    self.curves,_=read_curves(self.root);self.catalog=review(self.root);self.filter_candidates()
   except (ValueError,OSError) as exc:self.log.appendPlainText('候选结果需更新：'+str(exc))
  from .inversion import inversion_current
  models=readjson(self.root/'output/conditional_models.json',{'models':[]})['models'];frame=pd.DataFrame(models)
  if models and not inversion_current(self.root):self.log.appendPlainText('反演未完成或输入/选择/程序已改变：当前模型为历史或中间结果。')
  columns=['model_id','case','finite_layers','thickness_m','vs_m_s','phase_weighted_rms_fraction','maximum_group_median_fraction','independent_solver_pass']
  fill_table(self.modeltable,frame,columns) if len(frame) else fill_table(self.modeltable,pd.DataFrame(columns=['model_id','case','vs_m_s']))
  self.figurelist.clear()
  for p in sorted((self.root/'output/figures').glob('*.png')):self.figurelist.addItem(p.name)
  if self.figurelist.count():self.figurelist.setCurrentRow(0)
  self.hvsr.reload_result();self.geothermal.reload_result()
  if (self.root/'selection.json').exists():self.invjson.setPlainText(json.dumps(readjson(self.root/'selection.json'),ensure_ascii=False,indent=2))
 def open_report(self):
  if self.root:QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.root/'output/studio_report.html')))
 def open_output(self):
  if self.root:QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.root/'output')))
 def save_project_fields(self):
  self.require_root();spec=readjson(self.root/'project.json');spec.update(input_folder=self.folder.text().strip(),excluded_stations=[x.strip() for x in self.excluded.text().split(',') if x.strip()],automatic_geometry=self.autogeo.isChecked());write_json(self.root/'project.json',spec);self.save_dataset_checks();write_json(self.root/'inversion_settings.json',self.inv_settings())
 def start(self,job,auto=False):
  def act():
   self.require_root()
   if self.process and self.process.state()!=QProcess.ProcessState.NotRunning:raise ValueError('当前任务仍在运行')
   if job=='hvsr':self.hvsr.save()
   elif job in ('geothermal','geothermal-validate','tomography-export','tomography-invert'):self.geothermal.save()
   else:
    self.save_project_fields()
    self.persist_common()
   if job=='invert':self.persist_selection()
   if job in ('phase-export','phase-stack','all'):
    settings=self.phase_settings()
    if job=='phase-stack' and settings['time_mode']=='custom' and (not settings['start'] or not settings['end']):raise ValueError('自定义截取区间需要同时填写开始与结束UTC；完整记录模式无需填写')
    write_json(self.root/'phase_settings.json',settings)
   proc=QProcess(self);self.process=proc;self.current_job=job;self.buffer='';self.pending='';self.cancelled=False;self.decoder=codecs.getincrementaldecoder('utf-8')('replace')
   args=[str(APP_ROOT/'run.py'),job]
   if job=='inspect':args+=['--folder',self.folder.text(),'--profile',str(self.root/'mapping.json')]
   else:args+=['--project',str(self.root)]
   if job in ('phase-export','phase-stack'):args+=['--phase-settings',str(self.root/'phase_settings.json')]
   if job=='phase-velocity':args+=['--folder',str(self.root/'output/phase_guided'),'--output',str(self.root/'output/phase_guided/velocity_comparison')]
   if auto:args+=['--auto-invert','--max-cases',str(self.maxcases.value()),'--modes',self.automodes.text()]
   env=QProcessEnvironment.systemEnvironment();env.insert('PYTHONUNBUFFERED','1');env.insert('PYTHONIOENCODING','utf-8');proc.setProcessEnvironment(env);proc.setProgram(sys.executable);proc.setArguments(args);proc.setWorkingDirectory(str(APP_ROOT));proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
   stamp=datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f');logs=self.root/'logs';logs.mkdir(exist_ok=True);self.worker_log=(logs/(stamp+'_'+job+'.log')).open('w',encoding='utf-8')
   proc.readyReadStandardOutput.connect(self.consume);proc.finished.connect(self.finished);proc.errorOccurred.connect(self.process_error)
   self.tabs.setEnabled(False);self.newbutton.setEnabled(False);self.openbutton.setEnabled(False);self.runbutton.setEnabled(False);self.stopbutton.setEnabled(True);self.progress.setRange(0,0);self.status.setText('运行：'+job);self.log.appendPlainText('开始 '+job);proc.start()
  self.guard(act)
 def process_error(self,error):
  self.log.appendPlainText('进程状态: '+str(error))
  if error==QProcess.ProcessError.FailedToStart:
   if self.worker_log:self.worker_log.close();self.worker_log=None
   self.tabs.setEnabled(True);self.newbutton.setEnabled(True);self.openbutton.setEnabled(True);self.runbutton.setEnabled(True);self.stopbutton.setEnabled(False);self.progress.setRange(0,100);self.progress.setValue(0);self.status.setText('进程未启动，请检查Python环境')
 def consume(self):
  text=self.decoder.decode(bytes(self.process.readAllStandardOutput()));self.buffer+=text
  if self.worker_log:self.worker_log.write(text);self.worker_log.flush()
  if self.current_job=='inspect':return
  self.pending+=text
  while '\n' in self.pending:
   line,self.pending=self.pending.split('\n',1)
   if line.startswith('@@EVENT '):
    try:d=json.loads(line[8:]);self.status.setText(d['stage']+' · '+d['status'])
    except Exception:pass
   else:self.log.appendPlainText(line)
 def finished(self,code,status):
  self.consume()
  if self.worker_log:self.worker_log.close();self.worker_log=None
  self.tabs.setEnabled(True);self.newbutton.setEnabled(True);self.openbutton.setEnabled(True);self.runbutton.setEnabled(True);self.stopbutton.setEnabled(False);self.progress.setRange(0,100);self.progress.setValue(100 if code==0 else 0)
  if self.current_job=='inspect' and code in (0,2):
   try:self.populate_inventory(json.loads(self.buffer))
   except Exception as exc:self.log.appendPlainText('扫描结果无法解析: '+str(exc))
  else:self.refresh()
  from .phase_gui import update_phase_preview
  update_phase_preview(self)
  self.status.setText('已停止；当前阶段未完成' if self.cancelled else '完成' if code==0 else '未完成，请看日志')
  if code!=0 and self.current_job!='inspect' and not self.cancelled:QMessageBox.warning(self,'任务未完成',self.buffer[-2200:])
 def stop(self):
  proc=self.process
  if proc and proc.state()!=QProcess.ProcessState.NotRunning:
   self.cancelled=True;proc.terminate();QTimer.singleShot(1500,lambda p=proc:p.kill() if p.state()!=QProcess.ProcessState.NotRunning else None)
 def closeEvent(self,event):
  if self.process and self.process.state()!=QProcess.ProcessState.NotRunning:
   self.stop();event.ignore();self.status.setText('正在停止任务，稍后可关闭')
  else:event.accept()


def launch(project=None,screenshot=None,hvsr_tab=False):
 app=QApplication.instance() or QApplication(sys.argv[:1]);app.setApplicationName('Noise Dispersion Studio')
 font=APP_ROOT/'assets/NotoSansSC.otf'
 if font.exists():
  fid=QFontDatabase.addApplicationFont(str(font));families=QFontDatabase.applicationFontFamilies(fid)
  if families:app.setFont(QFont(families[0],10))
 window=MainWindow(project);window.show()
 if hvsr_tab:window.tabs.setCurrentWidget(window.hvsr)
 if screenshot:
  def capture():window.grab().save(str(screenshot));app.quit()
  QTimer.singleShot(500,capture)
 return app.exec()
