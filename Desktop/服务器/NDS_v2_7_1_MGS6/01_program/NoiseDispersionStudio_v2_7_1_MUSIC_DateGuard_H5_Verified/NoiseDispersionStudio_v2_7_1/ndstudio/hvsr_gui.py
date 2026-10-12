"""Single-station Qt panel; processing uses the shared cancellable child process."""
from pathlib import Path
import json
import pandas as pd
from PySide6.QtCore import Qt,QUrl
from PySide6.QtGui import QDesktopServices,QPixmap
from PySide6.QtWidgets import QWidget,QVBoxLayout,QHBoxLayout,QFormLayout,QLineEdit,QComboBox,QCheckBox,QLabel,QPushButton,QFileDialog,QSplitter,QTableWidget,QSizePolicy
from .hvsr import DEFAULTS,discover,settings,write_json


class HvsrPanel(QWidget):
    def __init__(self,owner):
        super().__init__();self.owner=owner
        from .gui import button,spin
        layout=QVBoxLayout(self);top=QHBoxLayout();self.folder=QLineEdit()
        top.addWidget(QLabel('单台三分量 H5 文件夹'));top.addWidget(self.folder,1)
        top.addWidget(button('选择',self.choose));top.addWidget(button('列出仪器 / 检查三分量',lambda:owner.guard(self.scan)));layout.addLayout(top)
        form=QFormLayout();self.instrument=QComboBox();form.addRow('同一仪器（台网.台站.位置.通道前缀）',self.instrument)
        self.start=QLineEdit();self.end=QLineEdit();self.start.setPlaceholderText('可空，采用三分量共同UTC');self.end.setPlaceholderText('可空；填写时结束UTC不包含端点')
        row=QHBoxLayout();row.addWidget(self.start);row.addWidget(self.end);form.addRow('开始 / 结束 UTC',row)
        self.fields={}
        entries=[('fmin_hz','频率下限 Hz',.001,100),('fmax_hz','频率上限 Hz',.002,1000),('window_s','窗长 s',1,7200),('overlap','重叠比例',0,.95),
                 ('ko_bandwidth','K–O 平滑带宽 b',5,150),('vs_min_m_s','等效覆盖层 Vs 下限 m/s',30,4000),('vs_max_m_s','等效覆盖层 Vs 上限 m/s',31,5000),('bedrock_vs_m_s','半空间 Vs 先验 m/s',100,6000)]
        # Pair fields into four compact rows to leave space for plots.
        for j in range(0,len(entries),2):
            row=QHBoxLayout()
            for k,label,lo,hi in entries[j:j+2]:
                w=spin(DEFAULTS[k],lo,hi,3,.01 if k in ('fmin_hz','fmax_hz','overlap') else 1.);self.fields[k]=w;row.addWidget(QLabel(label));row.addWidget(w)
            form.addRow(row)
        self.horizontal=QComboBox();self.horizontal.addItem('RMS（旋转不变）','rms');self.horizontal.addItem('几何平均','geometric');self.horizontal.addItem('总水平幅值','total');form.addRow('水平分量合成',self.horizontal)
        self.peak=QLineEdit();self.peak.setPlaceholderText('可空；有多峰时填写要分析的峰频率 Hz');form.addRow('指定峰（可选）',self.peak)
        self.transients=QCheckBox('剔除突发噪声（按原波形；不按H/V峰筛窗）');self.transients.setChecked(True)
        self.rayleigh=QCheckBox('启用条件 Rayleigh 椭圆率试算（须确认 Rayleigh 主导假设）')
        form.addRow(self.transients);form.addRow(self.rayleigh);layout.addLayout(form)
        row=QHBoxLayout();row.addWidget(button('计算 H/V 与浅层结构约束',lambda:owner.start('hvsr')));row.addWidget(button('高级 H/V 参数',lambda:owner.guard(self.advanced)));row.addWidget(button('打开 H/V 离线报告',self.open_report));row.addWidget(button('刷新结果',self.reload_result));layout.addLayout(row)
        self.info=QLabel('只需单台完整三分量，无须台对几何或频散预处理。低频需要足够长的窗与记录；Vs与层厚不能仅由H/V唯一确定。');self.info.setWordWrap(True);layout.addWidget(self.info)
        split=QSplitter();self.table=QTableWidget();split.addWidget(self.table);self.image=QLabel('计算后显示 H/V 曲线');self.image.setAlignment(Qt.AlignmentFlag.AlignCenter);self.image.setMinimumSize(160,160);self.image.setSizePolicy(QSizePolicy.Policy.Ignored,QSizePolicy.Policy.Ignored);split.addWidget(self.image);split.setSizes([420,700]);layout.addWidget(split,1)
        self.extra=dict(DEFAULTS)

    def choose(self):
        path=QFileDialog.getExistingDirectory(self,'选择下载H5文件夹',self.folder.text())
        if path:self.folder.setText(path);self.owner.guard(self.scan)

    def scan(self,selected=None):
        inv=discover(self.folder.text());current=selected or self.instrument.currentData();self.instrument.clear()
        for r in inv['instruments']:
            self.instrument.addItem(r['instrument_id']+' · '+','.join(r['channels'])+('' if r['ready'] else ' · '+','.join(r['issues'])),r['instrument_id'])
        i=self.instrument.findData(current)
        if i>=0:self.instrument.setCurrentIndex(i)
        errors='；'.join(x['error'] for x in inv['errors'][:3])
        self.info.setText(f"发现 {len(inv['instruments'])} 台仪器；完整三分量 {sum(x['ready'] for x in inv['instruments'])} 台。"+errors)

    def config(self):
        c=dict(self.extra)
        for k,w in self.fields.items():c[k]=w.value()
        c.update(start=self.start.text().strip(),end=self.end.text().strip(),horizontal=self.horizontal.currentData(),reject_transients=self.transients.isChecked(),rayleigh_assumption=self.rayleigh.isChecked(),peak_frequency_hz=float(self.peak.text()) if self.peak.text().strip() else None)
        return settings(c)

    def save(self):
        self.owner.require_root()
        if not self.folder.text().strip():raise ValueError('请选择三分量H5文件夹')
        if not self.instrument.currentData():self.scan()
        if not self.instrument.currentData():raise ValueError('未找到仪器；检查下载是否包含N/E/Z或有方向元数据的1/2/Z')
        write_json(self.owner.root/'hvsr_settings.json',dict(folder=str(Path(self.folder.text()).resolve()),instrument_id=self.instrument.currentData(),settings=self.config()))

    def load(self,root,default_folder=''):
        path=root/'hvsr_settings.json';data=json.loads(path.read_text(encoding='utf-8-sig')) if path.exists() else {}
        self.extra=dict(DEFAULTS)|data.get('settings',{});folder=data.get('folder') or default_folder
        if folder and not Path(folder).is_absolute():folder=str((root/folder).resolve())
        self.folder.setText(folder);self.instrument.clear()
        self.start.setText(self.extra['start']);self.end.setText(self.extra['end']);self.peak.setText(str(self.extra['peak_frequency_hz'] or ''))
        for k,w in self.fields.items():w.setValue(self.extra[k])
        self.horizontal.setCurrentIndex(max(0,self.horizontal.findData(self.extra['horizontal'])));self.transients.setChecked(self.extra['reject_transients']);self.rayleigh.setChecked(self.extra['rayleigh_assumption'])
        if data.get('instrument_id'):self.instrument.addItem(data['instrument_id'],data['instrument_id'])
        self.reload_result()

    def advanced(self):
        from .gui import JsonDialog
        d=JsonDialog('H/V 与条件浅层分析参数',self.config(),self)
        if d.exec():
            self.extra=settings(d.value)
            for k,w in self.fields.items():w.setValue(self.extra[k])
            self.horizontal.setCurrentIndex(self.horizontal.findData(self.extra['horizontal']));self.start.setText(self.extra['start']);self.end.setText(self.extra['end']);self.peak.setText(str(self.extra['peak_frequency_hz'] or ''));self.transients.setChecked(self.extra['reject_transients']);self.rayleigh.setChecked(self.extra['rayleigh_assumption'])

    def reload_result(self):
        if not self.owner.root:return
        out=self.owner.root/'output/hvsr'
        if (out/'hvsr_summary.json').exists():
            from .gui import fill_table
            s=json.loads((out/'hvsr_summary.json').read_text(encoding='utf-8-sig'))
            frame=pd.DataFrame(s['peaks']);cols=['f0_hz','amplitude','reliability_pass','clarity_pass_count','selected']
            fill_table(self.table,frame,cols) if len(frame) else fill_table(self.table,pd.DataFrame(columns=cols))
            status='完成' if (out/'hvsr_complete.json').exists() else '历史结果；当前运行未完成'
            state=s['structure']['status']
            text={'conditional':'条件约束','incomplete_peak_shape':'峰形不完整，暂停层厚估计',
                  'no_reliable_peak':'峰未通过门限，暂停层厚估计','ambiguous_multiple_peaks':'多峰待确认基频'}.get(state,state)
            self.info.setText(f"{status} · {s['instrument_id']} · 有效窗{s['accepted_windows']} · 独立窗{s['independent_windows']} · 浅层约束：{text}")
            self.display_image()
        else:self.table.setRowCount(0);self.image.clear()

    def display_image(self):
        if not self.owner.root:return
        p=QPixmap(str(self.owner.root/'output/hvsr/figures/hvsr_curve.png'))
        if not p.isNull():self.image.setPixmap(p.scaled(max(200,self.image.width()-10),max(150,self.image.height()-10),Qt.AspectRatioMode.KeepAspectRatio,Qt.TransformationMode.SmoothTransformation))

    def resizeEvent(self,event):
        super().resizeEvent(event);self.display_image()

    def open_report(self):
        if self.owner.root:QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.owner.root/'output/hvsr/hvsr_report.html')))
