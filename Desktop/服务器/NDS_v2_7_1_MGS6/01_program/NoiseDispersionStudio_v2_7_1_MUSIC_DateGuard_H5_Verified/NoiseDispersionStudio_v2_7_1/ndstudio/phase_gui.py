"""Frequency-adaptive phase workflow, shared by the desktop page and previews."""
import json
import numpy as np
import pandas as pd
from PySide6.QtCore import Qt,QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (QWidget,QVBoxLayout,QHBoxLayout,QGridLayout,QFormLayout,QGroupBox,
 QLabel,QLineEdit,QComboBox,QPushButton,QScrollArea,QTableWidget,QTableWidgetItem,QHeaderView,QAbstractItemView,QCheckBox,QPlainTextEdit)
from .phase_records import phase_window_geometry


def label(text,style=''):
    w=QLabel(text);w.setWordWrap(True)
    if style:w.setStyleSheet(style)
    return w


def build_phase_page(w,spin,button):
    scroll=QScrollArea();scroll.setWidgetResizable(True);scroll.setFrameShape(QScrollArea.Shape.NoFrame)
    page=QWidget();page.setStyleSheet('QGroupBox{background:#ffffff;border:1px solid #dce5ed;border-radius:6px;margin-top:12px;padding-top:8px;} QGroupBox::title{subcontrol-origin:margin;left:10px;color:#234e65;}');v=QVBoxLayout(page);v.setSpacing(10);scroll.setWidget(page);w.phasepage=scroll
    v.addWidget(label('整段波形 · 逐频筛选 · 相位叠加','font-size:20px;font-weight:600;color:#174963;'))
    steps=QGridLayout()
    for i,(title,text) in enumerate([
        ('01  全部共同记录','遍历全部记录；各台对只用共同有效窗，并顺延稳定时间。'),
        ('02  逐频 Goertzel','窗长 T=N/f，滑动 Δt=M/f；计算该频点的复数互谱。'),
        ('03  方向与45°筛选','共同UTC窗估计波矢；逐台对检查夹角、候选范围及歧义。'),
        ('04  筛相位与分组叠加','先确认方向，再筛相位；投影模式先统一各窗相位参考。'),
        ('05  连续频带与反演交接','重查最小连续带宽；保留源相位、周期与背景检验状态。')]):
        card=QGroupBox(title);box=QVBoxLayout(card);box.addWidget(label(text,'color:#536575;'));steps.addWidget(card,i//3,i%3)
    v.addLayout(steps)
    data=QGroupBox('1  处理数据');form=QFormLayout(data)
    w.phaseids=QLineEdit();w.phaseids.setPlaceholderText('留空处理全部候选；多个曲线编号用逗号分隔');form.addRow('候选曲线',w.phaseids)
    w.phasescope=label('全部已导入记录 · 自动确定UTC范围，无需填写开始/结束时间。')
    form.addRow('数据时长',w.phasescope)
    w.phasestabilization=spin(600.,0.,86400.,2,1.);w.phasestabilization.setToolTip('与附件一致：只把每对台站的共同起点顺延该秒数，默认600秒。不是每个窗口都丢弃600秒。')
    form.addRow('起点稳定时间 / s',w.phasestabilization);v.addWidget(data)
    rebuild=QGroupBox('低频实测重建与波矢投影');rf=QFormLayout(rebuild)
    w.phasereconstruct=QCheckBox('向低频外延先验，确认方向后用实测相位重建');rf.addRow(w.phasereconstruct)
    w.phaselow=spin(.2,.001,20.,3,.025);w.phaselowstep=spin(.05,.001,2.,3,.025)
    rf.addRow('低频下限 / Hz',w.phaselow);rf.addRow('低频步长 / Hz',w.phaselowstep)
    from .reconstruction import PROJECTION_NOTE
    rf.addRow(label(PROJECTION_NOTE));rf.addRow(label('须启用方向筛窗、精确Goertzel和曲线门限；先验预测不能直接反演。'));v.addWidget(rebuild)
    windows=QGroupBox('2  按频率设置时窗');grid=QGridLayout(windows)
    w.phaserepeat=spin(20.,.1,10000.,2,.1);w.phaseslide=spin(1.,.01,10000.,2,.01)
    w.phaserepeat.setToolTip('N个周期：每个频点的目标窗长 N/f 秒，取整到采样点。')
    w.phaseslide.setToolTip('M个周期：目标滑动步长 M/f 秒。遍历全部完整滑窗，与MATLAB一致。')
    form=QFormLayout();form.addRow('每窗周期数 N',w.phaserepeat);form.addRow('滑动周期数 M',w.phaseslide)
    w.phasetaper=QComboBox();w.phasetaper.addItem('矩形窗（none）','none');w.phasetaper.addItem('Hann 窗','hann');form.addRow('窗函数',w.phasetaper)
    w.phasefrequency=QComboBox();w.phasefrequency.addItem('最近 DFT 频点（与附件一致）','matlab_bin');w.phasefrequency.addItem('曲线精确频率','exact');form.addRow('Goertzel 频率',w.phasefrequency);grid.addLayout(form,0,0)
    w.phasepreview=QTableWidget();w.phasepreview.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    w.phasepreview.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection);w.phasepreview.verticalHeader().hide()
    w.phasepreview.setColumnCount(4);w.phasepreview.setHorizontalHeaderLabels(['曲线 f / Hz','变换 f / Hz','窗长 / s','滑动 / s'])
    w.phasepreview.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch);w.phasepreview.setFixedHeight(126);grid.addWidget(w.phasepreview,0,1)
    w.phasewindowhint=label('');grid.addWidget(w.phasewindowhint,1,0,1,2);grid.setColumnStretch(0,1);grid.setColumnStretch(1,1);v.addWidget(windows)
    from .direction import DEFAULTS
    direction=QGroupBox('3  MUSIC阵列方向约束');df=QGridLayout(direction)
    w.directionmode=QComboBox()
    for title,key in [('MUSIC逐频定向、投影筛窗 + 联合估计','filter'),('仅诊断方向','diagnostic'),('平面波联合估计 + 方向筛窗','joint'),('关闭方向（仅相位诊断）','off')]:w.directionmode.addItem(title,key)
    df.addWidget(QLabel('工作模式'),0,0);df.addWidget(w.directionmode,0,1,1,3)
    w.directionfields={}
    fields=[('最大台对夹角 / °','max_pair_angle_deg',45.,.1,89.,1.,1.),
        ('最大方向候选偏差 / °','max_direction_uncertainty_deg',20.,.1,89.,1.,.5),
        ('速度下限 / m/s','velocity_min_m_s',150.,1.,20000.,1.,10.),
        ('速度上限 / m/s','velocity_max_m_s',3500.,2.,50000.,1.,10.),
        ('方向步长 / °','direction_step_deg',5.,1.,30.,0.,1.),
        ('慢度网格点数','velocity_points',32,3.,256.,0.,1.),
        ('最低阵列相干度','min_coherence',.7,.01,.999,3.,.01),
        ('方向分组宽度 / °','direction_group_width_deg',15.,1.,90.,0.,1.),
        ('MUSIC频带半宽 / Hz','frequency_half_width_hz',0.,0.,5.,4.,.001),
        ('最少不重叠窗口','minimum_independent_windows',3,1.,10000.,0.,1.)]
    for i,(title,key,val,low,high,decimals,step) in enumerate(fields):
        field=spin(val,low,high,int(decimals),step);w.directionfields[key]=field
        df.addWidget(QLabel(title),1+i//2,2*(i%2));df.addWidget(field,1+i//2,2*(i%2)+1)
    w.directionconservative=QCheckBox('方向候选范围越过夹角阈值时排除（默认）');w.directionconservative.setChecked(True);df.addWidget(w.directionconservative,6,0,1,4)
    df.addWidget(label('角度按台对连线与波矢计算；反向传播单独分组。候选偏差是分辨率诊断，不是统计置信区间。MUSIC半宽与分析设置的最小连续带宽不同；半宽0表示逐频计算。前段不使用45°台对门槛，后段使用上方门槛。DPSS快照数/带宽可在高级JSON调整。'),7,0,1,4)
    advanced=QGroupBox('方向高级JSON（保存后同步上方参数）');advanced.setCheckable(True);advanced.setChecked(False);dv=QVBoxLayout(advanced)
    w.directioneditor=QPlainTextEdit();w.directioneditor.setFixedHeight(210);dv.addWidget(w.directioneditor)
    save=button('保存方向参数',lambda:w.guard(lambda:save_direction_settings(w)));dv.addWidget(save)
    advanced.toggled.connect(w.directioneditor.setVisible);advanced.toggled.connect(save.setVisible);w.directioneditor.hide();save.hide()
    df.addWidget(advanced,8,0,1,4);v.addWidget(direction)
    load_direction_settings(w,DEFAULTS)
    gate=QGroupBox('4  相位门限与分支叠加');form=QFormLayout(gate)
    w.phasescreening=QComboBox();w.phasescreening.addItem('频散曲线引导（使用逐频速度对应的相位范围）','curve');w.phasescreening.addItem('附件 MATLAB 原门限（分支初叠加相位 ±90°）','matlab');form.addRow('筛选模式',w.phasescreening)
    w.phaserange=QComboBox();w.phaserange.addItem('曲线预测相位 ± 固定容差','fixed');w.phaserange.addItem('原片段 bootstrap 相位范围 + 余量','measured')
    w.phasetol=spin(15.,0.,179.99,2,.1);w.phasetollabel=QLabel('固定相位半宽 / °')
    form.addRow('接受范围',w.phaserange);form.addRow(w.phasetollabel,w.phasetol)
    w.phaseamp=QCheckBox('剔除幅值与相邻窗口突变异常（5×MAD，与附件一致）');w.phaseamp.setChecked(True);form.addRow(w.phaseamp)
    form.addRow(label('互谱 C=G目标·conj(G参考)。旧模式按相位正负分支；新模式按已确认的传播方向筛窗，旋转到共同投影距离后叠加。'))
    form.addRow(label('包裹相位的正负与物理传播方向不同；投影模式允许同方向相位跨越±π，原始互谱与旋转角单独保存。'))
    v.addWidget(gate)
    row=QHBoxLayout();w.phaserun=button('处理全部记录并叠加',lambda:w.start('phase-stack'))
    w.phaserun.setStyleSheet('background:#176979;color:white;font-weight:600;padding:9px 18px;border-radius:4px;')
    row.addWidget(w.phaserun);row.addWidget(button('查看相位筛选报告',lambda:open_report(w,'phase_report.html')))
    row.addWidget(button('查看叠加曲线对比',lambda:open_report(w,'velocity_comparison/report.html')));row.addStretch();v.addLayout(row)
    row=QHBoxLayout();row.addWidget(button('查看方向与联合频散',lambda:open_report(w,'direction/report.html')))
    row.addWidget(button('查看逐阶段日志',lambda:open_report(w,'phase_processing.log')));row.addWidget(button('查看连续频带诊断',lambda:open_report(w,'phase_continuous_bands.csv')));row.addStretch();v.addLayout(row)
    w.phaselast=label('尚未生成相位叠加结果。','color:#536575;');v.addWidget(w.phaselast)
    extra=QGroupBox('辅助操作与可选区间');extra.setCheckable(True);extra.setChecked(False);ev=QVBoxLayout(extra)
    w.phaseadvanced=QWidget();av=QVBoxLayout(w.phaseadvanced);av.setContentsMargins(0,0,0,0)
    w.phasetimemode=QComboBox();w.phasetimemode.addItem('全部记录（默认）','full');w.phasetimemode.addItem('自定义UTC区间（仅调试/对照）','custom')
    af=QFormLayout();af.addRow('处理范围',w.phasetimemode);av.addLayout(af)
    w.phaseinterval=QWidget();tf=QFormLayout(w.phaseinterval);tf.setContentsMargins(0,0,0,0)
    w.phasestart=QLineEdit();w.phasestart.setPlaceholderText('例：2020-10-12T12:00:00Z')
    w.phaseend=QLineEdit();w.phaseend.setPlaceholderText('例：2020-10-12T12:30:00Z')
    tf.addRow('截取开始 UTC',w.phasestart);tf.addRow('截取结束 UTC（不含）',w.phaseend)
    tf.addRow(label('这里的时间仅用于截取整段记录，不是单窗长度；窗长始终按 N/f 计算。'));av.addWidget(w.phaseinterval)
    row=QHBoxLayout();row.addWidget(button('导出全部速度 / 相位表',lambda:w.start('phase-export')))
    row.addWidget(button('重新生成叠加曲线对比',lambda:w.start('phase-velocity')));av.addLayout(row)
    ev.addWidget(w.phaseadvanced);w.phaseadvanced.hide();w.phaseinterval.hide();extra.toggled.connect(w.phaseadvanced.setVisible);w.phaseextras=extra;v.addWidget(extra)
    v.addWidget(label('缺测、滤波边缘、无定义相位不参与叠加。曲线门限沿用源相位、时延分支与整数周期；相位正负分支与CCF正负时延分支不同。筛后集中不能独立验证输入曲线。','color:#667889;font-size:11px;'))
    row=QHBoxLayout();row.addWidget(button('按当前选择反演修正前后曲线',lambda:w.start('phase-compare')));row.addWidget(button('查看修正前后结构',lambda:open_report(w,'paired_comparison/report.html')));v.addLayout(row)
    v.addStretch();w.tabs.addTab(scroll,'6 曲线相位筛选')
    for control in (w.phaserepeat,w.phaseslide,w.phasestabilization):control.valueChanged.connect(lambda *_:update_phase_preview(w))
    w.phaseids.textChanged.connect(lambda *_:update_phase_preview(w))
    w.phasescreening.currentIndexChanged.connect(lambda *_:update_phase_preview(w))
    w.phasefrequency.currentIndexChanged.connect(lambda *_:update_phase_preview(w))
    w.phaserange.currentIndexChanged.connect(lambda *_:update_phase_preview(w))
    w.phasetimemode.currentIndexChanged.connect(lambda *_:update_phase_preview(w))
    update_phase_preview(w)
    def reconstruction_changed(enabled):
        if enabled:
            w.phasefrequency.setCurrentIndex(w.phasefrequency.findData('exact'))
            w.phasescreening.setCurrentIndex(w.phasescreening.findData('curve'))
            w.directionmode.setCurrentIndex(w.directionmode.findData('filter'))
    w.phasereconstruct.toggled.connect(reconstruction_changed)


def open_report(w,relative):
    def act():
        w.require_root();path=w.root/'output/phase_guided'/relative
        if not path.exists():raise ValueError('尚未生成该报告，请先完成相位筛选与叠加。')
        if relative.endswith('.html') and not (w.root/'output/phase_guided/phase_stack_complete.json').exists():raise ValueError('当前相位阶段未完成；请先查看逐阶段日志，不能把历史报告当作本次结果。')
        if relative.startswith('paired_comparison/'):
            from dispersion.common import sha
            from .inversion import selection_digest
            marker=w.root/'output/phase_guided/paired_comparison/comparison_complete.json'
            if not marker.exists() or json.loads(marker.read_text(encoding='utf-8-sig')).get('phase_completion_sha256')!=sha(w.root/'output/phase_guided/phase_stack_complete.json'):raise ValueError('配对反演未完成或相位结果已改变，请重新反演')
            selected=w.root/'selection.json'
            if selected.exists() and json.loads(marker.read_text(encoding='utf-8-sig')).get('selection_sha256')!=selection_digest(json.loads(selected.read_text(encoding='utf-8-sig'))):raise ValueError('曲线选择或反演参数已改变，请重新反演')
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))
    w.guard(act)


def direction_common(w):
    for field in w.directionfields.values():field.interpretText()
    return dict(mode=w.directionmode.currentData(),conservative_angle=w.directionconservative.isChecked(),
        **{key:int(field.value()) if key in ('velocity_points','minimum_independent_windows') else field.value() for key,field in w.directionfields.items()})


def load_direction_settings(w,c):
    from decimal import Decimal
    from .direction import settings
    c=settings(c);w.directionmode.setCurrentIndex(max(0,w.directionmode.findData(c['mode'])))
    for key,field in w.directionfields.items():
        blocked=field.blockSignals(True)
        field.setDecimals(max(field.decimals(),min(12,max(0,-Decimal(str(c[key])).as_tuple().exponent))))
        field.setRange(min(field.minimum(),c[key]),max(field.maximum(),c[key]));field.setValue(c[key]);field.blockSignals(blocked)
    w.directionconservative.setChecked(c['conservative_angle'])
    w.directioneditor.setPlainText(json.dumps(c,ensure_ascii=False,indent=2));w.direction_snapshot=direction_common(w)


def collect_direction_settings(w):
    from .direction import settings
    advanced=json.loads(w.directioneditor.toPlainText());current=direction_common(w)
    # Three-way merge: only UI fields changed since loading override the editor.
    for key,value in current.items():
        if value!=w.direction_snapshot.get(key):advanced[key]=value
    c=settings(advanced);load_direction_settings(w,c);return c


def save_direction_settings(w):
    from dispersion.common import write_json
    w.require_root();c=collect_direction_settings(w);phase=w.phase_settings();phase['direction']=c
    write_json(w.root/'phase_settings.json',phase)
    w.log.appendPlainText(f'方向参数已保存：mode={c["mode"]}，最大夹角={c["max_pair_angle_deg"]}°。重跑方向/相位阶段后生效；实际参数见phase_settings_used.json。')


def update_phase_preview(w):
    custom=w.phasetimemode.currentData()=='custom';w.phaseinterval.setVisible(custom)
    w.phaserun.setText('处理自定义区间并叠加' if custom else '处理全部记录并叠加')
    w.phasetollabel.setText('原范围外扩余量 / °' if w.phaserange.currentData()=='measured' else '固定相位半宽 / °')
    w.phasescope.setText('自定义截取区间（仅调试/对照）；仍按频率选择窗长。' if custom else f'每对台站的全部共同记录 · 共同起点顺延 {w.phasestabilization.value():g} s；无需手选开始/结束时间。')
    matlab=w.phasescreening.currentData()=='matlab';w.phaserange.setEnabled(not matlab);w.phasetol.setEnabled(not matlab)
    if matlab:w.phasetollabel.setText('分支初叠加相位 ±90°')
    frequencies=np.array([1.,2.,5.]);fs=20.;source='示例频点'
    try:
        if w.root:
            c=json.loads((w.root/'config.json').read_text(encoding='utf-8-sig'));fs=float(c.get('sample_rate_hz',20.))
            curves=pd.read_csv(w.root/'output/conditional_dispersion.csv',usecols=['curve_id','frequency_hz'])
            ids=[s.strip() for s in w.phaseids.text().split(',') if s.strip()]
            if ids:curves=curves[curves.curve_id.isin(ids)]
            frequencies=np.sort(curves.frequency_hz.unique());source='所选曲线频点'
    except (OSError,ValueError,KeyError):pass
    shown=frequencies[np.unique(np.linspace(0,len(frequencies)-1,min(3,len(frequencies)),dtype=int))] if len(frequencies) else []
    w.phasepreview.setRowCount(len(shown));geometry=[]
    for i,f in enumerate(shown):
        try:g=phase_window_geometry(float(f),fs,w.phaserepeat.value(),w.phaseslide.value());geometry.append(g);actual=f if w.phasefrequency.currentData()=='exact' else np.floor(f*g['window_samples']/fs+.5)*fs/g['window_samples'];values=[f,actual,g['window_duration_s'],g['slide_duration_s']]
        except ValueError:values=[f,'无效','无效','无效']
        for j,value in enumerate(values):
            item=QTableWidgetItem(f'{value:.6g}' if isinstance(value,(float,np.floating)) else str(value));item.setTextAlignment(Qt.AlignmentFlag.AlignCenter);w.phasepreview.setItem(i,j,item)
    warning=' 滑动步长大于窗长，会跳过部分样本；建议 M≤N。' if w.phaseslide.value()>w.phaserepeat.value() else ''
    w.phasewindowhint.setText(f'{source}预览 · T=N/f，Δt=M/f，按 {fs:g} Hz 取整到采样点；仅完整滑窗参与。'+warning)
    if w.root:
        path=w.root/'output/phase_guided/phase_stack_complete.json'
        try:
            r=json.loads(path.read_text(encoding='utf-8-sig'));scope='全部记录' if r.get('time_mode')=='full' else '自定义/旧版区间'
            w.phaselast.setText(f'上次输入范围：{scope} · {r["start_utc"]} → {r["end_utc_exclusive"]}（结束不含）\n{r["unique_pair_frequency_gates"]} 个逐频任务，{r["accepted_window_evaluations"]:,} 次窗口贡献参与最终叠加；起点稳定 {r.get("stabilization_s",0):g} s。')
            from dispersion.common import sha
            changed=r.get('settings_sha256')!=sha(w.root/'config.json')
            phasefile=w.root/'phase_settings.json'
            if r.get('phase_settings_sha256') and phasefile.exists():changed |= r['phase_settings_sha256']!=sha(phasefile)
            if changed:w.phaselast.setText(w.phaselast.text()+'\n参数已改变：以上为历史结果，请重跑方向与相位阶段。')
        except (OSError,ValueError,KeyError):w.phaselast.setText('尚未生成相位叠加结果；完成后显示实际UTC范围与窗口统计。')
