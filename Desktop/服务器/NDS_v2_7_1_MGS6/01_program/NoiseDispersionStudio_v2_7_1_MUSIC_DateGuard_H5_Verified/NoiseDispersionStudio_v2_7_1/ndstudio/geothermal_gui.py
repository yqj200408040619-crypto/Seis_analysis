"""Qt panel for geothermal/geology integration and tomography handoff."""
from __future__ import annotations
from pathlib import Path
import json
from PySide6.QtCore import QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (QWidget,QVBoxLayout,QHBoxLayout,QFormLayout,QLineEdit,QPushButton,QLabel,
 QComboBox,QCheckBox,QFileDialog,QMessageBox,QGroupBox,QPlainTextEdit)
from .geothermal import DEFAULT_SETTINGS,load_settings,save_settings,geothermal_status
from .tomography import tomography_status


def _button(text,fn):
    w=QPushButton(text);w.clicked.connect(fn);return w


class GeothermalPanel(QWidget):
    def __init__(self,parent=None):
        super().__init__(parent);self.main=parent;self.root=None
        v=QVBoxLayout(self)
        note=QLabel('下游解释层：Vs空间层析 + 地质/断裂/钻井/温泉。当前1-D条件结构不会被平铺成二维速度场。')
        note.setWordWrap(True);v.addWidget(note)
        form=QFormLayout();v.addLayout(form);self.fields={}
        for key,label,filters in [
            ('tomography_file','空间Vs层析（CSV/NC/GRD）','Spatial Vs (*.csv *.txt *.nc *.netcdf *.grd);;All files (*)'),
            ('geology','区域地质（SHP/GeoJSON/GPKG）','Vector (*.shp *.geojson *.json *.gpkg);;All files (*)'),
            ('faults','断裂（SHP/GeoJSON/GPKG）','Vector (*.shp *.geojson *.json *.gpkg);;All files (*)'),
            ('wells','钻井点（CSV/矢量）','Points (*.csv *.shp *.geojson *.gpkg);;All files (*)'),
            ('hot_springs','温泉点（CSV/矢量）','Points (*.csv *.shp *.geojson *.gpkg);;All files (*)')]:
            row=QHBoxLayout();edit=QLineEdit();self.fields[key]=edit;row.addWidget(edit,1);row.addWidget(_button('浏览',lambda _,k=key,f=filters:self.browse(k,f)));form.addRow(label,row)
        self.depths=QLineEdit('500,1000,2000,5000');form.addRow('Vs切片深度 / m',self.depths)
        self.depthtol=QLineEdit();self.depthtol.setPlaceholderText('留空=按深度层间距自动');form.addRow('深度匹配容差 / m',self.depthtol)
        self.backend=QComboBox();self.backend.addItems(['auto','pygmt','matplotlib']);form.addRow('绘图后端',self.backend)
        self.vsunit=QComboBox();self.vsunit.addItems(['auto','m/s','km/s']);form.addRow("含糊列名'vs'的单位",self.vsunit)
        self.griddepth=QLineEdit();self.griddepth.setPlaceholderText('仅2D Vs网格需要');form.addRow('二维网格代表深度 / m',self.griddepth)
        self.maxstd=QLineEdit();self.maxstd.setPlaceholderText('如 0.10；留空不筛');form.addRow('最大相对Vs标准差',self.maxstd)
        self.minres=QLineEdit();self.minres.setPlaceholderText('留空不筛');form.addRow('最小 resolution',self.minres)
        self.minrays=QLineEdit();self.minrays.setPlaceholderText('留空不筛');form.addRow('最少 ray_count',self.minrays)
        self.selected=QCheckBox('层析交接仅导出 selection.json 当前选定曲线');form.addRow(self.selected)
        self.labels=QCheckBox('地图标注台站名称');self.labels.setChecked(True);form.addRow(self.labels)
        self.title=QLineEdit('Geothermal interpretation');form.addRow('成果图标题',self.title)

        tomo=QGroupBox('二维相速度层析（直线路径）');tf=QFormLayout(tomo)
        self.tomocase=QLineEdit();self.tomocase.setPlaceholderText('可选：selection.json case 名称或序号');tf.addRow('层析选择 case',self.tomocase)
        self.tomoquality=QComboBox();self.tomoquality.addItems(['default_pass','strict_direction','all_valid']);tf.addRow('输入质量模式',self.tomoquality)
        self.tomodup=QComboBox();self.tomodup.addItems(['error','median']);tf.addRow('同台对/频率多候选',self.tomodup)
        self.tomogrid=QLineEdit('12,12');tf.addRow('网格 nx,ny',self.tomogrid)
        self.tomoangle=QLineEdit('15');tf.addRow('投影波最大台对夹角 / °',self.tomoangle)
        self.tomoreg=QLineEdit('0.15,1.0');tf.addRow('阻尼,平滑系数',self.tomoreg)
        self.tomominpaths=QLineEdit('6');tf.addRow('每频率最少路径',self.tomominpaths)
        self.tomomincell=QLineEdit('1');tf.addRow('单元最少穿越射线',self.tomomincell)
        self.tomocheck=QCheckBox('运行 checkerboard 分辨率测试');self.tomocheck.setChecked(True);tf.addRow(self.tomocheck)
        v.addWidget(tomo)

        row=QHBoxLayout();row.addWidget(_button('保存地热设置',self.save));row.addWidget(_button('导出层析射线路径',lambda:self.main.start('tomography-export')));row.addWidget(_button('运行二维相速度层析',lambda:self.main.start('tomography-invert')));row.addWidget(_button('验证空间Vs',lambda:self.main.start('geothermal-validate')));row.addWidget(_button('生成地热成果',lambda:self.main.start('geothermal')));row.addStretch();v.addLayout(row)
        row=QHBoxLayout();row.addWidget(_button('打开二维层析报告',self.open_tomography_report));row.addWidget(_button('打开地热报告',self.open_report));row.addWidget(_button('打开地热输出目录',self.open_output));row.addWidget(_button('安装/格式说明',self.open_guide));row.addStretch();v.addLayout(row)
        self.status=QPlainTextEdit();self.status.setReadOnly(True);self.status.setMaximumHeight(210);v.addWidget(self.status)
    def browse(self,key,filters):
        p,_=QFileDialog.getOpenFileName(self,'选择'+key,'',filters)
        if p:self.fields[key].setText(p)
    def _num(self,w):
        t=w.text().strip();return None if not t else float(t)
    def collect(self):
        grid=[int(x.strip()) for x in self.tomogrid.text().split(',') if x.strip()]
        reg=[float(x.strip()) for x in self.tomoreg.text().split(',') if x.strip()]
        if len(grid)!=2:raise ValueError('二维层析网格请输入 nx,ny，例如 12,12')
        if len(reg)!=2:raise ValueError('阻尼和平滑请输入两个数，例如 0.15,1.0')
        return dict(DEFAULT_SETTINGS,
            tomography_file=self.fields['tomography_file'].text().strip(),geology=self.fields['geology'].text().strip(),faults=self.fields['faults'].text().strip(),wells=self.fields['wells'].text().strip(),hot_springs=self.fields['hot_springs'].text().strip(),
            slice_depths_m=[float(x.strip()) for x in self.depths.text().split(',') if x.strip()],depth_tolerance_m=self._num(self.depthtol),backend=self.backend.currentText(),tomography_vs_unit=self.vsunit.currentText(),grid_depth_m=self._num(self.griddepth),max_relative_vs_std=self._num(self.maxstd),min_resolution=self._num(self.minres),min_ray_count=self._num(self.minrays),selected_paths_only=self.selected.isChecked(),label_points=self.labels.isChecked(),title=self.title.text().strip() or 'Geothermal interpretation',
            tomo_selection_case=self.tomocase.text().strip(),tomo_quality_mode=self.tomoquality.currentText(),tomo_duplicate_policy=self.tomodup.currentText(),tomo_nx=grid[0],tomo_ny=grid[1],tomo_max_direction_angle_deg=float(self.tomoangle.text().strip()),tomo_damping=reg[0],tomo_smoothing=reg[1],tomo_min_paths_per_frequency=int(self.tomominpaths.text().strip()),tomo_min_cell_rays=int(self.tomomincell.text().strip()),tomo_checkerboard=self.tomocheck.isChecked())
    def save(self):
        try:
            if self.root is None:raise ValueError('请先新建或打开项目')
            save_settings(self.root,self.collect());self.status.setPlainText('地热设置已保存：'+str(Path(self.root)/'geothermal_settings.json'))
        except Exception as exc:QMessageBox.warning(self,'地热设置未保存',str(exc))
    def load(self,root):
        self.root=Path(root);s=load_settings(self.root)
        for k in self.fields:self.fields[k].setText(str(s.get(k,'') or ''))
        self.depths.setText(','.join(f'{float(x):g}' for x in s['slice_depths_m']))
        self.depthtol.setText('' if s['depth_tolerance_m'] is None else str(s['depth_tolerance_m']))
        self.backend.setCurrentText(s['backend']);self.vsunit.setCurrentText(s['tomography_vs_unit']);self.griddepth.setText('' if s['grid_depth_m'] is None else str(s['grid_depth_m']))
        self.maxstd.setText('' if s['max_relative_vs_std'] is None else str(s['max_relative_vs_std']));self.minres.setText('' if s['min_resolution'] is None else str(s['min_resolution']));self.minrays.setText('' if s['min_ray_count'] is None else str(s['min_ray_count']))
        self.selected.setChecked(bool(s['selected_paths_only']));self.labels.setChecked(bool(s['label_points']));self.title.setText(s['title'])
        self.tomocase.setText(str(s.get('tomo_selection_case','') or ''));self.tomoquality.setCurrentText(s['tomo_quality_mode']);self.tomodup.setCurrentText(s['tomo_duplicate_policy']);self.tomogrid.setText(f"{int(s['tomo_nx'])},{int(s['tomo_ny'])}");self.tomoangle.setText(str(s['tomo_max_direction_angle_deg']));self.tomoreg.setText(f"{s['tomo_damping']},{s['tomo_smoothing']}");self.tomominpaths.setText(str(s['tomo_min_paths_per_frequency']));self.tomomincell.setText(str(s['tomo_min_cell_rays']));self.tomocheck.setChecked(bool(s['tomo_checkerboard']));self.reload_result()
    def reload_result(self):
        if self.root is None:return
        st=geothermal_status(self.root);ts=tomography_status(self.root)
        if not st.get('complete'):
            path=Path(self.root)/'output/geothermal/tomography_paths_manifest.json'
            summary={'geothermal':st,'tomography2d':ts}
            if path.exists():summary['tomography_paths_manifest']=json.loads(path.read_text(encoding='utf-8-sig'))
            self.status.setPlainText(json.dumps(summary,ensure_ascii=False,indent=2)[:7000]);return
        summary={'products':len(st.get('products',[])),'warnings':st.get('warnings',[]),'spatial_vs':st.get('spatial_vs'),'tomography_paths':st.get('tomography_paths'),'tomography2d':ts}
        self.status.setPlainText(json.dumps(summary,ensure_ascii=False,indent=2)[:7000])
    def open_tomography_report(self):
        if self.root:QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(self.root)/'output/geothermal/tomography2d/tomography_2d_report.html')))
    def open_report(self):
        if self.root:QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(self.root)/'output/geothermal/report.html')))
    def open_output(self):
        if self.root:QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(self.root)/'output/geothermal')))
    def open_guide(self):
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(__file__).resolve().parent.parent/'docs/Geothermal_Workflow_v260.md')))
