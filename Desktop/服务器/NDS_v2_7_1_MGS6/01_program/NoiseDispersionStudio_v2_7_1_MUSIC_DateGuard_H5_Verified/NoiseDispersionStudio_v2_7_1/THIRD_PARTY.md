# 数据、字体和依赖

- 新编写的源代码：MIT，见 LICENSE.txt。
- 示例合成波形和已知模型由本程序生成；只是数值测试，没有作为真实地球观测声明。
- DARE 示例观测由此前公开 Y7 台网数据处理得到：法国 2020–2021 Dense Array for site effect studies in the Rhône valley (DARE)，doi:10.14470/L27575187372，GEOFON/RESIF，CC BY 4.0。只附主段候选表与处理试验，未附原始波形。数据 DOI：https://doi.org/10.14470/L27575187372 。来源、观测日期与先前局限见 examples/dare_primary/README.md。
- Noto Sans CJK SC Regular 字体来自 https://github.com/notofonts/noto-cjk ，静态字体源 https://raw.githubusercontent.com/notofonts/noto-cjk/main/Sans/OTF/SimplifiedChinese/NotoSansCJKsc-Regular.otf ，SIL Open Font License，完整文本见 assets/OFL-NotoSansSC.txt。
- NumPy、SciPy、pandas、h5py、ObsPy、pyproj、Matplotlib、Pillow、PySide6/Qt、disba、Numba/llvmlite 通过 pip 独立安装，源码包未包含其运行库；各自许可证以项目发行包为准。可选 pysurf96 同样独立安装。
- 地热/空间制图可选依赖 GeoPandas、xarray、netCDF4、PyGMT 同样通过 pip 独立安装；PyGMT 还需要系统 GMT 6.5+。未安装这些可选包时，基础 Matplotlib 地热切片仍可工作（矢量地质文件叠加需要 GeoPandas）。
- Qt for Python 官方：https://doc.qt.io/qtforpython-6/ ，桌面耗时进程使用 QProcess。h5py 官方：https://docs.h5py.org/en/stable/ 。disba 官方：https://keurfonluu.github.io/disba/ ，单位km/km·s⁻¹与程序m/m·s⁻¹在边界转换，基阶mode=0。
- 方法参考 Bensen et al. (2007), Processing seismic ambient noise data to obtain reliable broad-band surface wave dispersion measurements. doi:10.1111/j.1365-246X.2007.03374.x 。本程序的阈值与检验属于可复现的工程实现，不声称完整实现论文所有方法或对所有数据保证可用。

HTML中嵌入的NDSDocs字体是上述Noto字体的重命名字符子集，SIL OFL；字体许可证也嵌入HTML注释，源全字体许可证随包提供。
