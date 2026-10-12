"""Format the completed Studio tests and their existing figures, without reanalysis."""
from pathlib import Path
import base64,html,json,re,shutil

SOURCE=Path(__file__).resolve().parents[2]
BASE=SOURCE.parents[1]
OUT=BASE/'output'

def embed(path,caption):
    data=base64.b64encode(Path(path).read_bytes()).decode()
    return f'<figure><img src="data:image/png;base64,{data}" alt="{html.escape(caption)}"><figcaption>{html.escape(caption)}</figcaption></figure>'

def update_document(path,title,section):
    old=path.read_text()
    old=re.sub(r'<section id="v22-real-test">.*?</section>','',old,flags=re.S)
    old=re.sub(r'<title>.*?</title>',f'<title>{title}</title>',old,count=1,flags=re.S)
    extra='<style>figure{margin:22px 0}figure img{width:100%;height:auto}pre{white-space:pre-wrap;overflow-wrap:anywhere} .v22-note{border-left:4px solid #196a7c;padding:12px 18px;background:#edf7f8}td,th{vertical-align:top}</style>'
    old=old.replace('</head>',extra+'</head>',1)
    old=re.sub(r'(<body\b[^>]*>)',lambda m:m.group(1)+'<section id="v22-real-test">'+section+'</section><hr><h2>此前版本说明与测试记录</h2>',old,count=1)
    path.write_text(old)

def main():
    roots=[BASE/'verification'/x for x in ['real_e2e','second_e2e']]
    summaries=[];rows=[]
    for root,label in zip(roots,['法国 DARE / Y7','爱尔兰都柏林 Geo_URBAN / 1V']):
        p=root/'project';d=root/'download'
        selection=json.loads((d/'analysis_array_selection.json').read_text())
        verify=json.loads((p/'output/verification.json').read_text())
        marker=json.loads((p/'output/inversion_complete.json').read_text())
        audit=json.loads((d/'acquisition_audit.json').read_text())
        assert verify['checks_passed'] and audit['checks_passed']
        models=json.loads((p/'output/conditional_models.json').read_text())['models']
        assert len(models)==marker['compatible_models'] and all(m['independent_solver_pass'] for m in models)
        entry=dict(array=label,network=selection['network'],requested_start_utc=selection['start_utc'],requested_end_utc=selection['end_utc'],
          stations=[r['station'] for r in selection['stations']],geometry=selection['geometry'],acquisition_audit=audit,
          verification=verify,inversion=marker,conditional_model_count=len(models))
        summaries.append(entry)
        g=selection['geometry']
        rows.append(f'<tr><td>{label}</td><td>{selection["start_utc"][:10]} 至 {selection["end_utc"][:10]}<br>96 h</td><td>{len(selection["stations"])}<br>{", ".join(entry["stations"])}</td><td>{g["mean_nearest_neighbor_m"]:.2f} m<br>孔径 {g["max_pair_distance_m"]:.2f} m</td><td>{audit["h5_files"]}</td><td>{verify["windows"]}</td><td>{verify["measured_segments"]}</td><td>{verify["conditional_curves"]}</td><td>{marker["cases"]}</td><td>{len(models)}</td></tr>')
    full_tests=(roots[0]/'full_tests_final.log').read_text()
    assert re.search(r'80 passed',full_tests)
    table='<table><thead><tr><th>阵列</th><th>请求 UTC</th><th>最终台站</th><th>最近邻平均间距 / 孔径</th><th>全部原始 H5</th><th>有效窗</th><th>稳健相位段</th><th>条件频散曲线</th><th>反演组</th><th>通过结构</th></tr></thead><tbody>'+''.join(rows)+'</tbody></table>'
    report='<h1>2.2：两套真实宽频带阵列端到端测试</h1><p>测试日期：2026-10-09。下载、选站、预处理、相位检测、背景对照、频散枚举及结构反演均通过本程序的模块和命令完成。报告汇总已完成的输出，不重新分析波形。第二轮为此前未测试的爱尔兰都柏林阵列，不复用 DARE 波形或曲线。</p>'+table
    report+='<div class="v22-note">这些是有明确整数周期、源相位和 Rayleigh 模态假设的条件频散及条件结构。通过独立数值求解器核对，说明数值实现一致；并不证明唯一地质结构、绝对周期或模态。层深度是拟合参数，搜索上限和图上深度不是已证实的探测深度。全部候选与拒绝尝试保留。</div>'
    report+='<h2>真实数据与选站记录</h2><p>两轮均先查询完整通道与响应元数据，再用程序按几何选出 10 台紧凑子阵列，要求每台最近邻距离在 200–2000 m、最大孔径不超过 6000 m，原生采样率均为 100 Hz。最终每轮保留 6 台，分别满足第一轮 5&lt;n&lt;50 和第二轮 5–100 台的要求。</p><p>DARE 使用 Trillium Compact 120；都柏林使用 CMG 6T（30 s–100 Hz）及 CMG 40T 宽频带仪器的 HHZ。都柏林 G012、G020 的请求被服务器明确返回无数据，其余 8 台共 32 个每日 H5 均保存。下载完成状态表示所有请求已处理，all_requested_chunks_have_h5=false 明确记录无数据，不能解释为 40 个文件全部下载成功。</p>'
    report+='<p>DARE 的 E02 出现极端直流偏置与很小动态范围，按明确的原始 counts DC/MAD 条件排除；这不等于已证明硬件故障。E01 在其后违反最近邻上限；C03 的时域有效比例不足 0.5。7 台交集的独立验证块不足后，仅按时间覆盖排序再排除 D01。都柏林 G015 的有效比例不足 0.5，7 台验证块不足后按时间覆盖再排除 G011。此过程未用相位或拟合结果选站，全部已下载原始记录保留，详细比例和几何见筛选 JSON 与预处理诊断表。</p>'
    report+='<h2>完整流程和原有门限</h2><p>原生 counts 和响应导入 → 对齐实际 UTC → 响应校正为 VEL → 带通与局部坏点/缺测掩码 → 20 Hz 处理、20 min 窗、共同起点后舍弃前 600 s → 多归一化/白化路线 → 训练与验证时间块 → 群到时和稳健相位段 → 40 次错时对照 → 固定每段整数周期和源相位条件的频散枚举 → 两两跨台对兼容的反演组 → 多起点层状 Rayleigh 反演 → 独立正演核对 → 分层深度表和图件。</p><p>保留原有质量门限：相位相对 RMS ≤3%、最大点误差 ≤5%、每带群速度相对误差中位数 ≤15%、独立相速度差 ≤0.5 m/s。使用 1/2 个有限层及半空间、8 起点、每起点最多 180 次评价，Vp/Vs=2、ρ=2 g/cm³；每段源相位固定为 π/4，模态假设为 0。缺失模态根不补值；高密度相邻频点不被当成独立样本。</p><p>DARE 检验默认前 12 个自动组。都柏林前 12 组没有通过模型，因此将自动组上限设为 1000，实际检验全部 219 个由程序生成的组，保持全部质量门限与先验不变。初始 12 组的拒绝记录单独归档，不能将新增检验误称为默认 12 组直接成功。</p>'
    report+='<h2>测试中修复和验证</h2><ul><li>并行下载的同一 StationXML 原子落盘临时文件冲突；下载任务串行汇总完整子阵列清单，断点完成状态核对实际 H5。</li><li>服务器重复、嵌套重叠记录的导入：只对 UTC 网格与数值逐点完全一致的前缀去重；原始记录仍完整保存在 H5。</li><li>重复响应查询仅因 StationXML Source/Created 头变化时允许导入；完整通道与响应语义必须一致，真实灵敏度或响应变化继续拒绝。</li><li>预处理失败前仍落盘站点、记录和局部质量诊断；时间覆盖筛选严格核对台数与每台最近邻。</li><li>独立求解器缺根异常改为记录该次拒绝并继续，禁止补根；新增 structure_layers.csv，明确有限层顶/底与半空间。</li></ul><p>完整回归：<strong>80 passed，0 skipped</strong>。实际 Qt/Tk 窗口运行检查通过；下图来自真实 GUI。两组原始 H5 均逐记录验证实际 UTC、原生采样率和 counts 与内嵌原服务器 MiniSEED 完全一致，并验证完整 StationXML 和原始字节 SHA256。程序验证确认原始 H5 在分析前后未改写。</p>'
    for root,label,shot in zip(roots,['DARE','Dublin'],['v22_dare_array.png','v22_dublin_array.png']):
        report+=f'<h2>{label} 真实输出图</h2>'
        report+=embed(SOURCE/'docs/screenshots'/shot,label+' 下载 GUI：六台子阵列、米制坐标与平均间距')
        for filename,caption in [('03_measured_coverage.png','经验证的实测相位覆盖'),('04_conditional_families.png','保留全部条件周期候选的频散曲线'),('10_structure_models.png','通过拟合条件的 Vs 分层深度及观测/预测分离')]:
            file=root/'project/output/figures'/filename
            if file.exists():report+=embed(file,label+'：'+caption)
    report+='<h2>交付与复现</h2><p>完整程序 ZIP 含源码、GUI、依赖列表、回归测试、复现驱动及本次报告。真实结果 ZIP 含两组全部 CSV/JSON/图件、分层深度、拒绝尝试、配置、原始数据审计和运行日志。每个阵列原始 H5 分为四个 UTC 日包，解压到同一目录即可合并。原服务器 MiniSEED 已完整内嵌在 H5，因此不再重复存放外置副本；未完成的 .part 临时文件不作为数据交付。缓存和 canonical 中间数组可由原始 H5 用程序重建，不作为原始数据。</p><p>需要准确路径时使用完整包中的 README 与复现驱动。原始文件摘要见原始H5_分包索引.json，分层结果在各项目 output/structure_layers.csv；观测点和模型预测分别保存，不能将插值预测当成新观测。</p>'
    manual='<h1>2.2：真实阵列下载、诊断与完整流程</h1><p>本版加入程序内置的几何选站、下载、完整 H5 审计、时域可用性诊断与子阵列筛选。GUI 下载前可查看阵列和米制最近邻平均台站间距；选中通道的几何重新计算，下载清单保存下载时的完整几何。详细实测结果见本版测试报告。</p>'
    manual+='<h2>先保存 H5，再做频散与深度反演</h2><ol><li>运行 python run.py download；选 Provider、阵列、日期、垂直宽频带通道，在 Arrays/Channels 查看台站分布、仪器与平均间距 m。</li><li>保持 H5 + full metadata 选中，下载后点击 Open H5 in dispersion。原始 H5 保留 counts、实际 UTC、NSLC、坐标/历元、完整 StationXML 和原 MiniSEED 字节，读取无需另找响应文件。</li><li>扫描/导入后检查记录质量与时间覆盖。台站排除和局部掩码分别记录；不因某组拟合好而选站。counts 必须通过真实响应校正。</li><li>逐步运行 prepare、scan、controls、families、review、verify、autoselect、invert、report。自动组合默认上限 12，可以在 GUI 修改；需穷尽自动组合时设 1000，仍保留未入选候选。</li><li>查看 measured_phase_segments.csv、conditional_dispersion.csv、inversion_attempts.csv 和 structure_layers.csv。无通过模型时完成报告并保留拒绝记录，不补造结构。</li></ol>'
    manual+='<h2>程序自身的命令与复现</h2><pre>python -m pip install -r requirements.txt\npython -m pip install -r requirements-optional.txt\n# 在线完整 DARE 测试\npython examples/real_dare_e2e/run_e2e.py --workdir DARE_真实端到端 --max-cases 12\n# 在线不同宽频带阵列：都柏林，检验全部自动组合\npython examples/real_dare_e2e/run_e2e.py --workdir Dublin_真实端到端 --network 1V --start 2019-04-01 --end 2019-04-05 --size 10 --max-cases 1000\n# 先将四个对应原始H5日包解压到同一目录，离线重建该阵列全部流程\npython examples/real_dare_e2e/run_e2e.py --offline --workdir Dublin_真实端到端 --network 1V --max-cases 1000\n# 只审计已保存原始H5，无需重新下载\npython run.py acquisition-audit --folder Dublin_真实端到端/download\n# 手动逐步计算：每条命令均使用本程序\npython run.py import --project Dublin_真实端到端/project\npython run.py prepare --project Dublin_真实端到端/project\npython run.py diagnose-prepare --project Dublin_真实端到端/project\npython run.py scan --project Dublin_真实端到端/project\npython run.py controls --project Dublin_真实端到端/project\npython run.py families --project Dublin_真实端到端/project\npython run.py review --project Dublin_真实端到端/project\npython run.py verify --project Dublin_真实端到端/project\npython run.py autoselect --project Dublin_真实端到端/project --max-cases 1000\npython run.py invert --project Dublin_真实端到端/project\npython run.py report --project Dublin_真实端到端/project</pre>'
    manual+='<p>复现驱动在默认门限下先检查原始时域质量、时间覆盖和几何，验证块不足时按可用时间选子阵列并重新导入/预处理。空候选或零个通过结构是记录明确的科学结果；complete 表示程序处理完成，structure_accepted 单独表示是否获得条件结构。可加 --require-model 将没有通过结构视作命令失败。网络/数据服务变化可能改变实际下载返回；已保存日包用于精确离线复现。</p>'
    manual+='<h2>阵列限制与原始数据保护</h2><p>array-select 默认台数限制 6–49；用 --min-stations 5 --max-stations 100 可配置第二轮所需范围。--size 为初始候选台数。平均间距采用每台到最近另一台的 WGS84 地表距离的平均值；全部台对平均距离另外显示。程序对每台最近邻上下限和孔径分别检查，不能只用平均数代替。screen-prepared --size 6 仅按时间覆盖选紧凑子阵列，排除记录原始 H5。</p><p>文件夹使用下载锁，重复启动不会相互覆盖；崩溃后操作系统释放锁。响应 XML 原子写入。续传既处理已下载 MiniSEED 的 H5 转换，也核对实际 H5 完成状态。明确无数据与网络错误分别记录。所有筛选只影响分析清单，原始信息和排除原因保留。</p>'
    manual+='<h2>读懂结构深度表</h2><p>structure_layers.csv 每行对应一个模型的一层：top_depth_m、bottom_depth_m、thickness_m、vs_m_s、vp_m_s 和 density_g_cm3；halfspace 的底深与厚度留空。模型 ID 和 case 对应明确的条件曲线、周期、源相位与模态。unique_structure_identified=false、depth_resolution_established=false 表明未证明唯一结构或深度分辨率。不能把多个候选模型合成有统计意义的置信区间。</p>'
    manual+=table
    update_document(OUT/'通用频散程序_测试报告.html','2.2 双阵列真实端到端 · 通用频散程序测试报告',report)
    update_document(OUT/'通用频散程序_使用手册.html','2.2 真实阵列完整流程 · 通用频散程序使用手册',manual)
    shutil.copyfile(OUT/'通用频散程序_测试报告.html',SOURCE/'docs/测试报告.html')
    shutil.copyfile(OUT/'通用频散程序_使用手册.html',SOURCE/'docs/使用手册.html')
    (OUT/'双阵列_端到端摘要.json').write_text(json.dumps({'program_version':'2.2.0','tests_passed':80,'tests_skipped':0,'arrays':summaries},ensure_ascii=False,indent=2)+'\n')
    print('Updated manual, test report and machine-readable summary.')

if __name__=='__main__':main()
