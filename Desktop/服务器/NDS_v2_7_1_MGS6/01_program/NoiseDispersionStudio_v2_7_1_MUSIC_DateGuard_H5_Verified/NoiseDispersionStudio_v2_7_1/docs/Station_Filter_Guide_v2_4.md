# Studio 2.4 台站范围与前缀筛选

主窗口点击“下载数据 / 阵列预览”，或双击 `run_download_windows.bat`。

下载窗口顶部新增 Station filters。先选择服务商和阵列，再填以下字段：

| 界面字段 | 含义 | 示例 |
|---|---|---|
| Station prefix | 台站代码的开头，不是台网代码。大小写不敏感 | `QS` |
| Latitude min | 纬度下限，南界，范围 −90 到 90 | `47` |
| Latitude max | 纬度上限，北界 | `49` |
| Longitude west | 经度西界，范围 −180 到 180 | `11` |
| Longitude east | 经度东界 | `11.04` |

经纬度边界留空表示该边不限。前缀留空表示不按前缀限制。范围和前缀取交集，边界上的台站也包含在内。
前缀写 `QS` 即可，不需要星号或省略号；也支持 `QS*`。多个前缀写成 `QS,QT`，它们之间取并集，再与地理范围取交集。
原来的 Stations 字段仍支持通配符，两种名称条件也取交集。
例如前缀 `QS`、Stations `QS0?` 表示只取 QS00 到 QS09。

## 从 100 个台站中取得 QS 开头的 20 个

1. 在 Arrays 页面选择阵列，等待台站预览出现。
2. 在 Station prefix 填入 `QS`，按需填经纬度边界。
3. 点击 Apply filters，或在筛选框中按 Enter。
4. 下方显示完整匹配名单和数量，例如 `20 matching stations`。
5. 选好要下载的通道，例如 HHZ，再点击 Download matched stations。
6. 用 Export stations 保存 CSV 名单，每个台站一行，包括坐标、通道、采样率和可用时间范围。

前缀匹配多少个就保留多少个，筛选不会为了凑够 20 个而自动裁剪台站。
下方预览表按通道历元列行，同一个台站有 HHZ、HHN、HHE 时会占三行，但匹配台站数按 NET.STA 去重。
Arrays 上方的台网台站数由服务商提供，可能是台网总数；实际筛选数量以下方 matching stations 为准。
只有窗口设定与通道历元相交的通道会形成下载任务，计划数量在下方 Selected 行和下载确认中显示。
每个子阵列建议使用独立输出目录，例如 `D:\data\QS_subset`，再从该目录打开频散项目。

## Channels 页面

同样填写顶部 Station filters，填好 Network、Station、Location、Channel 后点击 Browse 或 Apply filters。
列表显示全部匹配行，不再只显示前 5000 行；可多选要下载的通道。
Export stations 导出本次检索的全部匹配台站。Download selected 只下载列表中选中的通道。

## 边界与无匹配结果

经度按西界和东界填写。如果西界大于东界，如 `170` 到 `−170`，表示跨越 180° 经线。
启用地理范围时，坐标缺失或无效的记录不会匹配；未启用地理范围时，它们仍可按前缀列出。
无匹配结果时列表为空，不能继续用旧名单下载。修改条件会立即清空两页旧台站列表，过期查询返回的结果也会被忽略。
检索失败会显示错误，不会把失败当成完整的空结果。台站检索只读取元数据，列出通道不保证指定时段一定存在波形；下载仍记录缺测与服务商错误。

## 命令行

导出完整匹配名单，不执行紧凑子阵列裁剪：

```text
.venv\Scripts\python.exe run.py station-list --provider GFZ --network Y7 --channel "*" --station-prefix C --min-latitude 44.28 --max-latitude 44.32 --min-longitude 4.66 --max-longitude 4.68 --start 2020-10-12T12:00:00Z --end 2020-10-12T12:30:00Z --output "D:\station_catalog"
```

生成 `stations.csv`（台站去重）、`channels.csv`（全部通道历元）、`station_query.json`（完整查询条件与结果）。
这里的 Y7 / C 是 DARE 的示例；检索自己的 QS 阵列时替换台网代码和前缀。
`station-list` 默认通道 HHZ，列出所有通道需明确指定 `--channel "*"`。
指定起止时间时，只检索在该时段活动的通道元数据；不指定则查询服务商保存的匹配历元。

`array-select` 也接受同一组筛选参数，它在筛选后的候选台站中继续执行原来的间距与孔径约束。
`station-list` 的 JSON 是目录结果，不是 `fetch` 的阵列选择输入；命令行下载仍用 `array-select` 后接 `fetch`。

## 验证

100 台站离线测试阵列含 20 个 QS 台站、80 个 AA 台站，每站三个通道。
真实 Tk 界面筛选后显示 20 个匹配台站、60 个通道历元，选择 HHZ 后生成的 20 个下载任务全部属于 QS 台站。
覆盖了前缀与范围交集、包含边界、多个前缀、跨经线查询与去重、无匹配、无效输入、完整名单导出及旧结果失效。
使用既有 DARE StationXML：前缀 C、纬度 44.28–44.32、经度 4.66–4.68，正确得到 C00、C01、C02。
这次没有向服务商重新下载波形。检索和筛选未改变既有频散计算与相位筛选方法。
