import fnmatch
import os
import re
from io import BytesIO
import h5data_handoff as handoff
from dataclasses import dataclass, field
from pathlib import Path

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNNoDataException
import h5py
import numpy as np
from station_filters import (normalize_station_filters, geographic_queries, station_query_pattern,
                            filter_station_rows, unique_channel_rows, station_catalog, export_station_catalog)
from obspy import read

# --------------------------------------------------------------------------
# Constants and types used by the browse/download backend (and the GUI)
# --------------------------------------------------------------------------

DEFAULT_PROVIDER = "ETH"
COMMON_PROVIDERS = ("ETH", "IRIS", "GFZ", "ORFEUS", "INGV", "RESIF")
DEFAULT_CHUNK_SECONDS = 86400
BYTES_PER_SAMPLE = 4
MSEED_COMPRESSION_RATIO = 3.0
CONFIRM_RAW_BYTES = 500 * 1024 ** 2

TIME_FMT = "%Y-%m-%dT%H:%M:%S"
TIME_PRESETS = ("Entire available range", "Last hour", "Last 24 hours", "Last 7 days",
                "Last 30 days", "Yesterday (UTC)", "Custom")
PRESET_ENTIRE = TIME_PRESETS[0]
PRESET_CUSTOM = "Custom"
_RELATIVE_PRESETS = {"Last hour": 3600, "Last 24 hours": 86400,
                     "Last 7 days": 7 * 86400, "Last 30 days": 30 * 86400}
CHUNK_PRESETS = {"1 hour": 3600, "6 hours": 21600, "1 day": 86400, "7 days": 604800}
MIN_CHUNK_SECONDS = 60
MAX_CHUNK_SECONDS = 30 * 86400
CHANNEL_FILTER_PRESETS = ("all", "*", "HH?", "BH?", "HHZ", "BHZ")


class DataRequestError(RuntimeError):
    """Raised when an FDSN provider cannot be reached or rejects a request."""


@dataclass(frozen=True)
class ChannelRow:
    """One channel epoch returned by a station query."""
    network: str
    station: str
    location: str
    channel: str
    sample_rate: float | None
    latitude: float | None
    longitude: float | None
    start: UTCDateTime | None
    end: UTCDateTime | None
    elevation: float | None = None
    sensor: str = ""

    @property
    def seed_id(self) -> str:
        """SEED id such as ``CH.ZUR..HHZ``."""
        return f"{self.network}.{self.station}.{self.location}.{self.channel}"


@dataclass(frozen=True)
class DownloadJob:
    """A channel plus the time span to download for it."""
    row: ChannelRow
    start: UTCDateTime
    end: UTCDateTime


@dataclass
class DownloadResult:
    """Outcome of downloading one channel (possibly in many chunks)."""
    seed_id: str
    chunks_total: int
    chunks_done: int = 0
    files: list = field(default_factory=list)
    h5_files: list = field(default_factory=list)
    skipped: list = field(default_factory=list)
    errors: list = field(default_factory=list)
    cancelled: bool = False
    existing: list = field(default_factory=list)  # chunks already on disk (resume)
    aborted: str | None = None  # set when the consecutive-failure cutoff hit


@dataclass(frozen=True)
class NetworkRow:
    """One network ("array") returned by a network-level station query."""
    code: str
    description: str
    start: UTCDateTime | None
    end: UTCDateTime | None
    station_count: int | None


@dataclass(frozen=True)
class DownloadProgress:
    """Overall progress of a run_downloads() call."""
    chunks_done: int
    chunks_total: int
    channels_done: int
    channels_total: int
    stations_done: int
    stations_total: int
    seed_id: str


class FailureCounter:
    """Counts consecutive request failures; ``record_error`` is True at the limit."""

    def __init__(self, limit):
        self.limit = int(limit)
        self.count = 0

    def record_error(self):
        self.count += 1
        return self.count >= self.limit

    def reset(self):
        self.count = 0


def download_mseed_data():
    # 1. 初始化客户端 
    # "ETH" 代表瑞士地震局 (SED, ETH Zürich) 的 FDSN 服务器
    # 如果想换成全球数据，可以使用 "IRIS" 或 "GFZ" 等
    print("正在连接瑞士地震局(ETH)数据中心...")
    client = Client("ETH")

    # 2. 设置你想要抓取的台网、台站、通道和时间范围
    # 示例：CH台网（瑞士）下的 ZUR（苏黎世）台站
    network = "CH"
    station = "ZUR"
    location = "*"
    channel = "HHZ"  # 高采样率高灵敏度垂直向地动速度通道
    
    # 设置时间范围 (UTC时间)
    # 示例获取 2026年3月1日 某个小时内的数据
    start_time = UTCDateTime("2026-03-01T00:00:00")
    end_time = UTCDateTime("2026-03-01T01:00:00")

    print(f"正在请求数据: {network}.{station}.{location}.{channel}")
    print(f"时间范围: {start_time} 至 {end_time}")

    try:
        # 3. 抓取波形数据
        st = client.get_waveforms(network, station, location, channel, start_time, end_time)
        
        # 4. 打印抓取到的流(Stream)信息
        print("\n成功获取数据！")
        print(st)

        # 5. 保存为标准的 mseed 文件
        output_filename = f"{network}_{station}_{channel}.mseed"
        st.write(output_filename, format="MSEED")
        print(f"数据已成功保存至本地: {output_filename}")
        
        # 6. (可选) 绘制波形图查看
        # st.plot()

    except Exception as e:
        print(f"\n数据抓取失败，原因: {e}")
        print("请检查台站名是否正确，或者该时间段内该台站是否有公开数据。")


def convert_mseed_to_h5_strict(mseed_filename, h5_filename, verbose=False, **metadata):
    """Atomic, lossless archive; optional full metadata enables Studio import.

    Legacy two-argument calls retain counts and are marked incomplete. Downloads
    always request full metadata when HDF5 is enabled. See h5data_handoff.py.
    """
    return handoff.export_h5(mseed_filename, h5_filename, verbose=verbose, **metadata)


def convert_mseed_to_h5(mseed_filename, h5_filename):
    """Convert MiniSEED to HDF5, printing errors instead of raising.

    Thin wrapper around :func:`convert_mseed_to_h5_strict`. Behavior change
    from the original: if a failure happens after the output file was opened,
    the partially written ``.h5`` is now removed.
    """
    try:
        convert_mseed_to_h5_strict(mseed_filename, h5_filename, verbose=True)
    except Exception as e:
        print(f"转换失败: {e}")


# --------------------------------------------------------------------------
# Browsing available channels
# --------------------------------------------------------------------------

def make_client(provider=DEFAULT_PROVIDER, timeout=120):
    """Create an FDSN client (short name such as "ETH" or a base URL).

    Raises ValueError for a blank provider and DataRequestError when the
    client cannot be created (network error, unknown provider, ...).
    """
    provider = (provider or "").strip()
    if not provider:
        raise ValueError("Provider must not be empty.")
    try:
        return Client(provider, timeout=timeout)
    except Exception as e:
        raise DataRequestError(f"Cannot connect to FDSN provider '{provider}': {e}") from e


def _opt_float(value):
    return None if value is None else float(value)


def _sensor_text(cha):
    """Best available sensor/instrument description of a channel ("" if none)."""
    sensor = getattr(cha, "sensor", None)
    if sensor is None:
        return ""
    for attr in ("description", "type", "model"):
        v = getattr(sensor, attr, None)
        if v:
            return str(v)
    return ""


def _norm_filter(v):
    v = (v or "").strip()
    return v if v else "*"


def inventory_to_rows(inventory):
    """Flatten an obspy Inventory into a sorted list of ChannelRow (one per epoch)."""
    rows = []
    for net in inventory:
        for sta in net:
            for cha in sta:
                lat = cha.latitude if cha.latitude is not None else sta.latitude
                lon = cha.longitude if cha.longitude is not None else sta.longitude
                elev = cha.elevation if cha.elevation is not None else sta.elevation
                rows.append(ChannelRow(
                    network=net.code,
                    station=sta.code,
                    location=cha.location_code or "",
                    channel=cha.code,
                    sample_rate=_opt_float(cha.sample_rate),
                    latitude=_opt_float(lat),
                    longitude=_opt_float(lon),
                    start=cha.start_date,
                    end=cha.end_date,
                    elevation=_opt_float(elev),
                    sensor=_sensor_text(cha),
                ))

    def key(r):
        return (r.network, r.station, r.location, r.channel,
                r.start is not None, r.start.timestamp if r.start is not None else 0.0)

    rows.sort(key=key)
    return rows


def list_stations(provider=DEFAULT_PROVIDER, network="CH", station="*", location="*",
                  channel="HHZ", starttime=None, endtime=None, client=None, *,
                  station_prefix='', minlatitude=None, maxlatitude=None,
                  minlongitude=None, maxlongitude=None):
    """Query channel-level metadata and return it as ChannelRow objects.

    Blank filters become ``*``. Returns ``[]`` when the provider reports no
    data. Raises DataRequestError for any other provider/network failure.
    """
    def norm(v):
        v = (v or "").strip()
        return v if v else "*"

    filters=normalize_station_filters(station_prefix=station_prefix,minlatitude=minlatitude,
        maxlatitude=maxlatitude,minlongitude=minlongitude,maxlongitude=maxlongitude)
    if client is None:
        client = make_client(provider)
    kwargs = dict(network=norm(network), station=station_query_pattern(station,filters),
                  location=norm(location), channel=norm(channel), level="channel")
    if starttime is not None:
        kwargs["starttime"] = UTCDateTime(starttime)
    if endtime is not None:
        kwargs["endtime"] = UTCDateTime(endtime)
    rows=[]
    for bounds in geographic_queries(filters):
        try:
            inv = client.get_stations(**kwargs,**bounds)
        except FDSNNoDataException:
            continue
        except Exception as e:
            raise DataRequestError(f"Station query failed: {e}") from e
        rows.extend(inventory_to_rows(inv))
    return filter_station_rows(unique_channel_rows(rows),station=station,**filters)


# --------------------------------------------------------------------------
# Arrays (FDSN networks, real or virtual) and their preview
# --------------------------------------------------------------------------

def inventory_to_networks(inventory):
    """NetworkRow per network, sorted by code. Station count may be None."""
    rows = []
    for net in inventory:
        count = net.total_number_of_stations
        if count is None:
            count = net.selected_number_of_stations
        if count is None and len(net.stations):
            count = len(net.stations)
        rows.append(NetworkRow(code=net.code, description=net.description or "",
                               start=net.start_date, end=net.end_date,
                               station_count=None if count is None else int(count)))
    rows.sort(key=lambda r: r.code)
    return rows


def list_networks(provider=DEFAULT_PROVIDER, network="*", station="*", starttime=None,
                  endtime=None, client=None, *, station_prefix='', minlatitude=None,
                  maxlatitude=None, minlongitude=None, maxlongitude=None):
    """List networks ("arrays") with a ``level="network"`` query.

    Returns ``[]`` when the provider reports no data; raises DataRequestError
    for any other failure.
    """
    filters=normalize_station_filters(station_prefix=station_prefix,minlatitude=minlatitude,
        maxlatitude=maxlatitude,minlongitude=minlongitude,maxlongitude=maxlongitude)
    if client is None:
        client = make_client(provider)
    kwargs = dict(network=_norm_filter(network), station=station_query_pattern(station,filters),
                  level="network")
    if starttime is not None:
        kwargs["starttime"] = UTCDateTime(starttime)
    if endtime is not None:
        kwargs["endtime"] = UTCDateTime(endtime)
    networks={}
    for bounds in geographic_queries(filters):
        try:
            inv = client.get_stations(**kwargs,**bounds)
        except FDSNNoDataException:
            continue
        except Exception as e:
            raise DataRequestError(f"Network query failed: {e}") from e
        for row in inventory_to_networks(inv): networks.setdefault(row.code,row)
    return [networks[key] for key in sorted(networks)]


def list_array_channels(provider, network, station="*", client=None, **station_options):
    """All channel epochs of one array (network code), optionally a station subset."""
    code = (network or "").strip()
    if not code or any(ch in code for ch in "*?,[]"):
        raise ValueError("Select a single array (network code).")
    return list_stations(provider, code, station, "*", "*", client=client, **station_options)


def array_geometry(rows):
    return handoff.geometry(rows)


def station_points(rows):
    """Unique (network, station, lat, lon, elevation), skipping missing coordinates."""
    seen, points = set(), []
    for r in rows:
        key = (r.network, r.station)
        if key in seen or r.latitude is None or r.longitude is None:
            continue
        seen.add(key)
        points.append((r.network, r.station, r.latitude, r.longitude, r.elevation))
    return points


def summarize_array(rows, now=None):
    """Counts and overall span of an array's channel epochs."""
    rows = list(rows)
    stations = {(r.network, r.station) for r in rows}
    with_coords = {(r.network, r.station) for r in rows
                   if r.latitude is not None and r.longitude is not None}
    starts = [r.start for r in rows if r.start is not None]
    open_ended = any(r.end is None for r in rows)
    ends = [r.end for r in rows if r.end is not None]
    return {
        "stations": len(stations),
        "channels": len({r.seed_id for r in rows}),
        "epochs": len(rows),
        "start": min(starts) if starts else None,
        "end": None if open_ended or not ends else max(ends),
        "open_ended": open_ended,
        "missing_coords": len(stations - with_coords),
    }


# --------------------------------------------------------------------------
# Input choices and validation helpers (time presets, chunk size, filters)
# --------------------------------------------------------------------------

def format_utc(t):
    """``YYYY-MM-DDTHH:MM:SS`` for a UTCDateTime."""
    return UTCDateTime(t).strftime(TIME_FMT)


def _whole_seconds(t):
    t = UTCDateTime(t)
    return UTCDateTime(int(t.timestamp // 1))


def resolve_time_preset(name, now):
    """``(start, end)`` for a relative preset, ``None`` for Entire/Custom."""
    name = (name or "").strip()
    if name in (PRESET_ENTIRE, PRESET_CUSTOM):
        return None
    now = _whole_seconds(now)
    if name in _RELATIVE_PRESETS:
        return now - _RELATIVE_PRESETS[name], now
    if name == "Yesterday (UTC)":
        today = UTCDateTime(now.date)
        return today - 86400, today
    raise ValueError(f"Unknown time window preset '{name}'. "
                     f"Choose one of: {', '.join(TIME_PRESETS)}.")


def suggest_times(now, epoch_start=None, days=7):
    """Suggested Start/End strings: now, today 00:00, previous midnights, epoch start."""
    now = _whole_seconds(now)
    today = UTCDateTime(now.date)
    cands = [now, today] + [today - i * 86400 for i in range(1, days + 1)]
    if epoch_start is not None and UTCDateTime(epoch_start) < now:
        cands.append(UTCDateTime(epoch_start))
    out = []
    for t in cands:
        s = format_utc(t)
        if s not in out:
            out.append(s)
    return out


def parse_chunk_size(text):
    """Chunk size in seconds from a preset label, ``N`` or ``Ns``. Raises ValueError."""
    raw = "" if text is None else str(text)
    s = raw.strip().lower()
    for label, secs in CHUNK_PRESETS.items():
        if s == label:
            return secs
    m = re.fullmatch(r"(\d+)\s*s?", s)
    value = int(m.group(1)) if m else None
    if value is None or not MIN_CHUNK_SECONDS <= value <= MAX_CHUNK_SECONDS:
        raise ValueError(f"Invalid chunk size '{raw}'. Choose "
                         f"{', '.join(list(CHUNK_PRESETS)[:-1])} or {list(CHUNK_PRESETS)[-1]}"
                         f" or enter seconds ({MIN_CHUNK_SECONDS}–{MAX_CHUNK_SECONDS}).")
    return value


_CHANNEL_PATTERN = re.compile(r"[A-Z0-9?*\[\]]+")


def filter_channels(rows, pattern):
    """Keep rows whose channel matches comma-separated wildcard patterns.

    ``""``, ``all`` and ``*`` keep everything. Matching is case-insensitive.
    Raises ValueError for characters other than letters, digits, ``? * [ ]``.
    """
    rows = list(rows)
    text = (pattern or "").strip()
    if text.lower() in ("", "all", "*"):
        return rows
    pats = [p.strip().upper() for p in text.split(",")]
    if any(not _CHANNEL_PATTERN.fullmatch(p) for p in pats):
        raise ValueError(f"Invalid channel filter '{pattern}'. Use codes or wildcards "
                         "such as HHZ, HH?, BH? (comma-separated).")
    return [r for r in rows if any(fnmatch.fnmatchcase(r.channel.upper(), p) for p in pats)]


def _dedupe(items):
    out = []
    for i in items:
        if i not in out:
            out.append(i)
    return out


def channel_choices(rows):
    """Channel filter dropdown values: presets, band groups (``HH?``), exact codes."""
    codes = sorted({r.channel for r in rows})
    groups = sorted({c[:2] + "?" for c in codes if len(c) >= 2})
    return _dedupe(list(CHANNEL_FILTER_PRESETS) + groups + codes)


def code_choices(rows, attr, presets=("*",)):
    """Dropdown values for station/location codes (blank location shown as ``--``)."""
    codes = sorted({getattr(r, attr) or "--" for r in rows})
    return _dedupe(list(presets) + codes)


# --------------------------------------------------------------------------
# Time handling, job planning and size estimates
# --------------------------------------------------------------------------

def parse_utc(text, label="time"):
    """Parse a UTC time string (or pass through a UTCDateTime). Raises ValueError."""
    if isinstance(text, UTCDateTime):
        return text
    s = (text or "").strip() if isinstance(text, str) else text
    if s is None or s == "":
        raise ValueError(f"Invalid {label}: '{text}'. Use YYYY-MM-DDTHH:MM:SS (UTC).")
    try:
        return UTCDateTime(s)
    except Exception:
        raise ValueError(f"Invalid {label}: '{text}'. Use YYYY-MM-DDTHH:MM:SS (UTC).") from None


def validate_time_window(start, end, now=None):
    """Validate a user time window; returns (start, end) as UTCDateTime.

    Raises ValueError if either time is invalid, end <= start, or the
    window starts in the future. No clipping is done here.
    """
    now = UTCDateTime() if now is None else now
    s = parse_utc(start, "start time")
    e = parse_utc(end, "end time")
    if e <= s:
        raise ValueError(f"End time {e} must be after start time {s}.")
    if s >= now:
        raise ValueError(f"Start time {s} is in the future.")
    return s, e


def entire_range(row, now=None):
    """Full span of a channel epoch, with open/future ends clipped to now."""
    now = UTCDateTime() if now is None else now
    if row.start is None:
        raise ValueError(f"{row.seed_id}: epoch has no start date.")
    end = row.end if row.end is not None else now
    end = min(end, now)
    if end <= row.start:
        raise ValueError(f"{row.seed_id}: epoch starts in the future.")
    return row.start, end


def build_jobs(rows, mode, start=None, end=None, now=None):
    """Turn selected rows into DownloadJobs.

    ``mode`` is ``"entire"`` (each row's whole epoch) or ``"window"`` (the
    given start/end clipped to each epoch). Returns ``(jobs, skipped)`` where
    ``skipped`` is a list of ``(row, reason)``.
    """
    if mode not in ("entire", "window"):
        raise ValueError(f"Unknown mode '{mode}'.")
    rows = list(rows)
    if not rows:
        raise ValueError("No channels selected")
    now = UTCDateTime() if now is None else now
    jobs, skipped = [], []
    if mode == "entire":
        for r in rows:
            try:
                s, e = entire_range(r, now)
            except ValueError as exc:
                skipped.append((r, str(exc)))
                continue
            jobs.append(DownloadJob(r, s, e))
        return jobs, skipped

    ws, we = validate_time_window(start, end, now)
    for r in rows:
        s = max(ws, r.start) if r.start is not None else ws
        e = min(we, r.end if r.end is not None else now, now)
        if e <= s:
            skipped.append((r, "time window does not overlap the channel epoch"))
            continue
        jobs.append(DownloadJob(r, s, e))
    return jobs, skipped


def build_array_jobs(rows, mode, start=None, end=None, now=None, channel_filter="all"):
    """build_jobs() for a whole array after applying a channel filter."""
    filtered = filter_channels(rows, channel_filter)
    if not filtered:
        raise ValueError(f"No channels match filter '{channel_filter}'.")
    return build_jobs(filtered, mode, start, end, now)


def job_out_dir(out_dir, row, layout="flat"):
    """Output folder for a row: ``out`` (flat) or ``out/NET/STA`` (array)."""
    if layout == "flat":
        return Path(out_dir)
    if layout == "array":
        return Path(out_dir) / row.network / row.station
    raise ValueError(f"Unknown layout '{layout}'.")


def split_into_chunks(start, end, chunk_seconds=DEFAULT_CHUNK_SECONDS, align=True):
    """Split [start, end] into contiguous chunks of at most ``chunk_seconds``.

    With ``align`` the inner boundaries fall on multiples of ``chunk_seconds``
    since the epoch (UTC midnight for one day).
    """
    if chunk_seconds <= 0:
        raise ValueError("chunk_seconds must be positive.")
    start, end = UTCDateTime(start), UTCDateTime(end)
    if end <= start:
        raise ValueError("end must be after start.")
    chunks = []
    s = start
    while s < end:
        if align:
            nxt = UTCDateTime((int(s.timestamp // chunk_seconds) + 1) * chunk_seconds)
        else:
            nxt = s + chunk_seconds
        e = min(nxt, end)
        chunks.append((s, e))
        s = e
    return chunks


def count_chunks(start, end, chunk_seconds=DEFAULT_CHUNK_SECONDS, align=True):
    """``len(split_into_chunks(...))`` in O(1); same validation and errors."""
    if chunk_seconds <= 0:
        raise ValueError("chunk_seconds must be positive.")
    start, end = UTCDateTime(start), UTCDateTime(end)
    if end <= start:
        raise ValueError("end must be after start.")
    c_ns = int(round(chunk_seconds * 1_000_000_000))
    s_ns, e_ns = start.ns, end.ns
    if align:
        return -(-e_ns // c_ns) - (s_ns // c_ns)
    return -(-(e_ns - s_ns) // c_ns)


def estimate_bytes(jobs, chunk_seconds=DEFAULT_CHUNK_SECONDS):
    """Rough size estimate for a list of DownloadJobs."""
    seconds = 0.0
    raw = 0.0
    chunks = 0
    for j in jobs:
        span = j.end - j.start
        seconds += span
        raw += (j.row.sample_rate or 0.0) * span * BYTES_PER_SAMPLE
        chunks += count_chunks(j.start, j.end, chunk_seconds)
    return {"seconds": seconds, "raw_bytes": int(raw),
            "mseed_bytes": int(raw / MSEED_COMPRESSION_RATIO), "chunks": chunks}


def format_bytes(n):
    """Human-readable byte count, e.g. ``1.2 GB``."""
    n = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024


def chunk_filename(row, s, e, ext=".mseed"):
    """Windows-safe file name for one chunk of a channel."""
    fmt = "%Y%m%dT%H%M%S"
    def stamp(value):
        fraction = value.ns % 1000000000
        return value.strftime(fmt) + (f".{fraction:09d}" if fraction else "")
    return (f"{row.network}_{row.station}_{row.location or '--'}_{row.channel}_"
            f"{stamp(s)}_{stamp(e)}{ext}")


# --------------------------------------------------------------------------
# Chunked download
# --------------------------------------------------------------------------

def download_waveforms(client, network, station, location, channel, starttime, endtime,
                       out_dir, chunk_seconds=DEFAULT_CHUNK_SECONDS, convert_to_h5=False,
                       cancel_event=None, progress=None, row=None, resume=False,
                       failure_counter=None, metadata_root=None, metadata_cache=None,
                       provider="", array_geometry=None):
    """Download one channel in chunks, one MiniSEED file per chunk.

    Chunks with no data are recorded in ``skipped``; other per-chunk failures
    are recorded in ``errors`` and the download continues. ``cancel_event``
    is checked before each chunk. ``progress(done, total, message)`` is
    called after each chunk. Raises ValueError for invalid arguments.

    With ``resume`` stray ``.part`` files of this channel are deleted first,
    and chunks whose final ``.mseed`` (and ``.h5`` when converting) already
    exist are skipped without a request (recorded in ``existing``); a missing
    ``.h5`` next to an existing ``.mseed`` is only converted. A
    ``failure_counter`` (FailureCounter) counts consecutive request errors
    and stops the download (``aborted``) when its limit is reached.
    """
    chunks = split_into_chunks(starttime, endtime, chunk_seconds)
    if row is None:
        row = ChannelRow(network, station, location, channel, None, None, None, None, None)
    os.makedirs(out_dir, exist_ok=True)
    out_dir = Path(out_dir)
    result = DownloadResult(seed_id=row.seed_id, chunks_total=len(chunks))
    metadata_root = Path(metadata_root or out_dir).resolve()
    metadata_cache = metadata_cache if metadata_cache is not None else {}
    metadata_error = None

    def convert(final, target):
        nonlocal metadata_error
        if metadata_error is not None:
            raise ValueError(metadata_error)
        try:
            xml, relative = handoff.save_metadata(client, row, UTCDateTime(starttime),
                                                 UTCDateTime(endtime), metadata_root, metadata_cache)
        except Exception as exc:
            metadata_error = f"Full StationXML request failed for {row.seed_id}: {exc}"
            raise ValueError(metadata_error) from exc
        return convert_mseed_to_h5_strict(str(final), str(target), stationxml=xml,
             response_xml=relative, strict=True, provider=provider,
             request=dict(seed_id=row.seed_id, requested_start_utc=handoff.utc9(starttime),
                          requested_end_utc=handoff.utc9(endtime), chunk_start_utc=handoff.utc9(s),
                          chunk_end_utc=handoff.utc9(e)), array_geometry=array_geometry)


    if resume:
        prefix = f"{row.network}_{row.station}_{row.location or '--'}_{row.channel}_"
        for stray in list(out_dir.glob(glob_escape(prefix) + "*.mseed.part")) + list(out_dir.glob(glob_escape(prefix) + "*.h5.part")):
            try:
                os.remove(stray)
            except OSError:
                pass

    for s, e in chunks:
        if cancel_event is not None and cancel_event.is_set():
            result.cancelled = True
            break
        msg = f"{row.seed_id} {s} - {e}: "
        if resume:
            final = out_dir / chunk_filename(row, s, e)
            if final.exists():
                result.existing.append((s, e))
                h5 = final.with_suffix(".h5")
                if not convert_to_h5 or handoff.h5_is_complete(h5):
                    msg += "already downloaded, skipped"
                else:
                    try:
                        convert(final, h5)
                        result.h5_files.append(h5)
                        msg += f"already downloaded, converted {h5.name}"
                    except Exception as exc:
                        result.errors.append((s, e, f"HDF5 conversion failed: {exc}"))
                        msg += f"already downloaded (HDF5 conversion failed: {exc})"
                result.chunks_done += 1
                if progress is not None:
                    progress(result.chunks_done, result.chunks_total, msg)
                continue
        raw_mseed = None
        try:
            if convert_to_h5:
                buf = BytesIO()
                st = client.get_waveforms(network, station, location, channel, s, e, filename=buf)
                raw_mseed = buf.getvalue() or None
                if raw_mseed:
                    st = read(BytesIO(raw_mseed), format="MSEED")
            else:
                st = client.get_waveforms(network, station, location, channel, s, e)
        except FDSNNoDataException:
            st = None
        except Exception as exc:
            st = False
            result.errors.append((s, e, str(exc)))
            msg += f"error: {exc}"
            if failure_counter is not None and failure_counter.record_error():
                result.aborted = (f"Stopped after {failure_counter.limit} consecutive "
                                  f"request failures (last: {exc})")
                result.chunks_done += 1
                if progress is not None:
                    progress(result.chunks_done, result.chunks_total, msg)
                break
        if st is not False and failure_counter is not None:
            failure_counter.reset()
        if st is None or (st is not False and len(st) == 0):
            result.skipped.append((s, e))
            msg += "no data, skipped"
        elif st is not False:
            final = out_dir / chunk_filename(row, s, e)
            part = final.with_name(final.name + ".part")
            try:
                if raw_mseed is not None:
                    part.write_bytes(raw_mseed)
                else:
                    st.write(str(part), format="MSEED")
                os.replace(part, final)
                result.files.append(final)
                msg += f"saved {final.name}"
            except Exception as exc:
                result.errors.append((s, e, f"write failed: {exc}"))
                msg += f"write failed: {exc}"
                try:
                    os.remove(part)
                except OSError:
                    pass
                final = None
            if final is not None and convert_to_h5:
                h5 = final.with_suffix(".h5")
                try:
                    convert(final, h5)
                    result.h5_files.append(h5)
                    msg += f" + {h5.name}"
                except Exception as exc:
                    result.errors.append((s, e, f"HDF5 conversion failed: {exc}"))
                    msg += f" (HDF5 conversion failed: {exc})"
        result.chunks_done += 1
        if progress is not None:
            progress(result.chunks_done, result.chunks_total, msg)
    return result


def glob_escape(text):
    """Escape glob metacharacters (``[`` ``]`` ``*`` ``?``) in a literal prefix."""
    return re.sub(r"([\[\]*?])", r"[\1]", text)


def run_downloads(jobs, out_dir, client, chunk_seconds=DEFAULT_CHUNK_SECONDS,
                  convert_to_h5=False, cancel_event=None, progress=None, layout="flat",
                  resume=False, max_consecutive_errors=None, progress_detail=None, provider="",
                  write_manifest=True, array_geometry_override=None):

    """Run DownloadJobs sequentially with one client; stops after a cancel.

    ``progress(done_all, total_all, message)`` reports overall progress.
    ``progress_detail(DownloadProgress, message)`` additionally reports
    station/channel/chunk counts after every chunk and every job.
    ``layout`` is ``"flat"`` (all files in ``out_dir``) or ``"array"``
    (``out_dir/NET/STA``). ``max_consecutive_errors`` stops the whole run
    after that many consecutive request failures (across jobs).
    """
    jobs = list(jobs)
    total_all = sum(count_chunks(j.start, j.end, chunk_seconds) for j in jobs)
    counter = FailureCounter(max_consecutive_errors) if max_consecutive_errors else None
    station_jobs = {}
    for j in jobs:
        key = (j.row.network, j.row.station)
        station_jobs[key] = station_jobs.get(key, 0) + 1
    stations_total = len(station_jobs)
    remaining = dict(station_jobs)
    state = {"chunks": 0, "channels": 0, "stations": 0}

    def detail(seed_id, msg):
        if progress_detail is not None:
            progress_detail(DownloadProgress(state["chunks"], total_all, state["channels"],
                                             len(jobs), state["stations"], stations_total,
                                             seed_id), msg)

    results = []
    offset = 0
    metadata_cache = {}
    array_geometry = (array_geometry_override or handoff.geometry([j.row for j in jobs])) if convert_to_h5 else None
    for j in jobs:
        if cancel_event is not None and cancel_event.is_set():
            break
        r = j.row

        def cb(done, _total, msg, _off=offset, _sid=r.seed_id):
            state["chunks"] = _off + done
            if progress is not None:
                progress(_off + done, total_all, msg)
            detail(_sid, msg)

        target = job_out_dir(out_dir, r, layout)
        res = download_waveforms(client, r.network, r.station, r.location, r.channel,
                                 j.start, j.end, target, chunk_seconds=chunk_seconds,
                                 convert_to_h5=convert_to_h5, cancel_event=cancel_event,
                                 progress=cb, row=r, resume=resume,
                                 failure_counter=counter, metadata_root=Path(out_dir).resolve(),
                                 metadata_cache=metadata_cache, provider=provider,
                                 array_geometry=array_geometry)
        results.append(res)
        offset += res.chunks_total
        if res.cancelled or res.aborted:
            break
        state["channels"] += 1
        key = (r.network, r.station)
        remaining[key] -= 1
        if remaining[key] == 0:
            state["stations"] += 1
        detail(r.seed_id, f"{r.seed_id} finished")
    if convert_to_h5 and write_manifest:
        handoff.write_run_manifest(out_dir, jobs, results, provider=provider, layout=layout,
                                  chunk_seconds=chunk_seconds, resume=resume)
    return results


def summarize_results(results, jobs, chunk_seconds=DEFAULT_CHUNK_SECONDS,
                      cancel_requested=False):
    """Final state (Done / Cancelled / Stopped), counts and a status line."""
    total = sum(count_chunks(j.start, j.end, chunk_seconds) for j in jobs)
    done = sum(r.chunks_done for r in results)
    aborted = next((r.aborted for r in results if r.aborted), None)
    if aborted:
        state = "Stopped"
    elif any(r.cancelled for r in results) or (cancel_requested and done < total):
        state = "Cancelled"
    else:
        state = "Done"
    info = {
        "state": state,
        "files": sum(len(r.files) for r in results),
        "h5_files": sum(len(r.h5_files) for r in results),
        "existing": sum(len(r.existing) for r in results),
        "skipped": sum(len(r.skipped) for r in results),
        "errors": sum(len(r.errors) for r in results),
        "chunks_done": done,
        "chunks_total": total,
        "aborted": aborted,
    }
    text = f"{state}: {info['files']} file(s) written"
    if info["existing"]:
        text += f", {info['existing']} already present"
    if info["errors"]:
        text += f", {info['errors']} error(s)"
    info["text"] = text + "."
    return info


# 示例调用
if __name__ == "__main__":
    from h5data_gui import main
    main()
