"""Simple Tk GUI for browsing FDSN arrays/channels and downloading MiniSEED/HDF5.

Run:  python h5data_gui.py

Two tabs:
  * Arrays   - list networks ("arrays", real or virtual such as _GSN), preview
               the seismometers of the selected array (table + map) and
               download the whole array in one operation.
  * Channels - the original channel-level Browse / select / download workflow.

All network work runs on worker threads. Workers never touch widgets; they
post messages to a queue that the Tk main loop drains every 100 ms.
"""
import math
import os
import queue
import sys
import threading
import subprocess
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk
from tkinter import font as tkfont

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h5data  # noqa: E402
from obspy import UTCDateTime  # noqa: E402

MAX_TREE_ROWS = 5000
MAX_MAP_LABELS = 200
POLL_MS = 100
PREVIEW_DEBOUNCE_MS = 300
MAX_CONSECUTIVE_ERRORS = 5
DOWNLOAD_TIMEOUT = 60
TIME_FMT = h5data.TIME_FMT
COLUMNS = (
    ("network", "Network", 70), ("station", "Station", 80), ("location", "Location", 70),
    ("channel", "Channel", 70), ("rate", "Rate (Hz)", 80), ("lat", "Lat", 80),
    ("lon", "Lon", 80), ("start", "Start", 95), ("end", "End", 95),
)
ARRAY_COLUMNS = (
    ("code", "Code", 70), ("description", "Description", 330), ("start", "Start", 95),
    ("end", "End", 95), ("stations", "Stations", 70),
)
PREVIEW_COLUMNS = (
    ("station", "Station", 65), ("location", "Loc", 40), ("channel", "Channel", 60),
    ("lat", "Lat", 70), ("lon", "Lon", 70), ("elev", "Elev (m)", 65),
    ("rate", "Rate (Hz)", 65), ("sensor", "Sensor", 160), ("start", "Start", 85),
    ("end", "End", 85),
)
ARRAY_NET_PRESETS = ("*", "CH", "_GSN", "XG")
LOCATION_PRESETS = ("*", "--", "00", "01", "10")


def _fmt_date(t):
    return t.strftime("%Y-%m-%d") if t is not None else ""


def _fmt_num(v, fmt):
    return format(v, fmt) if v is not None else ""


def _fmt_span(seconds):
    days, rem = divmod(int(seconds), 86400)
    hours, rem = divmod(rem, 3600)
    minutes = rem // 60
    if days:
        return f"{days} d {hours} h"
    if hours:
        return f"{hours} h {minutes} min"
    return f"{minutes} min {int(seconds) % 60} s"


def _dedupe(items):
    out = []
    for i in items:
        if i not in out:
            out.append(i)
    return out


def project_points(points, width, height, pad=20):
    """Equirectangular bounding-box projection to canvas pixels.

    ``points`` are ``(net, sta, lat, lon, elev)``; returns ``[(sta, x, y)]``.
    A zero lat/lon span is padded by 0.5 degrees, so one station is centered.
    """
    if not points:
        return []
    lats = [p[2] for p in points]
    lons = [p[3] for p in points]
    lat0, lat1 = min(lats), max(lats)
    lon0, lon1 = min(lons), max(lons)
    if lat1 - lat0 == 0:
        lat0, lat1 = lat0 - 0.5, lat1 + 0.5
    if lon1 - lon0 == 0:
        lon0, lon1 = lon0 - 0.5, lon1 + 0.5
    w = max(width - 2 * pad, 1)
    h = max(height - 2 * pad, 1)
    return [(p[1], pad + (p[3] - lon0) / (lon1 - lon0) * w,
             pad + (lat1 - p[2]) / (lat1 - lat0) * h) for p in points]


def metric_points(points):
    """Local east/north in m, circular centre avoids dateline wrap."""
    if not points:
        return []
    from pyproj import CRS, Transformer
    lats = [p[2] for p in points]; lons = [p[3] for p in points]
    lon0 = math.degrees(math.atan2(sum(math.sin(math.radians(v)) for v in lons),
                                  sum(math.cos(math.radians(v)) for v in lons)))
    local = CRS.from_proj4(f"+proj=aeqd +lat_0={sum(lats)/len(lats)} +lon_0={lon0} +datum=WGS84 +units=m")
    x, y = Transformer.from_crs(4326, local, always_xy=True).transform(lons, lats)
    return [(p[0], p[1], float(north), float(east), p[4]) for p, east, north in zip(points, x, y)]


class StationMap(ttk.Frame):
    """Lon/lat scatter of stations (matplotlib if available, else a tk.Canvas)."""

    def __init__(self, master, use_matplotlib=None):
        super().__init__(master)
        self.points = []
        self.title = ""
        self.backend_name = "canvas"
        if use_matplotlib is not False:
            try:
                from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
                from matplotlib.figure import Figure
            except ImportError:
                if use_matplotlib:
                    raise
            else:
                self.backend_name = "matplotlib"
                self.figure = Figure(figsize=(3.8, 3.0), dpi=90)
                self.figure.subplots_adjust(left=0.18, right=0.96, top=0.9, bottom=0.15)
                self.ax = self.figure.add_subplot(111)
                self.mpl_canvas = FigureCanvasTkAgg(self.figure, master=self)
                self.mpl_canvas.get_tk_widget().pack(fill="both", expand=True)
        if self.backend_name == "canvas":
            self.canvas = tk.Canvas(self, width=340, height=270, background="white",
                                    highlightthickness=0)
            self.canvas.pack(fill="both", expand=True)
            self.canvas.bind("<Configure>", lambda _e: self._draw_canvas())
        self._draw()

    def plot(self, points, title=""):
        """Draw ``(net, sta, lat, lon, elev)`` points."""
        self.points = list(points)
        self.title = title
        self._draw()

    def _draw(self):
        if self.backend_name == "matplotlib":
            self._draw_mpl()
        else:
            self._draw_canvas()

    def _draw_mpl(self):
        ax = self.ax
        ax.clear()
        if not self.points:
            ax.text(0.5, 0.5, "No station coordinates", ha="center", va="center",
                    transform=ax.transAxes, color="gray")
            ax.set_xticks([])
            ax.set_yticks([])
        else:
            points = metric_points(self.points)
            x = [p[3] for p in points]
            y = [p[2] for p in points]
            ax.scatter(x, y, marker="^", s=40, color="tab:red", edgecolors="black",
                       linewidths=0.5, zorder=3)
            for i, p in enumerate(self.points[:MAX_MAP_LABELS]):
                ax.annotate(p[1], (x[i], y[i]), xytext=(3, 3), textcoords="offset points",
                            fontsize=7)
            ax.set_aspect("equal", adjustable="datalim")
            ax.margins(0.15)  # keep station labels inside the axes
            ax.set_xlabel("East (m)")
            ax.set_ylabel("North (m)")
            ax.grid(True, linewidth=0.3)
            ax.tick_params(labelsize=7)
        if self.title:
            ax.set_title(self.title, fontsize=9)
        self.mpl_canvas.draw_idle()

    def _draw_canvas(self):
        c = self.canvas
        c.delete("all")
        w = c.winfo_width() if c.winfo_width() > 1 else int(c["width"])
        h = c.winfo_height() if c.winfo_height() > 1 else int(c["height"])
        if not self.points:
            c.create_text(w / 2, h / 2, text="No station coordinates", fill="gray")
            return
        pad = 30
        c.create_rectangle(pad, pad, w - pad, h - pad, outline="#999")
        points = metric_points(self.points)
        lats = [p[2] for p in points]
        lons = [p[3] for p in points]
        c.create_text(pad, h - pad + 4, anchor="nw", text=f"{min(lons):.0f} m", fill="gray")
        c.create_text(w - pad, h - pad + 4, anchor="ne", text=f"{max(lons):.0f} m",
                      fill="gray")
        c.create_text(pad - 4, h - pad, anchor="se", text=f"{min(lats):.0f} m", fill="gray")
        c.create_text(pad - 4, pad, anchor="ne", text=f"{max(lats):.0f} m", fill="gray")
        if self.title:
            c.create_text(w / 2, 4, anchor="n", text=self.title)
        for i, (sta, x, y) in enumerate(project_points(points, w, h, pad)):
            c.create_polygon(x, y - 6, x - 5, y + 4, x + 5, y + 4, fill="red",
                             outline="black", tags=("station",))
            if i < MAX_MAP_LABELS:
                c.create_text(x + 6, y - 6, anchor="w", text=sta, font=("TkDefaultFont", 7),
                              tags=("label",))


class H5DataApp(ttk.Frame):
    """Browse arrays/channels from an FDSN provider and download them in chunks."""

    def __init__(self, master, backend=h5data, client_factory=None, ask_confirm=None,
                 show_error=None, now=None, use_matplotlib=None):
        super().__init__(master, padding=8)
        self.master = master
        self.backend = backend
        self._factory_injected = client_factory is not None
        self.client_factory = client_factory or backend.make_client
        self._ask_confirm = ask_confirm or messagebox.askyesno
        self._show_error = show_error or messagebox.showerror
        self.now = now or UTCDateTime
        self._use_mpl = use_matplotlib
        self.rows = []
        self.networks = []
        self.array_rows = []
        self.array_code = None
        self._busy = False
        self._queue = queue.Queue()
        self._cancel_event = threading.Event()
        self._preview_token = 0
        self._browse_token = 0
        self._network_token = 0
        self._preview_loading = False
        self._preview_after = None
        self._arraysta_after = None
        self._requested_array_code = None
        self._filling = False
        self._last_detail = ""
        self._geometry_tokens = {"array": 0, "channels": 0}
        self._geometry_keys = {}

        today = UTCDateTime(self.now().date)
        start = today - 86400
        self.provider_var = tk.StringVar(value=backend.DEFAULT_PROVIDER)
        self.network_var = tk.StringVar(value="CH")
        self.station_var = tk.StringVar(value="*")
        self.location_var = tk.StringVar(value="*")
        self.channel_var = tk.StringVar(value="HHZ")
        self.mode_var = tk.StringVar(value="window")
        self.window_var = tk.StringVar(value=backend.PRESET_CUSTOM)
        self.start_var = tk.StringVar(value=start.strftime(TIME_FMT))
        self.end_var = tk.StringVar(value=(start + 3600).strftime(TIME_FMT))
        self.chunk_var = tk.StringVar(value="1 day")
        self.outdir_var = tk.StringVar(value=os.path.join(os.getcwd(), "downloads"))
        self.convert_var = tk.BooleanVar(value=True)
        self.resume_var = tk.BooleanVar(value=True)
        self.status_var = tk.StringVar(value="Ready.")
        self.count_var = tk.StringVar(value="No channels loaded. Click Browse.")
        self.array_net_var = tk.StringVar(value="*")
        self.array_sta_var = tk.StringVar(value="*")
        self.array_chan_var = tk.StringVar(value="HHZ")
        self.arrays_count_var = tk.StringVar(value="No arrays loaded. Click List arrays.")
        self.array_summary_var = tk.StringVar(value="Select an array to preview it.")
        self.spacing_var = tk.StringVar(value="Mean station spacing (m): select an array")
        self.channel_spacing_var = tk.StringVar(value="Select channels to preview the download array")
        self.station_prefix_var = tk.StringVar(value='')
        self.geo_vars = {key:tk.StringVar(value='') for key in
                         ('minlatitude','maxlatitude','minlongitude','maxlongitude')}

        self._build()
        self._on_mode_change()
        self.mode_var.trace_add("write", self._on_mode_var)
        self.start_var.trace_add("write", self._on_time_edit)
        self.end_var.trace_add("write", self._on_time_edit)
        self.start_var.trace_add("write", self._update_channel_preview)
        self.end_var.trace_add("write", self._update_channel_preview)
        self.mode_var.trace_add("write", self._update_channel_preview)
        self.array_sta_var.trace_add("write", self._on_array_station_edit)
        for var in [self.station_prefix_var,*self.geo_vars.values()]:
            var.trace_add('write',self._on_station_filters_edit)
        self.provider_var.trace_add("write", self._on_provider_edit)
        for var in (self.window_var, self.start_var, self.end_var, self.chunk_var,
                    self.array_chan_var):
            var.trace_add("write", self._update_array_summary)
        try:
            master.protocol("WM_DELETE_WINDOW", self._on_close)
        except (AttributeError, tk.TclError):
            pass
        self._poll_id = self.after(POLL_MS, self._poll_queue)

    # ------------------------------------------------------------------ UI
    def _build(self):
        # Row height follows the font, so rows are not clipped on high-DPI screens.
        try:
            line = tkfont.nametofont("TkDefaultFont").metrics("linespace")
            ttk.Style(self).configure("Treeview", rowheight=line + 4)
        except tk.TclError:
            pass
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=3)
        self.rowconfigure(3, weight=1)

        src = ttk.Frame(self)
        src.grid(row=0, column=0, sticky="ew")
        provider_line=ttk.Frame(src);provider_line.pack(fill='x')
        ttk.Label(provider_line, text="Provider").pack(side="left")
        self.provider_box = ttk.Combobox(provider_line, textvariable=self.provider_var, width=14,
                                         values=self.backend.COMMON_PROVIDERS)
        self.provider_box.pack(side="left", padx=(4, 8))
        ttk.Label(provider_line, text="FDSN short name or base URL", foreground="gray"
                  ).pack(side="left")
        filters = ttk.LabelFrame(src, text='Station filters', padding=5)
        filters.pack(side='top',pady=(5,0),fill='x',expand=True)
        fields=[('Station prefix',self.station_prefix_var),('Latitude min',self.geo_vars['minlatitude']),
                ('Latitude max',self.geo_vars['maxlatitude']),('Longitude west',self.geo_vars['minlongitude']),
                ('Longitude east',self.geo_vars['maxlongitude'])]
        self.filter_entries=[]
        for col,(label,var) in enumerate(fields):
            ttk.Label(filters,text=label).grid(row=0,column=col,sticky='w',padx=3)
            entry=ttk.Entry(filters,textvariable=var,width=12)
            entry.grid(row=1,column=col,padx=3,sticky='ew');self.filter_entries.append(entry)
            entry.bind('<Return>',lambda _event:self.apply_station_filters())
        self.apply_filters_btn=ttk.Button(filters,text='Apply filters',command=self.apply_station_filters)
        self.apply_filters_btn.grid(row=1,column=5,padx=4)
        self.clear_filters_btn=ttk.Button(filters,text='Clear filters',command=self.clear_station_filters)
        self.clear_filters_btn.grid(row=1,column=6,padx=4)
        ttk.Label(filters,text='Prefix: QS or QS,QT. Blank bounds = unrestricted. West > east = cross dateline.',
                  foreground='gray').grid(row=2,column=0,columnspan=7,sticky='w')

        self.notebook = ttk.Notebook(self)
        self.notebook.grid(row=1, column=0, sticky="nsew", pady=(6, 0))
        self.arrays_tab = ttk.Frame(self.notebook, padding=6)
        self.channels_tab = ttk.Frame(self.notebook, padding=6)
        self.notebook.add(self.arrays_tab, text="Arrays")
        self.notebook.add(self.channels_tab, text="Channels")
        self._build_arrays_tab(self.arrays_tab)
        self._build_channels_tab(self.channels_tab)
        self._build_options()

        st = ttk.LabelFrame(self, text="Status", padding=6)
        st.grid(row=3, column=0, sticky="nsew", pady=(6, 0))
        st.columnconfigure(0, weight=1)
        st.rowconfigure(2, weight=1)
        self.progress = ttk.Progressbar(st, mode="determinate", maximum=1, value=0)
        self.progress.grid(row=0, column=0, sticky="ew")
        ttk.Label(st, textvariable=self.status_var).grid(row=1, column=0, sticky="w")
        self.log = scrolledtext.ScrolledText(st, height=5, state="disabled", wrap="word")
        self.log.grid(row=2, column=0, sticky="nsew")

    def _tree(self, parent, columns, height, selectmode):
        frame = ttk.Frame(parent)
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        tree = ttk.Treeview(frame, columns=[c[0] for c in columns], show="headings",
                            selectmode=selectmode, height=height)
        # Column widths are given at 96 dpi; scale them on high-DPI screens.
        scale = max(self.winfo_fpixels("1i") / 96.0, 1.0)
        for key, title, width in columns:
            tree.heading(key, text=title)
            tree.column(key, width=int(width * scale), minwidth=int(30 * scale),
                        anchor="w", stretch=key in ("description", "sensor"))
        vsb = ttk.Scrollbar(frame, orient="vertical", command=tree.yview)
        hsb = ttk.Scrollbar(frame,orient='horizontal',command=tree.xview)
        tree.configure(yscrollcommand=vsb.set,xscrollcommand=hsb.set)
        tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1,column=0,sticky='ew')
        return frame, tree

    def _build_arrays_tab(self, tab):
        tab.columnconfigure(0, weight=1)
        tab.rowconfigure(2, weight=1)
        top = ttk.Frame(tab)
        top.grid(row=0, column=0, sticky="ew")
        ttk.Label(top, text="Network / array").grid(row=0, column=0, sticky="w")
        self.array_net_box = ttk.Combobox(top, textvariable=self.array_net_var, width=12,
                                          values=ARRAY_NET_PRESETS)
        self.array_net_box.grid(row=1, column=0, padx=(0, 8))
        ttk.Label(top, text="Stations (sub-array)").grid(row=0, column=1, sticky="w")
        self.array_sta_box = ttk.Combobox(top, textvariable=self.array_sta_var, width=12,
                                          values=("*",))
        self.array_sta_box.grid(row=1, column=1, padx=(0, 8))
        self.arrays_btn = ttk.Button(top, text="List arrays", command=self.list_arrays)
        self.arrays_btn.grid(row=1, column=2)
        ttk.Label(top, text="Wildcards * ? allowed; virtual networks start with _ "
                  "(e.g. _GSN)", foreground="gray").grid(row=1, column=3, padx=8, sticky="w")

        frame, self.arrays_tree = self._tree(tab, ARRAY_COLUMNS, 5, "browse")
        frame.grid(row=1, column=0, sticky="nsew", pady=(6, 0))
        ttk.Label(frame, textvariable=self.arrays_count_var).grid(row=2, column=0,
                                                                  sticky="w")
        self.arrays_tree.bind("<<TreeviewSelect>>", self._on_array_select)

        pane = ttk.PanedWindow(tab, orient="horizontal")
        pane.grid(row=2, column=0, sticky="nsew", pady=(6, 0))
        left, self.preview_tree = self._tree(pane, PREVIEW_COLUMNS, 8, "browse")
        pane.add(left, weight=3)
        self.station_map = StationMap(pane, use_matplotlib=self._use_mpl)
        pane.add(self.station_map, weight=2)

        summary = ttk.Frame(tab)
        summary.grid(row=3, column=0, sticky="ew", pady=(4, 0))
        ttk.Label(summary, textvariable=self.array_summary_var, justify="left"
                  ).pack(anchor="w")
        ttk.Label(summary, textvariable=self.spacing_var, foreground="#125b73",
                  justify="left").pack(anchor="w", pady=(2, 0))
        bottom = ttk.Frame(tab)
        bottom.grid(row=4, column=0, sticky="ew", pady=(4, 0))
        ttk.Label(bottom, text="Channels").pack(side="left")
        self.array_chan_box = ttk.Combobox(bottom, textvariable=self.array_chan_var,
                                           width=14,
                                           values=self.backend.channel_choices([]))
        self.array_chan_box.pack(side="left", padx=(4, 8))
        ttk.Label(bottom, text="all, or codes/wildcards (HHZ, HH?, HHZ,BHZ)",
                  foreground="gray").pack(side="left")
        self.array_export_btn=ttk.Button(bottom,text='Export stations',
                                         command=lambda:self.export_stations('array'))
        self.array_export_btn.pack(side='right',padx=6)
        self.array_download_btn = ttk.Button(bottom, text="Download matched stations",
                                             command=self.download_array)
        self.array_download_btn.pack(side="right")

    def _build_channels_tab(self, tab):
        tab.columnconfigure(0, weight=1)
        tab.rowconfigure(1, weight=1)
        src = ttk.Frame(tab)
        src.grid(row=0, column=0, sticky="ew")
        boxes = []
        for col, (label, var, values) in enumerate((
                ("Network", self.network_var, ("*", "CH")),
                ("Station", self.station_var, ("*",)),
                ("Location", self.location_var, LOCATION_PRESETS),
                ("Channel", self.channel_var, self._chan_query_choices([])))):
            ttk.Label(src, text=label).grid(row=0, column=col, sticky="w")
            box = ttk.Combobox(src, textvariable=var, width=10, values=values)
            box.grid(row=1, column=col, padx=(0, 8))
            boxes.append(box)
        self.network_box, self.station_box, self.location_box, self.channel_box = boxes
        self.browse_btn = ttk.Button(src, text="Browse", command=self.browse)
        self.browse_btn.grid(row=1, column=4, padx=(4, 0))
        ttk.Label(src, text="Wildcards * ? allowed, e.g. HH? or Z*",
                  foreground="gray").grid(row=1, column=5, padx=8, sticky="w")

        pane = ttk.PanedWindow(tab, orient="horizontal")
        pane.grid(row=1, column=0, sticky="nsew", pady=(6, 0))
        frame, self.tree = self._tree(pane, COLUMNS, 10, "extended")
        pane.add(frame, weight=3)
        self.channel_map = StationMap(pane, use_matplotlib=self._use_mpl)
        pane.add(self.channel_map, weight=2)
        self.tree.bind("<<TreeviewSelect>>", self._update_channel_preview)
        bottom = ttk.Frame(tab)
        bottom.grid(row=2, column=0, sticky="ew")
        ttk.Label(bottom, textvariable=self.count_var).pack(side="left")
        self.download_btn = ttk.Button(bottom, text="Download selected",
                                       command=self.download)
        self.download_btn.pack(side="right", pady=(4, 0))
        self.channel_export_btn=ttk.Button(bottom,text='Export stations',
                                            command=lambda:self.export_stations('channels'))
        self.channel_export_btn.pack(side='right',padx=6)
        ttk.Label(tab, textvariable=self.channel_spacing_var, foreground="#125b73",
                  justify="left").grid(row=3, column=0, sticky="w", pady=(4, 0))

    def _build_options(self):
        dl = ttk.LabelFrame(self, text="Download options", padding=6)
        dl.grid(row=2, column=0, sticky="ew", pady=(6, 0))
        dl.columnconfigure(6, weight=1)
        ttk.Label(dl, text="Time window").grid(row=0, column=0, sticky="w")
        self.window_box = ttk.Combobox(dl, textvariable=self.window_var, width=22,
                                       state="readonly", values=self.backend.TIME_PRESETS)
        self.window_box.grid(row=1, column=0, padx=(0, 8))
        self.window_box.bind("<<ComboboxSelected>>", self._on_preset_change)
        suggestions = self.backend.suggest_times(self.now())
        ttk.Label(dl, text="Start (UTC)").grid(row=0, column=1, sticky="w")
        self.start_entry = ttk.Combobox(dl, textvariable=self.start_var, width=21,
                                        values=suggestions)
        self.start_entry.grid(row=1, column=1, padx=(0, 8))
        ttk.Label(dl, text="End (UTC)").grid(row=0, column=2, sticky="w")
        self.end_entry = ttk.Combobox(dl, textvariable=self.end_var, width=21,
                                      values=suggestions)
        self.end_entry.grid(row=1, column=2, padx=(0, 8))
        ttk.Label(dl, text="Chunk size").grid(row=0, column=3, sticky="w")
        self.chunk_box = ttk.Combobox(dl, textvariable=self.chunk_var, width=10,
                                      values=list(self.backend.CHUNK_PRESETS))
        self.chunk_box.grid(row=1, column=3, padx=(0, 8))
        ttk.Label(dl, text="Times: YYYY-MM-DDTHH:MM:SS; chunk: preset or seconds",
                  foreground="gray").grid(row=1, column=4, columnspan=3, sticky="w")

        ttk.Label(dl, text="Output folder").grid(row=2, column=0, sticky="w", pady=(6, 0))
        ttk.Entry(dl, textvariable=self.outdir_var).grid(row=2, column=1, columnspan=5,
                                                         sticky="ew", pady=(6, 0))
        ttk.Button(dl, text="Choose...", command=self._choose_dir).grid(row=2, column=6,
                                                                      sticky="w", padx=4,
                                                                      pady=(6, 0))
        ttk.Checkbutton(dl, text="H5 + full metadata", variable=self.convert_var
                        ).grid(row=3, column=0, sticky="w", pady=(6, 0))
        ttk.Checkbutton(dl, text="Resume: skip files already downloaded",
                        variable=self.resume_var
                        ).grid(row=3, column=1, columnspan=2, sticky="w", pady=(6, 0))
        self.cancel_btn = ttk.Button(dl, text="Cancel download", command=self.cancel,
                                     state="disabled")
        self.cancel_btn.grid(row=3, column=6, sticky="e", pady=(6, 0))
        self.open_dispersion_btn = ttk.Button(dl, text="Open H5 in dispersion",
                                             command=self.open_dispersion)
        self.open_dispersion_btn.grid(row=4, column=0, columnspan=3, sticky="w", pady=(6, 0))
        self.open_hvsr_btn = ttk.Button(dl, text="Open single-station H/V", command=self.open_hvsr)
        self.open_hvsr_btn.grid(row=5,column=0,columnspan=3,sticky="w",pady=(6,0))
        ttk.Label(dl,text="H/V requires Z + two horizontal channels (e.g. HH? / BH?); one station is enough",foreground="gray").grid(row=5,column=3,columnspan=4,sticky="w",padx=8)
        ttk.Label(dl, text="H5 retains raw samples, UTC, coordinates, full StationXML and original MiniSEED",
                  foreground="gray").grid(row=4, column=3, columnspan=4, sticky="w", padx=8)

    def _chan_query_choices(self, rows):
        # "all" is a download filter, not an FDSN query code.
        return [c for c in self.backend.channel_choices(rows) if c != "all"]

    # ------------------------------------------------- time window presets
    def _is_relative(self, name):
        return name in self.backend.TIME_PRESETS and name not in (
            self.backend.PRESET_ENTIRE, self.backend.PRESET_CUSTOM)

    def _set_times(self, start, end):
        self._filling = True
        try:
            self.start_var.set(self.backend.format_utc(start))
            self.end_var.set(self.backend.format_utc(end))
        finally:
            self._filling = False

    def _on_preset_change(self, _event=None):
        name = self.window_var.get()
        self._filling = True
        try:
            if name == self.backend.PRESET_ENTIRE:
                self.mode_var.set("entire")
            else:
                self.mode_var.set("window")
        finally:
            self._filling = False
        if self._is_relative(name):
            self._set_times(*self.backend.resolve_time_preset(name, self.now()))
        self._on_mode_change()
        self._update_array_summary()

    def _on_time_edit(self, *_):
        if not self._filling and self._is_relative(self.window_var.get()):
            self.window_var.set(self.backend.PRESET_CUSTOM)

    def _on_mode_var(self, *_):
        if not self._filling:
            if self.mode_var.get() == "entire":
                self.window_var.set(self.backend.PRESET_ENTIRE)
            elif self.window_var.get() == self.backend.PRESET_ENTIRE:
                self.window_var.set(self.backend.PRESET_CUSTOM)
        self._on_mode_change()

    def _on_mode_change(self):
        state = "normal" if self.mode_var.get() == "window" else "disabled"
        self.start_entry.configure(state=state)
        self.end_entry.configure(state=state)

    def _window_inputs(self, write_back=True):
        """``(mode, start, end)``; relative presets are re-resolved with now()."""
        mode = self.mode_var.get()
        name = self.window_var.get()
        if mode == "window" and self._is_relative(name):
            s, e = self.backend.resolve_time_preset(name, self.now())
            if write_back:
                self._set_times(s, e)
            return mode, s, e
        return mode, self.start_var.get(), self.end_var.get()

    def _chunk_seconds(self):
        return self.backend.parse_chunk_size(self.chunk_var.get())

    def _choose_dir(self):
        d = filedialog.askdirectory(parent=self.master, initialdir=self.outdir_var.get()
                                    or os.getcwd())
        if d:
            self.outdir_var.set(d)

    # ------------------------------------------------------------- helpers
    @property
    def busy(self):
        return self._busy

    def _set_busy(self, busy, downloading=False):
        self._busy = busy
        state = "disabled" if busy else "normal"
        for btn in (self.browse_btn, self.download_btn, self.arrays_btn,
                    self.array_download_btn, self.apply_filters_btn,self.clear_filters_btn,
                    self.array_export_btn,self.channel_export_btn,*self.filter_entries):
            btn.configure(state=state)
        self.cancel_btn.configure(state="normal" if busy and downloading else "disabled")
        self.open_dispersion_btn.configure(state="disabled" if busy else "normal")
        self.open_hvsr_btn.configure(state="disabled" if busy else "normal")

    def _log(self, text):
        self.log.configure(state="normal")
        self.log.insert("end", text + "\n")
        self.log.see("end")
        self.log.configure(state="disabled")

    def _error(self, title, msg):
        self._log(f"ERROR: {title}: {msg}")
        self.status_var.set(f"{title}.")
        self._show_error(title, msg, parent=self.master)

    def selected_rows(self):
        """ChannelRow objects for the selected tree items."""
        return [self.rows[int(i)] for i in self.tree.selection()]

    def _download_client(self, provider):
        if self._factory_injected:
            return self.client_factory(provider)
        return self.backend.make_client(provider, timeout=DOWNLOAD_TIMEOUT)

    def _station_filters(self):
        return self.backend.normalize_station_filters(station_prefix=self.station_prefix_var.get(),
                **{key:var.get() for key,var in self.geo_vars.items()})

    def _on_station_filters_edit(self,*_):
        # Invalidate pending queries as well as cached lists before any download.
        self._preview_token+=1;self._browse_token+=1;self._network_token+=1
        self.array_rows=[];self.rows=[];self._preview_loading=False
        self.preview_tree.delete(*self.preview_tree.get_children())
        self.tree.delete(*self.tree.get_children())
        self.count_var.set('Filters changed. Click Apply / Browse.')
        self.array_summary_var.set('Filters changed. Click Apply to refresh matching stations.')
        self._request_geometry([],'array');self._request_geometry([],'channels')

    def apply_station_filters(self):
        if self._busy:return
        try:self._station_filters()
        except ValueError as exc:self._error('Invalid station filters',str(exc));return
        if self.notebook.index('current')==0:
            if self._requested_array_code:self.preview_array(self._requested_array_code)
            else:self.list_arrays()
        else:self.browse()

    def clear_station_filters(self):
        if self._busy:return
        self.station_prefix_var.set('')
        for var in self.geo_vars.values():var.set('')
        self.apply_station_filters()

    def export_stations(self,surface):
        rows=self.array_rows if surface=='array' else self.rows
        if not rows:self._error('No stations','Query matching stations before exporting.');return
        filename=filedialog.asksaveasfilename(parent=self.master,title='Export all matching stations',
            defaultextension='.csv',filetypes=[('CSV station list','*.csv')],initialfile='matching_stations.csv')
        if filename:
            try:
                records=self.backend.export_station_catalog(rows,filename)
                self._log(f'Exported all {len(records)} matching stations: {filename}')
            except Exception as exc:self._error('Export failed',str(exc))

    # ------------------------------------------------------ actions: arrays
    def list_arrays(self):
        """List networks ("arrays") matching the filters (on a worker)."""
        if self._busy:
            return
        provider = self.provider_var.get().strip()
        net, sta = self.array_net_var.get(), self.array_sta_var.get()
        try:station_options=self._station_filters()
        except ValueError as exc:self._error('Invalid station filters',str(exc));return
        self._network_token+=1;token=self._network_token
        self.arrays_tree.delete(*self.arrays_tree.get_children())
        self.networks = []
        self.arrays_count_var.set("Searching...")
        self._set_busy(True)
        self.status_var.set(f"Listing arrays on {provider}...")
        self._log(f"List arrays {provider}: net={net} sta={sta}")

        def work():
            try:
                client = self.client_factory(provider)
                rows = self.backend.list_networks(provider, net, sta, client=client,**station_options)
                self._queue.put(("networks", token, rows))
            except Exception as e:
                self._queue.put(("error", "List arrays failed", str(e)))
            finally:
                self._queue.put(("idle",))

        threading.Thread(target=work, daemon=True).start()

    def _on_array_select(self, _event=None):
        sel = self.arrays_tree.selection()
        if not sel:
            return
        code = self.networks[int(sel[0])].code
        if self._preview_after is not None:
            self.after_cancel(self._preview_after)
        self._preview_after = self.after(PREVIEW_DEBOUNCE_MS,
                                         lambda: self._debounced_preview(code))

    def _debounced_preview(self, code):
        self._preview_after = None
        self.preview_array(code)

    def preview_array(self, code):
        """Fetch the stations/channels of one array (on a worker) and show them."""
        provider = self.provider_var.get().strip()
        try:station_options=self._station_filters()
        except ValueError as exc:self._error('Invalid station filters',str(exc));return
        self._requested_array_code = code
        sta = self.array_sta_var.get()
        self._preview_token += 1
        token = self._preview_token
        self._preview_loading = True
        self.array_rows = []
        self.array_code = None
        self._request_geometry([], "array")
        self.preview_tree.delete(*self.preview_tree.get_children())
        self.station_map.plot([], title=f"{code} (loading)")
        self.array_summary_var.set(f"Loading {code}...")
        self._log(f"Preview array {code} (stations {sta or '*'}) on {provider}")

        def work():
            try:
                client = self.client_factory(provider)
                rows = self.backend.list_array_channels(provider, code, sta, client=client,**station_options)
                self._queue.put(("preview", token, code, rows))
            except Exception as e:
                self._queue.put(("preview_error", token, code, str(e)))

        threading.Thread(target=work, daemon=True).start()

    def _on_array_station_edit(self, *_):
        code = self._requested_array_code
        if not code:
            return
        if self._arraysta_after is not None:
            self.after_cancel(self._arraysta_after)
        self._preview_token += 1
        self.array_rows = []
        self._preview_loading = True
        self._request_geometry([], "array")
        def refresh():
            self._arraysta_after = None
            self.preview_array(code)
        self._arraysta_after = self.after(PREVIEW_DEBOUNCE_MS, refresh)

    def _on_provider_edit(self, *_):
        self._preview_token += 1
        self._browse_token += 1
        self._network_token += 1
        self._requested_array_code = None
        self.array_code = None
        self.array_rows = []
        self.rows = []
        self.networks = []
        self._preview_loading = False
        for name in ('_arraysta_after','_preview_after'):
            token = getattr(self, name, None)
            if token is not None:
                self.after_cancel(token)
                setattr(self, name, None)
        for tree in (self.arrays_tree,self.preview_tree,self.tree):
            tree.delete(*tree.get_children())
        self._request_geometry([], "array")
        self._request_geometry([], "channels")
        self.array_summary_var.set("Provider changed. List arrays again.")
        self.count_var.set("Provider changed. Click Browse.")

    def download_array(self):
        """Download every (filtered) channel of the previewed array in one operation."""
        if self._busy:
            return
        if not self.array_rows:
            msg = ("Preview still loading." if self._preview_loading
                   else "Select an array and wait for its preview to load.")
            self._show_error("No array", msg, parent=self.master)
            return
        try:
            mode, s, e = self._window_inputs()
            chunk = self._chunk_seconds()
            jobs, skipped = self.backend.build_array_jobs(
                self.array_rows, mode, s, e, now=self.now(),
                channel_filter=self.array_chan_var.get())
        except ValueError as exc:
            self._error("Invalid input", str(exc))
            return
        for row, reason in skipped:
            self._log(f"Skipping {row.seed_id}: {reason}")
        if not jobs:
            self._error("Nothing to download",
                        "No channel of this array has data in the requested range.")
            return
        code = self.array_code
        out = self.outdir_var.get().strip() or os.getcwd()
        n_sta = len({(j.row.network, j.row.station) for j in jobs})
        title = [f"Array: {code}", f"Stations: {n_sta}",
                 f"Output: {os.path.join(out, code)}{os.sep}<STA>"]
        self._start_download(jobs, chunk, "array", title, True, mode)

    # ---------------------------------------------------- actions: channels
    def browse(self):
        """Query the provider for channels matching the filters (on a worker)."""
        if self._busy:
            return
        provider = self.provider_var.get().strip()
        filters = (self.network_var.get(), self.station_var.get(),
                   self.location_var.get(), self.channel_var.get())
        try:station_options=self._station_filters()
        except ValueError as exc:self._error('Invalid station filters',str(exc));return
        self._browse_token+=1;token=self._browse_token
        self.tree.delete(*self.tree.get_children())
        self.rows = []
        self.count_var.set("Searching...")
        self._set_busy(True)
        self.status_var.set(f"Browsing {provider}...")
        self._log(f"Browse {provider}: net={filters[0]} sta={filters[1]} "
                  f"loc={filters[2]} cha={filters[3]}")

        def work():
            try:
                client = self.client_factory(provider)
                rows = self.backend.list_stations(provider, *filters, client=client,**station_options)
                self._queue.put(("rows", token, rows))
            except Exception as e:
                self._queue.put(("error", "Browse failed", str(e)))
            finally:
                self._queue.put(("idle",))

        threading.Thread(target=work, daemon=True).start()

    def download(self):
        """Plan, confirm and start downloading the selected channels."""
        if self._busy:
            return
        sel = self.selected_rows()
        if not sel:
            self._show_error("No selection", "Select one or more channels.", parent=self.master)
            return
        try:
            mode, s, e = self._window_inputs()
            chunk = self._chunk_seconds()
            jobs, skipped = self.backend.build_jobs(sel, mode, s, e, now=self.now())
        except ValueError as exc:
            self._error("Invalid input", str(exc))
            return
        for row, reason in skipped:
            self._log(f"Skipping {row.seed_id}: {reason}")
        if not jobs:
            self._error("Nothing to download",
                        "None of the selected channels has data in the requested range.")
            return
        self._start_download(jobs, chunk, "flat", [], False, mode)

    def _start_download(self, jobs, chunk_seconds, layout, title_lines, always_confirm,
                        mode):
        est = self.backend.estimate_bytes(jobs, chunk_seconds)
        lines = list(title_lines) + [
            f"Channels: {len(jobs)}",
            f"Total span: {_fmt_span(est['seconds'])}",
            f"Chunks (1 file each): {est['chunks']}",
            f"Raw size: ~{self.backend.format_bytes(est['raw_bytes'])}",
            f"Estimated MiniSEED size: ~{self.backend.format_bytes(est['mseed_bytes'])}"]
        resume = bool(self.resume_var.get())
        if resume:
            lines.append("Files already on disk will be skipped (resume).")
        summary = "\n".join(lines) + "\n\nStart download?"
        self._log("; ".join(lines))
        if (always_confirm or mode == "entire"
                or est["raw_bytes"] > self.backend.CONFIRM_RAW_BYTES):
            if not self._ask_confirm("Confirm download", summary, parent=self.master):
                self._log("Download cancelled.")
                return

        out_dir = self.outdir_var.get().strip() or os.getcwd()
        try:
            os.makedirs(out_dir, exist_ok=True)
        except OSError as e:
            self._error("Output folder", str(e))
            return

        provider = self.provider_var.get().strip()
        convert = bool(self.convert_var.get())
        self._cancel_event = cancel_event = threading.Event()
        self.progress.configure(maximum=max(est["chunks"], 1), value=0)
        self._set_busy(True, downloading=True)
        self.status_var.set("Downloading...")
        self._log(f"Downloading to {out_dir}")

        def progress(done, total, msg):
            self._queue.put(("progress", done, total, msg))

        def detail(info, msg):
            self._queue.put(("detail", info, msg))

        def work():
            try:
                client = self._download_client(provider)
                results = self.backend.run_downloads(
                    jobs, out_dir, client, chunk_seconds=chunk_seconds,
                    convert_to_h5=convert, cancel_event=cancel_event, progress=progress,
                    layout=layout, resume=resume,
                    max_consecutive_errors=MAX_CONSECUTIVE_ERRORS, progress_detail=detail,
                    provider=provider)
                self._queue.put(("done", results, jobs, chunk_seconds,
                                 cancel_event.is_set()))
            except Exception as e:
                self._queue.put(("error", "Download failed", str(e)))
            finally:
                self._queue.put(("idle",))

        threading.Thread(target=work, daemon=True).start()

    def cancel(self):
        """Ask the running download to stop after the current chunk."""
        self._cancel_event.set()
        self._log("Cancelling after current chunk...")
        self.status_var.set("Cancelling...")

    def _on_close(self):
        self._cancel_event.set()
        self.master.destroy()

    def destroy(self):
        self._cancel_event.set()
        for name in ('_poll_id', '_preview_after', '_arraysta_after'):
            token = getattr(self, name, None)
            if token is not None:
                try:
                    self.after_cancel(token)
                except tk.TclError:
                    pass
                setattr(self, name, None)
        super().destroy()

    # ---------------------------------------------------------- main loop
    def _show_rows(self, rows):
        self.rows = self.backend.filter_station_rows(rows,**self._station_filters())
        self.tree.delete(*self.tree.get_children())
        for i, r in enumerate(self.rows):
            self.tree.insert("", "end", iid=str(i), values=(
                r.network, r.station, r.location or "--", r.channel,
                _fmt_num(r.sample_rate, "g"), _fmt_num(r.latitude, ".4f"),
                _fmt_num(r.longitude, ".4f"), _fmt_date(r.start),
                _fmt_date(r.end) if r.end is not None else "open"))
        n = len(self.rows)
        station_count=len({(r.network,r.station) for r in self.rows})
        self.count_var.set(f"{station_count} matching stations / {n} channel epochs (all listed)")
        self.station_box.configure(values=self.backend.code_choices(self.rows, "station"))
        self.location_box.configure(values=self.backend.code_choices(
            self.rows, "location", LOCATION_PRESETS))
        self.channel_box.configure(values=self._chan_query_choices(self.rows))
        self._log(f"{n} channels found.")
        self.status_var.set(f"{n} channels found.")

    def _show_networks(self, rows):
        self.networks = list(rows)
        for i, r in enumerate(self.networks[:MAX_TREE_ROWS]):
            self.arrays_tree.insert("", "end", iid=str(i), values=(
                r.code, r.description, _fmt_date(r.start),
                _fmt_date(r.end) if r.end is not None else "open",
                "?" if r.station_count is None else r.station_count))
        n = len(self.networks)
        self.arrays_count_var.set(f"{n} arrays found" + (
            f" (showing first {MAX_TREE_ROWS})" if n > MAX_TREE_ROWS else "")
            + ("; select one to preview it." if n else "."))
        codes = [r.code for r in self.networks]
        self.array_net_box.configure(values=_dedupe(list(ARRAY_NET_PRESETS) + codes))
        self.network_box.configure(values=_dedupe(["*", "CH"] + codes))
        self._log(f"{n} arrays found.")
        self.status_var.set(f"{n} arrays found.")

    def _show_preview(self, code, rows):
        self.array_rows = self.backend.filter_station_rows(rows,**self._station_filters())
        self.array_code = code
        tree = self.preview_tree
        tree.delete(*tree.get_children())
        for i, r in enumerate(self.array_rows):
            tree.insert("", "end", iid=str(i), values=(
                r.station, r.location or "--", r.channel, _fmt_num(r.latitude, ".4f"),
                _fmt_num(r.longitude, ".4f"), _fmt_num(r.elevation, ".0f"),
                _fmt_num(r.sample_rate, "g"), r.sensor, _fmt_date(r.start),
                _fmt_date(r.end) if r.end is not None else "open"))
        info = self.backend.summarize_array(self.array_rows)
        times = self.backend.suggest_times(self.now(), info["start"])
        self.start_entry.configure(values=times)
        self.end_entry.configure(values=times)
        self.array_chan_box.configure(values=self.backend.channel_choices(self.array_rows))
        self.array_sta_box.configure(values=self.backend.code_choices(self.array_rows,
                                                                      "station"))
        if self.array_rows:
            self._update_array_summary()
        else:
            self.array_summary_var.set(f"{code}: no matching stations/channels "
                                       "(check prefix, bounds and station filters).")
            self._request_geometry([],'array')
        self._log(f"{code}: {info['stations']} stations, {info['channels']} channels.")
        self.status_var.set(f"Preview of {code} loaded.")

    def _update_array_summary(self, *_):
        if not self.array_rows:
            return
        info = self.backend.summarize_array(self.array_rows)
        span_end = "open" if info["end"] is None else _fmt_date(info["end"])
        head = (f"{self.array_code}: {info['stations']} matching stations, {info['channels']} "
                f"channels ({info['epochs']} epochs); span {_fmt_date(info['start']) or '?'}"
                f" - {span_end}")
        if info["missing_coords"]:
            head += f"; {info['missing_coords']} station(s) without coordinates"
        head += '; all matching rows listed'
        try:
            mode, s, e = self._window_inputs(write_back=False)
            chunk = self._chunk_seconds()
            jobs, _ = self.backend.build_array_jobs(self.array_rows, mode, s, e,
                                                    now=self.now(),
                                                    channel_filter=self.array_chan_var.get())
            est = self.backend.estimate_bytes(jobs, chunk)
            fb = self.backend.format_bytes
            tail = (f"Selected: {len(jobs)} channels, {est['chunks']} chunks, "
                    f"~{fb(est['raw_bytes'])} raw / ~{fb(est['mseed_bytes'])} MiniSEED")
        except ValueError as exc:
            tail = f"Estimate unavailable: {exc}"
            jobs = []
        self.array_summary_var.set(head + "\n" + tail)
        self._request_geometry([j.row for j in jobs], "array")

    def _update_channel_preview(self, *_):
        try:
            mode, s, e = self._window_inputs(write_back=False)
            jobs, _ = self.backend.build_jobs(self.selected_rows(), mode, s, e, now=self.now())
            self._request_geometry([j.row for j in jobs], "channels")
        except ValueError as exc:
            self.channel_spacing_var.set(f"Preview unavailable: {exc}")

    def _request_geometry(self, rows, kind):
        # Geodesic all-pair statistics stay off the Tk event loop.
        key = tuple(sorted({(r.network, r.station, repr(r.latitude), repr(r.longitude), repr(r.elevation)) for r in rows}))
        if self._geometry_keys.get(kind) == key:
            return
        self._geometry_keys[kind] = key
        self._geometry_tokens[kind] += 1
        token = self._geometry_tokens[kind]
        var = self.spacing_var if kind == "array" else self.channel_spacing_var
        var.set("Mean station spacing (m): calculating...")
        def work():
            try:
                info = self.backend.array_geometry(rows)
                self._queue.put(("geometry", kind, token, info))
            except Exception as exc:
                self._queue.put(("geometry_error", kind, token, str(exc)))
        threading.Thread(target=work, daemon=True).start()

    def _show_geometry(self, kind, token, info):
        if token != self._geometry_tokens[kind]:
            return
        var = self.spacing_var if kind == "array" else self.channel_spacing_var
        view = self.station_map if kind == "array" else self.channel_map
        if info["mean_nearest_neighbor_m"] is None:
            text = "Mean station spacing (m): -- (needs 2 valid station coordinates)"
        else:
            text = (f"Mean station spacing: {info['mean_nearest_neighbor_m']:.2f} m (nearest neighbour)\n"
                    f"All-pair mean: {info['mean_pair_distance_m']:.2f} m; range: "
                    f"{info['min_pair_distance_m']:.2f}-{info['max_pair_distance_m']:.2f} m; "
                    f"{info['stations_with_coordinates']} unique stations")
        if info['missing_coordinates'] or info['coordinate_conflicts']:
            text += f"; missing coordinates: {info['missing_coordinates']}; conflicting stations: {len(info['coordinate_conflicts'])}"
        if info['coincident_pair_count']:
            text += f"; colocated pairs: {info['coincident_pair_count']}"
        var.set(text)
        view.plot(info['points'], title="Selected download array (m)")

    def open_hvsr(self):
        if self._busy:
            return
        try:
            from ndstudio.project import new_project
            from ndstudio.hvsr import discover,write_json,DEFAULTS
            import json
            folder=Path(self.outdir_var.get().strip() or os.getcwd()).resolve()
            inv=discover(folder);ready=[r for r in inv['instruments'] if r['ready']]
            if not ready: raise ValueError("No complete 3-component instrument found. Download Z + N/E (or oriented 1/2), same channel family and location; use HH? / BH? instead of HHZ only.")
            root=new_project(folder.parent/(folder.name+'_hvsr_project'))
            spec=json.loads((root/'project.json').read_text(encoding='utf-8-sig'));spec['input_folder']=str(folder);write_json(root/'project.json',spec)
            if not (root/'hvsr_settings.json').exists():write_json(root/'hvsr_settings.json',dict(folder=str(folder),instrument_id=ready[0]['instrument_id'],settings=DEFAULTS))
            subprocess.Popen([sys.executable,str(Path(__file__).with_name('run.py')),'gui','--project',str(root),'--hvsr-tab'],cwd=str(Path(__file__).parent))
            self._log(f"Opened single-station H/V: {root}; choose an instrument before analysis")
        except Exception as exc:
            self._error("Cannot open H/V",str(exc))

    def open_dispersion(self):
        """Create/reuse a portable Studio project and open its input page."""
        if self._busy:
            return
        try:
            from ndstudio.project import new_project
            from dispersion.common import write_json
            import json
            folder = Path(self.outdir_var.get().strip() or os.getcwd()).resolve()
            if not any(p.suffix.lower() in ('.h5','.hdf5','.hdf') for p in folder.rglob('*') if p.is_file()):
                raise ValueError("No H5 found in output folder. Download with H5 + full metadata first.")
            root = new_project(folder.parent/(folder.name + '_dispersion_project'))
            spec = json.loads((root/'project.json').read_text(encoding='utf-8-sig'))
            spec['input_folder'] = str(folder)
            write_json(root/'project.json', spec)
            subprocess.Popen([sys.executable, str(Path(__file__).with_name('run.py')), 'gui', '--project', str(root)],
                             cwd=str(Path(__file__).parent))
            self._log(f"Opened dispersion project: {root}; input: {folder}")
        except Exception as exc:
            self._error("Cannot open dispersion", str(exc))

    def _show_results(self, results, jobs, chunk_seconds, cancel_requested):
        for r in results:
            self._log(f"{r.seed_id}: {len(r.files)} MiniSEED, {len(r.h5_files)} HDF5, "
                      f"{len(r.existing)} already present, "
                      f"{len(r.skipped)} empty chunks skipped, {len(r.errors)} errors"
                      + (", cancelled" if r.cancelled else ""))
            for s, e, msg in r.errors:
                self._log(f"  error {s} - {e}: {msg}")
        info = self.backend.summarize_results(results, jobs, chunk_seconds,
                                              cancel_requested)
        if info["aborted"]:
            self._error("Download stopped", info["aborted"] + "\nCheck the provider/network "
                        "and try again; resume will skip files already downloaded.")
        self.status_var.set(info["text"])
        self._log(info["text"])

    def _handle(self, msg):
        kind = msg[0]
        if kind == "geometry":
            self._show_geometry(*msg[1:])
        elif kind == "geometry_error":
            _, surface, token, error = msg
            if token == self._geometry_tokens[surface]:
                var = self.spacing_var if surface == "array" else self.channel_spacing_var
                var.set(f"Spacing unavailable: {error}")
        elif kind == "rows":
            if msg[1]==self._browse_token:self._show_rows(msg[2])
        elif kind == "networks":
            if msg[1]==self._network_token:self._show_networks(msg[2])
        elif kind == "preview":
            _, token, code, rows = msg
            if token == self._preview_token:
                self._preview_loading = False
                self._show_preview(code, rows)
        elif kind == "preview_error":
            _, token, code, text = msg
            if token == self._preview_token:
                self._preview_loading = False
                self.array_summary_var.set(f"Preview of {code} failed.")
                self.station_map.plot([], title=code)
                self._error("Preview failed", text)
        elif kind == "progress":
            _, done, total, text = msg
            self.progress.configure(maximum=max(total, 1), value=done)
            self._log(text)
        elif kind == "detail":
            info = msg[1]
            self.progress.configure(maximum=max(info.chunks_total, 1),
                                    value=info.chunks_done)
            self._last_detail = (f"Stations {info.stations_done}/{info.stations_total}  |  "
                                 f"Channels {info.channels_done}/{info.channels_total}  |  "
                                 f"Chunks {info.chunks_done}/{info.chunks_total}")
            if not self._cancel_event.is_set():
                self.status_var.set(self._last_detail)
        elif kind == "log":
            self._log(msg[1])
        elif kind == "error":
            if msg[1] == "Browse failed" and not self.rows:
                self.count_var.set("No channels loaded.")
            elif msg[1] == "List arrays failed" and not self.networks:
                self.arrays_count_var.set("No arrays loaded.")
            self._error(msg[1], msg[2])
        elif kind == "done":
            self._show_results(*msg[1:])
        elif kind == "idle":
            self._set_busy(False)

    def _poll_queue(self):
        try:
            while True:
                msg = self._queue.get_nowait()
                try:
                    self._handle(msg)
                except Exception as e:  # keep the GUI alive whatever happens
                    try:
                        self._log(f"Internal error: {e}")
                    except tk.TclError:
                        pass
        except queue.Empty:
            pass
        try:
            self._poll_id = self.after(POLL_MS, self._poll_queue)
        except tk.TclError:
            pass  # window destroyed


def main():
    root = tk.Tk()
    root.title("FDSN station filters -> MiniSEED/HDF5 · Studio 2.5")
    root.minsize(1100, 800)
    H5DataApp(root).pack(fill="both", expand=True)
    root.mainloop()


if __name__ == "__main__":
    main()
