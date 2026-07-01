#!/usr/bin/env python3
"""
TimeSCape GUI  –  Circadian Analysis (Python v0.2)
===================================================
Python equivalent of TimeSCape_GUI.m (MATLAB v0.2) and app.R (R Shiny).

WORKFLOW:
  1. Load h5ad file  →  AnnData is read into memory.
  2. ① Define / Edit ZT Times  →  builds the tmeta DataFrame.
  3. ② Analysis Settings  →  cell type, normalization, period, workers.
  4. ③ Run Analysis  →  writes 6 CSV files per cell type.
  5. ④ Generate Heatmap  →  Z-score heatmap of confident genes.
  6. ⑤ Gene Explorer  →  batch PNG export or single-gene cosine plot.

Run from TimeSCape_py/:
    python app_gui.py

Or from anywhere:
    python /path/to/TimeSCape_py/app_gui.py
"""

from __future__ import annotations

import os
import sys
import re
import queue
import threading
import warnings
from datetime import datetime

# ── Set TkAgg BEFORE any other matplotlib import ──────────────────────────────
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext

import numpy as np
import pandas as pd

# ── TimeSCape package import ──────────────────────────────────────────────────
_PKG_DIR = os.path.dirname(os.path.abspath(__file__))
if _PKG_DIR not in sys.path:
    sys.path.insert(0, _PKG_DIR)

try:
    from timescape import (
        run_timescape,
        build_tmeta,
        plot_gene_single,
        generate_heatmap,
        save_batch_plots,
    )
    import timescape as _ts
    _TIMESCAPE_VERSION = _ts.__version__
    TIMESCAPE_OK = True
except ImportError as _e:
    TIMESCAPE_OK = False
    _TIMESCAPE_VERSION = "N/A"
    _IMPORT_ERR = str(_e)


# ── Colour palette (light theme — mirrors MATLAB constants) ──────────────────
_LIGHT = dict(
    bg        = "#FFFFFF",
    pan_bg    = "#F2F4FA",
    accent    = "#27529E",
    btn_run   = "#1A7A3D",
    btn_plot  = "#3870B0",
    btn_warn  = "#9E2E2E",
    btn_neu   = "#595959",
    fg        = "#111111",
    inp_bg    = "#FFFFFF",
    log_bg    = "#F8F8F8",
    ax_bg     = "#FFFFFF",
)

_DARK = dict(
    bg        = "#1E1E26",
    pan_bg    = "#2B2B38",
    accent    = "#7AB4FF",
    btn_run   = "#256B40",
    btn_plot  = "#3870B0",
    btn_warn  = "#9E2E2E",
    btn_neu   = "#595959",
    fg        = "#E0E0E6",
    inp_bg    = "#38384A",
    log_bg    = "#252530",
    ax_bg     = "#1A1A22",
)


def _safe_name(s: str) -> str:
    """Convert a cell type name to a filesystem-safe string (mirrors MATLAB safe_name)."""
    s = re.sub(r"[^a-zA-Z0-9_]", "_", str(s).strip())
    s = re.sub(r"_+", "_", s).strip("_")
    return s


# =============================================================================
#   MAIN GUI CLASS
# =============================================================================

class TimeSCapeGUI:
    """Main application window — mirrors MATLAB TimeSCape_GUI.m."""

    # ------------------------------------------------------------------
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("TimeSCape  –  Circadian Analysis  (v0.2)")
        self.root.geometry("1200x760")
        self.root.minsize(1000, 680)

        self._theme = "light"
        self._c     = _LIGHT.copy()

        # ── State variables ──────────────────────────────────────────────
        self.adata        = None        # loaded AnnData
        self.tmeta        = None        # ZT metadata DataFrame
        self.outdir       = ""          # last used output directory
        self._status_q: queue.Queue = queue.Queue()

        self._build_ui()
        self._apply_theme()
        self._check_queue()

        if not TIMESCAPE_OK:
            self._log(f"⚠  timescape package not found: {_IMPORT_ERR}", warn=True)
            self._log("   Run:  pip install -e .  (from TimeSCape_py/)", warn=True)
        else:
            self._log(f"TimeSCape v{_TIMESCAPE_VERSION} loaded.")

    # ------------------------------------------------------------------
    #   TOP-LEVEL LAYOUT
    # ------------------------------------------------------------------

    def _build_ui(self) -> None:
        self.root.columnconfigure(1, weight=1)
        self.root.rowconfigure(0, weight=1)

        # Left control panel (fixed width)
        self._left = tk.Frame(self.root, width=310, bd=0)
        self._left.grid(row=0, column=0, sticky="nsew", padx=(4, 2), pady=4)
        self._left.grid_propagate(False)
        self._left.columnconfigure(0, weight=1)

        # Right area (plot + toolbar + top buttons)
        self._right = tk.Frame(self.root)
        self._right.grid(row=0, column=1, sticky="nsew", padx=(2, 4), pady=4)
        self._right.rowconfigure(1, weight=1)
        self._right.columnconfigure(0, weight=1)

        self._build_left_panel()
        self._build_right_panel()

    # ------------------------------------------------------------------
    #   LEFT PANEL
    # ------------------------------------------------------------------

    def _build_left_panel(self) -> None:
        p = self._left

        # Scrollable canvas so the panel works on small screens
        self._lpane_canvas = tk.Canvas(p, highlightthickness=0)
        lsb = ttk.Scrollbar(p, orient="vertical",
                             command=self._lpane_canvas.yview)
        self._lpane_canvas.configure(yscrollcommand=lsb.set)

        lsb.pack(side="right", fill="y")
        self._lpane_canvas.pack(side="left", fill="both", expand=True)

        self._lpane = tk.Frame(self._lpane_canvas)
        self._lpane_window = self._lpane_canvas.create_window(
            (0, 0), window=self._lpane, anchor="nw"
        )

        self._lpane.bind("<Configure>", self._on_lpane_configure)
        self._lpane_canvas.bind("<Configure>", self._on_canvas_configure)

        # ── Sections ──────────────────────────────────────────────────
        self._build_load_section()
        self._build_metadata_section()
        self._build_settings_section()
        self._build_run_section()
        self._build_heatmap_section()
        self._build_gene_explorer_section()
        self._build_status_log()

    def _on_lpane_configure(self, _event) -> None:
        self._lpane_canvas.configure(
            scrollregion=self._lpane_canvas.bbox("all")
        )

    def _on_canvas_configure(self, event) -> None:
        self._lpane_canvas.itemconfig(
            self._lpane_window, width=event.width
        )

    # ------------------------------------------------------------------
    #   SECTION: LOAD DATA
    # ------------------------------------------------------------------

    def _build_load_section(self) -> None:
        f = self._section("Load h5ad Data")

        row = tk.Frame(f)
        row.pack(fill="x", pady=2)
        tk.Label(row, text="File:").pack(side="left")
        self._v_h5ad = tk.StringVar(value="")
        tk.Entry(row, textvariable=self._v_h5ad, width=22).pack(
            side="left", padx=4, fill="x", expand=True
        )
        self._btn(row, "Browse…", self._browse_h5ad, "neu").pack(side="left")

        row2 = tk.Frame(f)
        row2.pack(fill="x", pady=2)
        self._v_ct_col = tk.StringVar(value="cell_type")
        self._v_zt_col = tk.StringVar(value="ZT_time")
        tk.Label(row2, text="CellType col:").pack(side="left")
        tk.Entry(row2, textvariable=self._v_ct_col, width=12).pack(side="left", padx=2)
        tk.Label(row2, text="ZT col:").pack(side="left", padx=(6, 0))
        tk.Entry(row2, textvariable=self._v_zt_col, width=8).pack(side="left", padx=2)

        self._btn(f, "▶  Load Data", self._load_data, "run").pack(
            fill="x", pady=(4, 2)
        )
        self._lbl_load_status = tk.Label(
            f, text="No data loaded.", anchor="w", wraplength=280
        )
        self._lbl_load_status.pack(fill="x")

    # ------------------------------------------------------------------
    #   SECTION: ① ZT METADATA
    # ------------------------------------------------------------------

    def _build_metadata_section(self) -> None:
        f = self._section("① Define ZT Time Metadata")

        self._btn(f, "Define / Edit ZT Times", self._define_tmeta, "plot").pack(
            fill="x", pady=2
        )
        self._lbl_tmeta = tk.Label(
            f, text="⚠  Tmeta not defined yet",
            fg="#B35000", anchor="w", wraplength=280
        )
        self._lbl_tmeta.pack(fill="x")

    # ------------------------------------------------------------------
    #   SECTION: ② ANALYSIS SETTINGS
    # ------------------------------------------------------------------

    def _build_settings_section(self) -> None:
        f = self._section("② Analysis Settings")

        # Cell type
        row = tk.Frame(f)
        row.pack(fill="x", pady=2)
        tk.Label(row, text="Cell type:", width=12, anchor="w").pack(side="left")
        self._cb_celltype = ttk.Combobox(row, values=[], width=18, state="readonly")
        self._cb_celltype.pack(side="left", fill="x", expand=True, padx=2)

        # Normalization
        row2 = tk.Frame(f)
        row2.pack(fill="x", pady=2)
        tk.Label(row2, text="Normalization:", width=12, anchor="w").pack(side="left")
        self._cb_norm = ttk.Combobox(
            row2,
            values=["lib_size", "none", "logcounts"],
            width=18, state="readonly"
        )
        self._cb_norm.set("lib_size")
        self._cb_norm.pack(side="left", fill="x", expand=True, padx=2)

        # Period + workers
        row3 = tk.Frame(f)
        row3.pack(fill="x", pady=2)
        self._v_period12 = tk.BooleanVar(value=False)
        tk.Checkbutton(
            row3, text="Use 12-hr period",
            variable=self._v_period12
        ).pack(side="left")
        tk.Label(row3, text="Workers:").pack(side="left", padx=(10, 2))
        self._v_workers = tk.StringVar(value="-1")
        tk.Entry(row3, textvariable=self._v_workers, width=5).pack(side="left")

        # Test type
        row4 = tk.Frame(f)
        row4.pack(fill="x", pady=2)
        tk.Label(row4, text="Test type:", width=12, anchor="w").pack(side="left")
        self._cb_test = ttk.Combobox(
            row4, values=["Ftest", "LRT"], width=10, state="readonly"
        )
        self._cb_test.set("Ftest")
        self._cb_test.pack(side="left", padx=2)

        # Output directory
        row5 = tk.Frame(f)
        row5.pack(fill="x", pady=2)
        tk.Label(row5, text="Output dir:", width=12, anchor="w").pack(side="left")
        self._v_outdir = tk.StringVar(value=os.getcwd())
        tk.Entry(row5, textvariable=self._v_outdir, width=18).pack(
            side="left", fill="x", expand=True, padx=2
        )
        self._btn(row5, "…", self._browse_outdir, "neu").pack(side="left")

    # ------------------------------------------------------------------
    #   SECTION: ③ RUN ANALYSIS
    # ------------------------------------------------------------------

    def _build_run_section(self) -> None:
        f = self._section("③ Run Analysis  →  writes ALL + Confident files")

        info = (
            "Produces two result sets:\n"
            "  • ALL genes  (every gene tested)\n"
            "  • Confident  (F-test AND corr p < 0.05)"
        )
        tk.Label(f, text=info, anchor="w", justify="left",
                 wraplength=280).pack(fill="x", pady=(0, 4))

        row = tk.Frame(f)
        row.pack(fill="x", pady=2)
        self._btn(row, "▶  Run Selected", self._run_analysis, "run").pack(
            side="left", fill="x", expand=True, padx=(0, 2)
        )
        self._btn(row, "▶▶  Run All", self._run_all, "run").pack(
            side="left", fill="x", expand=True
        )

        # Progress bar
        self._pb = ttk.Progressbar(f, mode="indeterminate", length=280)
        self._pb.pack(fill="x", pady=2)

        self._lbl_run = tk.Label(f, text="", anchor="w", wraplength=280)
        self._lbl_run.pack(fill="x")

    # ------------------------------------------------------------------
    #   SECTION: ④ HEATMAP
    # ------------------------------------------------------------------

    def _build_heatmap_section(self) -> None:
        f = self._section("④ Generate Heatmap")

        self._v_strict = tk.BooleanVar(value=True)
        tk.Checkbutton(
            f, text="Confident genes only  (recommended)",
            variable=self._v_strict
        ).pack(anchor="w")

        self._v_circ_only = tk.BooleanVar(value=False)
        tk.Checkbutton(
            f, text="Restrict to core circadian gene set",
            variable=self._v_circ_only
        ).pack(anchor="w")

        row = tk.Frame(f)
        row.pack(fill="x", pady=2)
        tk.Label(row, text="Custom label:", width=13, anchor="w").pack(side="left")
        self._v_custom_name = tk.StringVar(value="")
        tk.Entry(row, textvariable=self._v_custom_name, width=16).pack(
            side="left", padx=2
        )

        self._btn(f, "Generate Heatmap", self._do_heatmap, "plot").pack(
            fill="x", pady=4
        )

    # ------------------------------------------------------------------
    #   SECTION: ⑤ GENE EXPLORER
    # ------------------------------------------------------------------

    def _build_gene_explorer_section(self) -> None:
        f = self._section("⑤ Gene Explorer")

        # Batch plots
        row = tk.Frame(f)
        row.pack(fill="x", pady=2)
        tk.Label(row, text="Batch plots:", width=12, anchor="w").pack(side="left")
        self._cb_plottype = ttk.Combobox(
            row,
            values=["Confident genes", "Non-confident genes", "Core circadian genes"],
            width=18, state="readonly"
        )
        self._cb_plottype.set("Confident genes")
        self._cb_plottype.pack(side="left", fill="x", expand=True, padx=2)

        self._btn(f, "Save Gene Plots (Batch)", self._batch_plots, "plot").pack(
            fill="x", pady=2
        )

        ttk.Separator(f, orient="horizontal").pack(fill="x", pady=4)

        # Single gene
        row2 = tk.Frame(f)
        row2.pack(fill="x", pady=2)
        tk.Label(row2, text="Single gene:", width=12, anchor="w").pack(side="left")
        self._cb_gene = ttk.Combobox(row2, values=[], width=18)
        self._cb_gene.pack(side="left", fill="x", expand=True, padx=2)

        row2b = tk.Frame(f)
        row2b.pack(fill="x", pady=0)
        tk.Label(row2b, text="Or type:", width=12, anchor="w").pack(side="left")
        self._v_gene_entry = tk.StringVar(value="")
        tk.Entry(row2b, textvariable=self._v_gene_entry, width=20).pack(
            side="left", padx=2
        )

        # Overlay options
        row3 = tk.Frame(f)
        row3.pack(fill="x", pady=2)
        self._v_show_cells = tk.BooleanVar(value=False)
        tk.Checkbutton(
            row3, text="Overlay single-cell data",
            variable=self._v_show_cells
        ).pack(side="left")

        row4 = tk.Frame(f)
        row4.pack(fill="x", pady=2)
        tk.Label(row4, text="SC style:", width=8, anchor="w").pack(side="left")
        self._v_violin = tk.BooleanVar(value=True)
        tk.Radiobutton(row4, text="Violin", variable=self._v_violin,
                       value=True).pack(side="left")
        tk.Radiobutton(row4, text="Dots", variable=self._v_violin,
                       value=False).pack(side="left")

        self._btn(f, "Plot Single Gene", self._plot_gene, "plot").pack(
            fill="x", pady=(4, 2)
        )

    # ------------------------------------------------------------------
    #   STATUS LOG
    # ------------------------------------------------------------------

    def _build_status_log(self) -> None:
        f = self._section("Status Log")
        self._log_widget = scrolledtext.ScrolledText(
            f, height=6, wrap="word", state="disabled", font=("Courier", 9)
        )
        self._log_widget.pack(fill="both", expand=True)

    # ------------------------------------------------------------------
    #   RIGHT PANEL: matplotlib canvas
    # ------------------------------------------------------------------

    def _build_right_panel(self) -> None:
        # Top-right buttons row
        btn_row = tk.Frame(self._right)
        btn_row.grid(row=0, column=0, sticky="ew", pady=(0, 2))
        btn_row.columnconfigure(0, weight=1)

        self._btn(btn_row, "💾  Save Figure…", self._save_figure, "plot").pack(
            side="right", padx=2
        )
        self._theme_btn = self._btn(
            btn_row, "🌙  Dark Mode", self._toggle_theme, "neu"
        )
        self._theme_btn.pack(side="right", padx=2)

        # Matplotlib figure embedded in canvas
        self._fig = Figure(figsize=(8, 5.5), dpi=96)
        self._ax  = self._fig.add_subplot(111)
        self._ax.set_xlabel("Zeitgeber Time (hrs)", fontsize=11)
        self._ax.set_ylabel("Expression (log-normalised)", fontsize=11)
        self._ax.set_title("← Use the Gene Explorer panel to plot a single gene",
                           fontsize=11)
        self._fig.tight_layout()

        self._canvas = FigureCanvasTkAgg(self._fig, master=self._right)
        self._canvas.get_tk_widget().grid(
            row=1, column=0, sticky="nsew"
        )

        # Navigation toolbar
        toolbar_frame = tk.Frame(self._right)
        toolbar_frame.grid(row=2, column=0, sticky="ew")
        self._toolbar = NavigationToolbar2Tk(self._canvas, toolbar_frame)
        self._toolbar.update()

    # ------------------------------------------------------------------
    #   SECTION HELPER
    # ------------------------------------------------------------------

    def _section(self, title: str) -> ttk.LabelFrame:
        lf = ttk.LabelFrame(self._lpane, text=f"  {title}  ",
                            padding=6)
        lf.pack(fill="x", padx=4, pady=3)
        return lf

    def _btn(self, parent, text, cmd, style="neu"):
        """Create a styled tk.Button (bg/fg from palette)."""
        colors = {
            "run":  (self._c["btn_run"],  "#FFFFFF"),
            "plot": (self._c["btn_plot"], "#FFFFFF"),
            "warn": (self._c["btn_warn"], "#FFFFFF"),
            "neu":  (self._c["btn_neu"],  "#FFFFFF"),
        }
        bg, fg = colors.get(style, (self._c["btn_neu"], "#FFFFFF"))
        b = tk.Button(parent, text=text, command=cmd,
                      bg=bg, fg=fg, activeforeground="#FFFFFF",
                      font=("TkDefaultFont", 9, "bold"),
                      relief="flat", padx=6, pady=3,
                      cursor="hand2")
        # Store style tag so theme updates can find it
        b._ts_style = style
        return b

    # ------------------------------------------------------------------
    #   CALLBACKS: LOAD DATA
    # ------------------------------------------------------------------

    def _browse_h5ad(self) -> None:
        path = filedialog.askopenfilename(
            title="Select h5ad file",
            filetypes=[("HDF5 AnnData", "*.h5ad"), ("All files", "*.*")]
        )
        if path:
            self._v_h5ad.set(path)
            # Set default outdir to the same directory
            self._v_outdir.set(os.path.dirname(path))

    def _browse_outdir(self) -> None:
        d = filedialog.askdirectory(title="Select output directory")
        if d:
            self._v_outdir.set(d)

    def _load_data(self) -> None:
        path = self._v_h5ad.get().strip()
        if not path:
            messagebox.showwarning("No file", "Please select an h5ad file first.")
            return
        if not os.path.isfile(path):
            messagebox.showerror("File not found", f"Not found:\n{path}")
            return
        ct_col = self._v_ct_col.get().strip()
        zt_col = self._v_zt_col.get().strip()

        self._lbl_load_status.config(text="⏳  Loading…")
        self.root.update_idletasks()

        def _worker():
            try:
                import scanpy as sc
                adata = sc.read_h5ad(path)
                self._status_q.put(("loaded", adata, ct_col, zt_col))
            except Exception as e:
                self._status_q.put(("load_err", str(e)))

        threading.Thread(target=_worker, daemon=True).start()

    # ------------------------------------------------------------------
    #   CALLBACKS: ZT METADATA EDITOR
    # ------------------------------------------------------------------

    def _define_tmeta(self) -> None:
        if self.adata is None:
            messagebox.showwarning("No data", "Load an h5ad file first (Load Data section).")
            return

        zt_col = self._v_zt_col.get().strip()
        if zt_col not in self.adata.obs.columns:
            messagebox.showerror("Column not found",
                                 f"ZT column '{zt_col}' not found.\n"
                                 f"Available: {list(self.adata.obs.columns)}")
            return

        batches = sorted(self.adata.obs[zt_col].dropna().unique().astype(str))

        # Count cells per batch
        counts = {b: int((self.adata.obs[zt_col].astype(str) == b).sum())
                  for b in batches}

        # Pre-fill ZT hours if labels look like ZTnn
        def _prefill(lbl):
            m = re.match(r"^ZT(\d+)$", lbl, re.IGNORECASE)
            return float(m.group(1)) if m else ""

        # Build dialog
        dlg = tk.Toplevel(self.root)
        dlg.title("Define ZT Tmeta")
        dlg.geometry("560x460")
        dlg.grab_set()
        dlg.resizable(False, True)
        dlg.configure(bg=self._c["bg"])

        tk.Label(
            dlg,
            text="Set the ZT hour (numeric) for each batch label.\n"
                 "Enter -1 to exclude a batch entirely.",
            bg=self._c["bg"], fg=self._c["fg"],
            anchor="w", justify="left", padx=12, pady=6
        ).pack(fill="x")

        # Scrollable table
        hdr = tk.Frame(dlg, bg=self._c["pan_bg"])
        hdr.pack(fill="x", padx=8)
        for txt, w in [("Original Label", 200), ("ZT Hour (numeric)", 130),
                       ("# Cells", 80)]:
            tk.Label(hdr, text=txt, width=w // 8, bg=self._c["pan_bg"],
                     fg=self._c["accent"], font=("TkDefaultFont", 9, "bold"),
                     anchor="w").pack(side="left", padx=2)

        canvas_dlg = tk.Canvas(dlg, bg=self._c["bg"], highlightthickness=0)
        sb_dlg = ttk.Scrollbar(dlg, orient="vertical",
                               command=canvas_dlg.yview)
        canvas_dlg.configure(yscrollcommand=sb_dlg.set)
        sb_dlg.pack(side="right", fill="y", padx=(0, 8))
        canvas_dlg.pack(fill="both", expand=True, padx=8)

        inner = tk.Frame(canvas_dlg, bg=self._c["bg"])
        canvas_dlg.create_window((0, 0), window=inner, anchor="nw")
        inner.bind("<Configure>",
                   lambda e: canvas_dlg.configure(
                       scrollregion=canvas_dlg.bbox("all")))

        entries = {}   # batch → Entry widget
        for i, b in enumerate(batches):
            row_bg = self._c["bg"] if i % 2 == 0 else self._c["pan_bg"]
            r = tk.Frame(inner, bg=row_bg)
            r.pack(fill="x")

            tk.Label(r, text=b, width=25, anchor="w",
                     bg=row_bg, fg=self._c["fg"]).pack(side="left", padx=2)
            e = tk.Entry(r, width=16, bg=self._c["inp_bg"],
                         fg=self._c["fg"], relief="solid", bd=1)
            pf = _prefill(b)
            if pf != "":
                e.insert(0, str(pf))
            e.pack(side="left", padx=4)
            entries[b] = e

            tk.Label(r, text=str(counts[b]), width=10, anchor="w",
                     bg=row_bg, fg=self._c["fg"]).pack(side="left")

        # Buttons
        btn_row = tk.Frame(dlg, bg=self._c["bg"])
        btn_row.pack(fill="x", padx=8, pady=8)

        def _save():
            rows = []
            for b in batches:
                txt = entries[b].get().strip()
                try:
                    zt = float(txt) if txt != "" else -1.0
                except ValueError:
                    messagebox.showerror("Invalid value",
                                         f"'{txt}' is not a number for batch '{b}'.")
                    return
                new_lbl = f"ZT{int(zt):02d}" if zt >= 0 else "-1"
                rows.append({"old_labels": b, "new_labels": new_lbl,
                              "ZT_times": zt})
            df = pd.DataFrame(rows)
            df = df[df["ZT_times"] >= 0].reset_index(drop=True)
            df = df.sort_values("ZT_times").reset_index(drop=True)
            self.tmeta = df
            n = len(df)
            self._lbl_tmeta.config(
                text=f"✓  Tmeta set: {n} time points",
                fg="#1A7A3D"
            )
            self._log(f"Tmeta saved: {n} time points — "
                      f"{list(df['ZT_times'].astype(int).values)}")
            dlg.destroy()

        def _auto_fill():
            """Pre-fill all entries using build_tmeta auto-parse."""
            if not TIMESCAPE_OK:
                messagebox.showerror("Package missing", "timescape not installed.")
                return
            try:
                tmp = build_tmeta(batches)
                lut = dict(zip(tmp["old_labels"].astype(str),
                               tmp["ZT_times"].astype(float)))
                for b, e in entries.items():
                    if b in lut:
                        e.delete(0, "end")
                        e.insert(0, str(lut[b]))
            except Exception as ex:
                messagebox.showwarning("Auto-fill failed", str(ex))

        tk.Button(btn_row, text="Auto-fill ZT",
                  command=_auto_fill,
                  bg=self._c["btn_neu"], fg="#FFFFFF",
                  font=("TkDefaultFont", 9, "bold"),
                  relief="flat", padx=6).pack(side="left", padx=2)
        tk.Button(btn_row, text="Save & Apply",
                  command=_save,
                  bg=self._c["accent"], fg="#FFFFFF",
                  font=("TkDefaultFont", 9, "bold"),
                  relief="flat", padx=6).pack(side="left", padx=2)
        tk.Button(btn_row, text="Close",
                  command=dlg.destroy,
                  bg=self._c["btn_warn"], fg="#FFFFFF",
                  font=("TkDefaultFont", 9, "bold"),
                  relief="flat", padx=6).pack(side="left", padx=2)

    # ------------------------------------------------------------------
    #   CALLBACKS: RUN ANALYSIS
    # ------------------------------------------------------------------

    def _get_period(self) -> float:
        return 12.0 if self._v_period12.get() else 24.0

    def _run_analysis(self) -> None:
        if not self._check_ready():
            return
        celltype = self._cb_celltype.get().strip()
        if not celltype:
            messagebox.showwarning("No cell type", "Select a cell type first.")
            return
        self._run_ts([celltype])

    def _run_all(self) -> None:
        if not self._check_ready():
            return
        ct_col = self._v_ct_col.get().strip()
        all_types = sorted(self.adata.obs[ct_col].dropna().unique().astype(str))
        self._run_ts(all_types)

    def _run_ts(self, celltypes: list) -> None:
        if not TIMESCAPE_OK:
            messagebox.showerror("Package missing", "timescape not installed.")
            return
        period   = self._get_period()
        norm_str = self._cb_norm.get()
        test_t   = self._cb_test.get()
        outdir   = self._v_outdir.get().strip() or os.getcwd()
        try:
            n_jobs = int(self._v_workers.get())
        except ValueError:
            n_jobs = -1
        zt_col  = self._v_zt_col.get().strip()
        ct_col  = self._v_ct_col.get().strip()

        self._pb.start(12)
        self._lbl_run.config(text=f"⏳  Running {len(celltypes)} cell type(s)…")
        self._log(f"Run started: {len(celltypes)} cell type(s), period={period}h")

        def _worker():
            try:
                custom = celltypes if len(celltypes) == 1 else None
                T1, T2 = run_timescape(
                    adata=self.adata,
                    tmeta=self.tmeta,
                    celltype_col=ct_col,
                    zt_col=zt_col,
                    custom_celltype=custom,
                    period=period,
                    norm_str=norm_str,
                    test_type=test_t,
                    rm_low_conf=True,
                    plot_heat=True,
                    outdir=outdir,
                    n_jobs=n_jobs,
                )
                n_all  = len(T1)
                n_conf = int(
                    ((T1["pvalue"] < 0.05) & (T1["pvalue_corr"] < 0.05)).sum()
                )
                self._status_q.put(("run_ok", n_all, n_conf, outdir))
            except Exception as e:
                self._status_q.put(("run_err", str(e)))

        threading.Thread(target=_worker, daemon=True).start()

    # ------------------------------------------------------------------
    #   CALLBACKS: HEATMAP
    # ------------------------------------------------------------------

    def _do_heatmap(self) -> None:
        if not self._check_ready():
            return
        if not TIMESCAPE_OK:
            messagebox.showerror("Package missing", "timescape not installed.")
            return
        celltype    = self._cb_celltype.get().strip()
        period      = self._get_period()
        strict      = self._v_strict.get()
        circ_only   = self._v_circ_only.get()
        custom_name = self._v_custom_name.get().strip()
        outdir      = self._v_outdir.get().strip() or os.getcwd()
        ct_safe     = _safe_name(celltype)
        ct_outdir   = os.path.join(outdir, ct_safe)

        self._log(f"Generating heatmap for {celltype}…")
        self._pb.start(12)

        def _worker():
            try:
                png = generate_heatmap(
                    celltype=celltype,
                    outdir=ct_outdir,
                    period=period,
                    strict=strict,
                    circ_only=circ_only,
                    custom_name=custom_name,
                )
                if png:
                    self._status_q.put(("heatmap_ok", png))
                else:
                    self._status_q.put(("heatmap_empty",))
            except Exception as e:
                self._status_q.put(("heatmap_err", str(e)))

        threading.Thread(target=_worker, daemon=True).start()

    # ------------------------------------------------------------------
    #   CALLBACKS: GENE EXPLORER
    # ------------------------------------------------------------------

    def _batch_plots(self) -> None:
        if not self._check_ready():
            return
        if not TIMESCAPE_OK:
            messagebox.showerror("Package missing", "timescape not installed.")
            return
        celltype = self._cb_celltype.get().strip()
        pt_map   = {
            "Confident genes": 1,
            "Non-confident genes": 2,
            "Core circadian genes": 3,
        }
        plot_type = pt_map.get(self._cb_plottype.get(), 1)
        period    = self._get_period()
        norm_str  = self._cb_norm.get()
        outdir    = self._v_outdir.get().strip() or os.getcwd()
        zt_col    = self._v_zt_col.get().strip()
        ct_col    = self._v_ct_col.get().strip()
        show_cells = self._v_show_cells.get()
        use_violin = self._v_violin.get()

        self._log(f"Batch plots: {celltype}  type={plot_type}")
        self._pb.start(12)

        def _worker():
            try:
                save_batch_plots(
                    adata=self.adata,
                    tmeta=self.tmeta,
                    celltype=celltype,
                    outdir=outdir,
                    celltype_col=ct_col,
                    zt_col=zt_col,
                    period=period,
                    norm_str=norm_str,
                    plot_type=plot_type,
                    show_cells=show_cells,
                    use_violin=use_violin,
                )
                ct_safe = _safe_name(celltype)
                subdir_map = {1: "plots_confident",
                              2: "plots_non_confident",
                              3: "plots_clock_genes"}
                subdir = os.path.join(outdir, ct_safe,
                                      subdir_map.get(plot_type, "plots"))
                self._status_q.put(("batch_ok", subdir))
            except Exception as e:
                self._status_q.put(("batch_err", str(e)))

        threading.Thread(target=_worker, daemon=True).start()

    def _plot_gene(self) -> None:
        if not self._check_ready():
            return
        if not TIMESCAPE_OK:
            messagebox.showerror("Package missing", "timescape not installed.")
            return

        # Typed entry takes priority over dropdown
        typed = self._v_gene_entry.get().strip()
        gene  = typed if typed else self._cb_gene.get().strip()
        if not gene:
            messagebox.showwarning("No gene", "Enter or select a gene name.")
            return

        if gene not in self.adata.var_names:
            messagebox.showerror("Gene not found",
                                 f"'{gene}' not found in dataset.\nCheck spelling.")
            return

        celltype   = self._cb_celltype.get().strip()
        period     = self._get_period()
        norm_str   = self._cb_norm.get()
        zt_col     = self._v_zt_col.get().strip()
        ct_col     = self._v_ct_col.get().strip()
        show_cells = self._v_show_cells.get()
        use_violin = self._v_violin.get()

        self._log(f"Plotting: {gene}  [{celltype}]")

        try:
            # Clear embedded axes and reuse the embedded figure
            self._ax.cla()
            self._ax.set_title("Plotting…")
            self._canvas.draw()
            self.root.update_idletasks()

            plot_gene_single(
                adata=self.adata,
                gene=gene,
                celltype=celltype,
                tmeta=self.tmeta,
                celltype_col=ct_col,
                zt_col=zt_col,
                period=period,
                norm_str=norm_str,
                show_cells=show_cells,
                use_violin=use_violin,
                ax=self._ax,
            )
            self._fig.tight_layout()
            self._canvas.draw()
            self._log(f"Plot complete: {gene}")
        except Exception as e:
            self._log(f"✗  Plot error: {e}", warn=True)
            messagebox.showerror("Plot error", str(e))

    # ------------------------------------------------------------------
    #   SAVE FIGURE
    # ------------------------------------------------------------------

    def _save_figure(self) -> None:
        path = filedialog.asksaveasfilename(
            title="Save Current Plot As",
            defaultextension=".png",
            initialfile="TimeSCape_plot.png",
            filetypes=[
                ("PNG image",  "*.png"),
                ("SVG vector", "*.svg"),
                ("PDF vector", "*.pdf"),
            ],
        )
        if not path:
            return
        try:
            self._fig.savefig(path, dpi=200, bbox_inches="tight",
                              facecolor=self._fig.get_facecolor())
            self._log(f"Figure saved: {path}")
            messagebox.showinfo("Saved", f"Figure saved:\n{path}")
        except Exception as e:
            messagebox.showerror("Save failed", str(e))

    # ------------------------------------------------------------------
    #   DARK / LIGHT THEME TOGGLE
    # ------------------------------------------------------------------

    def _toggle_theme(self) -> None:
        if self._theme == "light":
            self._theme = "dark"
            self._c = _DARK.copy()
            self._theme_btn.config(text="☀  Light Mode")
        else:
            self._theme = "light"
            self._c = _LIGHT.copy()
            self._theme_btn.config(text="🌙  Dark Mode")
        self._apply_theme()

    def _apply_theme(self) -> None:
        c = self._c
        self.root.configure(bg=c["bg"])
        self._left.configure(bg=c["bg"])
        self._lpane_canvas.configure(bg=c["bg"])
        self._lpane.configure(bg=c["bg"])

        # Style all LabelFrames + Frames recursively
        self._style_widget_tree(self._lpane)

        # Matplotlib axes
        self._fig.set_facecolor(c["ax_bg"])
        self._ax.set_facecolor(c["ax_bg"])
        txt_col = c["fg"]
        self._ax.title.set_color(txt_col)
        self._ax.xaxis.label.set_color(txt_col)
        self._ax.yaxis.label.set_color(txt_col)
        for spine in self._ax.spines.values():
            spine.set_edgecolor(txt_col)
        self._ax.tick_params(colors=txt_col)
        self._canvas.get_tk_widget().configure(bg=c["ax_bg"])
        self._canvas.draw()

    def _style_widget_tree(self, widget) -> None:
        c = self._c
        cls = type(widget).__name__

        if cls == "LabelFrame":
            widget.configure(bg=c["pan_bg"], fg=c["accent"])
        elif cls in ("Frame",):
            try:
                widget.configure(bg=c["bg"])
            except Exception:
                pass
        elif cls == "Label":
            try:
                widget.configure(bg=widget.master["bg"] if widget.master else c["bg"],
                                 fg=c["fg"])
            except Exception:
                pass
        elif cls == "Checkbutton":
            try:
                widget.configure(
                    bg=widget.master["bg"] if widget.master else c["bg"],
                    fg=c["fg"], selectcolor=c["bg"],
                    activebackground=c["bg"], activeforeground=c["fg"]
                )
            except Exception:
                pass
        elif cls == "Radiobutton":
            try:
                widget.configure(
                    bg=widget.master["bg"] if widget.master else c["bg"],
                    fg=c["fg"], selectcolor=c["bg"],
                    activebackground=c["bg"], activeforeground=c["fg"]
                )
            except Exception:
                pass
        elif cls == "Entry":
            try:
                widget.configure(bg=c["inp_bg"], fg=c["fg"],
                                 insertbackground=c["fg"])
            except Exception:
                pass
        elif cls == "Button":
            # Only update theme button — semantic buttons keep their colours
            if hasattr(widget, "_ts_style"):
                style_colors = {
                    "run":  (c["btn_run"],  "#FFFFFF"),
                    "plot": (c["btn_plot"], "#FFFFFF"),
                    "warn": (c["btn_warn"], "#FFFFFF"),
                    "neu":  (c["btn_neu"],  "#FFFFFF"),
                }
                bg, fg = style_colors.get(widget._ts_style, (c["btn_neu"], "#FFFFFF"))
                try:
                    widget.configure(bg=bg, fg=fg, activebackground=bg)
                except Exception:
                    pass
        elif cls == "ScrolledText":
            try:
                widget.configure(bg=c["log_bg"], fg=c["fg"],
                                 insertbackground=c["fg"])
            except Exception:
                pass

        for child in widget.winfo_children():
            self._style_widget_tree(child)

    # ------------------------------------------------------------------
    #   STATUS QUEUE POLLER
    # ------------------------------------------------------------------

    def _check_queue(self) -> None:
        while not self._status_q.empty():
            msg = self._status_q.get_nowait()
            kind = msg[0]

            if kind == "loaded":
                _, adata, ct_col, zt_col = msg
                self.adata = adata
                n_genes = adata.n_vars
                n_cells = adata.n_obs

                # Populate cell type combobox
                if ct_col in adata.obs.columns:
                    cts = sorted(adata.obs[ct_col].dropna().unique().astype(str))
                    self._cb_celltype.configure(values=cts)
                    if cts:
                        self._cb_celltype.set(cts[0])
                else:
                    self._log(f"⚠  Cell type col '{ct_col}' not found in obs.",
                              warn=True)

                # Populate gene dropdown (limit for performance)
                genes = sorted(adata.var_names.tolist())
                MAX_GENES_COMBO = 5000
                if len(genes) > MAX_GENES_COMBO:
                    self._cb_gene.configure(values=genes[:MAX_GENES_COMBO])
                    self._log(f"  Gene dropdown limited to {MAX_GENES_COMBO} "
                              f"(type name in 'Or type:' field for any gene)")
                else:
                    self._cb_gene.configure(values=genes)
                if genes:
                    self._cb_gene.set(genes[0])

                self._lbl_load_status.config(
                    text=f"✓  {n_genes:,} genes × {n_cells:,} cells",
                    fg="#1A7A3D"
                )
                self._log(
                    f"Loaded: {n_genes:,} genes × {n_cells:,} cells  "
                    f"[{os.path.basename(self._v_h5ad.get())}]"
                )
                # Auto-try tmeta from ZT column
                if TIMESCAPE_OK and zt_col in adata.obs.columns:
                    try:
                        zts = sorted(adata.obs[zt_col].dropna().unique().astype(str))
                        self.tmeta = build_tmeta(zts)
                        n_tp = len(self.tmeta)
                        self._lbl_tmeta.config(
                            text=f"✓  Auto-detected {n_tp} time points",
                            fg="#1A7A3D"
                        )
                        self._log(f"Tmeta auto-set from '{zt_col}': "
                                  f"{n_tp} time points "
                                  f"{list(self.tmeta['ZT_times'].astype(int))}")
                    except Exception:
                        pass  # Will be done manually

            elif kind == "load_err":
                _, err = msg
                self._lbl_load_status.config(text=f"✗  Load failed", fg="#9E2E2E")
                self._log(f"✗  Load error: {err}", warn=True)
                messagebox.showerror("Load failed", err)

            elif kind == "run_ok":
                _, n_all, n_conf, outdir = msg
                self._pb.stop()
                self.outdir = outdir
                self._lbl_run.config(
                    text=f"✓  {n_all:,} genes tested, {n_conf:,} confident"
                )
                self._log(f"Run done: {n_all:,} tested, {n_conf:,} confident → {outdir}")

            elif kind == "run_err":
                _, err = msg
                self._pb.stop()
                self._lbl_run.config(text="✗  Error — see log")
                self._log(f"✗  Run error: {err}", warn=True)
                messagebox.showerror("Analysis error", err)

            elif kind == "heatmap_ok":
                _, png = msg
                self._pb.stop()
                self._log(f"Heatmap saved: {png}")
                messagebox.showinfo("Done", f"Heatmap saved:\n{png}")

            elif kind == "heatmap_empty":
                self._pb.stop()
                self._log("⚠  No genes passed heatmap filters — skipped.", warn=True)
                messagebox.showwarning("No genes",
                                       "No genes passed the current heatmap filters.\n"
                                       "Run the analysis first (Step ③).")

            elif kind == "heatmap_err":
                _, err = msg
                self._pb.stop()
                self._log(f"✗  Heatmap error: {err}", warn=True)
                messagebox.showerror("Heatmap error", err)

            elif kind == "batch_ok":
                _, subdir = msg
                self._pb.stop()
                self._log(f"Batch plots saved → {subdir}")
                messagebox.showinfo("Done", f"Gene plots saved:\n{subdir}")

            elif kind == "batch_err":
                _, err = msg
                self._pb.stop()
                self._log(f"✗  Batch plot error: {err}", warn=True)
                messagebox.showerror("Batch plot error", err)

        self.root.after(150, self._check_queue)

    # ------------------------------------------------------------------
    #   STATUS LOG
    # ------------------------------------------------------------------

    def _log(self, msg: str, warn: bool = False) -> None:
        ts   = datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}]  {msg}\n"
        w    = self._log_widget
        w.configure(state="normal")
        w.insert("end", line)
        if warn:
            # Tag last line orange
            start = w.index("end - 2 lines linestart")
            end   = w.index("end - 1 lines lineend")
            w.tag_add("warn", start, end)
            w.tag_config("warn", foreground="#B35000")
        w.see("end")
        w.configure(state="disabled")

    # ------------------------------------------------------------------
    #   GUARD HELPERS
    # ------------------------------------------------------------------

    def _check_ready(self) -> bool:
        if self.adata is None:
            messagebox.showwarning("No data",
                                   "Load an h5ad file first (Load Data section).")
            return False
        if self.tmeta is None:
            messagebox.showwarning("Tmeta required",
                                   "Please define ZT times first (Step ①).")
            return False
        return True


# =============================================================================
#   ENTRY POINT
# =============================================================================

def main() -> None:
    if not TIMESCAPE_OK:
        # Still open the GUI but show a warning at startup
        pass

    root = tk.Tk()
    app = TimeSCapeGUI(root)

    # Center the window on screen
    root.update_idletasks()
    sw = root.winfo_screenwidth()
    sh = root.winfo_screenheight()
    w, h = 1200, 760
    root.geometry(f"{w}x{h}+{(sw - w) // 2}+{(sh - h) // 2}")

    root.mainloop()


if __name__ == "__main__":
    main()
