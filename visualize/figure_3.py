"""
figure_3.py — MERSI-II IR-enhanced cloud comparison figure  (6-panel)
======================================================================
Panel layout  (2 rows × 3 cols):
  (a) MERSI RGB true colour (L1B bands B1/B2/B3, or grey footprint)
  (b) IR-enhanced cloud image (L1B 10.8 µm, channel 24 of EV_1KM_Emissive)
  (c) MYD35 cloud mask  (resampled to MERSI grid)
  (d) MERSI recalibration CLM
  (e) Diff map: MYD35 − MERSI recal  (overlap region only)
  (f) Confusion matrix + validation stats  (recal vs MYD35)

Design notes
------------
* IR brightness-temperature (BT) in panel (b) is rendered with a
  "enhanced IR" palette: cold cloud tops (low BT) → white/light blue,
  warm surface (high BT) → grey/black.  Colourmap is inverted so bright
  = cold = cloudy, matching standard satellite imagery convention.
* Resampling (MYD35 → MERSI grid) is delegated to
  io_myd35.resample_to_mersi_grid(), matching Figure 2's approach.
* Diff panel (e) shows MYD35 − MERSI_recal (positive = MYD35 sees more
  cloud than MERSI recal) for intuitive sign convention.
* The confusion-matrix panel reuses figure_2._draw_confusion_panel().
* All 5 geo panels share the same Cartopy projection determined by
  choose_projection(lat), adapting to polar or mid-latitude swaths.

Usage
-----
# From pre-loaded arrays (typical call from run_validation.py):
from figure_3 import make_figure3

stats = make_figure3(
    mersi_lat       = lat,            # (H,W) float64
    mersi_lon       = lon,            # (H,W) float64
    mersi_rgb       = rgb,            # (H,W,3) float32  or  None
    mersi_bt108     = bt108,          # (H,W) float32  or  None
    recal_clm       = recal_clm,      # (H,W) int32
    myd35_data      = myd35_dict,     # from io_myd35.load_best_myd35_for_mersi()
    output          = "figure3.png",
    mersi_date_str  = "2023-06-06  14:40 UTC",
    step            = 4,
)

# From HDF5 / HDF files:
from figure_3 import make_figure3_from_files
make_figure3_from_files(
    recal_path   = "..._CLM_CLA_recal.h5",
    myd35_dirs   = ["/data/myd35/"],
    output       = "figure3.png",
    mersi_root   = "/data/Data_yuq/mersi",
)
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec

from plot_utils import (
    apply_nature_style, PANEL_WIDTH_IN, PANEL_HEIGHT_IN,
    FIGURE_FACE, TEXT_COLOR, MUTED_TEXT,
    make_geo_ax_with_caption, add_gridlines, panel_label, panel_title,
    panel_caption,
    plot_rgb, plot_rgb_placeholder, plot_clm, plot_diff,
    add_clm_colorbar, add_diff_colorbar,
    stats_caption_text, agreement_caption_text,
    get_extent, subsample, save_figure,
    CLM_LABEL_SHORT, choose_projection,
)
from io_mersi import (
    load_clm_hdf5, load_mersi_l1b, load_mersi_bt108,
    find_l1b_for_clm, parse_mersi_datetime, print_clm_distribution,
)
from io_myd35 import load_best_myd35_for_mersi

# Reuse the validation helpers from figure_2
from figure_2 import (
    compute_overlap_extent,
    overlap_mask_on_mersi_grid,
    compute_validation_stats,
    _draw_confusion_panel,
)


# ─────────────────────────────────────────────────────────────────────────────
# IR-enhanced BT rendering
# ─────────────────────────────────────────────────────────────────────────────

# Standard enhanced-IR colourmap limits (Kelvin)
# Cold cloud tops ≈ 200–230 K → white end; warm surface ≈ 300–310 K → dark end
BT_VMIN_K = 200.0
BT_VMAX_K = 310.0

# Colourmap: inverted grey-blue so cold = bright (cloud), warm = dark (surface).
# Uses matplotlib's gist_gray reversed + slight blue tint at the cold end.
IR_CMAP = "gray_r"   # simple inverted grey; swap to "Blues_r" for stronger contrast


def plot_ir_bt(
    ax: plt.Axes,
    lat: np.ndarray,
    lon: np.ndarray,
    bt: np.ndarray,                   # (H, W) float32/float64 in Kelvin
    step: int = 4,
    vmin: float = BT_VMIN_K,
    vmax: float = BT_VMAX_K,
):
    """
    Render MERSI 10.8 µm brightness temperature as an enhanced-IR image.

    Cold (low BT) → white/light;  warm (high BT) → dark grey.
    Uses pcolormesh on the subsampled 2-D grid — same approach as plot_clm.

    Returns
    -------
    ScalarMappable for the colourbar, or None if no valid data.
    """
    import cartopy.crs as ccrs

    la = subsample(lat, step)
    lo = subsample(lon, step)
    bt_s = subsample(bt, step)

    mask = np.isfinite(la) & np.isfinite(lo) & np.isfinite(bt_s)
    if not mask.any():
        ax.text(0.5, 0.5, "No BT data\navailable",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=8, color=MUTED_TEXT,
                bbox=dict(boxstyle="round,pad=0.35,rounding_size=0.08",
                          fc="white", ec="#D8D8D8", lw=0.5, alpha=0.9))
        return None

    bt_plot = np.where(mask, bt_s.astype(np.float64), np.nan)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    mesh = ax.pcolormesh(
        lo, la, bt_plot,
        cmap=IR_CMAP, norm=norm,
        transform=ccrs.PlateCarree(),
        shading="nearest", rasterized=True, zorder=3,
    )

    sm = plt.cm.ScalarMappable(cmap=plt.get_cmap(IR_CMAP), norm=norm)
    sm.set_array([])
    return sm


def add_bt_colorbar(
    fig: plt.Figure,
    ax: plt.Axes,
    sm,
    shrink: float = 0.74,
    pad: float = 0.055,
) -> None:
    """Add colourbar for BT panel (b)."""
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical",
                        shrink=shrink, pad=pad, aspect=24)
    cbar.set_label("BT (K)", fontsize=7.5, labelpad=4)
    cbar.ax.tick_params(labelsize=7, length=2.2, width=0.45,
                        colors=MUTED_TEXT, pad=2)
    cbar.outline.set_linewidth(0.45)
    cbar.outline.set_edgecolor("#AFAFAF")
    cbar.ax.yaxis.label.set_color(MUTED_TEXT)
    return cbar


def _bt_caption(bt: np.ndarray) -> str:
    """One-line brightness-temperature summary for panel_caption."""
    valid = bt[np.isfinite(bt)]
    if len(valid) == 0:
        return "No valid BT pixels"
    return (f"BT  min {valid.min():.1f} K  ·  "
            f"mean {valid.mean():.1f} K  ·  "
            f"max {valid.max():.1f} K")


# ─────────────────────────────────────────────────────────────────────────────
# Diff panel (MYD35 − MERSI recal)
# ─────────────────────────────────────────────────────────────────────────────

def _diff_panel_myd_minus_recal(
    ax,
    fig,
    mersi_lat:       np.ndarray,
    mersi_lon:       np.ndarray,
    recal_clm:       np.ndarray,
    myd35_resampled: np.ndarray,
    step: int = 4,
) -> dict:
    """
    Plot (MYD35_resampled − recal_clm) on MERSI grid, restricted to overlap.

    Positive values → MYD35 sees more cloud than MERSI recal.
    Returns validation stats dict.
    """
    ov_mask  = overlap_mask_on_mersi_grid(mersi_lat, mersi_lon, myd35_resampled)
    myd_ov   = np.where(ov_mask, myd35_resampled, -1)
    recal_ov = np.where(ov_mask, recal_clm,       -1)

    # plot_diff(a, b) renders (a − b); we want MYD35 − recal so pass myd first
    sm, _, _ = plot_diff(ax, mersi_lat, mersi_lon, myd_ov, recal_ov, step=step)
    add_diff_colorbar(fig, ax, sm, label="MYD35 − MERSI recal (class)")

    stats = compute_validation_stats(recal_clm, myd35_resampled, label="MERSI recal")
    if stats:
        cap = (f"Agreement {stats['agree_pct']:.1f}%  ·  "
               f"POD {stats['pod']:.1f}%  ·  FAR {stats['far']:.1f}%  ·  "
               f"HSS {stats['hss']:.3f}")
        panel_caption(ax, cap)
    return stats


# ─────────────────────────────────────────────────────────────────────────────
# Core figure builder
# ─────────────────────────────────────────────────────────────────────────────

def _build_figure3(
    mersi_lat:      np.ndarray,
    mersi_lon:      np.ndarray,
    mersi_rgb:      np.ndarray | None,   # (H,W,3) float32 or None
    mersi_bt108:    np.ndarray | None,   # (H,W) float32 in K or None
    recal_clm:      np.ndarray,          # (H,W) int32
    myd35_data:     dict,                # from load_best_myd35_for_mersi()
    mersi_date_str: str,
    output:         str,
    step:           int = 4,
) -> dict:
    """Render and save Figure 3.  Returns validation stats dict."""
    apply_nature_style()

    projection, is_polar = choose_projection(mersi_lat)

    myd35_resampled = myd35_data["clm_resampled"]
    dt_min = myd35_data.get("dt_diff_min", 0.0)

    # Extent: intersection of MERSI and MYD35 footprints
    overlap_extent = compute_overlap_extent(
        mersi_lat, mersi_lon, recal_clm,
        myd35_data["lat"], myd35_data["lon"], myd35_data["clm_native"],
        step=step)
    if overlap_extent is None:
        overlap_extent = get_extent(mersi_lat, mersi_lon, recal_clm, step=step)

    # ── Canvas ───────────────────────────────────────────────────────
    fig_w = PANEL_WIDTH_IN * 3 + 1.90
    fig_h = PANEL_HEIGHT_IN * 2 + 0.95
    fig   = plt.figure(figsize=(fig_w, fig_h), facecolor=FIGURE_FACE)

    gs = GridSpec(2, 3, figure=fig,
                  left=0.045, right=0.94,
                  top=0.875, bottom=0.055,
                  wspace=0.44, hspace=0.42)

    panel_specs = [gs[0, 0], gs[0, 1], gs[0, 2],
                   gs[1, 0], gs[1, 1], gs[1, 2]]
    letters = ["a", "b", "c", "d", "e", "f"]
    titles  = [
        "MERSI RGB true colour",
        "IR enhanced  (10.8 µm)",
        "MYD35 CLM (truth)",
        "MERSI recal CLM",
        "MYD35 − Recal diff",
        "Confusion matrix",
    ]

    # 5 geo axes + 1 plain axis for confusion matrix
    axs_geo = [
        make_geo_ax_with_caption(fig, panel_specs[i], projection=projection)
        for i in range(5)
    ]
    ax_conf = fig.add_subplot(panel_specs[5])

    # ── (a) RGB ─────────────────────────────────────────────────────
    ax = axs_geo[0]
    if mersi_rgb is not None:
        plot_rgb(ax, mersi_lat, mersi_lon, mersi_rgb, step=step)
    else:
        plot_rgb_placeholder(ax, mersi_lat, mersi_lon, recal_clm, step=step)
    if overlap_extent:
        ax.set_extent(overlap_extent)
    add_gridlines(ax)

    # ── (b) IR 10.8 µm BT ───────────────────────────────────────────
    ax = axs_geo[1]
    if mersi_bt108 is not None:
        sm_bt = plot_ir_bt(ax, mersi_lat, mersi_lon, mersi_bt108, step=step)
        if sm_bt is not None:
            add_bt_colorbar(fig, ax, sm_bt)
            panel_caption(ax, _bt_caption(mersi_bt108))
    else:
        ax.text(0.5, 0.5, "10.8 µm BT\nnot available",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=8, color=MUTED_TEXT,
                bbox=dict(boxstyle="round,pad=0.35,rounding_size=0.08",
                          fc="white", ec="#D8D8D8", lw=0.5, alpha=0.9))
    if overlap_extent:
        ax.set_extent(overlap_extent)
    add_gridlines(ax)

    # ── (c) MYD35 CLM on native grid ────────────────────────────────
    ax = axs_geo[2]
    plot_clm(ax, myd35_data["lat"], myd35_data["lon"],
             myd35_data["clm_native"], step=step)
    if overlap_extent:
        ax.set_extent(overlap_extent)
    add_gridlines(ax)
    add_clm_colorbar(fig, ax)
    dt_str = f"Δt = {dt_min:.1f} min" if dt_min else ""
    panel_caption(ax, stats_caption_text(myd35_data["clm_native"]) +
                  (f"   ({dt_str})" if dt_str else ""))

    # ── (d) MERSI recal CLM (overlap region) ────────────────────────
    ax = axs_geo[3]
    ov_mask  = overlap_mask_on_mersi_grid(mersi_lat, mersi_lon, myd35_resampled)
    recal_ov = np.where(ov_mask, recal_clm, -1)
    plot_clm(ax, mersi_lat, mersi_lon, recal_ov, step=step)
    if overlap_extent:
        ax.set_extent(overlap_extent)
    add_gridlines(ax)
    add_clm_colorbar(fig, ax)
    panel_caption(ax, stats_caption_text(recal_ov))

    # ── (e) MYD35 − Recal diff (overlap only) ───────────────────────
    ax = axs_geo[4]
    stats_recal = _diff_panel_myd_minus_recal(
        ax, fig, mersi_lat, mersi_lon, recal_clm, myd35_resampled, step=step)
    if overlap_extent:
        ax.set_extent(overlap_extent)
    add_gridlines(ax)

    # ── (f) Confusion matrix ─────────────────────────────────────────
    # Pass stats_onboard=None so only recal vs MYD35 is shown (one table)
    _draw_confusion_panel(ax_conf, stats_recal, stats_onboard=None)

    # ── Decorations ──────────────────────────────────────────────────
    all_axes = axs_geo + [ax_conf]
    for ax, ltr, ttl in zip(all_axes, letters, titles):
        panel_title(ax, ttl)
        panel_label(ax, ltr)

    lat_c = float(np.nanmedian(mersi_lat))
    lon_c = float(np.nanmedian(mersi_lon))
    myd_src = myd35_data.get("source", "")
    myd_basename = myd_src.split("/")[-1] if myd_src else ""
    fig.suptitle(
        f"FY-3D MERSI-II  IR-enhanced  vs  MYD35 (truth)   "
        f"MERSI: {mersi_date_str}   "
        f"centre {lat_c:.1f}°N {lon_c:.1f}°E"
        + (f"\nMYD35: {myd_basename}" if myd_basename else ""),
        fontsize=10.1, fontweight="semibold", color=TEXT_COLOR, y=0.975)

    save_figure(fig, output)
    return {"recal": stats_recal}


# ─────────────────────────────────────────────────────────────────────────────
# Public entry points
# ─────────────────────────────────────────────────────────────────────────────

def make_figure3(
    mersi_lat:      np.ndarray,
    mersi_lon:      np.ndarray,
    recal_clm:      np.ndarray,
    myd35_data:     dict,
    mersi_rgb:      np.ndarray | None = None,
    mersi_bt108:    np.ndarray | None = None,
    output:         str               = "figure3.png",
    mersi_date_str: str               = "",
    step:           int               = 4,
) -> dict:
    """
    Build Figure 3 from pre-loaded arrays.

    Parameters
    ----------
    mersi_lat, mersi_lon : (H,W) float64
        MERSI geolocation on 1-km grid.
    recal_clm : (H,W) int32
        MERSI recalibration CLM (0–3, -1 = invalid).
    myd35_data : dict
        Output of io_myd35.load_best_myd35_for_mersi().
        Required keys: clm_native, clm_resampled, lat, lon, dt, dt_diff_min.
    mersi_rgb : (H,W,3) float32 or None
        True-colour RGB.  None → grey footprint placeholder in panel (a).
    mersi_bt108 : (H,W) float32 or None
        Brightness temperature at 10.8 µm in Kelvin.
        None → "not available" text in panel (b).
        Load via io_mersi.load_mersi_bt108().
    output : str
        Output PNG path.
    mersi_date_str : str
        Human-readable date/time string for the suptitle.
    step : int
        Pixel subsampling stride (4 = every 4th pixel in each direction).

    Returns
    -------
    dict with key 'recal': validation stats between recal_clm and MYD35.
    """
    print(f"\n[FIG3] Building IR-enhanced validation figure → {output}")
    print_clm_distribution("MERSI recal CLM",     recal_clm)
    print_clm_distribution("MYD35 CLM (native)",  myd35_data["clm_native"])
    print_clm_distribution("MYD35 CLM (on MERSI grid)",
                            myd35_data["clm_resampled"])
    if mersi_bt108 is not None:
        valid_bt = mersi_bt108[np.isfinite(mersi_bt108)]
        if len(valid_bt):
            print(f"[FIG3] BT 10.8µm — "
                  f"min {valid_bt.min():.1f} K  mean {valid_bt.mean():.1f} K  "
                  f"max {valid_bt.max():.1f} K")

    return _build_figure3(
        mersi_lat=mersi_lat,
        mersi_lon=mersi_lon,
        mersi_rgb=mersi_rgb,
        mersi_bt108=mersi_bt108,
        recal_clm=recal_clm,
        myd35_data=myd35_data,
        mersi_date_str=mersi_date_str,
        output=output,
        step=step,
    )


def make_figure3_from_files(
    recal_path:      str,
    myd35_dirs:      list[str],
    output:          str   = "figure3.png",
    mersi_root:      str   = "/data/Data_yuq/mersi",
    step:            int   = 4,
    time_window_min: int   = 15,
    min_overlap:     float = 0.05,
) -> dict | None:
    """
    Build Figure 3 by reading all input files from disk.

    Parameters
    ----------
    recal_path : str
        Path to the recalibration CLM HDF5 file.
    myd35_dirs : list[str]
        Directories to search for matching MYD35 granules.
    output : str
        Output PNG path.
    mersi_root : str
        Root directory for MERSI L1B HDF files.
    step, time_window_min, min_overlap :
        Passed through to their respective loaders.

    Returns
    -------
    dict (validation stats) or None on failure.
    """
    import re
    import os

    # ── Load recal CLM ───────────────────────────────────────────────
    recal_data = load_clm_hdf5(recal_path)
    if recal_data is None:
        print("[ERROR] Could not load recal CLM — Figure 3 skipped.")
        return None

    lat = recal_data["lat"]
    lon = recal_data["lon"]
    recal_clm = recal_data["clm"]

    # ── Date string from filename ────────────────────────────────────
    m = re.search(r'(\d{8})_(\d{4})', os.path.basename(recal_path))
    date_str = ""
    if m:
        d, t = m.group(1), m.group(2)
        date_str = f"{d[:4]}-{d[4:6]}-{d[6:8]}  {t[:2]}:{t[2:]} UTC"

    # ── Load L1B (RGB + BT 10.8 µm) ─────────────────────────────────
    l1b_path = find_l1b_for_clm(recal_path, mersi_root)
    l1b_data = load_mersi_l1b(l1b_path) if l1b_path else None
    rgb = l1b_data["rgb"] if l1b_data else None

    bt108 = None
    if l1b_path:
        bt108 = load_mersi_bt108(l1b_path)

    # ── MYD35 matching ───────────────────────────────────────────────
    mersi_dt = parse_mersi_datetime(recal_path)
    if mersi_dt is None:
        print("[ERROR] Cannot parse MERSI datetime — Figure 3 skipped.")
        return None

    myd35_data = load_best_myd35_for_mersi(
        mersi_lat=lat,
        mersi_lon=lon,
        mersi_dt=mersi_dt,
        search_dirs=myd35_dirs,
        time_window_min=time_window_min,
        min_overlap=min_overlap,
    )
    if myd35_data is None:
        print("[WARN] No matching MYD35 granule — Figure 3 skipped.")
        return None

    return make_figure3(
        mersi_lat=lat,
        mersi_lon=lon,
        recal_clm=recal_clm,
        myd35_data=myd35_data,
        mersi_rgb=rgb,
        mersi_bt108=bt108,
        output=output,
        mersi_date_str=date_str,
        step=step,
    )
