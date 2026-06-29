"""
plot_utils.py — Shared Nature-journal style plotting primitives
===============================================================
Used by both figure_1.py and figure_2.py.

v2 changes (bug-fix pass)
--------------------------
- Figure canvas is now sized generously per-panel (not squeezed into a
  fixed 174 mm width) — small panels with long, thin satellite swaths
  were forcing tiny axes, which made stat boxes overflow and overlap
  with colorbars/titles.
- Per-class statistics are NO LONGER drawn as a floating text box on
  top of the map (this is what caused panel (b)/(c) numbers to get
  visually cut off and overlap the colorbar in earlier renders).
  They now render as a compact one-line caption directly below the
  panel, sized to the actual panel width.
- Basemap features are wrapped in `_SafeFeature`, required in
  network-restricted environments where Natural Earth shapefiles
  cannot be downloaded; without this wrapper cartopy raises during
  rendering and panels silently lose their land/ocean fill.

v2.1 — pcolormesh rendering
----------------------------
- plot_clm / plot_diff switched from layered scatter to single-pass
  pcolormesh(shading='nearest').  Invalid pixels are set to NaN
  (preserving 2-D grid structure) instead of being filtered out into
  1-D flat arrays.  This eliminates class-layering artefacts where
  cloudy classes (drawn last) would paint over clear classes, and
  guarantees every grid cell gets exactly one unambiguous colour.
- plot_rgb / plot_rgb_placeholder kept as scatter — single-pass
  rendering without class layering, where marker overlap is negligible
  at typical step=4 subsampling.

Exports
-------
apply_nature_style()
choose_projection(lat) → ccrs projection instance
make_geo_ax(fig, spec, projection) → ax
add_gridlines(ax)
panel_label(ax, letter)
panel_title(ax, text)
panel_caption(ax, text)
stats_caption_text(clm) -> str
agreement_caption_text(clm_a, clm_b) -> str
plot_rgb(ax, lat, lon, rgb, step)
plot_rgb_placeholder(ax, lat, lon, clm, step)
plot_clm(ax, lat, lon, clm, step)
plot_diff(ax, lat, lon, clm_a, clm_b, step)
add_clm_colorbar(fig, ax, ...)
add_diff_colorbar(fig, ax, sm, label)
get_extent(lat, lon, clm, step, pad)
subsample(arr, step)
"""

from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


# ─────────────────────────────────────────────────────────────────────────────
# Offline-safe basemap feature wrapper
# ─────────────────────────────────────────────────────────────────────────────
#
# In sandboxed / network-restricted environments, cartopy's NaturalEarth
# shapefiles cannot be downloaded on first use. cfeature.Feature defers the
# actual download until *render* time (inside fig.savefig), so wrapping the
# add_feature() CALL in try/except does NOT catch it. _SafeFeature fixes
# this by intercepting at geometry-fetch time and yielding nothing instead.

class _SafeFeature(cfeature.Feature):
    def __init__(self, feature, **kwargs):
        super().__init__(crs=feature.crs, **kwargs)
        self._feature = feature

    def geometries(self):
        try:
            yield from self._feature.geometries()
        except Exception:
            return  # shapefile unavailable — basemap layer simply omitted


def _safe(feat):
    return _SafeFeature(feat)


# ─────────────────────────────────────────────────────────────────────────────
# Style constants
# ─────────────────────────────────────────────────────────────────────────────

PANEL_WIDTH_IN  = 3.2
PANEL_HEIGHT_IN = 2.8
NATURE_DPI      = 300

FONT_FAMILY = ["Helvetica Neue", "Arial", "DejaVu Sans", "sans-serif"]

FIGURE_FACE = "#FBFAF7"
AX_FACE     = "#FFFFFF"
TEXT_COLOR  = "#252525"
MUTED_TEXT  = "#555555"

CLM_HEX = {
    0: "#2F5F9E",   # Cloudy             – deep blue
    1: "#8CC7DD",   # Prob. Cloudy       – light blue
    2: "#F2CF70",   # Prob. Clear        – soft amber
    3: "#D95F4A",   # Confident Clear    – muted red
}
CLM_LABEL = {
    0: "Cloudy",
    1: "Prob. Cloudy",
    2: "Prob. Clear",
    3: "Conf. Clear",
}
CLM_LABEL_SHORT = {0: "Cld", 1: "P.Cld", 2: "P.Clr", 3: "Clr"}

CLM_CMAP = mcolors.ListedColormap([CLM_HEX[v] for v in range(4)])
CLM_NORM = mcolors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], CLM_CMAP.N)

DELTA_CMAP = "RdBu_r"
DELTA_VMAX = 3.0

LAND_COLOR   = "#F1EEE8"
OCEAN_COLOR  = "#E7F0F6"
LAKE_COLOR   = "#DBE9F2"
COAST_COLOR  = "#5E5E5E"
BORDER_COLOR = "#AAAAAA"
GRID_COLOR   = "#BFC4C8"


# ─────────────────────────────────────────────────────────────────────────────
# rcParams
# ─────────────────────────────────────────────────────────────────────────────

def apply_nature_style() -> None:
    plt.rcParams.update({
        "font.family":          "sans-serif",
        "font.sans-serif":      FONT_FAMILY,
        "font.size":            9,
        "axes.titlesize":       9,
        "axes.labelsize":       8,
        "xtick.labelsize":      7,
        "ytick.labelsize":      7,
        "legend.fontsize":      8,
        "legend.framealpha":    0.9,
        "legend.edgecolor":     "#CCCCCC",
        "figure.facecolor":     FIGURE_FACE,
        "axes.facecolor":       AX_FACE,
        "axes.edgecolor":       "#BDBDBD",
        "axes.linewidth":       0.6,
        "axes.titlepad":        9,
        "text.color":           TEXT_COLOR,
        "xtick.major.width":    0.6,
        "ytick.major.width":    0.6,
        "xtick.major.size":     2.5,
        "ytick.major.size":     2.5,
        "figure.dpi":           150,
        "savefig.dpi":          NATURE_DPI,
        "savefig.bbox":         "tight",
        "savefig.facecolor":    FIGURE_FACE,
        "savefig.pad_inches":   0.05,
        "pdf.fonttype":         42,
        "ps.fonttype":          42,
        "lines.linewidth":      0.8,
        "patch.linewidth":      0.5,
    })


# ─────────────────────────────────────────────────────────────────────────────
# Array utilities
# ─────────────────────────────────────────────────────────────────────────────

def subsample(arr: np.ndarray, step: int) -> np.ndarray:
    return arr[::step, ::step]


def get_extent(
    lat: np.ndarray,
    lon: np.ndarray,
    clm: np.ndarray | None = None,
    step: int = 4,
    pad: float = 1.5,
) -> list | None:
    """Return [lon_min, lon_max, lat_min, lat_max] for map extent."""
    la = subsample(lat, step)
    lo = subsample(lon, step)
    if clm is not None:
        cl = subsample(clm, step)
        mask = np.isfinite(la) & np.isfinite(lo) & (cl >= 0)
    else:
        mask = np.isfinite(la) & np.isfinite(lo)
    if not mask.any():
        return None
    return [float(lo[mask].min()) - pad, float(lo[mask].max()) + pad,
            float(la[mask].min()) - pad, float(la[mask].max()) + pad]


def clamp_extent(extent: list, pad: float = 0.0) -> list:
    if extent is None:
        return extent
    return [
        max(extent[0], -180 + pad),
        min(extent[1],  180 - pad),
        max(extent[2],  -90 + pad),
        min(extent[3],   90 - pad),
    ]


# ─────────────────────────────────────────────────────────────────────────────
# Projection selector
# ─────────────────────────────────────────────────────────────────────────────

def choose_projection(lat: np.ndarray, polar_threshold: float = 0.5):
    """
    Choose an appropriate Cartopy projection based on the latitude distribution
    of the swath.

    Rules
    -----
    * If > ``polar_threshold`` (default 50 %) of valid pixels lie at
      |lat| ≥ 60°, use a polar-stereographic projection centred on
      the dominant pole (North or South).
    * Otherwise fall back to ``PlateCarree`` (standard equirectangular),
      which handles all mid-/low-latitude swaths without distortion.

    Parameters
    ----------
    lat : np.ndarray
        2-D (or 1-D) array of latitudes in decimal degrees.
    polar_threshold : float
        Fraction of valid pixels that must exceed |lat| ≥ 60° to trigger
        a polar projection.  Default 0.5 (50 %).

    Returns
    -------
    projection : cartopy.crs.CRS
        A Cartopy CRS instance ready to pass to ``add_subplot(projection=…)``.
    is_polar : bool
        True when a polar projection was selected.
    """
    valid = lat[np.isfinite(lat)]
    if len(valid) == 0:
        return ccrs.PlateCarree(), False

    polar_frac = np.mean(np.abs(valid) >= 60.0)
    if polar_frac <= polar_threshold:
        return ccrs.PlateCarree(), False

    # Dominant pole: whichever hemisphere contains more polar pixels
    n_north = np.sum(valid >= 60.0)
    n_south = np.sum(valid <= -60.0)
    central_lat = float(np.median(valid))

    if n_north >= n_south:
        return ccrs.NorthPolarStereo(central_longitude=0.0), True
    else:
        return ccrs.SouthPolarStereo(central_longitude=0.0), True


# ─────────────────────────────────────────────────────────────────────────────
# Axes / basemap helpers
# ─────────────────────────────────────────────────────────────────────────────

def _add_basemap_features(ax: plt.Axes) -> None:
    """Add land/ocean/coast/border features to any geo axis."""
    ax.set_facecolor(OCEAN_COLOR)
    ax.add_feature(_safe(cfeature.OCEAN),     facecolor=OCEAN_COLOR, zorder=0)
    ax.add_feature(_safe(cfeature.LAND),      facecolor=LAND_COLOR,  zorder=1)
    ax.add_feature(_safe(cfeature.LAKES),     facecolor=LAKE_COLOR,  zorder=1)
    ax.add_feature(_safe(cfeature.COASTLINE), linewidth=0.45,
                   edgecolor=COAST_COLOR, facecolor="none", zorder=4)
    ax.add_feature(_safe(cfeature.BORDERS),   linewidth=0.28,
                   edgecolor=BORDER_COLOR, facecolor="none", linestyle="--", zorder=4)


def make_geo_ax(fig: plt.Figure, spec, projection=None) -> plt.Axes:
    """Create a Cartopy geo axis with Nature-style basemap.

    Parameters
    ----------
    projection : cartopy.crs.CRS or None
        Pass the result of ``choose_projection(lat)`` to use a polar
        stereographic projection for high-latitude swaths.  Defaults to
        ``PlateCarree`` when None.
    """
    if projection is None:
        projection = ccrs.PlateCarree()
    ax = fig.add_subplot(spec, projection=projection)
    _add_basemap_features(ax)
    _style_map_frame(ax)
    # Tag so gridline helper knows which formatter to use
    ax._is_polar = not isinstance(projection, ccrs.PlateCarree)
    return ax


def make_geo_ax_with_caption(
    fig: plt.Figure,
    spec,
    caption_frac: float = 0.12,
    projection=None,
) -> plt.Axes:
    """
    Create a Cartopy map axis PLUS a dedicated caption sub-axes stacked
    directly beneath it, using a nested GridSpec inside `spec`.

    Parameters
    ----------
    projection : cartopy.crs.CRS or None
        Pass the result of ``choose_projection(lat)`` to activate polar
        stereographic projection for high-latitude swaths.  Defaults to
        ``PlateCarree`` when None.

    This guarantees the caption has its own reserved vertical band that
    never overlaps the map, its gridline tick labels, or the panel
    title — unlike a fixed `ax.text(y=-0.3, ...)` offset, which breaks
    down whenever panel aspect ratio or figure size changes.

    The created caption axes is stashed on the returned map axes as
    `ax._caption_ax`, so `panel_caption(ax, text)` can find it.
    """
    from matplotlib.gridspec import GridSpecFromSubplotSpec

    if projection is None:
        projection = ccrs.PlateCarree()

    inner = GridSpecFromSubplotSpec(
        2, 1, subplot_spec=spec,
        height_ratios=[1 - caption_frac, caption_frac],
        hspace=0.0,
    )
    ax     = fig.add_subplot(inner[0, 0], projection=projection)
    cap_ax = fig.add_subplot(inner[1, 0])
    cap_ax.set_axis_off()
    ax._caption_ax = cap_ax

    _add_basemap_features(ax)
    _style_map_frame(ax)
    # Tag so gridline helper knows which formatter to use
    ax._is_polar = not isinstance(projection, ccrs.PlateCarree)
    return ax


def _style_map_frame(ax: plt.Axes) -> None:
    """Give map panels a quiet frame without touching the plotted data."""
    for spine in ax.spines.values():
        spine.set_linewidth(0.55)
        spine.set_edgecolor("#B9B9B9")
    if "geo" in ax.spines:
        ax.spines["geo"].set_linewidth(0.55)
        ax.spines["geo"].set_edgecolor("#B9B9B9")


def add_gridlines(ax: plt.Axes) -> None:
    """Add gridlines with degree labels, adapting to polar vs. PlateCarree."""
    is_polar = getattr(ax, "_is_polar", False)
    gl = ax.gridlines(draw_labels=True, linewidth=0.25,
                      color=GRID_COLOR, alpha=0.65, linestyle=":")
    gl.top_labels   = False
    gl.right_labels = False
    if is_polar:
        # Polar-stereo projections do not support LONGITUDE_FORMATTER /
        # LATITUDE_FORMATTER — leave formatters as the Cartopy default
        # (plain degree numbers), which renders correctly.
        gl.xlabel_style = {"size": 6.4, "color": MUTED_TEXT}
        gl.ylabel_style = {"size": 6.4, "color": MUTED_TEXT}
    else:
        gl.xformatter   = LONGITUDE_FORMATTER
        gl.yformatter   = LATITUDE_FORMATTER
        gl.xlabel_style = {"size": 6.4, "color": MUTED_TEXT}
        gl.ylabel_style = {"size": 6.4, "color": MUTED_TEXT}


def panel_label(ax: plt.Axes, letter: str, fontsize: int = 10) -> None:
    """Bold panel label (a), (b), … outside top-left corner."""
    ax.text(-0.045, 1.075, f"({letter})",
            transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold",
            va="bottom", ha="right",
            color=TEXT_COLOR, fontfamily="sans-serif",
            bbox=dict(boxstyle="round,pad=0.16,rounding_size=0.03",
                      fc=FIGURE_FACE, ec="none", alpha=0.95))


def panel_title(ax: plt.Axes, text: str) -> None:
    """
    Panel title with headroom so the (letter) label above it doesn't
    collide. Left-aligned (loc='left'), so titles must be kept short
    enough to fit within the axes width — long titles will visually
    run into an adjacent colorbar's tick labels since the colorbar is
    drawn outside the axes bounding box. Keep titles under ~22 chars
    for 3-column layouts, ~28 chars for 2-column layouts.
    """
    ax.set_title(text, fontsize=9, fontweight="semibold",
                 color=TEXT_COLOR, pad=9, loc="left")


def panel_caption(ax: plt.Axes, text: str, fontsize: float = 7.5) -> None:
    """
    Render a one-line caption in a DEDICATED caption axes that was
    pre-attached to `ax` via make_geo_ax_with_caption(). Falls back to a
    below-axes text offset if no caption axes exists (e.g. for plain
    non-map axes), but the dedicated-axes path is what avoids overlap
    with cartopy's own tick labels and prevents right-edge clipping on
    long strings.
    """
    cap_ax = getattr(ax, "_caption_ax", None)
    if cap_ax is not None:
        cap_ax.clear()
        cap_ax.set_axis_off()
        cap_ax.text(0.0, 0.95, text, transform=cap_ax.transAxes,
                    fontsize=fontsize, color=MUTED_TEXT,
                    va="top", ha="left", linespacing=1.3,
                    wrap=True)
        return

    # Fallback: no caption axes attached — use a conservative offset
    ax.text(0.0, -0.34, text,
            transform=ax.transAxes,
            fontsize=fontsize, color=MUTED_TEXT,
            va="top", ha="left", linespacing=1.4)


def stats_caption_text(clm: np.ndarray) -> str:
    """Compact one-line class-distribution string for panel_caption."""
    total = int((clm >= 0).sum())
    if total == 0:
        return "No valid pixels"
    parts = [f"{CLM_LABEL_SHORT[v]} {100.0*(clm==v).sum()/total:.0f}%"
             for v in range(4)]
    return "  ·  ".join(parts)


def agreement_caption_text(clm_a: np.ndarray, clm_b: np.ndarray) -> str:
    """Compact one-line agreement/POD/FAR string for panel_caption."""
    mask = (clm_a >= 0) & (clm_b >= 0)
    if not mask.any():
        return "No overlap"
    agree_pct = 100 * np.mean(clm_a[mask] == clm_b[mask])
    a_cloud = clm_a[mask] <= 1
    b_cloud = clm_b[mask] <= 1
    pod = (a_cloud & b_cloud).sum() / (b_cloud.sum() + 1e-9) * 100
    far = (a_cloud & ~b_cloud).sum() / (a_cloud.sum() + 1e-9) * 100
    return f"Agreement {agree_pct:.1f}%  ·  POD {pod:.1f}%  ·  FAR {far:.1f}%"


# ─────────────────────────────────────────────────────────────────────────────
# Data renderers
# ─────────────────────────────────────────────────────────────────────────────

def plot_rgb(
    ax: plt.Axes,
    lat: np.ndarray,
    lon: np.ndarray,
    rgb: np.ndarray,
    step: int = 4,
    max_pts: int = 250_000,
) -> None:
    """
    Scatter RGB true-colour pixels onto geo axis.

    rgb MUST share the same (rows, cols) axis order as lat/lon, i.e.
    rgb.shape[:2] == lat.shape. A shape mismatch here is the #1 cause
    of "RGB looks like random noise" — fix the mismatch at the loader
    (io_mersi.py), not by reshaping here.
    """
    if rgb.shape[:2] != lat.shape:
        _no_data_text(
            ax, f"RGB/geo shape mismatch\n{rgb.shape[:2]} vs {lat.shape}")
        return

    la = subsample(lat, step)
    lo = subsample(lon, step)
    rg = subsample(rgb, step)
    mask = np.isfinite(la) & np.isfinite(lo)
    if not mask.any():
        _no_data_text(ax, "No L1b\navailable")
        return
    laf, lof, rgf = la[mask], lo[mask], rg[mask]
    if len(laf) > max_pts:
        rng = np.random.default_rng(0)
        idx = rng.choice(len(laf), max_pts, replace=False)
        laf, lof, rgf = laf[idx], lof[idx], rgf[idx]
    ax.scatter(lof, laf, c=np.clip(rgf, 0, 1), s=1.2, linewidths=0,
               transform=ccrs.PlateCarree(), zorder=3,
               edgecolors="none", rasterized=True)


def plot_rgb_placeholder(
    ax: plt.Axes,
    lat: np.ndarray,
    lon: np.ndarray,
    clm: np.ndarray,
    step: int = 4,
) -> None:
    """Grey footprint when L1B is unavailable."""
    la = subsample(lat, step)
    lo = subsample(lon, step)
    cl = subsample(clm, step)
    mask = np.isfinite(la) & np.isfinite(lo) & (cl >= 0)
    ax.scatter(lo[mask], la[mask], c="#CCCCCC", s=1.2,
               linewidths=0, transform=ccrs.PlateCarree(),
               zorder=3, alpha=0.6, rasterized=True)
    _no_data_text(ax, "No L1b file\n(CLM footprint)")


def plot_clm(
    ax: plt.Axes,
    lat: np.ndarray,
    lon: np.ndarray,
    clm: np.ndarray,
    step: int = 4,
) -> None:
    """
    Render CLM via pcolormesh — one-shot, no class-layering artefacts.

    Keeps the 2-D grid intact (invalid pixels → NaN) so pcolormesh
    can assign every cell a unique colour without any "later class
    paints over earlier class" scatter-ordering problem.
    """
    la = subsample(lat, step)
    lo = subsample(lon, step)
    cl = subsample(clm, step)
    mask = np.isfinite(la) & np.isfinite(lo) & (cl >= 0)
    if not mask.any():
        _no_data_text(ax, "No valid CLM pixels")
        return
    cl_masked = np.where(mask, cl.astype(np.float64), np.nan)
    ax.pcolormesh(lo, la, cl_masked, cmap=CLM_CMAP, norm=CLM_NORM,
                  transform=ccrs.PlateCarree(), shading="nearest",
                  rasterized=True, zorder=3)


def plot_diff(
    ax: plt.Axes,
    lat: np.ndarray,
    lon: np.ndarray,
    clm_a: np.ndarray,
    clm_b: np.ndarray,
    step: int = 4,
    vmax: float = DELTA_VMAX,
):
    """
    Render class difference (clm_a − clm_b) via pcolormesh.

    Keeps the 2-D grid intact with NaN for invalid / non-overlapping
    pixels so every cell gets ONE unambiguous colour.  Returns
    (ScalarMappable, mask, diff) for colorbar + stats reuse.
    """
    la = subsample(lat, step)
    lo = subsample(lon, step)
    ca = subsample(clm_a, step)
    cb = subsample(clm_b, step)
    mask = np.isfinite(la) & np.isfinite(lo) & (ca >= 0) & (cb >= 0)

    diff   = np.where(mask, ca.astype(np.float64) - cb.astype(np.float64), np.nan)
    d_cmap = plt.get_cmap(DELTA_CMAP)
    d_norm = mcolors.Normalize(vmin=-vmax, vmax=vmax)

    if mask.any():
        ax.pcolormesh(lo, la, diff, cmap=d_cmap, norm=d_norm,
                      transform=ccrs.PlateCarree(), shading="nearest",
                      rasterized=True, zorder=3)
    else:
        _no_data_text(ax, "No overlap")

    sm = plt.cm.ScalarMappable(cmap=d_cmap, norm=d_norm)
    sm.set_array([])
    return sm, mask, diff


# ─────────────────────────────────────────────────────────────────────────────
# Colorbars
# ─────────────────────────────────────────────────────────────────────────────

def add_clm_colorbar(
    fig: plt.Figure,
    ax: plt.Axes,
    shrink: float = 0.74,
    pad: float = 0.055,
) -> None:
    sm = plt.cm.ScalarMappable(cmap=CLM_CMAP, norm=CLM_NORM)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical",
                        shrink=shrink, pad=pad, aspect=24,
                        ticks=[0, 1, 2, 3])
    cbar.ax.set_yticklabels(
        [CLM_LABEL[v] for v in range(4)], fontsize=7.5)
    _style_colorbar(cbar)
    return cbar


def add_diff_colorbar(
    fig: plt.Figure,
    ax: plt.Axes,
    sm,
    label: str = "Δ class",
    shrink: float = 0.74,
    pad: float = 0.055,
) -> None:
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical",
                        shrink=shrink, pad=pad, aspect=24)
    cbar.set_label(label, fontsize=7.5, labelpad=4)
    _style_colorbar(cbar)
    return cbar


def _style_colorbar(cbar) -> None:
    cbar.ax.tick_params(labelsize=7, length=2.2, width=0.45,
                        colors=MUTED_TEXT, pad=2)
    cbar.outline.set_linewidth(0.45)
    cbar.outline.set_edgecolor("#AFAFAF")
    cbar.ax.yaxis.label.set_color(MUTED_TEXT)


def source_label(ax: plt.Axes, text: str, fontsize: float = 6.5) -> None:
    """Small grey source annotation inside the bottom-right of the map."""
    ax.text(0.985, 0.02, text,
            transform=ax.transAxes, va="bottom", ha="right",
            fontsize=fontsize, color="#777777",
            bbox=dict(fc="white", ec="none", pad=2, alpha=0.75))


def _no_data_text(ax: plt.Axes, msg: str) -> None:
    ax.text(0.5, 0.5, msg, transform=ax.transAxes,
            ha="center", va="center", fontsize=8, color=MUTED_TEXT,
            bbox=dict(boxstyle="round,pad=0.35,rounding_size=0.08",
                      fc="white", ec="#D8D8D8", lw=0.5, alpha=0.9))


# ─────────────────────────────────────────────────────────────────────────────
# Figure saving
# ─────────────────────────────────────────────────────────────────────────────

def save_figure(fig: plt.Figure, output: str) -> None:
    import os
    fig.savefig(output, dpi=NATURE_DPI, bbox_inches="tight",
                facecolor=FIGURE_FACE, pad_inches=0.05,
                pil_kwargs={"compression": 6})
    plt.close(fig)
    print(f"[SAVE] {output}  ({os.path.getsize(output)/1024:.0f} KB)")
