"""
Publication figure style (IEEE).

Two layouts of the same four metrics:
  * two_column_figure: 2x2 grid, 7.16 in wide (full page width, figure*)
  * one_column_figure: 4x1 stack with a shared x axis, 3.5 in wide (one column)
Include the files at their native width (\\textwidth / \\columnwidth) so that the
font sizes are preserved.

Output: PDF (TrueType fonts embedded) + SVG (text kept as text) + PNG 300 dpi preview.
Colors: blue / orange (colorblind-safe). Series also differ by marker and line
style, so they stay readable in black-and-white print.
Model = lines, measurement = markers without lines.
"""

import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

DOUBLE_COL_IN = 7.16
SINGLE_COL_IN = 3.5

R_COLORS = {1: "#2a78d6", 3: "#eb6834"}
R_MARKERS = {1: "o", 3: "s"}
R_LINESTYLES = {1: "-", 3: "--"}

INK = "#1a1a1a"
GRID = "#e1e0d9"
AXIS = "#8a897f"

# (metric key, y label in the 2x2 layout, y label in the one-column layout)
PANELS = [
    ("mean_delay_ms", "Mean delay [ms]", "Mean delay [ms]"),
    ("jitter_ms", "Jitter [ms]", "Jitter [ms]"),
    ("loss_prob", "Packet loss probability", "Packet loss prob."),
    ("pctl_delay_ms", None, None),  # label depends on the percentile
]


def apply(base_size=8):
    """base_size: axis labels; tick labels and legend are 1 pt smaller (0.5 pt
    when base_size >= 9)."""
    small = base_size - (0.5 if base_size >= 9 else 1)
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": base_size,
        "axes.labelsize": base_size,
        "axes.titlesize": base_size,
        "xtick.labelsize": small,
        "ytick.labelsize": small,
        "legend.fontsize": small,
        "axes.linewidth": 0.6,
        "axes.edgecolor": AXIS,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": AXIS,
        "ytick.color": AXIS,
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.4,
        "lines.linewidth": 1.2,
        "lines.markersize": 4.5,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    })


def pctl_label(pctl, short=False):
    if short:
        return f"P{pctl:g} delay [ms]"
    n = int(pctl) if float(pctl).is_integer() else None
    if n is None:
        return f"{pctl:g}th percentile delay [ms]"
    suf = "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{n}{suf} percentile delay [ms]"


def two_column_figure(xlabel, pctl, model=None, meas=None, extra=None, xticks=None,
                      model_markers=False, xscale="linear", xtick_fmt=None):
    """2x2 layout: (a) mean delay, (b) jitter, (c) packet loss, (d) delay percentile.

    model : {R: (xs, [metrics])}  -> lines "Analytical, R=..."
    meas  : {R: (xs, [metrics])}  -> markers "Measured, R=..."
    extra : [(label, color, xs, [metrics])] -> additional (dotted) model lines
    metrics: dicts with keys mean_delay_ms, jitter_ms, loss_prob, pctl_delay_ms.
    xscale: "linear" or "log2" (x values doubling -> equal spacing).
    """
    apply(8)
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COL_IN, 4.6))
    for k, (ax, (key, ylabel, _)) in enumerate(zip(axes.flat, PANELS)):
        _draw_panel(ax, key, ylabel or pctl_label(pctl), xlabel, model, meas, extra,
                    xticks, model_markers, xscale, xtick_fmt)
        ax.set_title(f"({'abcd'[k]})", loc="left", pad=3)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=min(len(labels), 4),
               bbox_to_anchor=(0.5, 1.0), handlelength=2.4, columnspacing=1.6)
    fig.tight_layout(rect=(0, 0, 1, 0.95), h_pad=1.2, w_pad=1.6)
    return fig


def one_column_figure(xlabel, pctl, model=None, meas=None, extra=None, xticks=None,
                      model_markers=False, xscale="linear", xtick_fmt=None):
    """4x1 layout, one column wide (3.5 in), shared x axis, larger fonts
    (9 pt axis labels, 8.5 pt tick labels and legend)."""
    apply(9)
    fig, axes = plt.subplots(4, 1, figsize=(SINGLE_COL_IN, 7.0), sharex=True)
    for k, (ax, (key, _, ylabel)) in enumerate(zip(axes, PANELS)):
        _draw_panel(ax, key, ylabel or pctl_label(pctl, short=True), xlabel, model, meas,
                    extra, xticks, model_markers, xscale, xtick_fmt)
        ax.set_title(f"({'abcd'[k]})", loc="right", pad=2)
        if k < 3:
            ax.set_xlabel("")
        ax.yaxis.set_label_coords(-0.13, 0.5)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2,
               bbox_to_anchor=(0.5, 1.0), handlelength=2.2, columnspacing=1.0,
               handletextpad=0.5)
    fig.tight_layout(rect=(0, 0, 1, 0.925), h_pad=0.3)
    return fig


def _draw_panel(ax, key, ylabel, xlabel, model, meas, extra, xticks,
                model_markers, xscale, xtick_fmt):
    for R, (xs, ms) in (model or {}).items():
        ax.plot(xs, [m[key] for m in ms], color=R_COLORS.get(R, INK),
                linestyle=R_LINESTYLES.get(R, "-"),
                marker=R_MARKERS.get(R, "o") if model_markers else None,
                markerfacecolor="white", markersize=3.5,
                label=f"Analytical, $R={R}$", zorder=2)
    for label, color, xs, ms in (extra or []):
        ax.plot(xs, [m[key] for m in ms], color=color, linestyle=":",
                linewidth=1.2, label=label, zorder=2)
    for R, (xs, ms) in (meas or {}).items():
        ax.plot(xs, [m[key] for m in ms], linestyle="none",
                marker=R_MARKERS.get(R, "o"), color=R_COLORS.get(R, INK),
                markeredgecolor="white", markeredgewidth=0.6, markersize=5.5,
                label=f"Measured, $R={R}$", zorder=3)
    ax.set_ylim(bottom=0)
    if xscale == "log2":
        ax.set_xscale("log", base=2)
        ax.minorticks_off()
    if xticks is not None:
        ax.set_xticks(xticks)
        if xtick_fmt:
            ax.set_xticklabels([xtick_fmt.format(x) for x in xticks])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_axisbelow(True)


def save(fig, base_path):
    """Save base_path.pdf / .svg (vector) + .png (300 dpi preview).
    Retries if the file is temporarily locked (Windows: preview / antivirus)."""
    out = []
    for ext, kw in (("pdf", {}), ("svg", {}), ("png", {"dpi": 300})):
        path = f"{base_path}.{ext}"
        for attempt in range(10):
            try:
                fig.savefig(path, **kw)
                break
            except OSError:
                if attempt == 9:
                    raise
                time.sleep(0.5)
        out.append(path)
    plt.close(fig)
    return out
