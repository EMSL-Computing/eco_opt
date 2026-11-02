"""Shared matplotlib style for the manuscript figures.

Not part of the eco_opt package -- a manuscript build-time helper.
Call ``apply_style()`` once at the top of each figure script.

Figures are drawn at their printed size (the elsarticle ``3p`` text
width, 468 pt = 6.5 in) and included with ``width=\\linewidth``, so font
sizes here are the sizes on the page.

Colours carry one meaning across all figures:

- ``STATE_COLOR``   soil pools (states)
- ``CONTROL_COLOR`` the investment (control)
- ``PRICE_COLOR``   shadow prices (costates)
- ``PAYOFF_COLOR``  payoff (objective, growth)
- ``REF_COLOR``     analytical or reference solution
- simple rules compared with the optimum are drawn in greys

Layout rules:

- One shared legend above the panels (``shared_legend``) instead of
  per-panel legends; line styles mean the same in every panel.
- No legend or annotation may cover a plotted line. ``save`` runs
  ``check_overlaps`` and prints every text box that a line crosses;
  fix the placement until it reports none.
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl

# Okabe-Ito, colour-blind safe.
STATE_COLOR = "#0072B2"      # blue
CONTROL_COLOR = "#D55E00"    # vermillion
PRICE_COLOR = "#CC79A7"      # reddish purple
PAYOFF_COLOR = "#009E73"     # bluish green
REF_COLOR = "#222222"
SECONDARY_COLOR = PAYOFF_COLOR
TERTIARY_COLOR = PRICE_COLOR
GRID_COLOR = "#E6E6E6"

PALETTE = [STATE_COLOR, CONTROL_COLOR, PAYOFF_COLOR, PRICE_COLOR,
           "#E69F00", "#56B4E9"]

# Line conventions: reference solutions solid and thin, eco_opt dotted
# and thick on top, simple rules dashed.
REF_LW = 1.2
OPT_LW = 2.4
OPT_LS = (0, (1, 1.2))
RULE_LS = (0, (4, 2))

# Widths in inches.
TEXT_W = 468 / 72.27  # elsarticle final,3p \textwidth (468 pt)
HALF_W = TEXT_W / 2
# Legacy names used by older scripts.
ONE_COL_IN = HALF_W
TWO_COL_IN = TEXT_W

FIG_DIR = Path(os.environ.get("ECO_OPT_FIGURE_DIR", Path(__file__).resolve().parent / "figures")).expanduser()

FONT = "DejaVu Sans"  # TrueType outlines embed consistently in PDF viewers.


def apply_style(*, large: bool = False) -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [FONT, "Helvetica", "Arial", "DejaVu Sans"],
        "mathtext.fontset": "custom",
        "mathtext.rm": FONT,
        "mathtext.it": f"{FONT}:italic",
        "mathtext.bf": f"{FONT}:bold",
        "font.size": 7.5,
        "axes.titlesize": 7.5,
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelsize": 7.5,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5,
        "legend.frameon": False,
        "legend.handlelength": 2.0,
        "axes.linewidth": 0.6,
        "axes.edgecolor": "#333333",
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.prop_cycle": mpl.cycler(color=PALETTE),
        "lines.linewidth": REF_LW,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    if large:
        mpl.rcParams.update({"font.size": 8.5, "axes.labelsize": 8.5,
                             "axes.titlesize": 8.5})


def panel_label(ax, letter, x=-0.02, y=1.02):
    """Bold lower-case panel letter at the top left, outside the axes."""
    ax.text(x, y, letter, transform=ax.transAxes, fontweight="bold",
            fontsize=8.5, ha="right", va="bottom")


def shared_legend(fig, entries, y=1.0):
    """One legend row above all panels.

    ``entries`` is a list of (label, line kwargs). Line styles mean the
    same in every panel, so per-panel legends are not needed.
    """
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], **kw) for _, kw in entries]
    fig.legend(handles, [lab for lab, _ in entries], loc="lower center",
               bbox_to_anchor=(0.5, y), ncol=len(entries), fontsize=6.5,
               handlelength=2.4, columnspacing=1.6, borderaxespad=0.0)


def _line_points(line, n=400):
    """Display coordinates of a line, densified along each segment."""
    import numpy as np
    xy = line.get_transform().transform(np.column_stack(line.get_data()))
    xy = xy[np.all(np.isfinite(xy), axis=1)]
    if len(xy) < 2 or line.get_linestyle() in ("None", "", " "):
        return xy
    seg = np.linalg.norm(np.diff(xy, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    t = np.linspace(0.0, s[-1], max(n, len(xy)))
    return np.column_stack([np.interp(t, s, xy[:, 0]), np.interp(t, s, xy[:, 1])])


def check_overlaps(fig, pad=1.0):
    """List legends and annotations that a plotted line passes through.

    Checks, in each axes, every legend and every text artist placed by
    the script (not tick labels, axis labels, or titles) against the
    lines drawn in that axes. Returns a list of messages.
    """
    import numpy as np
    from matplotlib.legend import Legend
    from matplotlib.text import Text
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    problems = []
    for ax in fig.axes:
        lines = [ln for ln in ax.get_lines() if ln.get_visible()]
        if not lines:
            continue
        boxes = [(f"legend", a.get_window_extent(renderer))
                 for a in ax.get_children() if isinstance(a, Legend)]
        # Text.get_window_extent: the text only, not an annotation's leader.
        boxes += [(repr(t.get_text()), Text.get_window_extent(t, renderer))
                  for t in ax.texts if t.get_visible() and t.get_text().strip()]
        clip = ax.get_window_extent(renderer)
        for label, bb in boxes:
            bb = bb.expanded(1.0, 1.0).padded(pad)
            for ln in lines:
                pts = _line_points(ln)
                inside = ((pts[:, 0] > bb.x0) & (pts[:, 0] < bb.x1)
                          & (pts[:, 1] > bb.y0) & (pts[:, 1] < bb.y1)
                          & (pts[:, 0] >= clip.x0) & (pts[:, 0] <= clip.x1)
                          & (pts[:, 1] >= clip.y0) & (pts[:, 1] <= clip.y1))
                if np.any(inside):
                    name = ln.get_label()
                    name = "" if name.startswith("_") else f" ({name})"
                    problems.append(f"{label} is crossed by a line{name}")
                    break
    return problems


def save(fig, name):
    """Write ``figures/<name>.pdf`` (for LaTeX) and a local PNG preview."""
    for msg in check_overlaps(fig):
        print("overlap:", msg)
    preview = FIG_DIR / "preview"
    preview.mkdir(parents=True, exist_ok=True)
    pdf = FIG_DIR / f"{name}.pdf"
    fig.savefig(pdf)
    fig.savefig(preview / f"{name}.png", dpi=200)
    print(f"wrote {pdf}")
    return pdf
