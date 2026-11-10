"""Shared colour palette for all CREST-GV paper figures.

Values are the validated default palette from the dataviz skill
(references/palette.md), used as-is: fixed categorical hue order (never
cycled), one hue for sequential/magnitude encodings, a blue<->red diverging
pair for signed enrichment scores, and a reserved status palette for the
discrete QC classes (Good/Decent/Bad). These are static/print (matplotlib)
figures for a PDF manuscript, so only the colour-formula and mark-spec parts
of the skill apply here -- no hover/dark-mode/interactive layer.
"""

import matplotlib as mpl

# Chart chrome & ink (light mode only -- print figures)
SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"

# Categorical palette, fixed order (slots 1-8)
CATEGORICAL = {
    "blue": "#2a78d6",
    "orange": "#eb6834",
    "aqua": "#1baf7a",
    "yellow": "#eda100",
    "magenta": "#e87ba4",
    "green": "#008300",
    "violet": "#4a3aa7",
    "red": "#e34948",
}

# Two-series slots used throughout (e.g. raw vs standardised, atlas vs upload)
SERIES_1 = CATEGORICAL["blue"]
SERIES_2 = CATEGORICAL["orange"]

# Diverging pair for signed enrichment scores (Z / Z-tilde)
DIVERGING_POS = CATEGORICAL["blue"]
DIVERGING_NEG = CATEGORICAL["red"]
DIVERGING_MID = "#f0efec"

# Sequential single-hue ramp (blue), light -> dark, for magnitude-only encodings
SEQUENTIAL_BLUE = [
    "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec",
    "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab",
    "#184f95", "#104281", "#0d366b",
]

# Status palette (fixed, never themed) -- used only for the QC classes
STATUS_GOOD = "#0ca30c"
STATUS_WARNING = "#fab219"
STATUS_CRITICAL = "#d03b3b"

QC_LABEL_COLOR = {"Good": STATUS_GOOD, "Decent": STATUS_WARNING, "Bad": STATUS_CRITICAL}


def apply_base_style():
    """Thin marks, white background, no gridlines, system sans -- applied once per script."""
    mpl.rcParams.update({
        "figure.facecolor": "#ffffff",
        "axes.facecolor": "#ffffff",
        "savefig.facecolor": "#ffffff",
        "axes.edgecolor": BASELINE,
        "axes.labelcolor": INK_PRIMARY,
        "xtick.color": INK_MUTED,
        "ytick.color": INK_MUTED,
        "text.color": INK_PRIMARY,
        "axes.grid": False,
        "axes.axisbelow": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.family": "sans-serif",
        "font.size": 10,
        "figure.dpi": 150,
        "savefig.dpi": 300,
    })
