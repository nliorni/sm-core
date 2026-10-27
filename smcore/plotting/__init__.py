"""sm-core figure style: one look for every plot the suite draws.

The style follows the visual language of Chen Liu's figures4papers
(https://github.com/ChenLiu-1996/figures4papers): sans-serif type, no top or
right spines, heavy axis lines, frameless legends, a small semantic palette
(deep blue, red-pink, green, neutral grey) and values written on the figure
rather than read off a grid. This module is an independent implementation of
that style for sm-core; no code is copied from that repository, which is
licensed CC BY-NC 4.0.

Pipeline scripts run inside their own conda environments, which only need
matplotlib. They import this module by putting the sm-core root on sys.path
(rules pass it as a param). Notebooks can do the same, or use the matplotlib
style file next to this module:

    import matplotlib.pyplot as plt
    plt.style.use("<sm-core>/smcore/plotting/smcore.mplstyle")

Colour roles are fixed across the suite, so a figure reads the same in every
pipeline:

    up / increased / positive    UP      deep red
    down / decreased / negative  DOWN    deep blue
    not significant / background NS      light grey
    categories (conditions...)   CATEGORICAL, in order
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")  # headless: pipelines run on cluster nodes

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, to_rgb  # noqa: E402

STYLE_FILE = Path(__file__).with_name("smcore.mplstyle")

PALETTE = {
    "blue": "#0F4D92",
    "blue_light": "#3775BA",
    "blue_pale": "#C9D9EE",
    "red": "#B64342",
    "red_light": "#E9A6A1",
    "red_pale": "#F6CFCB",
    "green": "#8BCF8B",
    "green_light": "#AADCA9",
    "green_pale": "#DDF3DE",
    "teal": "#42949E",
    "violet": "#9A4D8E",
    "gold": "#E0B000",
    "grey": "#CFCECE",
    "grey_dark": "#767676",
    "ink": "#272727",
}

UP = PALETTE["red"]
DOWN = PALETTE["blue"]
NS = PALETTE["grey"]

# Categories in a fixed order: the first two are the suite's anchor colours, so
# a two-group comparison is always blue vs red.
CATEGORICAL = [
    PALETTE["blue"],
    PALETTE["red"],
    PALETTE["green"],
    PALETTE["teal"],
    PALETTE["violet"],
    PALETTE["gold"],
    PALETTE["blue_light"],
    PALETTE["red_light"],
    PALETTE["grey_dark"],
]

DIVERGING = LinearSegmentedColormap.from_list(
    "smcore_diverging", [PALETTE["blue"], PALETTE["blue_light"], "#FFFFFF", PALETTE["red_light"], PALETTE["red"]]
)
SEQUENTIAL_RED = LinearSegmentedColormap.from_list("smcore_red", ["#FFFFFF", PALETTE["red_pale"], PALETTE["red"]])
SEQUENTIAL_BLUE = LinearSegmentedColormap.from_list("smcore_blue", ["#FFFFFF", PALETTE["blue_pale"], PALETTE["blue"]])


_FONT_DIRS = ["/usr/share/fonts", "/usr/local/share/fonts", "~/.fonts", "~/.local/share/fonts",
              "/Library/Fonts", "/System/Library/Fonts", "C:/Windows/Fonts"]
_FONT_FILES = ("Arial*.ttf", "arial*.ttf", "Helvetica*.tt[fc]", "LiberationSans-*.ttf", "Arimo*.ttf")


def _register_system_fonts() -> None:
    """Make Arial / Helvetica / Liberation Sans visible to matplotlib. Inside a
    conda environment matplotlib often only sees its bundled DejaVu, and the
    style asks for an Arial-like face (Liberation Sans and Arimo have Arial's
    metrics). Missing fonts are fine: the style falls back to DejaVu Sans."""
    from matplotlib import font_manager

    known = {Path(f.fname).name for f in font_manager.fontManager.ttflist}
    for directory in _FONT_DIRS:
        root = Path(directory).expanduser()
        if not root.is_dir():
            continue
        for pattern in _FONT_FILES:
            for path in root.rglob(pattern):
                if path.name not in known:
                    try:
                        font_manager.fontManager.addfont(str(path))
                        known.add(path.name)
                    except Exception:  # an unreadable font is not worth failing a plot
                        pass


def apply_style(font_size: float = 14) -> None:
    """Load the sm-core rcParams; `font_size` scales every text element."""
    _register_system_fonts()
    plt.style.use(STYLE_FILE)
    plt.rcParams.update(
        {
            "font.size": font_size,
            "axes.labelsize": font_size,
            "axes.titlesize": font_size + 1,
            "xtick.labelsize": font_size - 2,
            "ytick.labelsize": font_size - 2,
            "legend.fontsize": font_size - 2,
        }
    )


def colors_for(levels: Iterable[str]) -> dict:
    """A stable colour per category level, in the order given."""
    return {level: CATEGORICAL[i % len(CATEGORICAL)] for i, level in enumerate(levels)}


def lighten(color: str, amount: float) -> tuple:
    """Blend `color` towards white; amount 0 = unchanged, 1 = white."""
    r, g, b = to_rgb(color)
    return (r + (1 - r) * amount, g + (1 - g) * amount, b + (1 - b) * amount)


def save(fig, path_without_suffix: str | Path, formats: Sequence[str] = ("png", "pdf")) -> list[Path]:
    """Write the figure once per format (PNG at 300 dpi for viewing, PDF with
    editable text for figures), then close it."""
    base = Path(path_without_suffix)
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(pad=1.5)
    written = []
    for fmt in formats:
        out = base.with_name(f"{base.name}.{fmt}")
        fig.savefig(out, dpi=300, facecolor="white")
        written.append(out)
    plt.close(fig)
    return written


def label_bars(ax, bars, labels: Sequence[str], inside: bool = False, fontsize: float | None = None) -> None:
    """Write a value next to each horizontal bar: the style reads numbers off
    the figure, not off a grid."""
    fontsize = fontsize or plt.rcParams["font.size"] - 3
    for bar, text in zip(bars, labels):
        width = bar.get_width()
        y = bar.get_y() + bar.get_height() / 2
        ha = "left" if width >= 0 else "right"
        pad = 0.01 * (ax.get_xlim()[1] - ax.get_xlim()[0])
        x = width + (pad if width >= 0 else -pad)
        if inside:
            ha = "right" if width >= 0 else "left"
            x = width - (pad if width >= 0 else -pad)
        ax.text(x, y, text, va="center", ha=ha, fontsize=fontsize, color=PALETTE["ink"])
