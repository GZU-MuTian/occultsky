"""
Single-axis timeline visualization for night periods and rise/set events.
"""

import logging
from datetime import datetime, timedelta
from typing import Optional

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from ..main import Observer

logger = logging.getLogger(__name__)

# Body colors: Moon fixed silver, planets from matplotlib default palette
_MOON_COLOR = "#808080"
_PLANET_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2"]
_NIGHT_BOX_COLOR = "#333333"
_SOLAR_NOON_COLOR = "#CC6600"

# Vertical geometry (data units, axis line at y=0)
_NIGHT_BOX_HALF = 0.4
_MARKER_SIZE = 6
# Two staggered label levels per side to avoid overlap
_RISE_LEVELS = (0.75, 1.15)
_SET_LEVELS = (-0.75, -1.15)
_LABEL_GAP = 0.18


def draw_timeline(
    observer: Observer,
    date_str: str,
    planets: Optional[list[str]] = None,
    output: Optional[str] = None,
) -> None:
    """Draw a single-axis astronomical timeline chart.

    Night periods are drawn as hollow boxes around the time axis; rise and
    set events are drawn as milestone markers (colored dot + bold time and
    body label), with rise labels above the axis and set labels below.

    Args:
        observer: An Observer instance.
        date_str: Local date in ``YYYY-MM-DD`` format.
        planets: Optional list of planet names (e.g. ``['mars', 'jupiter']``).
        output: If provided, save to this file path instead of showing.
    """
    if planets is None:
        planets = []

    # Compute time window: solar noon ± 1 day
    solar_noon = observer.get_solarnoon(date_str)
    window_start = solar_noon - timedelta(days=1)
    window_end = solar_noon + timedelta(days=1)
    x_min = mdates.date2num(window_start)
    x_max = mdates.date2num(window_end)

    # Gather data: night periods + rise/set events per body
    previous_dusk, dawn, dusk, next_dawn = observer.get_twilight(date_str)
    night_bands = [(previous_dusk, dawn), (dusk, next_dawn)]

    bodies: list[tuple[str, str, list[datetime], list[datetime]]] = [
        ("Moon", _MOON_COLOR, observer.get_moonrise(date_str), observer.get_moonset(date_str)),
    ]
    for i, planet in enumerate(planets):
        color = _PLANET_COLORS[i % len(_PLANET_COLORS)]
        bodies.append((
            planet.capitalize(),
            color,
            observer.get_planet_rise(planet, date_str),
            observer.get_planet_set(planet, date_str),
        ))

    # Single-axis figure. Datetimes are LMT-aware; date2num converts them
    # to absolute (UTC) numbers, so ticks are aligned to LMT multiples below.
    fig, ax = plt.subplots(figsize=(14, 3.5), layout="constrained")
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-1.55, 1.6)
    ax.set_yticks([])
    ax.axhline(0, color="black", linewidth=1, zorder=1)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Night periods: hollow boxes around the axis
    for start, end in night_bands:
        x0 = mdates.date2num(start)
        width = mdates.date2num(end) - x0
        ax.add_patch(Rectangle(
            (x0, -_NIGHT_BOX_HALF), width, 2 * _NIGHT_BOX_HALF,
            fill=False, edgecolor=_NIGHT_BOX_COLOR, linewidth=1.5, zorder=2,
        ))

    # Collect in-window events, sorted chronologically per side
    rise_events: list[tuple[datetime, str, str]] = []
    set_events: list[tuple[datetime, str, str]] = []
    for name, color, rises, sets in bodies:
        for t in rises:
            if window_start <= t <= window_end:
                rise_events.append((t, name, color))
        for t in sets:
            if window_start <= t <= window_end:
                set_events.append((t, name, color))
    rise_events.sort(key=lambda e: e[0])
    set_events.sort(key=lambda e: e[0])

    # Milestone markers: colored dot + bold time + body label
    for i, (t, name, color) in enumerate(rise_events):
        x = mdates.date2num(t)
        level = _RISE_LEVELS[i % len(_RISE_LEVELS)]
        ax.plot(x, 0, "o", color=color, markersize=_MARKER_SIZE, zorder=3)
        ax.text(
            x, level, t.strftime("%H:%M"),
            color=color, fontweight="bold", fontsize=8,
            ha="center", va="bottom", zorder=4,
        )
        ax.text(
            x, level - _LABEL_GAP, f"{name} rise",
            color=color, fontsize=7,
            ha="center", va="bottom", zorder=4,
        )
    for i, (t, name, color) in enumerate(set_events):
        x = mdates.date2num(t)
        level = _SET_LEVELS[i % len(_SET_LEVELS)]
        ax.plot(x, 0, "o", color=color, markersize=_MARKER_SIZE, zorder=3)
        ax.text(
            x, level, t.strftime("%H:%M"),
            color=color, fontweight="bold", fontsize=8,
            ha="center", va="top", zorder=4,
        )
        ax.text(
            x, level + _LABEL_GAP, f"{name} set",
            color=color, fontsize=7,
            ha="center", va="top", zorder=4,
        )

    # True solar noon marker: vertical line + label on the axis
    x_noon = mdates.date2num(solar_noon)
    ax.vlines(x_noon, -1.4, 1.4, color=_SOLAR_NOON_COLOR, linewidth=1,
              linestyle="--", zorder=2)
    ax.text(x_noon, 1.48, f"Solar noon {solar_noon.strftime('%H:%M')}",
            color=_SOLAR_NOON_COLOR, fontsize=8,
            ha="center", va="bottom", zorder=4)

    # Timezone (LMT) annotation in the upper-right corner
    offset = observer.tz.utcoffset(None)
    hours, remainder = divmod(abs(int(offset.total_seconds())), 3600)
    tz_str = f"UTC{'+' if offset >= timedelta(0) else '-'}{hours:02d}:{remainder // 60:02d} (LMT)"
    ax.annotate(tz_str, xy=(1.0, 1.0), xycoords="axes fraction",
                xytext=(0, -4), textcoords="offset points",
                ha="right", va="top", fontsize=8, color="gray")

    # Color legend for bodies, placed below the timezone annotation
    handles = [Line2D([0], [0], color=color, lw=2, label=name) for name, color, _, _ in bodies]
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.0, 0.86),
              fontsize=8, framealpha=0.9)

    # X-axis format: 24-hour clock in LMT; ticks aligned to LMT multiples
    base = datetime.strptime(date_str, "%Y-%m-%d").replace(tzinfo=observer.tz) - timedelta(days=1)
    major_ticks, minor_ticks = [], []
    tick = base
    while tick <= window_end:
        hours = int((tick - base) / timedelta(hours=1))
        if hours % 6 == 0:
            major_ticks.append(mdates.date2num(tick))
        else:
            minor_ticks.append(mdates.date2num(tick))
        tick += timedelta(hours=2)
    ax.set_xticks(major_ticks)
    ax.set_xticks(minor_ticks, minor=True)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d %H:%M", tz=observer.tz))
    for label in ax.get_xticklabels():
        label.set_fontsize(8)

    # Title
    lat = observer.topos.latitude.degrees
    lon = observer.topos.longitude.degrees
    lat_str = f"{abs(lat):.1f}°{'N' if lat >= 0 else 'S'}"
    lon_str = f"{abs(lon):.1f}°{'E' if lon >= 0 else 'W'}"
    ax.set_title(f"{date_str} ({lat_str}, {lon_str})", fontsize=12)

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
        logger.info("Timeline saved to %s", output)
    else:
        plt.show()
