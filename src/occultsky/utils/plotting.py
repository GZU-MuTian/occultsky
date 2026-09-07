"""
Chart drawing helpers: the night timeline and Stellarium constellation charts.
"""

import logging
from datetime import timedelta
from typing import Any, Optional

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from ..main import Observer

logger = logging.getLogger(__name__)

# Body colors: Moon fixed silver, planets from matplotlib default palette
_MOON_COLOR = "#808080"
_PLANET_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2"]
_NIGHT_BAND_COLOR = "midnightblue"
_NIGHT_BAND_ALPHA = 0.12
# Axes-fraction margin left empty above and below the night band
_NIGHT_BAND_MARGIN = 0.35
_SOLAR_NOON_COLOR = "#CC6600"

# Vertical geometry (data units, axis line at y=0)
_MARKER_SIZE = 6
# Two staggered label levels per side to avoid overlap
_RISE_LEVELS = (0.75, 1.15)
_SET_LEVELS = (-0.75, -1.15)
_LABEL_GAP = 0.18


def draw_timeline(
    lat: float,
    lon: float,
    elevation: float,
    date_str: str,
    planets: list[str] = [],
    output: Optional[str] = None,
    **observer_kwargs: Any,
) -> None:
    """Draw a single-axis astronomical timeline chart.

    Night periods are drawn as semi-transparent bands around the time axis; rise and set events are drawn as milestone markers (colored dot + bold time and body label), with rise labels above the axis and set labels below.

    Args:
        lat: Latitude in degrees (north positive).
        lon: Longitude in degrees (east positive).
        elevation: Elevation above sea level in meters.
        date_str: Local date in ``YYYY-MM-DD`` format.
        planets: Optional list of planet names (e.g. ``['mars', 'jupiter']``).
        output: If provided, save to this file path instead of showing.
        **observer_kwargs: Extra keyword arguments forwarded to ``Observer`` (e.g. ``skyfield_data``, ``ephemeris``).
    """
    observer = Observer(
        lat=lat, lon=lon, elevation=elevation, **observer_kwargs
    )

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 3.5), layout="constrained")

    # Compute time window: solar noon ± 1 day
    solar_noon = observer.get_solarnoon(date_str)
    window_start = solar_noon - timedelta(days=1)
    window_end = solar_noon + timedelta(days=1)
    # Convert to matplotlib date numbers
    x_min = mdates.date2num(window_start)
    x_max = mdates.date2num(window_end)

    # Set axis limits and remove y-axis
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-1.55, 1.6)
    ax.set_yticks([])

    # Single-axis figure
    ax.axhline(0, color="black", linewidth=1, zorder=1)
    ax.spines[:].set_visible(False)

    # X-axis format: 24-hour clock in LMT
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=6, tz=observer.tz))
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=2, tz=observer.tz))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d %H:%M", tz=observer.tz))
    ax.tick_params(axis="x", which="major", labelsize=12)

    # Night periods: semi-transparent bands centered on the axis line.
    previous_dusk, dawn, dusk, next_dawn = observer.get_twilight(date_str)
    night_bands = [(previous_dusk, dawn), (dusk, next_dawn)]
    for start, end in night_bands:
        ax.axvspan(
            start, end,
            ymin=_NIGHT_BAND_MARGIN, ymax=1 - _NIGHT_BAND_MARGIN,
            color=_NIGHT_BAND_COLOR, alpha=_NIGHT_BAND_ALPHA, zorder=0,
        )

    # True solar noon marker: vertical line + label on the axis
    x_noon = mdates.date2num(solar_noon)
    ax.vlines(x_noon, -1.4, 1.4, color=_SOLAR_NOON_COLOR, linewidth=1, linestyle="--", zorder=2)
    ax.text(x_noon, 1.48, f"Solar noon {solar_noon.strftime('%H:%M')}", color=_SOLAR_NOON_COLOR, fontsize=12, ha="center", va="bottom", zorder=4)

    # bodies: list[tuple[str, str, list[datetime], list[datetime]]] = [
    #     ("Moon", _MOON_COLOR, observer.get_moonrise(date_str), observer.get_moonset(date_str)),
    # ]
    # for i, planet in enumerate(planets):
    #     color = _PLANET_COLORS[i % len(_PLANET_COLORS)]
    #     bodies.append((
    #         planet.capitalize(),
    #         color,
    #         observer.get_planet_rise(planet, date_str),
    #         observer.get_planet_set(planet, date_str),
    #     ))




    # Collect in-window events, sorted chronologically per side
    # rise_events: list[tuple[datetime, str, str]] = []
    # set_events: list[tuple[datetime, str, str]] = []
    # for name, color, rises, sets in bodies:
    #     for t in rises:
    #         if window_start <= t <= window_end:
    #             rise_events.append((t, name, color))
    #     for t in sets:
    #         if window_start <= t <= window_end:
    #             set_events.append((t, name, color))
    # rise_events.sort(key=lambda e: e[0])
    # set_events.sort(key=lambda e: e[0])

    # Milestone markers: colored dot + bold time + body label
    # for i, (t, name, color) in enumerate(rise_events):
    #     x = mdates.date2num(t)
    #     level = _RISE_LEVELS[i % len(_RISE_LEVELS)]
    #     ax.plot(x, 0, "o", color=color, markersize=_MARKER_SIZE, zorder=3)
    #     ax.text(
    #         x, level, t.strftime("%H:%M"),
    #         color=color, fontweight="bold", fontsize=8,
    #         ha="center", va="bottom", zorder=4,
    #     )
    #     ax.text(
    #         x, level - _LABEL_GAP, f"{name} rise",
    #         color=color, fontsize=7,
    #         ha="center", va="bottom", zorder=4,
    #     )
    # for i, (t, name, color) in enumerate(set_events):
    #     x = mdates.date2num(t)
    #     level = _SET_LEVELS[i % len(_SET_LEVELS)]
    #     ax.plot(x, 0, "o", color=color, markersize=_MARKER_SIZE, zorder=3)
    #     ax.text(
    #         x, level, t.strftime("%H:%M"),
    #         color=color, fontweight="bold", fontsize=8,
    #         ha="center", va="top", zorder=4,
    #     )
    #     ax.text(
    #         x, level + _LABEL_GAP, f"{name} set",
    #         color=color, fontsize=7,
    #         ha="center", va="top", zorder=4,
    #     )



    # Timezone (LMT) annotation in the upper-right corner
    # offset = observer.tz.utcoffset(None)
    # hours, remainder = divmod(abs(int(offset.total_seconds())), 3600)
    # tz_str = f"UTC{'+' if offset >= timedelta(0) else '-'}{hours:02d}:{remainder // 60:02d} (LMT)"
    # ax.annotate(tz_str, xy=(1.0, 1.0), xycoords="axes fraction",
    #             xytext=(0, -4), textcoords="offset points",
    #             ha="right", va="top", fontsize=8, color="gray")

    # Color legend for bodies, placed below the timezone annotation
    # handles = [Line2D([0], [0], color=color, lw=2, label=name) for name, color, _, _ in bodies]
    # ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.0, 0.86),
    #           fontsize=8, framealpha=0.9)



    # Title
    # lat = observer.topos.latitude.degrees
    # lon = observer.topos.longitude.degrees
    # lat_str = f"{abs(lat):.1f}°{'N' if lat >= 0 else 'S'}"
    # lon_str = f"{abs(lon):.1f}°{'E' if lon >= 0 else 'W'}"
    # ax.set_title(f"{date_str} ({lat_str}, {lon_str})", fontsize=12)

    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
        logger.info("Timeline saved to %s", output)
    else:
        plt.show()



