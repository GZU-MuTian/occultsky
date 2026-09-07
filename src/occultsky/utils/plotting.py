"""
Chart drawing helpers: the night timeline and Stellarium constellation charts.
"""

import logging
from datetime import datetime, timedelta
from typing import Any, Optional

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from ..main import Observer

logger = logging.getLogger(__name__)

# Body colors: Moon fixed silver, planets from matplotlib default palette
_BODY_COLORS = [
    "#808080",
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2", "#17becf",
]
_NIGHT_BAND_COLOR = "midnightblue"
_NIGHT_BAND_ALPHA = 0.12
# Axes-fraction margin left empty above and below the night band
_NIGHT_BAND_MARGIN = 0.37
_SOLAR_NOON_COLOR = "#CC6600"
# Axes-fraction blank left above and below the solar-noon line
_SOLAR_NOON_MARGIN = 0.1

# Vertical geometry (axis line at y=0)
_MARKER_SIZE = 8
# Four staggered label levels per side to avoid overlap, as axes
# fractions of the axes height (0.5 is the axis line)
_RISE_LEVELS = (0.58, 0.68, 0.78, 0.88)
_SET_LEVELS = (0.42, 0.32, 0.22, 0.12)
# Axes-fraction gap between the time line and the body line of a label
_LABEL_GAP = 0.04
# Label typography (pt)
_TIME_FONTSIZE = 14
_BODY_FONTSIZE = 14


def add_timestamp(
    ax: plt.Axes,
    time: datetime,
    level: float,
    time_str: str,
    body_str: str,
    color: str,
    va: str,
) -> None:
    """Add a milestone timestamp: a colored dot on the axis plus a two-line label.

    ``level`` is an axes fraction of the axes height (0.5 = the axis line); the body label sits one ``_LABEL_GAP`` away, above for ``va="bottom"`` and below for ``va="top"``.
    """
    x = mdates.date2num(time)
    trans = ax.get_xaxis_transform()
    ax.plot(x, 0, "o", color=color, markersize=_MARKER_SIZE, zorder=3)
    ax.text(
        x, level, time_str,
        color=color, fontweight="bold", fontsize=_TIME_FONTSIZE,
        ha="center", va=va, transform=trans, zorder=4,
    )
    offset = -_LABEL_GAP if va == "bottom" else _LABEL_GAP
    ax.text(
        x, level + offset, body_str,
        color=color, fontsize=_BODY_FONTSIZE,
        ha="center", va=va, transform=trans, zorder=4,
    )


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
    fig, ax = plt.subplots(figsize=(14, 5), layout="constrained")

    # Compute time window: solar noon ± 1 day
    solar_noon = observer.get_solarnoon(date_str)
    window_start = solar_noon - timedelta(days=1)
    window_end = solar_noon + timedelta(days=1)
    # Convert to matplotlib date numbers
    x_min = mdates.date2num(window_start)
    x_max = mdates.date2num(window_end)

    # Set axis limits and remove y-axis
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-2.8, 2.8)
    ax.set_yticks([])

    # Single-axis figure
    ax.axhline(0, color="black", linewidth=1, zorder=1)
    ax.spines[:].set_visible(False)

    # X-axis format: 24-hour clock in LMT
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=6, tz=observer.tz))
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=2, tz=observer.tz))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d %H:%M", tz=observer.tz))
    ax.tick_params(axis="x", which="major", labelsize=12)

    # Timezone note as the time-axis label.
    z = solar_noon.strftime('%z')
    tz_str = f"UTC{z[:3]}:{z[3:5]}"
    ax.set_xlabel(f"Local Mean Time ({tz_str})", loc="right", fontsize=10, color="gray")

    # Title
    lat = observer.topos.latitude.degrees
    lon = observer.topos.longitude.degrees
    lat_str = f"{abs(lat):.2f}°{'N' if lat >= 0 else 'S'}"
    lon_str = f"{abs(lon):.2f}°{'E' if lon >= 0 else 'W'}"
    ax.set_title(f"{date_str} ({lat_str}, {lon_str})", fontsize=16)

    # Night periods: semi-transparent bands centered on the axis line.
    previous_dusk, dawn, dusk, next_dawn = observer.get_twilight(date_str)
    night_bands = [(previous_dusk, dawn), (dusk, next_dawn)]
    for start, end in night_bands:
        ax.axvspan(
            start, end,
            ymin=_NIGHT_BAND_MARGIN, ymax=1 - _NIGHT_BAND_MARGIN,
            color=_NIGHT_BAND_COLOR, alpha=_NIGHT_BAND_ALPHA, zorder=0,
        )

    # True solar noon marker: dashed line with blank margins + label above it
    x_noon = mdates.date2num(solar_noon)
    ax.axvline(x_noon, ymin=_SOLAR_NOON_MARGIN, ymax=1 - _SOLAR_NOON_MARGIN, color=_SOLAR_NOON_COLOR, linewidth=1, linestyle="--", zorder=2)
    # Label lives in the blank band above the line (x in data, y in axes fraction)
    ax.text(x_noon, 1 - _SOLAR_NOON_MARGIN, f"Solar noon {solar_noon.strftime('%H:%M')}", color=_SOLAR_NOON_COLOR, fontsize=14, ha="center", va="bottom", transform=ax.get_xaxis_transform(), zorder=4)

    # Bodies: Moon and optionally planets; colors derive from list order
    bodies: list[tuple[str, list[datetime], list[datetime]]] = [
        ("Moon", observer.get_moonrise(date_str), observer.get_moonset(date_str)),
    ]
    for planet in planets:
        bodies.append((
            planet.capitalize(),
            observer.get_planet_rise(planet, date_str),
            observer.get_planet_set(planet, date_str),
        ))

    # Milestone markers per body: colored dot + bold time + body label.
    for i, (name, rises, sets) in enumerate(bodies):
        color = _BODY_COLORS[i]
        rise_level = _RISE_LEVELS[i % len(_RISE_LEVELS)]
        set_level = _SET_LEVELS[i % len(_SET_LEVELS)]
        for t in sorted(rises):
            add_timestamp(
                ax, t, rise_level,
                t.strftime("%H:%M"), f"{name} rise", color, "bottom",
            )
        for t in sorted(sets):
            add_timestamp(
                ax, t, set_level,
                t.strftime("%H:%M"), f"{name} set", color, "top",
            )

    # Color legend for bodies
    handles = [Line2D([0], [0], color=_BODY_COLORS[i], lw=2, label=name) for i, (name, _, _) in enumerate(bodies)]
    ax.legend(handles=handles, loc="best", bbox_to_anchor=(1.0, 0.98), fontsize=10, framealpha=0.9)

    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
        logger.info("Timeline saved to %s", output)
    else:
        plt.show()



