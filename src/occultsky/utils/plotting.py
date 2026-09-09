"""
Chart drawing helpers: the night timeline and Stellarium constellation charts.
"""

import logging
from datetime import datetime, timedelta
from typing import Any, Optional

from astropy.visualization.wcsaxes import WCSAxes
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from shapely import affinity
from shapely.geometry import LineString
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

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
# Skyview chart geometry: square canvas field of view (deg) and pixel size
_SKYVIEW_FOV = 6
_SKYVIEW_SIZE = 800


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

    # Color legend: night band plus bodies
    handles = [
        Patch(
            facecolor=_NIGHT_BAND_COLOR, alpha=_NIGHT_BAND_ALPHA,
            label="Night",
        ),
    ]
    handles += [Line2D([0], [0], color=_BODY_COLORS[i], lw=2, label=name) for i, (name, _, _) in enumerate(bodies)]
    ax.legend(handles=handles, loc="best", bbox_to_anchor=(1.0, 0.98), fontsize=10, framealpha=0.9)

    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
        logger.info("Timeline saved to %s", output)
    else:
        plt.show()


def add_moon_track(
    ax: WCSAxes,
    coordinates: list[tuple[float, float]],
    angular_radius: float,
) -> None:
    """Draw the sky region the Moon's disk sweeps along its nightly track.

    The region is a filled band of width ``2 * angular_radius`` around the sampled ICRS path.

    Args:
        ax: WCS axes to draw on.
        coordinates: ``(ra, dec)`` tuples in degrees, as returned by ``Observer.get_radec``.
        angular_radius: Apparent angular radius of the Moon in degrees, as returned by ``Observer.get_moon_angular_radius``.
    """
    track_ra = np.degrees(
        np.unwrap(np.radians([item[0] for item in coordinates]))
    )
    track_dec = np.array([item[1] for item in coordinates], dtype=float)

    # A degree of right ascension spans only cos(dec) degrees of arc on the sky, so buffering raw (ra, dec) coordinates would squash the disk horizontally by cos(dec). Compress the path in RA first (median dec gives a local plane approximation, good to <1% over one night), buffer, then stretch back: the band keeps the disk's true angular size on the celestial sphere. The stretch must be anchored at the origin, or the path itself would shift.
    cos_dec = np.cos(np.radians(float(np.median(track_dec))))
    envelope = LineString(zip(track_ra * cos_dec, track_dec)).buffer(
        angular_radius, quad_segs=32
    )
    envelope = affinity.scale(
        envelope, xfact=1.0 / cos_dec, yfact=1.0, origin=(0.0, 0.0)
    )
    xs, ys = envelope.exterior.xy
    ax.fill(
        xs, ys,
        color=_BODY_COLORS[0], alpha=0.3,
        transform=ax.get_transform("world"),
    )


def add_planet_track(
    ax: WCSAxes,
    coordinates: list[tuple[float, float]],
    color: str,
) -> None:
    """Draw one planet's nightly positions as overlapping dots.

    Args:
        ax: WCS axes to draw on.
        coordinates: ``(ra, dec)`` tuples in degrees, as returned by ``Observer.get_radec``.
        color: Color of the dots.
    """
    track_ra = np.array([item[0] for item in coordinates], dtype=float)
    track_dec = np.array([item[1] for item in coordinates], dtype=float)
    ax.plot(
        track_ra, track_dec, "o",
        color=color, markersize=4,
        transform=ax.get_transform("world"),
    )


def draw_skyview(
    lat: float,
    lon: float,
    elevation: float,
    date_str: str,
    planets: list[str] = [],
    interval: int = 15,
    projection: str = "STG",
    output: Optional[str] = None,
    stars: bool = False,
    star_names: bool = False,
    constellation_lines: bool = False,
    **observer_kwargs: Any,
) -> None:
    """Draw a two-subplot sky view chart of the Moon's nightly tracks.

    One square WCS-projected subplot per night window of ``date_str``: the
    previous night (``previous_dusk`` to ``dawn``) on the left, the night of
    ``date_str`` (``dusk`` to ``next_dawn``) on the right. Each subplot centers
    on the midpoint of its own Moon track, so the whole track fits in the
    6° field of view. A night with no visible Moon is left blank.

    Args:
        lat: Latitude in degrees (north positive).
        lon: Longitude in degrees (east positive).
        elevation: Elevation above sea level in meters.
        date_str: Local date in ``YYYY-MM-DD`` format.
        planets: Optional planet names drawn alongside the Moon (e.g. ``['mars', 'jupiter']``).
        interval: Sampling step in minutes (default: 15).
        projection: WCS projection code (default: ``STG``).
        output: If provided, save to this file path instead of showing.
        stars: Not implemented here; passing ``True`` raises ``ValueError``.
        star_names: Not implemented here; passing ``True`` raises ``ValueError``.
        constellation_lines: Not implemented here; passing ``True`` raises
            ``ValueError``.
        **observer_kwargs: Extra keyword arguments forwarded to ``Observer``
            (e.g. ``skyfield_data``, ``ephemeris``).

    Raises:
        ValueError: If ``stars``, ``star_names`` or ``constellation_lines`` is
            enabled.
    """
    # Star and constellation layers are drawn by utils.skyculture; they are
    # not part of this chart yet, so fail loudly instead of dropping them.
    if stars or star_names or constellation_lines:
        raise ValueError(
            "stars/star_names/constellation_lines are not implemented in the "
            "skyview chart yet; see utils.skyculture for constellation charts"
        )

    from .wcs import build_wcs

    observer = Observer(
        lat=lat, lon=lon, elevation=elevation, **observer_kwargs
    )

    # One subplot per night window. Each subplot centers on the midpoint of its own Moon track.
    nights = observer.get_visible_times(
        date_str, target="moon", interval=interval
    )

    fig = plt.figure(figsize=(18, 9))

    # Figure title
    lat_deg = observer.topos.latitude.degrees
    lon_deg = observer.topos.longitude.degrees
    lat_str = f"{abs(lat_deg):.1f}°{'N' if lat_deg >= 0 else 'S'}"
    lon_str = f"{abs(lon_deg):.1f}°{'E' if lon_deg >= 0 else 'W'}"
    fig.suptitle(f"{date_str} ({lat_str}, {lon_str})", fontsize=14)
    # The Moon's disk radius, drawn to scale along the track
    moon_radius = observer.get_moon_angular_radius(date_str)

    # Sample each planet once per night window, bundled with its color
    planet_tracks = [
        (
            planet,
            observer.get_visible_times(
                date_str, target=planet, interval=interval
            ),
            color,
        )
        for planet, color in zip(planets, _BODY_COLORS[1:])
    ]

    for index, times in enumerate(nights):
        # No visible Moon that night: leave the subplot blank.
        if not times:
            ax = fig.add_subplot(1, 2, index + 1)
            ax.set_axis_off()
            continue

        coordinates = observer.get_radec(times, target="moon")

        # Find the approximate center of the Moon track
        center = coordinates[len(times) // 2]
        wcs = build_wcs(projection, center, _SKYVIEW_FOV, _SKYVIEW_SIZE)

        # Add the WCS axes to the figure
        ax = fig.add_subplot(1, 2, index + 1, projection=wcs)
        ax.set_xlim(-0.5, _SKYVIEW_SIZE - 0.5)
        ax.set_ylim(-0.5, _SKYVIEW_SIZE - 0.5)

        # Grid and axis labels
        for coord in ax.coords:
            coord.grid(True, color="#9a9a9a", alpha=0.35, linestyle=":", linewidth=0.6)
        ax.coords[0].set_axislabel("Right ascension (ICRS)")
        ax.coords[1].set_axislabel("Declination (ICRS)")

        # Add the Moon's track as a filled band
        add_moon_track(ax, coordinates, angular_radius=moon_radius)

        # Add the planet tracks as overlapping dots
        for planet, planet_nights, color in planet_tracks:
            planet_times = planet_nights[index]
            if not planet_times:
                continue
            add_planet_track(
                ax,
                observer.get_radec(planet_times, target=planet),
                color=color,
            )

    # One shared legend for both subplots
    handles = [
        Patch(color=_BODY_COLORS[0], alpha=0.3, label="Moon"),
    ]
    for planet, _, color in planet_tracks:
        handles.append(
            Line2D(
                [0], [0],
                color=color,
                marker="o", markersize=4,
                label=planet.capitalize(),
            )
        )
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.99, 0.97))

    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
        logger.info("Skyview saved to %s", output)
    else:
        plt.show()
