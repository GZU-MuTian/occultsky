"""
Night-sky observation calculator.

Reconstruct the sky as ancient observers would have seen it, using
Skyfield high-precision ephemerides (DE441).

Supported range: 13 200 BCE – 17 191 CE.
Julian calendar is used automatically before 1582-10-15.
"""

from __future__ import annotations

import logging
import math
import os
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

from dotenv import load_dotenv

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

from skyfield import almanac
from skyfield.api import Loader, wgs84
from skyfield.constellationlib import load_constellation_map

load_dotenv()

logger = logging.getLogger(__name__)

# ── Constants ────────────────────────────────────────────────────────

_JULIAN_CUTOFF_JD: int = 2299160.5  # 1582-10-15

_PLANET_ALIASES: dict[str, list[str]] = {
    "mercury": ["mercury", "mercury barycenter"],
    "venus": ["venus", "venus barycenter"],
    "mars": ["mars barycenter"],
    "jupiter": ["jupiter barycenter"],
    "saturn": ["saturn barycenter"],
    "uranus": ["uranus barycenter"],
    "neptune": ["neptune barycenter"],
}

DEFAULT_PLANETS: list[str] = [
    "mercury",
    "venus",
    "mars",
    "jupiter",
    "saturn",
]

_MOON_PHASE_NAMES: list[str] = [
    "New Moon",
    "Waxing Crescent",
    "First Quarter",
    "Waxing Gibbous",
    "Full Moon",
    "Waning Gibbous",
    "Last Quarter",
    "Waning Crescent",
]

# HIP number → common name for naked-eye bright stars (mag ≲ 2.0).
_STAR_NAMES: dict[int, str] = {
    # Canis Major
    32349: "Sirius",
    30438: "Mirzam",
    34444: "Wezen",
    35037: "Adhara",
    # Carina
    45238: "Miaplacidus",
    # Centaurus
    68702: "Hadar",
    60718: "Rigil Kent",
    # Gemini
    37826: "Pollux",
    36850: "Castor",
    30324: "Alhena",
    # Orion
    27989: "Betelgeuse",
    24436: "Rigel",
    27366: "Bellatrix",
    26311: "Mintaka",
    26207: "Alnilam",
    25930: "Alnitak",
    28614: "Saiph",
    # Taurus
    21421: "Aldebaran",
    17702: "Elnath",
    20889: "Alcyone",
    # Scorpius
    80763: "Antares",
    80331: "Shaula",
    78401: "Dschubba",
    78265: "Acrab",
    # Leo
    49669: "Regulus",
    57632: "Denebola",
    50583: "Algieba",
    # Ursa Major (Big Dipper)
    41704: "Merak",
    42911: "Dubhe",
    46853: "Mizar",
    45571: "Alioth",
    45394: "Megrez",
    43587: "Phecda",
    37764: "Alkaid",
    # Ursa Minor
    11767: "Polaris",
    # Lyra
    91262: "Vega",
    # Cygnus
    102098: "Deneb",
    # Aquila
    97649: "Altair",
    # Auriga
    24608: "Capella",
    # Bootes
    69673: "Arcturus",
    # Virgo
    65474: "Spica",
    # Crux
    62434: "Gacrux",
    # Corona Borealis
    76267: "Alphecca",
    # Triangulum Australe
    77952: "Atria",
    # Ophiuchus
    86032: "Rasalhague",
    # Pegasus
    113963: "Markab",
    113881: "Scheat",
    109268: "Enif",
    # Andromeda
    31775: "Alpheratz",
    # Cassiopeia
    3179: "Schedar",
}

# 古代中国都城坐标字典（经度, 纬度），单位：度
# 供skyfield使用
capitals: dict[str, tuple[float, float]] = {
    'luoyang':   (112.4, 34.7),
    'nanjing':   (118.8, 32.0),
    'linzhang':  (114.6, 36.3),
    'xian':      (108.9, 34.2),
    'kaifeng':   (114.6, 34.5),
    'hangzhou':  (120.2, 30.2),
    'beijing':   (116.4, 39.6)
}

# ── Plot style constants ─────────────────────────────────────────────
_STYLE_BG = "#06061a"      # figure / axes background
_STYLE_TEXT = "#ffffff"     # label text colour
_STYLE_SPINE = "#333355"   # axes spine colour


# ── Helpers ──────────────────────────────────────────────────────────


def _resolve_body(eph: Any, name: str) -> Any:
    """Resolve a body from the ephemeris, trying aliases."""
    aliases = _PLANET_ALIASES.get(name.lower().strip(), [name.lower().strip()])
    for alias in aliases:
        try:
            return eph[alias]
        except (KeyError, ValueError):
            continue
    logger.warning("body not resolvable in ephemeris name=%s", name)
    return None


def _parse_date(date_str: str) -> datetime:
    """Parse a local date string ``"YYYY-MM-DD"`` into a UTC datetime.

    The input is treated as a local calendar date.  No timezone
    conversion is performed — sunrise/sunset calculations use the
    observer's lat/lon to determine local solar time, so the date
    components are what matter for Skyfield.

    The time is set to UTC midnight (00:00); actual observation
    times are determined by sunrise/sunset calculations, not by
    this function.

    Returns:
        UTC datetime at midnight.

    Raises:
        ValueError: If the string is not a valid ISO date.
    """
    dt = datetime.fromisoformat(date_str.strip())
    return dt.replace(tzinfo=timezone.utc)


def _angular_sep_deg(
    ra1_hours: float, dec1_deg: float, ra2_hours: float, dec2_deg: float
) -> float:
    """Angular separation (degrees) between two equatorial positions.

    Uses the spherical law of cosines.  RA is in hours, Dec in degrees.
    """
    r1 = math.radians(ra1_hours * 15.0)
    d1 = math.radians(dec1_deg)
    r2 = math.radians(ra2_hours * 15.0)
    d2 = math.radians(dec2_deg)
    cos_d = (
        math.sin(d1) * math.sin(d2)
        + math.cos(d1) * math.cos(d2) * math.cos(r1 - r2)
    )
    return math.degrees(math.acos(min(1.0, max(-1.0, cos_d))))


def _radec_to_altaz(
    ra_deg: Any,
    dec_deg: Any,
    lat: float,
    lon: float,
    t_mid: Any,
) -> tuple[Any, Any]:
    """Convert equatorial (RA, Dec) to horizontal (alt, az).

    Manual geometry — avoids Skyfield ``.apparent()`` which can
    trigger deflector computation out of ephemeris range for
    many stars.  The astrometric vs apparent difference is
    negligible (a few arc-seconds).

    Parameters:
        ra_deg: Right ascension array (degrees).
        dec_deg: Declination array (degrees).
        lat: Observer latitude (degrees).
        lon: Observer longitude (degrees).
        t_mid: Skyfield ``Time`` object (for GAST).

    Returns:
        ``(alt_deg, az_deg)`` arrays (degrees).
    """
    import numpy as np

    ra_deg = np.atleast_1d(ra_deg)
    dec_deg = np.atleast_1d(dec_deg)
    lst = t_mid.gast * 15.0 + lon
    ha_rad = np.radians((lst - ra_deg) % 360)
    dec_rad = np.radians(dec_deg)
    lat_rad = math.radians(lat)

    sin_alt = (
        np.sin(dec_rad) * math.sin(lat_rad)
        + np.cos(dec_rad) * np.cos(lat_rad) * np.cos(ha_rad)
    )
    alt_deg = np.degrees(np.arcsin(np.clip(sin_alt, -1, 1)))

    cos_az = (np.sin(dec_rad) - math.sin(lat_rad) * sin_alt) / (
        math.cos(lat_rad) * np.cos(np.radians(alt_deg))
    )
    cos_az = np.clip(cos_az, -1, 1)
    az_deg = np.degrees(np.arccos(cos_az))
    sin_az = -np.cos(dec_rad) * np.cos(lat_rad) * np.sin(ha_rad)
    az_deg = np.where(sin_az < 0, 360 - az_deg, az_deg)
    return alt_deg, az_deg


# ── Star / constellation data loaders ──────────────────────────────


def _load_stars(
    data_dir: str,
    max_mag: float = 6.0,
) -> Any:
    """Load bright stars from the Hipparcos catalog.

    Parameters:
        data_dir: Directory containing ``hip_main.dat``.
        max_mag: Maximum visual magnitude (default 6.0 = naked
            eye limit).

    Returns:
        Skyfield ``Star`` object with all matching stars, and
        the filtered DataFrame (hip, magnitude, bv_color).
    """
    from io import BytesIO

    from skyfield.data import hipparcos
    from skyfield.starlib import Star

    hip_path = os.path.join(data_dir, "hip_main.dat")
    with open(hip_path, "rb") as f:
        raw = f.read()

    # Parse main catalog (HIP, mag, RA, Dec, …)
    df = hipparcos.load_dataframe(BytesIO(raw))

    # Parse B-V colour index from the same file
    compression: str | None = "gzip" if raw[:2] == b"\x1f\x8b" else None

    from pandas import read_csv

    bv_df = read_csv(
        BytesIO(raw),
        sep="|",
        names=hipparcos._COLUMN_NAMES,
        compression=compression,
        usecols=["HIP", "B-V"],
        na_values=["     ", "       ", "        ", "            "],
    )
    bv_df = bv_df.dropna(subset=["B-V"]).set_index("HIP")
    bv_df.columns = ["bv_color"]
    from pandas import to_numeric

    bv_df["bv_color"] = to_numeric(bv_df["bv_color"], errors="coerce")
    bv_df = bv_df.dropna()

    bright = df[df["magnitude"] <= max_mag].copy()
    bright = bright.join(bv_df, how="left")
    bright["bv_color"] = to_numeric(bright["bv_color"], errors="coerce").fillna(0.0)

    star = Star.from_dataframe(bright)
    return star, bright


def _load_constellation_lines(
    data_dir: str,
) -> list[tuple[str, list[tuple[int, int]]]]:
    """Parse Stellarium ``constellationship.fab``.

    Parameters:
        data_dir: Directory containing the file.

    Returns:
        List of ``(name, [(hip1, hip2), ...])`` tuples.
    """
    from skyfield.data import stellarium

    fab_path = os.path.join(data_dir, "constellationship.fab")
    with open(fab_path, "rb") as f:
        return stellarium.parse_constellations(f.readlines())


# B-V colour index → (R, G, B) in [0, 1].
# Piecewise-linear table from Harre & Heller (2021), Table A.1.
_BV_TABLE: list[tuple[float, tuple[float, float, float]]] = [
    (-0.40, (0.62, 0.69, 1.00)),
    (-0.30, (0.67, 0.73, 1.00)),
    (-0.20, (0.73, 0.78, 1.00)),
    (-0.10, (0.79, 0.84, 1.00)),
    ( 0.00, (0.86, 0.89, 1.00)),
    ( 0.10, (0.93, 0.93, 1.00)),
    ( 0.20, (1.00, 0.97, 0.98)),
    ( 0.30, (1.00, 0.98, 0.91)),
    ( 0.40, (1.00, 0.95, 0.81)),
    ( 0.50, (1.00, 0.91, 0.71)),
    ( 0.60, (1.00, 0.86, 0.62)),
    ( 0.70, (1.00, 0.81, 0.54)),
    ( 0.80, (1.00, 0.75, 0.46)),
    ( 0.90, (1.00, 0.69, 0.39)),
    ( 1.00, (1.00, 0.63, 0.32)),
    ( 1.10, (1.00, 0.57, 0.27)),
    ( 1.20, (1.00, 0.52, 0.22)),
    ( 1.30, (1.00, 0.47, 0.18)),
    ( 1.40, (1.00, 0.43, 0.15)),
    ( 1.50, (1.00, 0.39, 0.12)),
    ( 1.60, (1.00, 0.35, 0.10)),
    ( 1.70, (1.00, 0.32, 0.08)),
    ( 1.80, (1.00, 0.29, 0.07)),
    ( 1.90, (1.00, 0.27, 0.06)),
    ( 2.00, (1.00, 0.25, 0.05)),
]


def _bv_to_rgba_array(bv_arr):  # np.ndarray -> np.ndarray
    """Vectorised B-V → RGBA array (N, 4) with full opacity."""
    import numpy as np

    bv_t = np.array([row[0] for row in _BV_TABLE])
    r_t = np.array([row[1][0] for row in _BV_TABLE])
    g_t = np.array([row[1][1] for row in _BV_TABLE])
    b_t = np.array([row[1][2] for row in _BV_TABLE])

    bv_clamped = np.clip(bv_arr, bv_t[0], bv_t[-1])
    r = np.interp(bv_clamped, bv_t, r_t)
    g = np.interp(bv_clamped, bv_t, g_t)
    b = np.interp(bv_clamped, bv_t, b_t)
    return np.column_stack([r, g, b, np.ones_like(r)])


# ── CelestialBody base class ────────────────────────────────────────


class CelestialBody:
    """Base class for a celestial body at a given time and location.

    Subclasses must implement ``_skyfield_body`` and ``name``.
    Positions are computed eagerly at construction time; coordinate
    properties are lazily cached.

    Parameters:
        observer: The Observer instance.
        date: ``"YYYY-MM-DD"`` date string.  Defaults to the
            Observer's bound date if available.
        dt: Timezone-aware datetime (alternative to *date*).
    """

    def __init__(
        self,
        observer: Observer,
        date: str | None = None,
        dt: datetime | None = None,
    ) -> None:
        if date is not None:
            dt_val = _parse_date(date)
            has_explicit_time = False
        elif dt is not None:
            if dt.tzinfo is None:
                raise ValueError("datetime must be timezone-aware")
            dt_val = dt.astimezone(timezone.utc)
            has_explicit_time = True
        elif observer._dt_val is not None:
            dt_val = observer._dt_val
            has_explicit_time = False
        else:
            raise ValueError("either date or dt is required")

        year, month, day = dt_val.year, dt_val.month, dt_val.day
        hour, minute = dt_val.hour, dt_val.minute
        second = dt_val.second + dt_val.microsecond / 1e6

        # Default to 20:00 local time when only a date was given
        if not has_explicit_time:
            tz_delta = timedelta(hours=observer.tz)
            local_dt = datetime(
                year, month, day, 20, 0, 0,
                tzinfo=timezone(tz_delta),
            )
            utc_dt = local_dt.astimezone(timezone.utc)
            year, month, day = utc_dt.year, utc_dt.month, utc_dt.day
            hour, minute = utc_dt.hour, utc_dt.minute
            second = 0.0

        self._observer = observer
        self._year = year
        self._month = month
        self._day = day
        self._hour = hour
        self._minute = minute
        self._second = second

        # Reuse Observer's Loader and Timescale, and prefer cached eph
        self._loader = observer._loader
        self._eph = observer._get_eph(year)
        self._ts = observer._ts

        self._t = self._ts.utc(year, month, day, hour, minute, second)

        obs_vector = observer._to_earth_topos(self._eph)
        self._apparent = obs_vector.at(self._t).observe(self._skyfield_body).apparent()

        # Lazy caches
        self._altaz: Any = None
        self._radec: Any = None
        self._constellation: str | None = None
        self._constellation_map: Any = None

    @property
    def _skyfield_body(self) -> Any:
        """Skyfield body object.  Must be overridden by subclasses."""
        raise NotImplementedError

    @property
    def name(self) -> str:
        """Short English name (e.g. ``"Moon"``, ``"mars"``)."""
        raise NotImplementedError

    # ── Horizontal coordinates ───────────────────────────────────────

    @property
    def alt(self) -> float:
        """Altitude above horizon (degrees)."""
        if self._altaz is None:
            self._altaz = self._apparent.altaz()
        return self._altaz[0].degrees

    @property
    def az(self) -> float:
        """Azimuth (degrees, north=0, east=90)."""
        if self._altaz is None:
            self._altaz = self._apparent.altaz()
        return self._altaz[1].degrees

    # ── Equatorial coordinates ───────────────────────────────────────

    @property
    def ra(self) -> float:
        """Right ascension (hours)."""
        if self._radec is None:
            self._radec = self._apparent.radec()
        return self._radec[0].hours

    @property
    def dec(self) -> float:
        """Declination (degrees)."""
        if self._radec is None:
            self._radec = self._apparent.radec()
        return self._radec[1].degrees

    # ── Constellation ────────────────────────────────────────────────

    @property
    def constellation(self) -> str:
        """IAU constellation (abbreviation, full name)."""
        if self._constellation is None:
            if self._constellation_map is None:
                self._constellation_map = load_constellation_map()
            self._constellation = self._constellation_map(self._apparent)
        return self._constellation

    # ── Visibility ───────────────────────────────────────────────────

    @property
    def is_above_horizon(self) -> bool:
        """Whether the body is above the horizon."""
        return self.alt > 0

    # ── Snapshot ─────────────────────────────────────────────────────

    def to_position(self) -> dict[str, Any]:
        """Return a JSON-serialisable position snapshot."""
        return {
            "name": self.name,
            "altitude": self.alt,
            "azimuth": self.az,
            "ra_hours": self.ra,
            "dec_deg": self.dec,
            "constellation": self.constellation,
            "above_horizon": self.is_above_horizon,
        }

    @property
    def observer(self) -> Observer:
        """The Observer that created this body."""
        return self._observer

    @property
    def t(self) -> Any:
        """Skyfield Time object."""
        return self._t

    # ── Trajectory ───────────────────────────────────────────────────

    def trajectory(self, n_points: int = 100) -> list[dict[str, Any]]:
        """Sample the body's path from sunset to sunrise.

        Parameters:
            n_points: Number of sample points.

        Returns:
            List of dicts with keys:
            ``time``, ``altitude``, ``azimuth``, ``ra``,
            ``dec``, ``constellation``.
        """
        sunset_t, sunrise_t = self._observer._find_sun_times(
            self._year,
            self._month,
            self._day,
        )
        if sunset_t is None or sunrise_t is None:
            logger.debug("sunrise/sunset not found, returning empty trajectory")
            return []

        jd_start = sunset_t.tt
        jd_end = sunrise_t.tt

        if self._constellation_map is None:
            self._constellation_map = load_constellation_map()

        results: list[dict[str, Any]] = []
        obs_vector = self._observer._to_earth_topos(self._eph)

        for i in range(n_points):
            frac = i / (n_points - 1) if n_points > 1 else 0
            jd = jd_start + frac * (jd_end - jd_start)
            t_i = self._ts.tt(jd=jd)

            app = obs_vector.at(t_i).observe(self._skyfield_body).apparent()
            altaz = app.altaz()
            radec = app.radec()
            const = self._constellation_map(app)

            results.append(
                {
                    "time": t_i,
                    "altitude": altaz[0].degrees,
                    "azimuth": altaz[1].degrees,
                    "ra": radec[0].hours,
                    "dec": radec[1].degrees,
                    "constellation": const,
                }
            )

        return results


# ── Moon ─────────────────────────────────────────────────────────────


class Moon(CelestialBody):
    """The Moon — computes position at construction time.

    Parameters:
        observer: The Observer instance.
        date: ``"YYYY-MM-DD"`` date string.  Defaults to the
            Observer's bound date if available.
        dt: Timezone-aware datetime (alternative to *date*).

    Example::

        obs = Observer(date="1054-07-04", lat=34.8, lon=114.3)
        moon = Moon(obs)
        print(moon.alt, moon.phase_name)
    """

    @property
    def _skyfield_body(self) -> Any:
        return self._eph["moon"]

    @property
    def name(self) -> str:
        return "Moon"

    # ── Moon-specific properties ─────────────────────────────────────

    @property
    def phase_angle(self) -> float:
        """Phase angle (degrees).  0° = new, 180° = full."""
        return almanac.moon_phase(self._eph, self._t).degrees

    @property
    def illumination(self) -> float:
        """Illuminated fraction (0–1)."""
        moon_app = (
            self._observer._to_earth_topos(self._eph)
            .at(self._t)
            .observe(self._eph["moon"])
            .apparent()
        )
        return moon_app.fraction_illuminated(self._eph["sun"])

    @property
    def phase_name(self) -> str:
        """Moon phase name (English)."""
        return _MOON_PHASE_NAMES[int((self.phase_angle % 360) / 45) % 8]

    @property
    def mag(self) -> float:
        """Apparent magnitude (approximate).

        Uses the Astronomical Almanac formula where *phi* is the
        Sun–Moon–Earth angle (0° = full, 180° = new).
        """
        phi = abs(self.phase_angle - 180.0)
        return -12.73 + 0.026 * phi + 4.0e-9 * phi**4

    @property
    def size(self) -> float:
        """Angular diameter (arc-seconds)."""
        if self._radec is None:
            self._radec = self._apparent.radec()
        dist_km = self._radec[2].km
        return 2 * math.degrees(math.atan(1737.4 / dist_km)) * 3600


# ── Planet ───────────────────────────────────────────────────────────


class Planet(CelestialBody):
    """A planet — computes position at construction time.

    Parameters:
        observer: The Observer instance.
        name: Planet name
            (mercury, venus, mars, jupiter, saturn, uranus, neptune).
        date: ``"YYYY-MM-DD"`` date string.  Defaults to the
            Observer's bound date if available.
        dt: Timezone-aware datetime (alternative to *date*).

    Example::

        obs = Observer(date="1054-07-04", lat=34.8, lon=114.3)
        mars = Planet(obs, "mars")
        print(mars.alt, mars.az)
    """

    def __init__(
        self,
        observer: Observer,
        name: str,
        date: str | None = None,
        dt: datetime | None = None,
    ) -> None:
        self._name = name.lower().strip()
        super().__init__(observer, date=date, dt=dt)
        if self._sf_body is None:
            raise ValueError(f"body not resolvable in ephemeris: {name}")

    @property
    def _skyfield_body(self) -> Any:
        if not hasattr(self, "_sf_body"):
            self._sf_body = _resolve_body(self._eph, self._name)
        return self._sf_body

    @property
    def name(self) -> str:
        return self._name

# --------------------------------
# Observer
# --------------------------------

class Observer:
    """An observer at a geographic location on a specific date.

    A geographic point on Earth's surface used as the origin for
    astronomical calculations.  A *date* is required at construction
    time so the Observer is always time-bound.

    Parameters:
        date: Observation date in ``"YYYY-MM-DD"`` format
            (e.g. ``"1054-07-04"``).
        lat: Latitude (degrees, north positive). Defaults to Guiyang.
        lon: Longitude (degrees, east positive). Defaults to Guiyang.
        elevation: Elevation above sea level (metres). Defaults to Guiyang.
        tz: UTC offset in hours (positive = east).
            Defaults to 8.0 (China Standard Time).
        directory: Skyfield ephemeris directory.
            Defaults to ``GSTARPY_DATA_DIR`` env var if ``None``.

    Example::

        obs = Observer(date="1054-07-04", lat=34.8, lon=114.3)
        snap = obs.snapshot()
        moon = obs.Moon()
    """

    def __init__(
        self,
        date: str,
        lat: float = 26.65,
        lon: float = 106.63,
        elevation: float = 1100.0,
        tz: float = 8.0,
        directory: str | None = None,
    ) -> None:
        self.lat = lat
        self.lon = lon
        self.elevation = elevation
        self.tz = tz
        directory = directory or os.environ.get("GSTARPY_DATA_DIR")
        if directory is None:
            raise ValueError("directory or GSTARPY_DATA_DIR env var must be provided")
        self._loader = Loader(directory)
        self._ts = self._loader.timescale()
        self._ts.julian_calendar_cutoff = _JULIAN_CUTOFF_JD

        # Parse date and preload ephemeris
        self._dt_val: datetime = _parse_date(date)
        self._eph: Any = self._load_eph(self._dt_val.year)

    def _load_eph(self, year: int) -> Any:
        """Open the ephemeris file for *year*, with caching."""
        if year <= 1969:
            return self._loader("de441_part-1.bsp")
        return self._loader("de421.bsp")

    def _get_eph(self, year: int) -> Any:
        """Return ephemeris for *year*, reusing preloaded copy when possible."""
        if self._eph is not None and year == self._dt_val.year:  # type: ignore[union-attr]
            return self._eph
        return self._load_eph(year)

    def _to_earth_topos(self, eph: Any) -> Any:
        """Create a Skyfield observer vector (earth + topos)."""
        topos = wgs84.latlon(self.lat, self.lon, elevation_m=self.elevation)
        return eph["earth"] + topos

    def _night_midpoint(self, dt_val: datetime | None = None) -> Any:
        """Return the Skyfield ``Time`` at the midpoint of tonight's night.

        Defaults to the Observer's bound date if *dt_val* is ``None``.

        Falls back to midnight UTC when sunset/sunrise data is
        unavailable (e.g. polar regions in summer/winter).
        """
        if dt_val is None:
            dt_val = self._dt_val
        sun_info = self.find_sun_times(dt_val.year, dt_val.month, dt_val.day)
        sunset_t = sun_info["sunset"]
        sunrise_t = sun_info["sunrise"]
        if sunset_t is not None and sunrise_t is not None:
            return self._ts.tt_jd((sunset_t.tt + sunrise_t.tt) / 2.0)
        return self._ts.utc(dt_val.year, dt_val.month, dt_val.day, 0)

    # ── Body factory methods ─────────────────────────────────────────

    def Moon(self, date: str | None = None) -> Moon:
        """Create a Moon object, defaulting to the Observer's date.

        Example::

            obs = Observer(date="1054-07-04", lat=34.8, lon=114.3)
            moon = obs.Moon()
            print(moon.phase_name)

        You may override the date explicitly::

            moon = obs.Moon("1066-12-25 00:00")
        """
        return Moon(self, date=date)

    def Planet(self, name: str, date: str | None = None) -> Planet:
        """Create a Planet object, defaulting to the Observer's date.

        Example::

            obs = Observer(date="1054-07-04", lat=34.8, lon=114.3)
            mars = obs.Planet("mars")
            print(mars.alt, mars.az)

        You may override the date explicitly::

            mars = obs.Planet("mars", "1066-12-25 00:00")
        """
        return Planet(self, name, date=date)

    # ── Core API ─────────────────────────────────────────────────────

    def snapshot(
        self,
        date: str | None = None,
        planets: list[str] | None = None,
        span_days: int = 0,
    ) -> dict[str, Any]:
        """Generate a night-sky position snapshot.

        When the Observer was created with a *date*, that date is
        used by default and the *date* argument may be omitted.

        Parameters:
            span_days: When > 0, scan ±*span_days* around the given
                date for the best Moon–planet closest approaches.
                This corrects for UTC calendar mismatch: an event
                recorded on local date D may actually occur on UTC
                date D−1 (if in the early morning hours).  Defaults
                to 0 (single-night snapshot).

        Returns:
            Dict with keys: ``year``, ``month``, ``day``,
            ``observer_lat``, ``observer_lon``, ``moon``,
            ``planets``, ``visible_planets``, ``sunset_utc``,
            ``sunrise_utc``, ``night_hours``, ``moon_phase_deg``,
            ``moon_illumination``, ``moon_phase_name``,
            ``moon_mag``, ``moon_size_arcsec``,
            ``visible_constellations``.
            When *span_days* > 0, ``moon_closest_approaches``
            includes ``date_offset`` and ``search_date`` keys.

        Example::

            obs = Observer(date="1054-07-04", lat=34.8, lon=114.3)
            snap = obs.snapshot()
            print(snap["moon_phase_name"])

            # With UTC calendar mismatch correction:
            snap = obs.snapshot(date=(1499, 8, 22), span_days=1)
            ca = snap["moon_closest_approaches"]["mars"]
            print(ca["separation_deg"], ca["date_offset"], ca["search_date"])
        """
        if date is not None:
            dt_val = _parse_date(date)
        else:
            dt_val = self._dt_val
        moon = Moon(self, dt=dt_val)

        planet_bodies: dict[str, Planet] = {}
        planet_positions: dict[str, dict[str, Any]] = {}
        for name in planets or DEFAULT_PLANETS:
            try:
                p = Planet(self, name, dt=dt_val)
                planet_bodies[name] = p
                planet_positions[name] = p.to_position()
            except (ValueError, KeyError) as exc:
                logger.debug(
                    "planet position failed name=%s error=%s",
                    name,
                    exc,
                )

        sunset_t, sunrise_t = self._find_sun_times(
            dt_val.year,
            dt_val.month,
            dt_val.day,
        )
        night_h = None
        if sunset_t is not None and sunrise_t is not None:
            night_h = (sunrise_t.tt - sunset_t.tt) * 24.0

        # Angular separation between Moon and each visible planet
        separations: dict[str, float] = {}
        for name, pos in planet_positions.items():
            if not pos["above_horizon"]:
                continue
            sep = _angular_sep_deg(
                moon.ra, moon.dec, pos["ra_hours"], pos["dec_deg"]
            )
            separations[name] = round(sep, 2)

        # Closest approach between Moon and each planet over the night(s)
        if span_days > 0:
            closest = self._closest_approaches_multi_day(
                dt_val, planets or DEFAULT_PLANETS, span_days,
            )
        else:
            closest = self._closest_approaches_single_night(
                moon, planet_positions, planet_bodies,
            )

        return {
            "year": dt_val.year,
            "month": dt_val.month,
            "day": dt_val.day,
            "observer_lat": self.lat,
            "observer_lon": self.lon,
            "moon": moon.to_position(),
            "planets": planet_positions,
            "visible_planets": [
                p for p in planet_positions.values() if p["above_horizon"]
            ],
            "sunset_utc": sunset_t,
            "sunrise_utc": sunrise_t,
            "night_hours": night_h,
            "moon_phase_deg": moon.phase_angle,
            "moon_illumination": moon.illumination,
            "moon_phase_name": moon.phase_name,
            "moon_mag": round(moon.mag, 2),
            "moon_size_arcsec": round(moon.size, 1),
            "moon_planet_separations_deg": separations,
            "moon_closest_approaches": closest,
            "visible_constellations": self._visible_constellations(dt_val),
        }

    # ── Closest-approach helpers ─────────────────────────────────────

    @staticmethod
    def _closest_approaches_single_night(
        moon: Moon,
        planet_positions: dict[str, dict[str, Any]],
        planet_bodies: dict[str, Planet],
    ) -> dict[str, dict[str, Any]]:
        """Compute closest Moon–planet approaches within a single night."""
        n_traj = 200
        moon_traj = moon.trajectory(n_points=n_traj)
        closest: dict[str, dict[str, Any]] = {}
        for name in planet_positions:
            planet_body = planet_bodies.get(name)
            if planet_body is None:
                continue
            try:
                planet_traj = planet_body.trajectory(n_points=n_traj)
            except (ValueError, KeyError):
                continue
            if not planet_traj or not moon_traj:
                continue

            best_sep = float("inf")
            best_idx = 0
            for j in range(min(len(moon_traj), len(planet_traj))):
                m_pt = moon_traj[j]
                p_pt = planet_traj[j]
                if m_pt["altitude"] <= 0 and p_pt["altitude"] <= 0:
                    continue
                sep_val = _angular_sep_deg(
                    m_pt["ra"], m_pt["dec"], p_pt["ra"], p_pt["dec"]
                )
                if sep_val < best_sep:
                    best_sep = sep_val
                    best_idx = j

            if best_sep < float("inf"):
                m_best = moon_traj[best_idx]
                p_best = planet_traj[best_idx]
                t_best = m_best["time"]
                closest[name] = {
                    "separation_deg": round(best_sep, 2),
                    "time_utc": (
                        f"{int(t_best.utc[3]):02d}:{int(t_best.utc[4]):02d}"
                    ),
                    "moon_alt": round(m_best["altitude"], 2),
                    "planet_alt": round(p_best["altitude"], 2),
                    "moon_ra": round(m_best["ra"], 4),
                    "moon_dec": round(m_best["dec"], 2),
                    "planet_ra": round(p_best["ra"], 4),
                    "planet_dec": round(p_best["dec"], 2),
                }
        return closest

    def _closest_approaches_multi_day(
        self,
        dt_val: datetime,
        planets: list[str],
        span_days: int,
    ) -> dict[str, dict[str, Any]]:
        """Find best Moon–planet approach across ±*span_days* nights.

        Corrects for UTC calendar mismatch: an event on local date D
        may occur during the night of UTC date D−1 (early morning).
        """
        from datetime import timedelta as td  # noqa: PLC0415

        merged: dict[str, dict[str, Any]] = {}
        for delta in range(-span_days, span_days + 1):
            if delta == 0:
                d = dt_val
            else:
                d = dt_val + td(days=delta)
            moon = Moon(self, date=d)
            planet_bodies: dict[str, Planet] = {}
            planet_positions: dict[str, dict[str, Any]] = {}
            for name in planets:
                try:
                    p = Planet(self, name, date=d)
                    planet_bodies[name] = p
                    planet_positions[name] = p.to_position()
                except (ValueError, KeyError):
                    pass
            night_closest = self._closest_approaches_single_night(
                moon, planet_positions, planet_bodies,
            )
            for name, ca in night_closest.items():
                sep = ca["separation_deg"]
                prev = merged.get(name)
                if prev is None or sep < prev["separation_deg"]:
                    ca["date_offset"] = delta
                    ca["search_date"] = (
                        f"{d.year}.{d.month}.{d.day}"
                    )
                    merged[name] = ca

        return merged

    def _visible_constellations(
        self, dt_val: datetime
    ) -> list[dict[str, str]]:
        """Return constellations with ≥1 star above the horizon at *dt_val*.

        Returns:
            List of ``{"abbrev": "...", "name": "..."}`` dicts.
        """
        data_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "data", "Skyfield")
        data_dir = os.path.normpath(data_dir)

        centers = _compute_all_constellation_centers(data_dir)
        star_info = _compute_star_positions(self, dt_val)
        ra_arr = star_info["ra_deg"]
        dec_arr = star_info["dec_deg"]
        alt_arr = star_info["alt_deg"]

        visible: list[dict[str, str]] = []
        for abbr, en_name in _CONSTELLATIONS:
            if abbr not in centers:
                continue
            ra_c, dec_c, ra_h, dec_h = centers[abbr]
            for i in range(len(ra_arr)):
                if alt_arr[i] <= 0:
                    continue
                dra = abs(ra_arr[i] - ra_c)
                if dra > 180:
                    dra = 360 - dra
                if dra <= ra_h and abs(dec_arr[i] - dec_c) <= dec_h:
                    visible.append({"abbrev": abbr, "name": en_name})
                    break
        return visible

    def positions_at(
        self,
        dt: datetime | None = None,
        planets: list[str] | None = None,
    ) -> dict[str, dict[str, Any]]:
        """Compute all body positions at a specific time.

        When the Observer was created with a *date*, *dt* defaults
        to that date.

        Parameters:
            dt: Time (must be timezone-aware).  Defaults to the
                Observer's bound date.
            planets: Planet names to compute.

        Returns:
            ``{name: position_dict}`` mapping.
        """
        if dt is not None:
            if dt.tzinfo is None:
                raise ValueError("datetime must be timezone-aware")
            utc_dt = dt.astimezone(timezone.utc)
        elif self._dt_val is not None:
            utc_dt = self._dt_val
        else:
            raise ValueError("dt is required when Observer has no date")
        moon = Moon(self, dt=utc_dt)
        result: dict[str, dict[str, Any]] = {"Moon": moon.to_position()}

        for name in planets or DEFAULT_PLANETS:
            try:
                p = Planet(self, name, dt=utc_dt)
                result[name] = p.to_position()
            except (ValueError, KeyError) as exc:
                logger.debug(
                    "planet position failed name=%s error=%s",
                    name,
                    exc,
                )
        return result

    def moon_phase_at(
        self,
        year: int | None = None,
        month: int | None = None,
        day: int | None = None,
        hour: int = 0,
        minute: int = 0,
    ) -> dict[str, float]:
        """Query the moon phase at a specific time.

        When the Observer was created with a *date*, all arguments
        default to that date.

        Returns:
            ``{"phase_deg": float, "illumination": float}``
        """
        if year is None:
            if self._dt_val is None:
                raise ValueError("date arguments are required when Observer has no date")
            year, month, day = self._dt_val.year, self._dt_val.month, self._dt_val.day
        moon = Moon(self, date=(year, month, day, hour, minute))
        return {
            "phase_deg": moon.phase_angle,
            "illumination": moon.illumination,
        }

    def moon_phases_in_range(
        self,
        year: int,
        month_start: int,
        month_end: int,
    ) -> list[dict[str, Any]]:
        """Query discrete moon-phase events in a date range.

        Returns:
            List of dicts with keys:
            ``phase``, ``time_utc``, ``phase_code``.
        """
        eph = self._get_eph(year)
        ts = self._ts
        t0 = ts.utc(year, month_start, 1)
        t1 = ts.utc(
            year + (1 if month_end >= 12 else 0),
            1 if month_end >= 12 else month_end + 1,
            1,
        )
        f = almanac.moon_phases(eph)
        times, codes = almanac.find_discrete(t0, t1, f)
        return [
            {
                "phase": almanac.MOON_PHASES[c],
                "time_utc": t,
                "phase_code": c,
            }
            for t, c in zip(times, codes)
        ]

    def find_sun_times(
        self,
        year: int | None = None,
        month: int | None = None,
        day: int | None = None,
    ) -> dict[str, Any | None]:
        """Query sunset and sunrise for a given date.

        When the Observer was created with a *date*, all arguments
        default to that date.

        Returns:
            ``{"sunset": Time|None, "sunrise": Time|None,
            "night_hours": float|None}``
        """
        if year is None:
            if self._dt_val is None:
                raise ValueError("date arguments are required when Observer has no date")
            year, month, day = self._dt_val.year, self._dt_val.month, self._dt_val.day
        sunset_t, sunrise_t = self._find_sun_times(year, month, day)
        night_h = None
        if sunset_t is not None and sunrise_t is not None:
            night_h = (sunrise_t.tt - sunset_t.tt) * 24.0
        return {
            "sunset": sunset_t,
            "sunrise": sunrise_t,
            "night_hours": night_h,
        }

    # ── Internal ─────────────────────────────────────────────────────

    def _find_sun_times(
        self,
        year: int,
        month: int,
        day: int,
    ) -> tuple[Any | None, Any | None]:
        """Find sunset and next-day sunrise.

        Returns:
            ``(sunset_time, sunrise_time)`` — ``None`` if not found.
        """
        eph = self._get_eph(year)
        ts = self._ts
        topos = wgs84.latlon(self.lat, self.lon, elevation_m=self.elevation)
        sun = eph["sun"]

        t = ts.utc(year, month, day, 0)
        f = almanac.risings_and_settings(eph, sun, topos, horizon_degrees=0.0)

        sunset_t: Any | None = None
        sunrise_t: Any | None = None
        try:
            times, events = almanac.find_discrete(t, ts.tt(jd=t.tt + 2.0), f)
            for time, event in zip(times, events):
                if not event and sunset_t is None:
                    # False = setting
                    sunset_t = time
                elif event and sunset_t is not None:
                    # True = rising after sunset
                    sunrise_t = time
                    break
        except (ValueError, IndexError) as exc:
            logger.debug("sun times lookup failed error=%s", exc)

        return sunset_t, sunrise_t


# ── Plotting ─────────────────────────────────────────────────────────


def _draw_stars_and_constellations(
    ax: Any,
    observer: Observer,
    dt_val: datetime,
    wcs: Any,
    ra_center: float,
    dec_center: float,
    fov: float,
    dec_span: float,
    naxis1: int,
    naxis2: int,
    *,
    show_stars: bool,
    show_constellations: bool,
    label_stars: bool,
    legend_handles: list[Any],
) -> None:
    """Draw naked-eye stars and constellation lines on *ax*.

    Loads the Hipparcos catalog and (optionally) the Stellarium
    constellation line data, computes apparent positions at the
    observation epoch, and draws those that fall within the WCS
    field of view.
    """
    import numpy as np
    from skyfield.api import wgs84 as _wgs84

    data_dir = observer._loader.directory

    # Load bright stars
    stars, star_df = _load_stars(data_dir, max_mag=6.0)

    # Compute apparent position at night midpoint
    t_mid = observer._night_midpoint(dt_val)

    eph = observer._get_eph(dt_val.year)
    earth = eph["earth"]
    topos = _wgs84.latlon(
        observer.lat,
        observer.lon,
        elevation_m=observer.elevation,
    )
    observer_vec = earth + topos

    app = observer_vec.at(t_mid).observe(stars)
    # Skip .apparent() — it triggers deflector computation
    # that can go out of ephemeris range for many stars.
    # For stars the astrometric vs apparent difference is
    # negligible (a few arc-seconds).
    ra_obj, dec_obj, _ = app.radec()

    ra_deg = np.atleast_1d(ra_obj._degrees) % 360
    dec_deg = np.atleast_1d(dec_obj._degrees)
    alt_deg, az_deg = _radec_to_altaz(
        ra_deg, dec_deg, observer.lat, observer.lon, t_mid,
    )

    # Field-of-view filter
    half_ra = fov / 2.0
    half_dec = dec_span / 2.0
    dra = (ra_deg - ra_center + 180) % 360 - 180
    ddec = dec_deg - dec_center
    in_fov = (alt_deg > 0) & (np.abs(dra) < half_ra) & (np.abs(ddec) < half_dec)

    if not in_fov.any():
        return

    mag_arr = star_df["magnitude"].values

    # ── Draw stars ──
    if show_stars:
        sizes = 60.0 * 10.0 ** (-0.4 * (mag_arr[in_fov] - 1.0))
        sizes = np.clip(sizes, 3.0, 120.0)

        bv_arr = star_df["bv_color"].values.astype(float)
        star_colors = _bv_to_rgba_array(bv_arr[in_fov])
        star_colors[:, 3] = 0.92  # global alpha

        ax.scatter(
            ra_deg[in_fov],
            dec_deg[in_fov],
            transform=ax.get_transform("world"),
            s=sizes,
            facecolor=star_colors,
            edgecolor="none",
            zorder=2,
        )

    # ── Draw constellation lines ──
    if show_constellations:
        from matplotlib.collections import LineCollection

        hip_visible = set(
            star_df.index.values[in_fov].tolist(),
        )
        hip_to_radec: dict[int, tuple[float, float]] = {}
        for idx in np.where(in_fov)[0]:
            hip_id = int(star_df.index.values[idx])
            hip_to_radec[hip_id] = (
                float(ra_deg[idx]),
                float(dec_deg[idx]),
            )

        const_data = _load_constellation_lines(data_dir)
        segments: list[list[tuple[float, float]]] = []
        for _name, edges in const_data:
            for hip1, hip2 in edges:
                if hip1 in hip_visible and hip2 in hip_visible:
                    ra1, dec1 = hip_to_radec[hip1]
                    ra2, dec2 = hip_to_radec[hip2]
                    if abs(ra1 - ra2) > 180:
                        continue
                    segments.append([(ra1, dec1), (ra2, dec2)])

        if segments:
            lc = LineCollection(
                segments,
                colors="#4499dd",
                linewidths=1.8,
                alpha=0.8,
                zorder=1,
                transform=ax.get_transform("world"),
            )
            ax.add_collection(lc)

        legend_handles.append(
            ax.scatter(
                [],
                [],
                marker=".",
                s=30,
                color="white",
                label="Stars (mag < 6)",
            ),
        )

    # ── Label named stars ──
    if label_stars and show_stars:
        fov_indices = np.where(in_fov)[0]
        for idx in fov_indices:
            hip_id = int(star_df.index.values[idx])
            star_name = _STAR_NAMES.get(hip_id)
            if star_name is None:
                continue
            ra_val = float(ra_deg[idx])
            dec_val = float(dec_deg[idx])
            mag_val = float(mag_arr[idx])
            # Brighter stars get larger font
            fontsize = max(8.0, 12.0 - 1.2 * mag_val)
            ax.text(
                ra_val,
                dec_val,
                f" {star_name}",
                transform=ax.get_transform("world"),
                color=_STYLE_TEXT,
                fontsize=fontsize,
                alpha=0.9,
                ha="left",
                va="bottom",
                zorder=5,
            )


# ── Standalone star/constellation plots ────────────────────────────


def _compute_star_positions(
    observer: Observer,
    dt_val: datetime,
) -> dict[str, Any]:
    """Compute star RA/Dec/alt/az arrays at the night midpoint."""
    import numpy as np
    from skyfield.api import wgs84 as _wgs84

    data_dir = observer._loader.directory
    stars, star_df = _load_stars(data_dir, max_mag=6.0)

    t_mid = observer._night_midpoint(dt_val)

    eph = observer._get_eph(dt_val.year)
    earth = eph["earth"]
    topos = _wgs84.latlon(
        observer.lat, observer.lon, elevation_m=observer.elevation,
    )
    app = (earth + topos).at(t_mid).observe(stars)
    ra_obj, dec_obj, _ = app.radec()

    ra_deg = np.atleast_1d(ra_obj._degrees) % 360
    dec_deg = np.atleast_1d(dec_obj._degrees)
    alt_deg, az_deg = _radec_to_altaz(
        ra_deg, dec_deg, observer.lat, observer.lon, t_mid,
    )

    return {
        "star_df": star_df,
        "ra_deg": ra_deg,
        "dec_deg": dec_deg,
        "alt_deg": alt_deg,
        "az_deg": az_deg,
        "mag": star_df["magnitude"].values,
        "bv_color": star_df["bv_color"].values.astype(float),
    }


def plot_stars(
    observer: Observer,
    date: str,
    *,
    fov: float = 20.0,
    center: tuple[float, float] | None = None,
    show_constellations: bool = True,
    label_stars: bool = True,
    save_path: str | Path | None = None,
) -> Any:
    """Plot visible naked-eye stars on a plain RA/Dec chart.

    Only stars above the horizon at the night midpoint are drawn.
    Constellation lines and star-name labels are optional.

    Parameters:
        observer: Observer instance.
        date: Date string in ``"YYYY-MM-DD"`` format.
        fov: Field of view diameter in degrees (default 20).
        center: ``(ra_hours, dec_deg)`` field centre.  Defaults to
            the zenith at the night midpoint.
        show_constellations: Draw Stellarium constellation lines.
        label_stars: Annotate named bright stars.
        save_path: If given, save the figure to this path.

    Returns:
        matplotlib Figure.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.collections import LineCollection

    dt_val = _parse_date(date)
    data_dir = observer._loader.directory
    sp = _compute_star_positions(observer, dt_val)

    ra_all = sp["ra_deg"]
    dec_all = sp["dec_deg"]
    alt_all = sp["alt_deg"]
    mag_all = sp["mag"]
    bv_all = sp["bv_color"]
    star_df = sp["star_df"]

    # Default centre: zenith at night midpoint
    if center is not None:
        ra_c = center[0] * 15.0
        dec_c = center[1]
    else:
        t_mid = observer._night_midpoint(dt_val)
        lst = t_mid.gast * 15.0 + observer.lon
        ra_c = lst % 360
        dec_c = observer.lat

    half = fov / 2.0
    dra = (ra_all - ra_c + 180) % 360 - 180
    visible = (alt_all > 0) & (np.abs(dra) < half) & (np.abs(dec_all - dec_c) < half)

    if not visible.any():
        raise ValueError("No visible stars in the specified field of view.")

    ra_v = dra[visible]
    dec_v = dec_all[visible]
    mag_v = mag_all[visible]
    bv_v = bv_all[visible]
    idx_v = np.where(visible)[0]

    # ── Figure ──
    fig, ax = plt.subplots(figsize=(10, 10), facecolor=_STYLE_BG)
    ax.set_facecolor(_STYLE_BG)

    sizes = 60.0 * 10.0 ** (-0.4 * (mag_v - 1.0))
    sizes = np.clip(sizes, 3.0, 120.0)
    colors = _bv_to_rgba_array(bv_v)
    colors[:, 3] = 0.92

    ax.scatter(ra_v, dec_v, s=sizes, facecolor=colors, edgecolor="none", zorder=2)

    # ── Constellation lines ──
    if show_constellations:
        hip_vis = set(star_df.index.values[visible].tolist())
        hip_pos: dict[int, tuple[float, float]] = {}
        for i, idx in enumerate(idx_v):
            hip_pos[int(star_df.index.values[idx])] = (float(dra[idx]), float(dec_all[idx]))

        segments: list[list[tuple[float, float]]] = []
        for _name, edges in _load_constellation_lines(data_dir):
            for h1, h2 in edges:
                if h1 in hip_vis and h2 in hip_vis:
                    segments.append([hip_pos[h1], hip_pos[h2]])

        if segments:
            lc = LineCollection(segments, colors="#4499dd",
                                linewidths=1.8, alpha=0.8, zorder=1)
            ax.add_collection(lc)

    # ── Star labels ──
    if label_stars:
        for i, idx in enumerate(idx_v):
            hip_id = int(star_df.index.values[idx])
            name = _STAR_NAMES.get(hip_id)
            if name is None or mag_v[i] > 3.5:
                continue
            fs = max(8.0, 12.0 - 1.2 * mag_v[i])
            ax.text(float(dra[idx]), float(dec_all[idx]),
                    f"  {name}", color=_STYLE_TEXT, fontsize=fs,
                    alpha=0.9, zorder=3, va="center")

    ax.set_xlim(half, -half)
    ax.set_ylim(dec_c - half, dec_c + half)
    ax.set_xlabel("RA offset (deg)", color="white", fontsize=10)
    ax.set_ylabel("Dec (deg)", color="white", fontsize=10)
    ax.set_title(
        f"Visible Stars  {dt_val.strftime('%Y-%m-%d')}",
        color="white", fontsize=16, pad=12,
    )
    ax.tick_params(colors="white", labelsize=8)
    ax.grid(True, color=_STYLE_SPINE, alpha=0.35, linewidth=0.3)
    for spine in ax.spines.values():
        spine.set_edgecolor(_STYLE_SPINE)

    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight", facecolor=_STYLE_BG)
    return fig


# Built-in constellation data with pre-computed center/spans.
_CONSTELLATIONS: list[tuple[str, str]] = [
    # (abbreviation, english_name)
    ("UMi", "Ursa Minor"), ("Dra", "Draco"), ("Cep", "Cepheus"),
    ("Cas", "Cassiopeia"), ("Cam", "Camelopardalis"),
    ("UMa", "Ursa Major"), ("Lyn", "Lynx"), ("CVn", "Canes Venatici"),
    ("LMi", "Leo Minor"), ("Per", "Perseus"), ("Aur", "Auriga"),
    ("Lyr", "Lyra"), ("Cyg", "Cygnus"), ("Vul", "Vulpecula"),
    ("Sge", "Sagitta"), ("Del", "Delphinus"), ("Equ", "Equuleus"),
    ("Lac", "Lacerta"), ("Her", "Hercules"), ("CrB", "Corona Borealis"),
    ("Boo", "Boötes"), ("Com", "Coma Berenices"),
    ("And", "Andromeda"), ("Tri", "Triangulum"), ("Peg", "Pegasus"),
    ("Ari", "Aries"), ("Tau", "Taurus"), ("Gem", "Gemini"),
    ("Cnc", "Cancer"), ("Leo", "Leo"), ("Vir", "Virgo"),
    ("Lib", "Libra"), ("Sco", "Scorpius"), ("Sgr", "Sagittarius"),
    ("Cap", "Capricornus"), ("Aqr", "Aquarius"), ("Psc", "Pisces"),
    ("Ori", "Orion"), ("CMi", "Canis Minor"), ("CMa", "Canis Major"),
    ("Mon", "Monoceros"), ("Aql", "Aquila"), ("Oph", "Ophiuchus"),
    ("Ser", "Serpens"), ("Sct", "Scutum"), ("Hya", "Hydra"),
    ("Crt", "Crater"), ("Crv", "Corvus"), ("Sex", "Sextans"),
    ("Lep", "Lepus"), ("Col", "Columba"), ("Eri", "Eridanus"),
    ("For", "Fornax"), ("Scl", "Sculptor"), ("PsA", "Piscis Austrinus"),
    ("Cet", "Cetus"), ("Ant", "Antlia"), ("Pyx", "Pyxis"),
    ("Cru", "Crux"), ("Car", "Carina"), ("Vel", "Vela"),
    ("Pup", "Puppis"), ("Cen", "Centaurus"), ("Lup", "Lupus"),
    ("CrA", "Corona Australis"), ("TrA", "Triangulum Australe"),
    ("Mus", "Musca"), ("Nor", "Norma"), ("Ara", "Ara"),
    ("Cir", "Circinus"), ("Aps", "Apus"), ("Cha", "Chamaeleon"),
    ("Pav", "Pavo"), ("Gru", "Grus"), ("Phe", "Phoenix"),
    ("Tuc", "Tucana"), ("Hyi", "Hydrus"), ("Dor", "Dorado"),
    ("Vol", "Volans"), ("Men", "Mensa"), ("Mic", "Microscopium"),
    ("Tel", "Telescopium"), ("Ind", "Indus"), ("Oct", "Octans"),
    ("Ret", "Reticulum"), ("Hor", "Horologium"), ("Cae", "Caelum"),
    ("Pic", "Pictor"),
]
_CONST_ABBREVS: dict[str, str] = {c[0]: c[1] for c in _CONSTELLATIONS}
_CONST_ABBREVS_LOWER: dict[str, str] = {c[0].lower(): c[0] for c in _CONSTELLATIONS}


def plot_constellation(
    abbrev: str,
    *,
    observer: Observer | None = None,
    save_path: str | Path | None = None,
) -> Any:
    """Plot a single constellation with all its stars and lines.

    Draws all Hipparcos stars (mag ≤ 6) in the constellation area,
    regardless of horizon visibility, using Stellarium constellation
    lines and B-V star colours.

    Parameters:
        abbrev: IAU three-letter abbreviation (e.g. ``"Ori"``,
            ``"Sco"``, ``"cru"`` — case-insensitive).
        observer: If given, use its data directory.  Defaults to
            the default data path.
        save_path: If given, save the figure to this path.

    Returns:
        matplotlib Figure.

    Raises:
        ValueError: Unknown constellation abbreviation.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.collections import LineCollection

    # Resolve abbreviation (case-insensitive)
    key = abbrev.strip()
    if key in _CONST_ABBREVS:
        canon = key
    elif key.lower() in _CONST_ABBREVS_LOWER:
        canon = _CONST_ABBREVS_LOWER[key.lower()]
    else:
        raise ValueError(
            f"Unknown constellation '{abbrev}'. "
            f"Valid: {', '.join(sorted(_CONST_ABBREVS))}"
        )
    english_name = _CONST_ABBREVS[canon]

    data_dir = (
        observer._loader.directory if observer is not None
        else os.environ.get("GSTARPY_DATA_DIR", "data/Skyfield")
    )

    centers = _compute_all_constellation_centers(data_dir)
    if canon not in centers:
        raise ValueError(f"No center data for {canon} ({english_name}).")
    ra_h, dec_c, ra_half, dec_half = centers[canon]
    ra_c = ra_h * 15.0

    # Load all stars
    stars, star_df = _load_stars(data_dir, max_mag=6.0)
    loader = Loader(data_dir)
    ts = loader.timescale()
    earth = loader("de421.bsp")["earth"]
    t = ts.utc(2026, 1, 1, 0, 0)
    app = earth.at(t).observe(stars)
    ra_obj, dec_obj, _ = app.radec()
    ra_all = np.asarray(ra_obj._degrees) % 360
    dec_all = np.asarray(dec_obj._degrees)
    mag_all = star_df["magnitude"].values
    bv_all = star_df["bv_color"].values.astype(float)

    dra = (ra_all - ra_c + 180) % 360 - 180
    in_box = (np.abs(dra) < ra_half) & (np.abs(dec_all - dec_c) < dec_half)

    if not in_box.any():
        raise ValueError(f"No stars found for {canon} ({english_name}).")

    ra_v = dra[in_box]
    dec_v = dec_all[in_box]
    mag_v = mag_all[in_box]
    bv_v = bv_all[in_box]
    idx_v = np.where(in_box)[0]

    # ── Figure ──
    aspect = ra_half / dec_half if dec_half > 0 else 1.0
    fig_h = 8.0
    fig_w = max(6.0, min(16.0, fig_h * aspect))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor=_STYLE_BG)
    ax.set_facecolor(_STYLE_BG)

    sizes = 60.0 * 10.0 ** (-0.4 * (mag_v - 1.0))
    sizes = np.clip(sizes, 3.0, 120.0)
    colors = _bv_to_rgba_array(bv_v)
    colors[:, 3] = 0.92
    ax.scatter(ra_v, dec_v, s=sizes, facecolor=colors, edgecolor="none", zorder=2)

    # Constellation lines
    hip_in = set(star_df.index.values[in_box].tolist())
    hip_pos: dict[int, tuple[float, float]] = {}
    for i, idx in enumerate(idx_v):
        hip_pos[int(star_df.index.values[idx])] = (float(dra[idx]), float(dec_all[idx]))

    segments: list[list[tuple[float, float]]] = []
    for const_name, edges in _load_constellation_lines(data_dir):
        if const_name != canon:
            continue
        for h1, h2 in edges:
            if h1 in hip_in and h2 in hip_in:
                segments.append([hip_pos[h1], hip_pos[h2]])

    if segments:
        lc = LineCollection(segments, colors="#4499dd",
                            linewidths=1.8, alpha=0.8, zorder=1)
        ax.add_collection(lc)

    # Star labels
    for i, idx in enumerate(idx_v):
        hip_id = int(star_df.index.values[idx])
        name = _STAR_NAMES.get(hip_id)
        if name is None or mag_v[i] > 3.5:
            continue
        fs = max(8.0, 12.0 - 1.2 * mag_v[i])
        ax.text(float(dra[idx]), float(dec_all[idx]),
                f"  {name}", color=_STYLE_TEXT, fontsize=fs,
                alpha=0.9, zorder=3, va="center")

    ax.set_xlim(ra_half, -ra_half)
    ax.set_ylim(dec_c - dec_half, dec_c + dec_half)
    ax.set_xlabel("RA offset (deg)", color="white", fontsize=10)
    ax.set_ylabel("Dec (deg)", color="white", fontsize=10)
    ax.set_title(english_name, color="white", fontsize=16, pad=12)
    ax.tick_params(colors="white", labelsize=8)
    ax.grid(True, color=_STYLE_SPINE, alpha=0.35, linewidth=0.3)
    for spine in ax.spines.values():
        spine.set_edgecolor(_STYLE_SPINE)

    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight", facecolor=_STYLE_BG)
    return fig


def plot_night(
    observer: Observer,
    date: str,
    *,
    bodies: list[str] | None = None,
    n_points: int = 100,
    fov: float = 20.0,
    center: tuple[float, float] | None = None,
    show_stars: bool = False,
    show_constellations: bool = False,
    label_stars: bool = False,
    save_path: str | Path | None = None,
) -> plt.Figure:
    """Plot night-sky trajectories.

    By default the field of view follows the Moon's path from
    sunset to sunrise.  Pass *center=(ra_hours, dec_deg)* to
    point at a specific region instead.

    The Moon is drawn as a filled disc with a shadow overlay
    showing the phase.  Other bodies are plotted as simple
    markers if they fall within the field.

    Parameters:
        observer: Observer instance.
        date: Date string in ``"YYYY-MM-DD"`` format
            (e.g. ``"1054-07-04"``).
        bodies: Body names to plot.  Defaults to Moon + five
            classical planets.
        n_points: Trajectory sample points.
        fov: Field of view width in degrees (default 20).
        center: ``(ra_hours, dec_deg)`` to override the default
            Moon-based field centre.
        show_stars: Overlay naked-eye stars from the Hipparcos
            catalog.
        show_constellations: Draw constellation stick-figure
            lines (requires Stellarium data file).
        label_stars: Annotate bright stars with common names
            (requires ``show_stars=True``).
        save_path: If given, save the figure to this path.

    Returns:
        matplotlib Figure with WCS Tan axes.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.wcs import WCS
    from matplotlib.patches import Circle

    dt_val = _parse_date(date)

    # Build body objects
    names = bodies if bodies is not None else ["moon", *DEFAULT_PLANETS]
    body_objects: list[CelestialBody] = []
    for name in names:
        try:
            if name.lower() == "moon":
                body_objects.append(Moon(observer, dt=dt_val))
            else:
                body_objects.append(Planet(observer, name, dt=dt_val))
        except (ValueError, KeyError):
            logger.debug("skipping body, not resolvable name=%s", name)

    if not body_objects:
        raise ValueError("no valid celestial bodies to plot")

    # Gather trajectories
    body_traj: list[tuple[CelestialBody, list[dict[str, Any]]]] = []
    for body in body_objects:
        pts = body.trajectory(n_points=n_points)
        # Keep only bodies with at least some visible points
        if pts and any(p["altitude"] > 0 for p in pts):
            body_traj.append((body, pts))

    if not body_traj:
        raise ValueError("no trajectories computed")

    # Find the moon object and its trajectory
    moon_obj: Moon | None = None
    moon_pts: list[dict[str, Any]] = []
    for body, pts in body_traj:
        if isinstance(body, Moon):
            moon_obj = body
            moon_pts = pts
            break

    # Field centre — explicit override or Moon trajectory default
    if center is not None:
        ra_center = (center[0] * 15.0) % 360  # hours → degrees
        dec_center = center[1]
    elif moon_pts:
        ra_vis = [p["ra"] * 15.0 for p in moon_pts if p["altitude"] > 0]
        dec_vis = [p["dec"] for p in moon_pts if p["altitude"] > 0]
        if ra_vis:
            ra_arr = np.array(ra_vis)
            dec_arr = np.array(dec_vis)
        else:
            mid = len(moon_pts) // 2
            ra_arr = np.array([moon_pts[mid]["ra"] * 15.0])
            dec_arr = np.array([moon_pts[mid]["dec"]])
        # Circular mean for RA (handles wrap-around)
        ra_sin = float(np.mean(np.sin(np.radians(ra_arr))))
        ra_cos = float(np.mean(np.cos(np.radians(ra_arr))))
        ra_center = np.degrees(np.arctan2(ra_sin, ra_cos)) % 360
        dec_center = float(np.mean(dec_arr))
    else:
        ra_center = body_traj[0][0].ra * 15.0
        dec_center = body_traj[0][0].dec

    dec_span = fov * 0.75
    naxis1 = naxis2 = 800
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [naxis1 / 2, naxis2 / 2]
    wcs.wcs.cdelt = [fov / naxis1, dec_span / naxis2]
    wcs.wcs.crval = [ra_center, dec_center]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    BG = _STYLE_BG

    fig = plt.figure(figsize=(10, 9), facecolor=BG)
    ax = fig.add_subplot(111, projection=wcs)
    ax.set_facecolor(BG)
    ax.set_xlim(0.5, naxis1 + 0.5)
    ax.set_ylim(0.5, naxis2 + 0.5)

    legend_handles: list[Any] = []

    # ── Stars and constellations ──
    if show_stars or show_constellations:
        _draw_stars_and_constellations(
            ax,
            observer,
            dt_val,
            wcs,
            ra_center,
            dec_center,
            fov,
            dec_span,
            naxis1,
            naxis2,
            show_stars=show_stars,
            show_constellations=show_constellations,
            label_stars=label_stars,
            legend_handles=legend_handles,
        )

    # ── Plot each body ──
    n_bodies = max(1, len(body_traj))
    color_cycle = plt.cm.tab10(np.linspace(0, 1, n_bodies))

    for i, (body, pts) in enumerate(body_traj):
        ra_vis = [p["ra"] * 15.0 for p in pts if p["altitude"] > 0]
        dec_vis = [p["dec"] for p in pts if p["altitude"] > 0]
        ra_all = [p["ra"] * 15.0 for p in pts]
        dec_all = [p["dec"] for p in pts]

        # Check if any visible point is within the FOV
        half_ra = fov / 2.0
        half_dec = dec_span / 2.0
        in_fov = any(
            (abs((ra - ra_center + 180) % 360 - 180) < half_ra
             and abs(dec - dec_center) < half_dec)
            for ra, dec in zip(ra_vis, dec_vis)
        )

        # Trajectory line + shadow area for Moon
        if len(ra_vis) > 1:
            if isinstance(body, Moon) and moon_obj is not None:
                # Shadow: overlapping circles along trajectory
                moon_r = moon_obj.size / 3600.0 / 2.0  # arcsec → deg (radius)
                for ra, dec in zip(ra_vis, dec_vis):
                    circ = Circle(
                        (ra, dec),
                        radius=moon_r,
                        transform=ax.get_transform("world"),
                        facecolor="#334455",
                        edgecolor="none",
                        alpha=0.5,
                        zorder=1,
                    )
                    ax.add_patch(circ)

                ax.plot(
                    ra_vis,
                    dec_vis,
                    transform=ax.get_transform("world"),
                    color="#667788",
                    linewidth=1.0,
                    alpha=0.9,
                )
            else:
                ax.plot(
                    ra_vis,
                    dec_vis,
                    transform=ax.get_transform("world"),
                    color=color_cycle[i],
                    linewidth=1.5,
                    alpha=0.7,
                )

        if isinstance(body, Moon) and moon_obj is not None:
            # ── Moon: shadow area only ──
            ra_deg = body.ra * 15.0
            dec_deg = body.dec
            moon_px = ax.wcs.world_to_pixel_values(ra_deg, dec_deg)
            moon_r_world = moon_obj.size / 3600.0 / 2.0  # arcsec → deg (radius)
            moon_px2 = ax.wcs.world_to_pixel_values(
                ra_deg + moon_r_world, dec_deg,
            )
            moon_r_px = abs(moon_px2[0] - moon_px[0])

            phase = moon_obj.phase_angle % 360
            illum = moon_obj.illumination

            # Draw the shadow portion of the moon disc
            if illum > 0.02:
                cos_ph = abs(math.cos(math.radians(phase)))
                shadow_r = moon_r_px * cos_ph
                shadow_shift = moon_r_px * (2.0 * illum - 1.0)

                if phase < 180:
                    offset_x = moon_px[0] - shadow_shift
                else:
                    offset_x = moon_px[0] + shadow_shift
                offset_y = moon_px[1]

                shadow = Circle(
                    (offset_x, offset_y),
                    radius=shadow_r,
                    transform=ax.transData,
                    facecolor="#0a0a20",
                    edgecolor="none",
                    alpha=0.9,
                    zorder=11,
                    clip_on=True,
                )
                ax.add_patch(shadow)

            # Legend proxy (invisible scatter)
            if in_fov:
                legend_handles.append(
                    ax.scatter(
                        [],
                        [],
                        marker="o",
                        s=100,
                        facecolor="#ffffcc",
                        edgecolor="#aabb88",
                        linewidth=1.0,
                        label=(f"Moon ({body.constellation})"),
                    )
                )
        else:
            # ── Other bodies: simple marker ──
            ax.scatter(
                ra_all,
                dec_all,
                transform=ax.get_transform("world"),
                marker="o",
                s=20,
                color=color_cycle[i],
                alpha=0.5,
                zorder=4,
            )
            if in_fov:
                legend_handles.append(
                    ax.scatter(
                        [],
                        [],
                        marker="o",
                        s=80,
                        color=color_cycle[i],
                        edgecolors="#cccccc",
                        linewidth=0.8,
                        label=(f"{body.name} ({body.constellation})"),
                    )
                )

    ax.coords.grid(color="#222244", alpha=0.6, linestyle="--")
    ax.coords[0].set_axislabel("RA", color="#8888aa")
    ax.coords[1].set_axislabel("Dec", color="#8888aa")
    ax.coords[0].tick_params(colors="#555577")
    ax.coords[1].tick_params(colors="#555577")
    ax.legend(
        handles=legend_handles,
        loc="lower right",
        fontsize="small",
        facecolor=BG,
        edgecolor="#333366",
        labelcolor="#ccccdd",
    )

    date_str = f"{dt_val.year}-{dt_val.month:02d}-{dt_val.day:02d}"
    names_str = ", ".join(b.name for b, _ in body_traj)
    ax.set_title(f"{date_str}  ({names_str})", color="#aabbdd", fontweight="bold")

    if save_path:
        fig.savefig(str(save_path), dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())

    return fig


def _compute_all_constellation_centers(
    data_dir: str,
) -> dict[str, tuple[float, float, float, float]]:
    """Compute (ra_hours, dec_deg, ra_half, dec_half) for all constellations.

    Loads the star catalog once and reuses it for every constellation,
    avoiding repeated file I/O and Skyfield computations.
    """
    import numpy as np

    stars, star_df = _load_stars(data_dir, max_mag=6.0)
    const_data = _load_constellation_lines(data_dir)

    loader = Loader(data_dir)
    ts = loader.timescale()
    earth = loader("de421.bsp")["earth"]
    t = ts.utc(2026, 1, 1, 0, 0)
    app = earth.at(t).observe(stars)
    ra_obj, dec_obj, _ = app.radec()
    ra_all = np.asarray(ra_obj._degrees) % 360
    dec_all = np.asarray(dec_obj._degrees)

    hip_map: dict[int, tuple[float, float]] = {}
    for i, hip in enumerate(star_df.index.values):
        hip_map[int(hip)] = (float(ra_all[i]), float(dec_all[i]))

    # Group edges by constellation
    edges_by_const: dict[str, list[tuple[int, int]]] = {}
    for const_name, edges in const_data:
        edges_by_const.setdefault(const_name, []).extend(edges)

    result: dict[str, tuple[float, float, float, float]] = {}
    for canon, _english in _CONSTELLATIONS:
        edges = edges_by_const.get(canon, [])
        all_hips: set[int] = set()
        for h1, h2 in edges:
            all_hips.add(h1)
            all_hips.add(h2)

        if not all_hips:
            result[canon] = (0.0, 0.0, 15.0, 15.0)
            continue

        positions = [hip_map[h] for h in all_hips if h in hip_map]
        ra_arr = np.array([p[0] for p in positions])
        dec_arr = np.array([p[1] for p in positions])
        ra_mean = (
            np.degrees(
                np.arctan2(np.mean(np.sin(np.radians(ra_arr))),
                           np.mean(np.cos(np.radians(ra_arr))))
            )
            % 360
        )
        dec_mean = float(np.mean(dec_arr))
        dra = (ra_arr - ra_mean + 180) % 360 - 180
        ra_half = float(np.max(np.abs(dra)) + 3.0)
        dec_half = float(np.max(np.abs(dec_arr - dec_mean)) + 3.0)
        result[canon] = (ra_mean / 15.0, dec_mean, ra_half, dec_half)

    return result


def plot_visible_constellations(
    observer: Observer,
    date: str,
    *,
    save_path: str | Path | None = None,
) -> Any:
    """Plot all constellations visible at the given date and location.

    Computes star positions at the night midpoint, identifies which
    constellations have stars above the horizon, then draws their
    stars and constellation lines on a single RA/Dec chart.

    Parameters:
        observer: Observer instance.
        date: Date string in ``"YYYY-MM-DD"`` format.
        save_path: If given, save the figure to this path.

    Returns:
        matplotlib Figure.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.collections import LineCollection

    dt_val = _parse_date(date)
    data_dir = observer._loader.directory

    # Compute star positions at the night midpoint (once)
    sp = _compute_star_positions(observer, dt_val)
    ra_all = sp["ra_deg"]
    dec_all = sp["dec_deg"]
    alt_all = sp["alt_deg"]
    mag_all = sp["mag"]
    bv_all = sp["bv_color"]
    star_df = sp["star_df"]

    const_data = _load_constellation_lines(data_dir)

    # Pre-compute all constellation centers (once, single star load)
    all_centers = _compute_all_constellation_centers(data_dir)

    # Find visible constellations: any star in box is above horizon
    visible_abbrevs: list[str] = []
    for canon, _english in _CONSTELLATIONS:
        ra_h, dec_c, ra_half, dec_half = all_centers[canon]
        ra_c = ra_h * 15.0

        dra_center = (ra_all - ra_c + 180) % 360 - 180
        near_center = (
            (np.abs(dra_center) < ra_half)
            & (np.abs(dec_all - dec_c) < dec_half)
            & (alt_all > 0)
        )
        if near_center.any():
            visible_abbrevs.append(canon)

    if not visible_abbrevs:
        raise ValueError(
            "No constellations visible at the given date and location."
        )

    # Collect all visible constellation stars and segments
    seen_hips: set[int] = set()
    all_ra: list[float] = []
    all_dec: list[float] = []
    all_mag: list[float] = []
    all_bv: list[np.ndarray] = []
    all_idx: list[int] = []
    segments: list[list[tuple[float, float]]] = []

    for canon in visible_abbrevs:
        ra_h, dec_c, ra_half, dec_half = all_centers[canon]
        ra_c = ra_h * 15.0

        dra = (ra_all - ra_c + 180) % 360 - 180
        in_box = (
            (np.abs(dra) < ra_half)
            & (np.abs(dec_all - dec_c) < dec_half)
        )

        for idx in np.where(in_box)[0]:
            i = int(idx)
            hip_id = int(star_df.index.values[i])
            if hip_id not in seen_hips:
                seen_hips.add(hip_id)
                all_ra.append(float(dra[idx]))
                all_dec.append(float(dec_all[i]))
                all_mag.append(float(mag_all[i]))
                all_bv.append(bv_all[i])
                all_idx.append(i)

        # Constellation lines for this constellation
        hip_in_box = set(star_df.index.values[in_box].tolist())
        hip_pos: dict[int, tuple[float, float]] = {}
        for idx in np.where(in_box)[0]:
            i = int(idx)
            hip_id = int(star_df.index.values[i])
            hip_pos[hip_id] = (float(dra[idx]), float(dec_all[i]))

        for const_name, edges in const_data:
            if const_name != canon:
                continue
            for h1, h2 in edges:
                if h1 in hip_in_box and h2 in hip_in_box:
                    segments.append([hip_pos[h1], hip_pos[h2]])

    ra_arr = np.array(all_ra)
    dec_arr = np.array(all_dec)
    mag_arr = np.array(all_mag)
    bv_arr = np.array(all_bv)

    # Compute plot bounds from all visible stars
    ra_span = (float(np.max(ra_arr)) - float(np.min(ra_arr))) / 2.0 + 3.0
    dec_span = (float(np.max(dec_arr)) - float(np.min(dec_arr))) / 2.0 + 3.0

    # ── Figure ──
    aspect = ra_span / dec_span if dec_span > 0 else 1.0
    fig_h = 8.0
    fig_w = max(6.0, min(16.0, fig_h * aspect))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor=_STYLE_BG)
    ax.set_facecolor(_STYLE_BG)

    sizes = 60.0 * 10.0 ** (-0.4 * (mag_arr - 1.0))
    sizes = np.clip(sizes, 3.0, 120.0)
    colors = _bv_to_rgba_array(bv_arr)
    colors[:, 3] = 0.92

    ax.scatter(ra_arr, dec_arr, s=sizes, facecolor=colors,
               edgecolor="none", zorder=2)

    if segments:
        lc = LineCollection(segments, colors="#4499dd",
                            linewidths=1.8, alpha=0.8, zorder=1)
        ax.add_collection(lc)

    # Star labels
    for i, idx in enumerate(all_idx):
        hip_id = int(star_df.index.values[idx])
        name = _STAR_NAMES.get(hip_id)
        if name is None or mag_arr[i] > 3.5:
            continue
        fs = max(8.0, 12.0 - 1.2 * mag_arr[i])
        ax.text(float(ra_arr[i]), float(dec_arr[i]),
                f"  {name}", color=_STYLE_TEXT, fontsize=fs,
                alpha=0.9, zorder=3, va="center")

    ax.set_xlim(ra_span, -ra_span)
    ax.set_ylim(-dec_span, dec_span)
    ax.set_xlabel("RA offset (deg)", color="white", fontsize=10)
    ax.set_ylabel("Dec (deg)", color="white", fontsize=10)
    ax.set_title(
        f"Visible Constellations  {dt_val.strftime('%Y-%m-%d')}",
        color="white", fontsize=16, pad=12,
    )
    ax.tick_params(colors="white", labelsize=8)
    ax.grid(True, color=_STYLE_SPINE, alpha=0.35, linewidth=0.3)
    for spine in ax.spines.values():
        spine.set_edgecolor(_STYLE_SPINE)

    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight",
                    facecolor=_STYLE_BG)
    return fig
