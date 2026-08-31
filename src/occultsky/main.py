"""
Reproducing Ancient Astronomical Records with Skyfield.
"""

import logging
import os
from datetime import datetime, timedelta, timezone
from typing import Optional

from dotenv import load_dotenv
from skyfield import almanac
from skyfield.api import Loader, wgs84
from skyfield.api import GREGORIAN_START
from skyfield.framelib import ecliptic_frame

load_dotenv()

logger = logging.getLogger(__name__)


class Observer:
    """
    Represents an astronomical observation site on Earth.
    """

    def __init__(
        self, 
        lat: float,
        lon: float,
        elevation: float,
        skyfield_data: Optional[str] = None,
        ephemeris: str = 'de441.bsp'
    ):
        """
        Initialize an observer at a given geographic location.

        For ancient dates, consider using DE441 which has longer coverage.

        Args:
            lat: Latitude in degrees (north positive).
            lon: Longitude in degrees (east positive). Also used to calculate the Local Mean Time (LMT) offset.
            elevation: Elevation above sea level in meters.
            skyfield_data: Path to Skyfield data directory. Falls back to SKYFIELD_DATA environment variable if not provided.
            ephemeris: JPL ephemeris filename (default: de441.bsp).

        Raises:
            ValueError: If skyfield_data directory is not specified.
        """
        # Set Skyfield data directory.
        skyfield_data = skyfield_data or os.getenv('SKYFIELD_DATA')
        if skyfield_data is None:
            raise ValueError("Skyfield data directory not specified. Pass it as a parameter or set SKYFIELD_DATA environment variable")
        load = Loader(skyfield_data)

        # Load the timescale.
        self.ts = load.timescale()
        # Dates before 1582-10-15 use the Julian calendar
        self.ts.julian_calendar_cutoff = GREGORIAN_START

        # Load the JPL ephemeris.
        self.eph = load(ephemeris)

        # Geographic position (from Earth's center to location)
        self.topos = wgs84.latlon(lat, lon, elevation_m=elevation)
        # Observer vector (from Solar System center to location)
        self.observer = self.eph['earth'] + self.topos

        # timezone offset: each degree of longitude corresponds to 4 minutes
        self.tz = timezone(timedelta(seconds=round(lon * 240)))


    def get_solarnoon(self, date_str: str) -> str:
        """
        Return the time of true solar noon (Sun's upper meridian transit) for the given local date.

        Args:
            date_str: Local date in ``YYYY-MM-DD`` format.

        Returns:
            LMT ISO string (``YYYY-MM-DDTHH:MM:SS±HH:MM``).

        Raises:
            ValueError: If ``date_str`` is not ``YYYY-MM-DD`` format, or if no solar meridian transit is found.
        """
        try:
            date = datetime.strptime(date_str, "%Y-%m-%d")
        except ValueError:
            raise ValueError(
                f"date_str must be 'YYYY-MM-DD', got: {date_str!r}"
            ) from None
        # 24h window from local midnight covers the full local day
        start = datetime(date.year, date.month, date.day, tzinfo=self.tz)
        end = start + timedelta(days=1)

        t0 = self.ts.from_datetime(start)
        t1 = self.ts.from_datetime(end)
        times = almanac.find_transits(self.observer, self.eph["Sun"], t0, t1)
        if len(times) == 0:
            raise ValueError(
                f"No solar meridian transit found on {date_str}"
            )
        solar_noon = times[0].astimezone(self.tz)
        return solar_noon.isoformat()


    def get_twilight(self, date_str: str) -> list:
        """Return astronomical twilight times for the given local date.

        Searches within solar noon ± 1 day for astronomical dawn and dusk.

        Args:
            date_str: Local date in ``YYYY-MM-DD`` format.

        Returns:
            List of 4 local ISO time strings in fixed order:
            [previous_dusk, dawn, dusk, next_dawn]
        """
        # Get solar noon time
        solar_noon_str = self.get_solarnoon(date_str)
        solar_noon = datetime.fromisoformat(solar_noon_str)

        # Search window: solar noon ± 1 day
        start = solar_noon - timedelta(days=1)
        end = solar_noon + timedelta(days=1)

        f = almanac.dark_twilight_day(self.eph, self.topos)
        t0 = self.ts.from_datetime(start)
        t1 = self.ts.from_datetime(end)
        times, events = almanac.find_discrete(t0, t1, f)

        # Collect dusk and dawn times separately
        dusks = []
        dawns = []
        previous_e = f(t0).item()
        for ti, ei in zip(times, events):
            local_iso = ti.astimezone(self.tz).isoformat()
            if previous_e == 0 and ei == 1:
                dawns.append(local_iso)
            elif previous_e == 1 and ei == 0:
                dusks.append(local_iso)
            previous_e = ei

        # Return in fixed order: [previous_dusk, dawn, dusk, next_dawn]
        return [dusks[0], dawns[0], dusks[1], dawns[1]]


    def get_moon_phase(self, date_str: str) -> dict:
        """Return moon phase and illumination at solar noon + 12h for the given local date.

        Args:
            date_str: Local date in ``YYYY-MM-DD`` format.

        Returns:
            Dict with keys: ``phase_angle`` (0–360°), ``illumination`` (0–1), ``phase_name`` (string).
        """
        # Get solar noon and add 12 hours for midnight observation
        solar_noon_str = self.get_solarnoon(date_str)
        solar_noon = datetime.fromisoformat(solar_noon_str)
        midnight = solar_noon + timedelta(hours=12)

        t = self.ts.from_datetime(midnight)
        sun, moon, earth = self.eph['sun'], self.eph['moon'], self.eph['earth']
        
        e = earth.at(t)
        s = e.observe(sun).apparent()
        m = e.observe(moon).apparent()

        _, slon, _ = s.frame_latlon(ecliptic_frame)
        _, mlon, _ = m.frame_latlon(ecliptic_frame)
        phase_angle = (mlon.degrees - slon.degrees) % 360.0

        illumination = float(m.fraction_illuminated(sun))

        idx = round(phase_angle / 90.0) % 4
        phase_name = almanac.MOON_PHASES[idx]

        return {
            "phase_angle": round(phase_angle, 2),
            "illumination": round(illumination, 4),
            "phase_name": phase_name,
        }


    # def moon(self, dt: datetime):
    #     astrometric = self.observer.at(t).observe(self.eph['moon'])
    #     ra, dec, _ = astrometric.apparent().radec()
    #     return {
    #         "right_ascension": float(ra.hours),
    #         "declination": float(dec.degrees),
    #     }

    # def planet(self, planet, dt: datetime) -> dict:

    #     astrometric = self.observer.at(t).observe(self.eph[planet])
    #     return astrometric.apparent().radec()

    # def get_angular_distance(
    #     self, 
    #     body1: str, 
    #     body2: str, 
    #     dt: datetime
    # ) -> float:
    #     """Return angular separation (degrees) between two celestial bodies.

    #     Args:
    #         body1: Name of the first body (e.g. ``"Moon"``, ``"Sun"``).
    #         body2: Name of the second body (e.g. ``"Jupiter"``).
    #         dt: Timezone-aware datetime (any offset is accepted).

    #     Returns:
    #         Angular distance in degrees.
    #     """
    #     try:
    #         b1 = self.eph[body1]
    #         b2 = self.eph[body2]
    #     except KeyError as e:
    #         raise KeyError(
    #             f"Unknown body {e.args[0]!r}; available: "
    #             f"{sorted(self.eph.keys())}"
    #         ) from None
    #     pos1 = self.observer.at(t).observe(b1)
    #     pos2 = self.observer.at(t).observe(b2)
    #     return float(pos1.apparent().separation_from(pos2.apparent()).degrees)