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

from .utils.calc import moon_angular_radius
from .utils.constants import PLANET_MAP

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


    def get_solarnoon(self, date_str: str) -> datetime:
        """
        Return the time of true solar noon (Sun's upper meridian transit) for the given local date.

        Args:
            date_str: Local date in ``YYYY-MM-DD`` format.

        Returns:
            Timezone-aware datetime object.

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
        return solar_noon


    def get_twilight(self, date_str: str) -> list[datetime]:
        """Return astronomical twilight times for the given local date.

        Searches within solar noon ± 1 day for astronomical dawn and dusk.

        Args:
            date_str: Local date in ``YYYY-MM-DD`` format.

        Returns:
            List of 4 timezone-aware datetime objects in fixed order:
            [previous_dusk, dawn, dusk, next_dawn]
        """
        # Get solar noon time
        solar_noon = self.get_solarnoon(date_str)

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
            local_dt = ti.astimezone(self.tz)
            if previous_e == 0 and ei == 1:
                dawns.append(local_dt)
            elif previous_e == 1 and ei == 0:
                dusks.append(local_dt)
            previous_e = ei

        # Return in fixed order: [previous_dusk, dawn, dusk, next_dawn]
        return [dusks[0], dawns[0], dusks[1], dawns[1]]


    def get_moonrise(self, date_str: str) -> list[datetime]:
        """Return the time of moonrise for the given local date.

        Searches within solar noon ± 1 day for moonrise.

        Args:
            date_str: Local date in ``YYYY-MM-DD`` format.

        Returns:
            List of timezone-aware datetime objects.
        """
        # Get solar noon time
        solar_noon = self.get_solarnoon(date_str)

        # Search window: solar noon ± 1 day
        start = solar_noon - timedelta(days=1)
        end = solar_noon + timedelta(days=1)

        moon = self.eph['moon']
        t0 = self.ts.from_datetime(start)
        t1 = self.ts.from_datetime(end)

        t, y = almanac.find_risings(self.observer, moon, t0, t1)
        return [ti.astimezone(self.tz) for ti in t]


    def get_moonset(self, date_str: str) -> list[datetime]:
        """Return the time of moonset for the given local date.

        Searches within solar noon ± 1 day for moonset.

        Args:
            date_str: Local date in ``YYYY-MM-DD`` format.

        Returns:
            List of timezone-aware datetime objects.
        """
        # Get solar noon time
        solar_noon = self.get_solarnoon(date_str)

        # Search window: solar noon ± 1 day
        start = solar_noon - timedelta(days=1)
        end = solar_noon + timedelta(days=1)

        moon = self.eph['moon']
        t0 = self.ts.from_datetime(start)
        t1 = self.ts.from_datetime(end)

        t, y = almanac.find_settings(self.observer, moon, t0, t1)
        return [ti.astimezone(self.tz) for ti in t]


    def get_planet_rise(self, planet: str, date_str: str) -> list[datetime]:
        """Return the time of rising for a planet on the given local date.

        Searches within solar noon ± 1 day for planet rise.

        Args:
            planet: Planet name (e.g. ``'mars'``, ``'jupiter'``, ``'saturn'``).
            date_str: Local date in ``YYYY-MM-DD`` format.

        Returns:
            List of timezone-aware datetime objects.

        Raises:
            KeyError: If planet name is not found in ephemeris.
        """
        # Get solar noon time
        solar_noon = self.get_solarnoon(date_str)

        # Search window: solar noon ± 1 day
        start = solar_noon - timedelta(days=1)
        end = solar_noon + timedelta(days=1)

        # Map planet name to ephemeris key
        planet_key = PLANET_MAP.get(planet.lower(), planet.lower())
        try:
            body = self.eph[planet_key]
        except (KeyError, ValueError):
            raise KeyError(f"Unknown planet {planet!r}") from None
        t0 = self.ts.from_datetime(start)
        t1 = self.ts.from_datetime(end)

        t, y = almanac.find_risings(self.observer, body, t0, t1)
        return [ti.astimezone(self.tz) for ti in t]


    def get_planet_set(self, planet: str, date_str: str) -> list[datetime]:
        """Return the time of setting for a planet on the given local date.

        Searches within solar noon ± 1 day for planet set.

        Args:
            planet: Planet name (e.g. ``'mars'``, ``'jupiter'``, ``'saturn'``).
            date_str: Local date in ``YYYY-MM-DD`` format.

        Returns:
            List of timezone-aware datetime objects.

        Raises:
            KeyError: If planet name is not found in ephemeris.
        """
        # Get solar noon time
        solar_noon = self.get_solarnoon(date_str)

        # Search window: solar noon ± 1 day
        start = solar_noon - timedelta(days=1)
        end = solar_noon + timedelta(days=1)

        # Map planet name to ephemeris key
        planet_key = PLANET_MAP.get(planet.lower(), planet.lower())
        try:
            body = self.eph[planet_key]
        except (KeyError, ValueError):
            raise KeyError(f"Unknown planet {planet!r}") from None
        t0 = self.ts.from_datetime(start)
        t1 = self.ts.from_datetime(end)

        t, y = almanac.find_settings(self.observer, body, t0, t1)
        return [ti.astimezone(self.tz) for ti in t]


    def get_moon_phase(self, date_str: str) -> dict:
        """Return moon phase and illumination at solar noon + 12h for the given local date.

        Args:
            date_str: Local date in ``YYYY-MM-DD`` format.

        Returns:
            Dict with keys: ``phase_angle`` (0–360°), ``illumination`` (0–1), ``phase_name`` (string).
        """
        # Get solar noon and add 12 hours for midnight observation
        solar_noon = self.get_solarnoon(date_str)
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


    def get_moon_angular_radius(self, date_str: str) -> float:
        """Return the Moon's apparent angular radius in arcseconds.

        Evaluated at solar noon + 12h (local midnight) using the topocentric distance from the observer to the Moon's center, so the result matches what a ground-based observer actually sees.

        Args:
            date_str: Local date in ``YYYY-MM-DD`` format.

        Returns:
            Apparent angular radius of the Moon in arcseconds.
        """
        # Get solar noon and add 12 hours for midnight observation
        solar_noon = self.get_solarnoon(date_str)
        midnight = solar_noon + timedelta(hours=12)

        t = self.ts.from_datetime(midnight)
        moon = self.eph['moon']

        # Topocentric distance from the observer to the Moon's center
        distance_km = self.observer.at(t).observe(moon).apparent().distance().km
        return moon_angular_radius(distance_km)


    def get_visible_times(
        self,
        date_str: str,
        target: str = 'moon',
        interval: int = 15
    ) -> list[datetime]:
        """Return the times when the target is visible on the given local date.

        Samples the [previous_dusk, dawn] and [dusk, next_dawn] windows at the given interval and keeps the points where the target is above the horizon.

        Args:
            date_str: Local date in ``YYYY-MM-DD`` format.
            target: Target body name (e.g. ``'moon'``, ``'mars'``, ``'jupiter'``).
            interval: Time interval for sampling in minutes (default: 15).

        Returns:
            List of timezone-aware datetimes when the target is visible. Returns empty list if the target is not visible in either window.

        Raises:
            KeyError: If target name is not found in ephemeris.
        """
        # Get twilight times
        previous_dusk, dawn, dusk, next_dawn = self.get_twilight(date_str)

        # Get target body from ephemeris
        target_key = PLANET_MAP.get(target.lower(), target.lower())
        try:
            body = self.eph[target_key]
        except (KeyError, ValueError):
            raise KeyError(f"Unknown target {target!r}") from None

        # "Is the target above the horizon" function, evaluated on the sample grid
        f = almanac.risings_and_settings(self.eph, body, self.topos)

        # Sample each night window at the given interval and keep the points where the target is up.
        step = timedelta(minutes=interval)
        result = []
        for window_start, window_end in ((previous_dusk, dawn), (dusk, next_dawn)):
            n = (window_end - window_start) // step  # number of steps in the window
            samples = [window_start + k * step for k in range(n + 1)]
            t = self.ts.from_datetimes(samples)
            above_horizon = f(t)
            result.extend(dt for dt, visible in zip(samples, above_horizon) if visible)

        return result


    def get_altaz(
        self,
        datetimes: list[datetime],
        target: str = 'moon'
    ) -> list[tuple[float, float]]:
        """Return the altitude and azimuth of the target at the given datetimes.

        Args:
            datetimes: List of timezone-aware datetimes.
            target: Target body name (e.g. ``'moon'``, ``'mars'``, ``'jupiter'``).

        Returns:
            List of tuples containing the altitude and azimuth in degrees.

        Raises:
            KeyError: If target name is not found in ephemeris.
            ValueError: If datetimes list is empty.
        """
        if len(datetimes) == 0:
            raise ValueError("datetimes list cannot be empty")

        # Map target name to ephemeris key
        target_key = PLANET_MAP.get(target.lower(), target.lower())
        try:
            body = self.eph[target_key]
        except (KeyError, ValueError):
            raise KeyError(f"Unknown target {target!r}") from None

        # Compute altitude and azimuth for all times
        t = self.ts.from_datetimes(datetimes)
        pos = self.observer.at(t).observe(body).apparent()
        alt, az, _ = pos.altaz()
        return list(zip(alt.degrees.tolist(), az.degrees.tolist()))


    def get_radec(
        self,
        datetimes: list[datetime],
        target: str = 'moon'
    ) -> list[tuple[float, float]]:
        """Return the right ascension and declination of the target at the given datetimes.

        Args:
            datetimes: List of timezone-aware datetimes.
            target: Target body name (e.g. ``'moon'``, ``'mars'``, ``'jupiter'``).

        Returns:
            List of tuples containing the right ascension and declination in degrees.

        Raises:
            KeyError: If target name is not found in ephemeris.
            ValueError: If datetimes list is empty.
        """
        if len(datetimes) == 0:
            raise ValueError("datetimes list cannot be empty")

        # Map target name to ephemeris key
        target_key = PLANET_MAP.get(target.lower(), target.lower())
        try:
            body = self.eph[target_key]
        except (KeyError, ValueError):
            raise KeyError(f"Unknown target {target!r}") from None

        t = self.ts.from_datetimes(datetimes)
        pos = self.observer.at(t).observe(body).apparent()
        ra, dec, _ = pos.radec()
        return list(zip(ra.degrees.tolist(), dec.degrees.tolist()))


    def get_angular_distance(
        self, 
        body1: str, 
        body2: str, 
        datetimes: list[datetime]
    ) -> list[float]:
        """Return angular separation (degrees) between two celestial bodies.

        Args:
            body1: Name of the first body (e.g. ``"Moon"``, ``"sun"``).
            body2: Name of the second body (e.g. ``"Jupiter"``, ``"mars"``).
            datetimes: List of timezone-aware datetimes.

        Returns:
            List of angular distances in degrees.

        Raises:
            KeyError: If body name is not found in ephemeris.
            ValueError: If datetimes list is empty.
        """
        if len(datetimes) == 0:
            raise ValueError("datetimes list cannot be empty")

        # Map body names to ephemeris keys
        key1 = PLANET_MAP.get(body1.lower(), body1.lower())
        key2 = PLANET_MAP.get(body2.lower(), body2.lower())
        try:
            b1 = self.eph[key1]
            b2 = self.eph[key2]
        except (KeyError, ValueError):
            raise KeyError(f"Unknown body {body1!r} or {body2!r}") from None

        # Compute separation for all times
        t = self.ts.from_datetimes(datetimes)
        pos1 = self.observer.at(t).observe(b1).apparent()
        pos2 = self.observer.at(t).observe(b2).apparent()
        sep = pos1.separation_from(pos2)
        return [float(d) for d in sep.degrees]
