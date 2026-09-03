"""Astronomical calculation helpers."""

import numpy as np

from .constants import MOON_RADIUS_KM

# Radians to arcseconds conversion factor (180 * 3600 / pi)
_RAD_TO_ARCSEC = 180.0 * 3600.0 / np.pi


def moon_angular_radius(distance_km: float) -> float:
    """Return the Moon's apparent angular radius in arcseconds.

    The angular radius is derived from the Moon's physical radius and its distance via ``arcsin(radius / distance)``.

    Args:
        distance_km: Topocentric distance from the ground-based observer to the Moon's center in kilometers.

    Returns:
        Apparent angular radius in arcseconds.

    Raises:
        ValueError: If ``distance_km`` is not greater than the Moon's radius.
    """
    if distance_km <= MOON_RADIUS_KM:
        raise ValueError(f"distance_km must exceed {MOON_RADIUS_KM} km, got: {distance_km}")
    radius_rad = np.arcsin(MOON_RADIUS_KM / distance_km)
    return float(radius_rad * _RAD_TO_ARCSEC)
