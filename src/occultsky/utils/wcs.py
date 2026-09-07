"""
Set the WCS information manually by setting properties of the WCS object.
"""

from astropy.wcs import WCS


DEFAULT_PROJECTION = "STG"
DEFAULT_FOV = 30.0
DEFAULT_SIZE = 800


def build_wcs(
    projection: str = DEFAULT_PROJECTION,
    center: tuple[float, float] = (0.0, 0.0),
    fov: float = DEFAULT_FOV,
    size: int = DEFAULT_SIZE,
) -> WCS:
    """
    Build a WCS for a square, zenithal sky chart.

    Args:
        projection: WCS projection code in uppercase FITS form. Use a zenithal code such as ``STG``, ``TAN``, ``SIN`` or ``ZEA``. Do **not** pass a full-sky code (``AIT``, ``CAR``, ``CEA``, ``MOL``): a square grid cannot frame the whole sky and the result comes out distorted.
        center: ``(ra_degrees, dec_degrees)`` placed at the middle of the chart, ICRS. Defaults to ``(0, 0)``.
        fov: Width of the chart in degrees.
        size: Side length of the square chart in pixels.

    Returns:
        The WCS object.
    """
    ra, dec = center
    n = int(size)
    # One pixel covers this many degrees, so the chart spans exactly fov
    scale = fov / (n - 1)

    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [(n + 1) / 2.0, (n + 1) / 2.0]
    wcs.wcs.crval = [ra % 360.0, dec]
    wcs.wcs.ctype = [f"RA---{projection}", f"DEC--{projection}"]
    wcs.wcs.cunit = ["deg", "deg"]

    # Right ascension grows eastward, which is to the left on a sky chart
    wcs.wcs.cdelt = [-scale, scale]
    wcs.wcs.radesys = "ICRS"
    wcs.pixel_shape = (n, n)

    return wcs
