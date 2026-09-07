"""
OccultSky CLI
"""

import logging
import sys
from typing import Optional

import typer
from dotenv import load_dotenv

load_dotenv()


app = typer.Typer(
    no_args_is_help=True,
    help="OccultSky - Night sky observation calculator.",
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger(__name__)


@app.callback()
def main(
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging"),
) -> None:
    """OccultSky CLI."""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)


@app.command()
def timeline(
    date: str = typer.Argument(help="Local date in YYYY-MM-DD format"),
    lat: float = typer.Option(
        ..., help="Latitude in degrees (north positive)"
    ),
    lon: float = typer.Option(
        ..., help="Longitude in degrees (east positive)"
    ),
    elevation: float = typer.Option(
        ..., help="Elevation above sea level in meters"
    ),
    planets: Optional[str] = typer.Option(
        None, "--planets", "-p", help="Comma-separated planet names (e.g. mars,jupiter)"
    ),
    skyfield_data: Optional[str] = typer.Option(
        None,
        envvar="SKYFIELD_DATA",
        help="Path to Skyfield data directory (falls back to SKYFIELD_DATA env var)",
    ),
    ephemeris: str = typer.Option(
        "de441.bsp", help="JPL ephemeris filename"
    ),
    output: Optional[str] = typer.Option(
        None, "--output", "-o", help="Save to file instead of showing window"
    ),
) -> None:
    """Draw an astronomical timeline chart."""
    from .utils.plotting import draw_timeline

    planet_list = [p.strip() for p in planets.split(",") if p.strip()] if planets else []

    try:
        draw_timeline(
            lat=lat,
            lon=lon,
            elevation=elevation,
            date_str=date,
            planets=planet_list,
            output=output,
            skyfield_data=skyfield_data,
            ephemeris=ephemeris,
        )
    except (ValueError, KeyError) as e:
        logger.error("%s", e)
        raise typer.Exit(code=1)


# @app.command()
# def skyview(
#     date: str = typer.Argument(help="Local date in YYYY-MM-DD format"),
#     lat: float = typer.Option(..., help="Latitude in degrees (north positive)"),
#     lon: float = typer.Option(..., help="Longitude in degrees (east positive)"),
#     elevation: float = typer.Option(..., help="Elevation above sea level in meters"),
#     planets: Optional[str] = typer.Option(
#         None, "--planets", "-p", help="Comma-separated planet names (e.g. mars,jupiter)"
#     ),
#     interval: int = typer.Option(15, help="Sampling step in minutes"),
#     skyfield_data: Optional[str] = typer.Option(
#         None,
#         envvar="SKYFIELD_DATA",
#         help="Path to Skyfield data directory (falls back to SKYFIELD_DATA env var)",
#     ),
#     ephemeris: str = typer.Option(
#         "de441.bsp", help="JPL ephemeris filename"
#     ),
#     projection: Optional[str] = typer.Option(
#         None,
#         help="WCS projection code (default: STG, ignored when --fits-file is given)",
#     ),
#     ra: Optional[float] = typer.Option(
#         None, help="Chart center right ascension in degrees (default: fitted to the tracks)"
#     ),
#     dec: Optional[float] = typer.Option(
#         None, help="Chart center declination in degrees (default: fitted to the tracks)"
#     ),
#     fov: Optional[float] = typer.Option(
#         None, help="Chart width in degrees (default: fitted to the tracks)"
#     ),
#     size: Optional[int] = typer.Option(
#         None, help="Chart width in pixels (default: 800)"
#     ),
#     output: Optional[str] = typer.Option(
#         None, "--output", "-o", help="Save to file instead of showing window"
#     ),
# ):
#     """Draw a sky view chart of the Moon's nightly track plus any --planets."""
#     from .utils.plotting import draw_skyview

#     planet_list = [p.strip() for p in planets.split(",") if p.strip()] if planets else []

#     try:
#         observer = Observer(
#             lat=lat,
#             lon=lon,
#             elevation=elevation,
#             skyfield_data=skyfield_data,
#             ephemeris=ephemeris,
#         )
#         draw_skyview(
#             observer=observer,
#             date_str=date,
#             planets=planet_list,
#             interval=interval,
#             projection=projection,
#             ra=ra,
#             dec=dec,
#             fov=fov,
#             size=size,
#             output=output,
#         )
#     except (ValueError, KeyError) as e:
#         logger.error("%s", e)
#         raise typer.Exit(code=1)
