# OccultSky

OccultSky is a Python toolkit for searching historical lunar occultation events of planets, powered by [Skyfield](https://rhodesmill.org/skyfield/). Given a date and an observing site, it computes the rise/set times, positions, and mutual angular distances of the Moon and planets, and renders night-sky charts.

## Features

- **High-precision ephemeris** — JPL DE441 SPICE kernel covering 13,200 BCE – 17,191 CE, suitable for historical dates.
- **Historical-friendly time model** — dates before 1582-10-15 follow the Julian calendar; local times use Local Mean Time (LMT) derived directly from the site's longitude.
- **Two chart types**:
  - `timeline` — a single-axis chart of the night with Moon/planet rise and set milestones and shaded night/twilight bands;
  - `preview` — a two-panel sky chart in a WCS projection showing the Moon's nightly track (as a buffered band), planet positions, and sky-culture layers (constellation lines, star names, asterisms, the 28 lunar mansions).
- **CLI + Python API** — a Typer command-line interface and an `Observer` class for programmatic use.

## Requirements

- Python >= 3.13 (pinned by `.python-version`)
- [uv](https://docs.astral.sh/uv/) (recommended) or pip

## Installation

```powershell
uv sync                 # create .venv and install dependencies from uv.lock
uv run occultsky --help # verify the CLI entry point
```

## Configuration

OccultSky needs a local data directory for Skyfield caches. Set it in one of these ways (highest priority first):

1. CLI option: `--skyfield-data <dir>`
2. Environment variable: `SKYFIELD_DATA=<dir>`
3. `.env` file in the working directory (loaded via python-dotenv):

```ini
SKYFIELD_DATA=C:\path\to\occultsky-data
```

```powershell
# PowerShell example
$env:SKYFIELD_DATA = "C:\path\to\occultsky-data"
```

On first run, Skyfield downloads the required data into this directory automatically:

- `de441.bsp` — JPL ephemeris (also downloadable manually from [NASA JPL NAIF](https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/));
- `hip_main.dat` — Hipparcos star catalogue;
- `skyculture_<culture>.json` — Stellarium sky-culture data.

If `SKYFIELD_DATA` is not configured, the toolkit raises a `ValueError` with instructions.

## Quick Start

Dates are **local** dates in `YYYY-MM-DD`; zero-pad years shorter than four digits (e.g. `0303-12-05`).

### Timeline chart

```powershell
uv run occultsky timeline 0303-12-05 --lat 39.93 --lon 116.37 --elevation 52 --planets venus -o timeline.png
```

### Sky preview chart

```powershell
uv run occultsky preview 0303-12-05 --lat 39.93 --lon 116.37 --elevation 52 --planets venus --stars --constellation --mansions -o preview.png
```

Without `-o`, charts are shown in an interactive window instead of being saved.

## CLI Options

Common options (both commands):

| Option | Description |
| --- | --- |
| `DATE` (argument) | Local date in `YYYY-MM-DD` |
| `--lat`, `--lon`, `--elevation` | Observing site: latitude/longitude in degrees (north/east positive), elevation in meters (all required) |
| `--planets`, `-p` | Comma-separated planet names, e.g. `mars,jupiter` (lowercase English) |
| `--skyfield-data` | Skyfield data directory (falls back to the `SKYFIELD_DATA` env var) |
| `--ephemeris` | JPL ephemeris filename (default: `de441.bsp`) |
| `--output`, `-o` | Save the chart to a file instead of showing a window |
| `--verbose`, `-v` | Enable verbose logging |

`preview` additionally accepts:

| Option | Description |
| --- | --- |
| `--interval` | Sampling step in minutes (default: 15) |
| `--projection` | WCS projection code (default: `STG`) |
| `--fov` | Field of view in degrees (default: 6) |
| `--stars`, `--constellation`, `--asterism`, `--mansions` | Toggle sky-culture layers |
| `--skyculture` | Stellarium sky culture (default: `chinese`) |

Run `uv run occultsky <command> --help` for the full option list.

## External Data

| Data | Source | Purpose |
| --- | --- | --- |
| `de441.bsp` | [NASA JPL NAIF SPICE](https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/) | Lunar and planetary ephemeris, 13,200 BCE – 17,191 CE |
| `hip_main.dat` | Hipparcos catalogue (via Skyfield) | Star positions and magnitudes |
| `skyculture_<culture>.json` | [Stellarium](https://github.com/Stellarium/stellarium) sky cultures | Constellation lines, star names, lunar mansions |

## Star Name Annotations

Star names are passed through verbatim from Stellarium's sky-culture data. In the `chinese` culture, suffixes mark stars added from Pan Nai's *The History of Stellar Observation in China* (中国恒星观测史): `*` = certified, `?` = uncertain identification (e.g. `天籥二?`). See Stellarium's [chinese sky-culture description](https://github.com/Stellarium/stellarium/blob/master/skycultures/chinese/description.md).

These suffixes are kept in charts on purpose: for historical-record verification, an uncertain identification is itself evidence and must not be silently dropped.
