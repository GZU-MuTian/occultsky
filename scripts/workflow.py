import bisect
import json
import logging
import sys
import time
from pathlib import Path

from dotenv import load_dotenv

from occultsky.main import Observer
from occultsky.utils.constants import CAPITALS
from occultsky.utils.plotting import draw_preview, draw_timeline

PROJECT_ROOT = Path(__file__).resolve().parents[1]

# Batch configuration: edit these defaults.
DEFAULT_INPUT = PROJECT_ROOT / "tmp" / "成功校验汇总.json"
OUTPUT_DIR = PROJECT_ROOT / "output"
# Output is always the output directory (JSONL goes to events.jsonl inside it).

# Load SKYFIELD_DATA (and other config) from the project .env regardless of CWD.
load_dotenv(PROJECT_ROOT / ".env")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger(__name__)

# Chinese body names used in the dataset mapped to Skyfield ephemeris keys.
BODY_NAME_MAP: dict[str, str] = {
    "月": "moon",
    "水星": "mercury",
    "金星": "venus",
    "火星": "mars",
    "木星": "jupiter",
    "土星": "saturn",
}

DYNASTY_CAPITALS: dict[str, tuple[tuple[str, ...], tuple[int, ...]]] = {
    # Relocating dynasties (capital-table notes a-d): seats, then each move year.
    "西晋":   (("luoyang", "xian"),    (313,)),   # Luoyang 265-313, Xi'an 313-316
    "北魏":   (("datong", "luoyang"),  (494,)),   # Datong 398-494, Luoyang 494-534
    "金朝":   (("beijing", "kaifeng"), (1214,)),  # Beijing 1153-1214, Kaifeng 1214-1234
    "明朝":   (("nanjing", "beijing"), (1420,)),  # Nanjing 1368-1420, Beijing 1420-1644
    # Single-seat dynasties, grouped by city.
    "唐朝":   (("xian",), ()),
    "曹魏":   (("luoyang",), ()),
    "后唐":   (("luoyang",), ()),
    "东晋":   (("nanjing",), ()),
    "南朝宋": (("nanjing",), ()),
    "南朝齐": (("nanjing",), ()),
    "后晋":   (("kaifeng",), ()),
    "北宋":   (("kaifeng",), ()),
    "南宋":   (("hangzhou",), ()),
    "元朝":   (("beijing",), ()),
    "清朝":   (("beijing",), ()),
    "东魏":   (("linzhang",), ()),
    "北齐":   (("linzhang",), ()),
}


def parse_western_date(raw_date: str) -> tuple[str, int]:
    """Convert a ``Y.M.D`` western-calendar string to ``(YYYY-MM-DD, year)``.

    Args:
        raw_date: Date string such as ``"234.11.22"``.

    Returns:
        Tuple of the zero-padded ``YYYY-MM-DD`` local date and the integer year.
    """
    year_str, month_str, day_str = raw_date.split(".")
    year = int(year_str)
    date_str = f"{year:04d}-{int(month_str):02d}-{int(day_str):02d}"
    return date_str, year


def resolve_capital(dynasty: str, year: int) -> str:
    """Return the CAPITALS key that was the dynasty's capital in ``year``.

    Args:
        dynasty: Dynasty name as it appears in the dataset (Chinese).
        year: Western-calendar year of the record.

    Returns:
        Capital city key present in ``CAPITALS``.

    Raises:
        KeyError: If the dynasty is unknown.
    """
    seats, moves = DYNASTY_CAPITALS[dynasty]
    # Moves already made by ``year`` = index of the seat then in effect; a year
    # before the first move (or a dynasty that never moved) stays on ``seats[0]``.
    return seats[bisect.bisect_right(moves, year)]


def get_observer(capital_key: str) -> Observer:
    """Build a fresh Observer for ``capital_key``.

    Each call constructs a new Observer (no caching), so state can never leak
    between records that happen to share a capital.

    Args:
        capital_key: Key into ``CAPITALS`` (stored as ``(lon, lat, elevation)``).

    Returns:
        The Observer located at that capital.
    """
    longitude, latitude, elevation = CAPITALS[capital_key]
    # CAPITALS stores (lon, lat, elev); Observer expects (lat, lon, elevation).
    return Observer(lat=latitude, lon=longitude, elevation=elevation)


def process_record(record: dict) -> dict:
    """Reproduce one dataset record and return a JSON-serializable result.

    Args:
        record: A single entry from the dataset.

    Returns:
        Flat dict with the record identifiers.

    Raises:
        KeyError: If a body name, dynasty, or ephemeris target is unknown.
        ValueError: If the date string is malformed.
    """
    date_str, year = parse_western_date(record["western_calendar_validated"])
    dynasty = record["dynasties"]
    capital_key = resolve_capital(dynasty, year)

    # Extract the planet name (the non-Moon body)
    planet = record["object"]["name"] if record["subject"]["name"] == "月" else record["subject"]["name"]
    planet_en = BODY_NAME_MAP[planet]

    # Generate output path for the preview chart
    output_path = OUTPUT_DIR / f"{record['id']}_{date_str}.png"

    result: dict = {
        "id": record["id"],
        "date": date_str,
        "dynasty": dynasty,
        "capital": capital_key,
        "planet": planet,
        "planet_en": planet_en,
        "output": str(output_path),
    }

    return result


def main(
    input_path: Path = DEFAULT_INPUT,
    output_dir: Path = OUTPUT_DIR,
) -> None:
    """Reproduce every record in the dataset and emit JSONL results.

    Args:
        input_path: JSON dataset to read (one record per element).
        output_dir: Directory to write charts and JSONL results to.
    """
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")

    with input_path.open(encoding="utf-8") as input_file:
        records = json.load(input_file)
    logger.info(
        "batch start input=%s records=%d output_dir=%s",
        input_path, len(records), output_dir,
    )

    # Ensure the output directory exists.
    output_dir.mkdir(parents=True, exist_ok=True)
    events_path = output_dir / "events.jsonl"

    start = time.monotonic()
    succeeded = 0
    with events_path.open("w", encoding="utf-8") as out:
        for record in records:
            try:
                result = process_record(record)
            except (KeyError, ValueError) as exc:
                logger.error("record failed id=%s error=%s", record.get("id"), exc)
                out.write(json.dumps({"id": record.get("id"), "error": str(exc)},
                                     ensure_ascii=False) + "\n")
                continue

            # Generate preview and timeline charts
            chart_path = Path(result["output"])
            chart_path.parent.mkdir(parents=True, exist_ok=True)
            timeline_path = chart_path.with_stem(chart_path.stem + "_timeline")
            try:
                lon, lat, elev = CAPITALS[result["capital"]]
                draw_preview(
                    lat=lat, lon=lon, elevation=elev,
                    date_str=result["date"],
                    planets=[result["planet_en"]],
                    output=str(chart_path),
                    stars=True, constellation=True,
                    mansions=True, asterism=True,
                    fov=10
                )
                draw_timeline(
                    lat=lat, lon=lon, elevation=elev,
                    date_str=result["date"],
                    planets=[result["planet_en"]],
                    output=str(timeline_path),
                )
                logger.info("charts saved %s %s", chart_path, timeline_path)
            except Exception as exc:
                logger.error("chart failed id=%s error=%s", result["id"], exc)

            succeeded += 1
            logger.debug("record done id=%s", result["id"])
            out.write(json.dumps(result, ensure_ascii=False) + "\n")

    elapsed = time.monotonic() - start
    logger.info(
        "batch done succeeded=%d total=%d elapsed=%.2fs",
        succeeded, len(records), elapsed,
    )


if __name__ == "__main__":
    main()
