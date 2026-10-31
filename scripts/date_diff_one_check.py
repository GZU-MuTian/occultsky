"""Verify erratum records whose validation mentions "日期差一天" (off by one day).

For every record of the erratum YAML report whose ``validation`` text contains
"日期差一天", this script recomputes, for each of its dates, four consecutive
night windows (local dusk to dawn, 1-minute grid):

- ``prev2_night``: dusk of date-2 to dawn of date-1 (the "previous day morning");
- ``prev_night`` : dusk of date-1 to dawn of date  (previous night);
- ``curr_night`` : dusk of date   to dawn of date+1 (the night of the date);
- ``next_night`` : dusk of date+1 to dawn of date+2 (the "next day evening").

For each window it records the Moon's visible time range and, per classical
planet, the closest Moon-planet angular separation among samples where the
planet is also above the geometric horizon, with its LMT.

Results are written as UTF-8 JSONL: the full filtered records and one computed
record per line.
"""

import argparse
import json
import logging
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

import yaml
from dotenv import load_dotenv

from occultsky.main import Observer
from occultsky.utils.constants import CAPITALS

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "inbox" / "月掩行星勘误筛选报告.yaml"
DEFAULT_RECORDS_OUTPUT = PROJECT_ROOT / "output" / "date_diff_one_records.jsonl"
DEFAULT_WINDOWS_OUTPUT = PROJECT_ROOT / "output" / "date_diff_one_windows.jsonl"

# Load SKYFIELD_DATA (and other config) from the project .env regardless of CWD.
load_dotenv(PROJECT_ROOT / ".env")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger(__name__)

# The five classical planets whose closest approach to the Moon is recorded.
PLANETS: tuple[str, ...] = ("mercury", "venus", "mars", "jupiter", "saturn")

# The planet each record's erratum actually refers to, keyed by record line.
FOCUS_PLANET: dict[int, str] = {
    5979: "venus",    # text says 荧惑 (mars) but the erratum identifies Venus
    6749: "venus",    # 月犯太白
    6754: "saturn",   # 犯镇星
    6779: "venus",    # 月犯太白
    6859: "jupiter",  # 月复掩岁星
    6869: "venus",    # 月掩太白
    7474: "mars",     # 荧惑掩月
    7711: "venus",    # 犯太白
    7716: "jupiter",  # 又犯岁星
    7721: "saturn",   # 又犯填星
    7726: "saturn",   # 犯填星
    10945: "saturn",  # 太阴犯填星
    11120: "mercury", # 犯辰星
    11132: "mars",    # 犯荧惑
    11137: "venus",   # 犯太白
    11142: "venus",   # 犯太白
    11147: "saturn",  # 犯填星
}

# Marker string identifying the erratum type to filter on.
TARGET_MARKER: str = "日期差一天"

# Observer cache keyed by capital so the large ephemeris is parsed once per city.
_observer_cache: dict[str, Observer] = {}


def parse_western_date(raw_date: str) -> str:
    """Convert a ``Y.M.D`` western-calendar string to zero-padded ``YYYY-MM-DD``.

    Args:
        raw_date: Date string such as ``"303.12.5"``.

    Returns:
        Zero-padded ``YYYY-MM-DD`` local date string.

    Raises:
        ValueError: If ``raw_date`` is not in ``Y.M.D`` form.
    """
    parts = raw_date.split(".")
    if len(parts) != 3:
        raise ValueError(f"date must be 'Y.M.D', got: {raw_date!r}")
    year_str, month_str, day_str = parts
    return f"{int(year_str):04d}-{int(month_str):02d}-{int(day_str):02d}"


def shift_date(date_str: str, days: int) -> str:
    """Return ``date_str`` shifted by ``days`` calendar days.

    Args:
        date_str: Date in ``YYYY-MM-DD`` format.
        days: Day offset (may be negative).

    Returns:
        Shifted date in ``YYYY-MM-DD`` format.
    """
    shifted = datetime.strptime(date_str, "%Y-%m-%d") + timedelta(days=days)
    return shifted.strftime("%Y-%m-%d")


def get_observer(capital_key: str) -> Observer:
    """Return a cached Observer for ``capital_key``.

    Args:
        capital_key: Key into ``CAPITALS`` (stored as ``(lon, lat, elevation)``).

    Returns:
        The Observer located at that capital (one instance per capital).

    Raises:
        KeyError: If the capital key is unknown.
    """
    if capital_key not in _observer_cache:
        longitude, latitude, elevation = CAPITALS[capital_key]
        # CAPITALS stores (lon, lat, elev); Observer expects (lat, lon, elevation).
        _observer_cache[capital_key] = Observer(
            lat=latitude, lon=longitude, elevation=elevation
        )
    return _observer_cache[capital_key]


def compute_window(observer: Observer, moon_samples: list[datetime]) -> dict:
    """Compute the Moon's visible range and closest planet approaches in a window.

    Args:
        observer: Observer located at the record's capital.
        moon_samples: 1-minute LMT samples where the Moon is above the horizon
            (from ``get_visible_times``); empty when the Moon never rose.

    Returns:
        JSON-serializable dict with ``moon_visible`` (first/last sample) and one
        ``closest`` entry per planet; fields stay ``None`` when nothing is
        visible.
    """
    window: dict = {
        "moon_visible": None,
        "closest": {planet: None for planet in PLANETS},
    }
    if not moon_samples:
        # The Moon never rose in this window: every field stays empty.
        return window

    window["moon_visible"] = {
        "start": moon_samples[0].isoformat(),
        "end": moon_samples[-1].isoformat(),
    }

    for planet in PLANETS:
        # Keep only samples where the planet is above the geometric horizon.
        altaz = observer.get_altaz(moon_samples, target=planet)
        up = [dt for dt, (alt, _az) in zip(moon_samples, altaz) if alt > 0.0]
        if not up:
            # Planet below the horizon all night: keep the key with an empty value.
            continue
        distances = observer.get_angular_distance("moon", planet, up)
        index = min(range(len(distances)), key=distances.__getitem__)
        window["closest"][planet] = {
            "min_dist_deg": distances[index],
            "lmt": up[index].isoformat(),
        }
    return window


def process_record(record: dict) -> dict:
    """Recompute one date off-by-one record across all of its dates.

    Args:
        record: A single entry from the erratum YAML report (must carry
            ``line``, ``capital`` and ``date``).

    Returns:
        JSON-serializable dict keyed by ``line`` with one event per date, each
        holding the four consecutive night windows.

    Raises:
        KeyError: If a capital or ephemeris target is unknown.
        ValueError: If a date string is malformed.
    """
    line = record["line"]
    capital = record["capital"]
    observer = get_observer(capital)

    raw_dates = record["date"] if isinstance(record["date"], list) else [record["date"]]
    events = []
    for raw_date in raw_dates:
        date_str = parse_western_date(raw_date)
        # Main windows: the night of the date and the previous night.
        prev_visible, curr_visible = observer.get_visible_times(
            date_str, target="moon", interval=1
        )
        # Expansion windows to locate the true close approach: the previous
        # day's own night (early morning of date-1) and the next day's night
        # (evening of date+1), so all four windows tile the timeline.
        prev2_visible, _ = observer.get_visible_times(
            shift_date(date_str, -1), target="moon", interval=1
        )
        _, next_visible = observer.get_visible_times(
            shift_date(date_str, 1), target="moon", interval=1
        )
        events.append({
            "date": raw_date,
            "windows": {
                "prev2_night": compute_window(observer, prev2_visible),
                "prev_night": compute_window(observer, prev_visible),
                "curr_night": compute_window(observer, curr_visible),
                "next_night": compute_window(observer, next_visible),
            },
        })
        logger.debug("event done line=%d date=%s", line, raw_date)

    return {
        "line": line,
        "capital": capital,
        "focus_planet": FOCUS_PLANET.get(line),
        "status": record.get("status"),
        "events": events,
    }


def main(input_path: Path, records_output: Path, windows_output: Path) -> None:
    """Filter the erratum report, recompute the records and write JSONL results.

    Args:
        input_path: YAML erratum report to read.
        records_output: JSONL file for the full filtered records (UTF-8).
        windows_output: JSONL file for the computed window results (UTF-8).
    """
    start = time.monotonic()
    with input_path.open(encoding="utf-8") as input_file:
        records = yaml.safe_load(input_file)

    # Keep records whose validation text mentions the off-by-one-day erratum.
    selected = [
        r for r in records if TARGET_MARKER in (r.get("validation") or "")
    ]
    logger.info(
        "filter done input=%s total=%d selected=%d",
        input_path, len(records), len(selected),
    )

    records_output.parent.mkdir(parents=True, exist_ok=True)
    with records_output.open("w", encoding="utf-8") as out:
        for record in selected:
            out.write(json.dumps(record, ensure_ascii=False, default=str) + "\n")

    succeeded = 0
    with windows_output.open("w", encoding="utf-8") as out:
        for record in selected:
            line = record.get("line")
            try:
                result = process_record(record)
            except (KeyError, ValueError) as exc:
                logger.error("record failed line=%s error=%s", line, exc)
                out.write(
                    json.dumps({"line": line, "error": str(exc)}, ensure_ascii=False)
                    + "\n"
                )
                continue
            succeeded += 1
            out.write(json.dumps(result, ensure_ascii=False) + "\n")

    elapsed = time.monotonic() - start
    logger.info(
        "batch done succeeded=%d/%d records=%s windows=%s elapsed=%.2fs",
        succeeded, len(selected), records_output, windows_output, elapsed,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Verify "日期差一天" erratum records with night windows.'
    )
    parser.add_argument(
        "--input", type=Path, default=DEFAULT_INPUT,
        help="Path to the erratum YAML report",
    )
    parser.add_argument(
        "--records-output", type=Path, default=DEFAULT_RECORDS_OUTPUT,
        help="Path of the UTF-8 JSONL output for the filtered records",
    )
    parser.add_argument(
        "--windows-output", type=Path, default=DEFAULT_WINDOWS_OUTPUT,
        help="Path of the UTF-8 JSONL output for the computed windows",
    )
    args = parser.parse_args()
    main(args.input, args.records_output, args.windows_output)
