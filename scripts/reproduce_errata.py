"""Reproduce lunar-occultation errata records and emit JSONL verification data.

For every in-scope record of the erratum YAML report this script computes, per
date and per night (previous night and the night of the date):
- the Moon's (ra, dec) at the centre of its visible window;
- the closest Moon-planet angular separation and its LMT for each of the five
  classical planets, searching only where the Moon is up and the planet is
  above the geometric horizon.

Results are written as JSONL, one record per line.
"""

import argparse
import json
import logging
import sys
import time
from datetime import datetime
from pathlib import Path

import yaml
from dotenv import load_dotenv

from occultsky.main import Observer
from occultsky.utils.constants import CAPITALS

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "inbox" / "月掩行星勘误筛选报告.yaml"
DEFAULT_OUTPUT = PROJECT_ROOT / "output" / "erratum_repro.jsonl"

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

# Records with this status are out of scope; only 未核验/已核验 are processed.
EXCLUDED_STATUS: str = "已排除"

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


def compute_night(
    observer: Observer,
    moon_samples: list[datetime],
) -> dict:
    """Compute the Moon position at its visible-window centre and closest planet approaches.

    Args:
        observer: Observer located at the record's capital.
        moon_samples: 1-minute LMT samples where the Moon is above the horizon
            (from ``get_visible_times``); empty when the Moon never rose.

    Returns:
        JSON-serializable dict with ``moon_radec_mid`` and one ``closest`` entry
        per planet; fields stay ``None`` when nothing is visible.
    """
    night: dict = {
        "moon_radec_mid": None,
        "closest": {planet: None for planet in PLANETS},
    }
    if not moon_samples:
        # The Moon never rose that night: every field stays empty.
        return night

    # Centre of the Moon's visible window: the middle sample of the
    # get_visible_times result (1-minute grid), always a moment the Moon is up.
    mid = moon_samples[len(moon_samples) // 2]
    ra, dec = observer.get_radec([mid], target="moon")[0]
    night["moon_radec_mid"] = {"ra": ra, "dec": dec}

    for planet in PLANETS:
        # Keep only samples where the planet is above the geometric horizon.
        altaz = observer.get_altaz(moon_samples, target=planet)
        up = [dt for dt, (alt, _az) in zip(moon_samples, altaz) if alt > 0.0]
        if not up:
            # Planet below the horizon all night: keep the key with an empty value.
            continue
        distances = observer.get_angular_distance("moon", planet, up)
        index = min(range(len(distances)), key=distances.__getitem__)
        night["closest"][planet] = {
            "min_dist_deg": distances[index],
            "lmt": up[index].isoformat(),
        }
    return night


def process_record(record: dict) -> dict:
    """Reproduce one erratum record across all of its dates.

    Args:
        record: A single entry from the erratum YAML report (must carry
            ``line``, ``capital`` and ``date``).

    Returns:
        JSON-serializable dict keyed by ``line`` with one event per date.

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
        # 1-minute samples: the Moon-visible window is also the search grid.
        previous_visible, current_visible = observer.get_visible_times(
            date_str, target="moon", interval=1
        )
        events.append({
            "date": raw_date,
            "prev_night": compute_night(observer, previous_visible),
            "curr_night": compute_night(observer, current_visible),
        })
        logger.debug("event done line=%d date=%s", line, raw_date)

    return {"line": line, "capital": capital, "events": events}


def main(input_path: Path, output_path: Path) -> None:
    """Read the erratum report, compute reproductions and write JSONL results.

    Args:
        input_path: YAML erratum report to read.
        output_path: JSONL file to write, one record per line (UTF-8, written
            directly to avoid shell redirection encoding issues).
    """
    start = time.monotonic()
    with input_path.open(encoding="utf-8") as input_file:
        records = yaml.safe_load(input_file)

    # Only unreviewed (未核验) and verified (已核验) records are in scope.
    in_scope = [r for r in records if r.get("status") != EXCLUDED_STATUS]
    logger.info(
        "batch start input=%s total=%d in_scope=%d",
        input_path, len(records), len(in_scope),
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    succeeded = 0
    events_total = 0
    with output_path.open("w", encoding="utf-8") as out:
        for record in in_scope:
            line = record.get("line")
            if not record.get("capital"):
                logger.warning("record skipped line=%s reason=missing capital", line)
                continue
            if not record.get("date"):
                logger.warning("record skipped line=%s reason=missing date", line)
                continue
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
            events_total += len(result["events"])
            out.write(json.dumps(result, ensure_ascii=False) + "\n")

    elapsed = time.monotonic() - start
    logger.info(
        "batch done succeeded=%d events=%d output=%s elapsed=%.2fs",
        succeeded, events_total, output_path, elapsed,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Reproduce lunar-occultation errata records into JSONL."
    )
    parser.add_argument(
        "--input", type=Path, default=DEFAULT_INPUT,
        help="Path to the erratum YAML report",
    )
    parser.add_argument(
        "--output", type=Path, default=DEFAULT_OUTPUT,
        help="Path of the UTF-8 JSONL output file",
    )
    args = parser.parse_args()
    main(args.input, args.output)
