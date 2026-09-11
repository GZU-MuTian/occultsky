import json
import logging
import os
from dataclasses import dataclass

import numpy as np

from skyfield.data import hipparcos
from skyfield.api import Loader

logger = logging.getLogger(__name__)


# URLs and filenames for Stellarium skyculture data.
_SKYCULTURE_URL = (
    "https://raw.githubusercontent.com/Stellarium/stellarium/master/skycultures/{culture}/index.json"
)
_SKYCULTURE_FILE = "skyculture_{culture}.json"


@dataclass(frozen=True)
class Constellation:
    """A constellation figure with its lines as tuples of HIP numbers."""
    name: str
    lines: tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class NamedStars:
    """A collection of HIP numbers paired with display names."""
    hip: np.ndarray   # HIP number
    name: list[str]   # display name per entry

    def __len__(self):
        return len(self.hip)


class Skyculture:
    """
    Loads the Stellarium JSON and Hipparcos catalog for drawing stars and constellation lines.
    """

    def __init__(self, culture: str, skyfield_data: str = None):
        # Set Skyfield data directory.
        skyfield_data = skyfield_data or os.environ.get("SKYFIELD_DATA")
        if not skyfield_data:
            raise ValueError("Skyfield data directory not specified. Pass it as a parameter or set SKYFIELD_DATA environment variable")
        load = Loader(skyfield_data)

        # Load Hipparcos catalog using skyfield Loader
        with load.open(hipparcos.URL, filename="hip_main.dat") as handle:
            self.hips = hipparcos.load_dataframe(handle)
        
        # Load Stellarium JSON file for the specified culture.
        url = _SKYCULTURE_URL.format(culture=culture)
        filename = _SKYCULTURE_FILE.format(culture=culture)
        with load.open(url, filename=filename) as f:
            self.skyculture = json.load(f)

    def within_validity(self, year: int) -> bool:
        """Check if a year falls within the skyculture's validity range.
        
        Args:
            year: Calendar year to check.
            
        Returns:
            True if year is within beginTime..endTime, or if range is not specified.
        """
        begin = self.skyculture.get('beginTime')
        end = self.skyculture.get('endTime')
        
        if begin is not None and year < begin:
            return False
        if end is not None and year > end:
            return False
        return True

    def constellations(self, name_field: str = 'native') -> list[Constellation]:
        """Return the constellation figures from the skyculture data.

        Args:
            name_field: Key in each ``common_name`` dict to use as the display name (e.g. ``'native'``, ``'english'``, ``'pronounce'``).
        """
        return [
            Constellation(
                name=(e.get('common_name') or {}).get(name_field, ''),
                lines=tuple(tuple(r) for r in e.get('lines', [])),
            )
            for e in self.skyculture.get('constellations', [])
        ]

    def asterisms(self, name_field: str = 'native') -> list[Constellation]:
        """Return the asterism figures from the skyculture data.

        Args:
            name_field: Key in each ``common_name`` dict to use as the display name (e.g. ``'native'``, ``'english'``, ``'pronounce'``).
        """
        return [
            Constellation(
                name=(e.get('common_name') or {}).get(name_field, ''),
                lines=tuple(tuple(r) for r in e.get('lines', [])),
            )
            for e in self.skyculture.get('asterisms', [])
        ]

    def mansions(self, name_field: str = 'native') -> NamedStars:
        """Return the lunar mansions (e.g. Chinese Xiu) from the skyculture data.

        Args:
            name_field: Key in each name dict to use as the display name (e.g. ``'native'``, ``'english'``, ``'pronounce'``).

        Returns:
            A NamedStars with arrays of HIP numbers and names. Empty arrays if the skyculture has no lunar_system section.

        Raises:
            ValueError: If the lunar_system has no ``defining_stars``, or if
                ``defining_stars`` and ``names`` have different lengths.
        """
        ls = self.skyculture.get('lunar_system')
        if not ls:
            return NamedStars(hip=np.array([], dtype=int), name=[])

        defining_stars = ls.get('defining_stars')
        if not defining_stars:
            raise ValueError(f"lunar_system has no defining_stars for '{self.skyculture.get('id', '?')}'")

        names_list = ls.get('names', [])
        if len(defining_stars) != len(names_list):
            raise ValueError(
                f"lunar_system defining_stars ({len(defining_stars)}) and names ({len(names_list)}) length mismatch"
            )

        return NamedStars(
            hip=np.array(defining_stars),
            name=[entry.get(name_field, '') for entry in names_list],
        )

    def star_names(self, name_field: str = 'native') -> NamedStars:
        """Return the common star names from the skyculture data.

        Args:
            name_field: Key in each name entry to use as the display name (e.g. ``'native'``, ``'english'``, ``'pronounce'``).

        Returns:
            A NamedStars with arrays of HIP numbers and names.
        """
        raw = self.skyculture.get('common_names', {})
        hips: list[int] = []
        names: list[str] = []
        for key, entries in raw.items():
            if not key.startswith('HIP '):
                continue
            hip = int(key[4:])
            if entries:
                hips.append(hip)
                names.append(entries[0].get(name_field, ''))
        return NamedStars(hip=np.array(hips), name=names)
