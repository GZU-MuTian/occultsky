from typing import Dict, Tuple

# Moon mean radius in kilometers (IAU standard)
MOON_RADIUS_KM: float = 1737.4

# Ancient Chinese capitals (lon, lat) in degrees
CAPITALS: Dict[str, Tuple[float, float]] = {
    'luoyang':   (112.4, 34.7),
    'nanjing':   (118.8, 32.0),
    'linzhang':  (114.6, 36.3),
    'xian':      (108.9, 34.2),
    'kaifeng':   (114.6, 34.5),
    'hangzhou':  (120.2, 30.2),
    'beijing':   (116.4, 39.6),
}

# Outer planets only support barycenter names in DE441
PLANET_MAP: Dict[str, str] = {
    'mercury': 'mercury',
    'venus': 'venus',
    'mars': 'mars barycenter',
    'jupiter': 'jupiter barycenter',
    'saturn': 'saturn barycenter',
    'uranus': 'uranus barycenter',
    'neptune': 'neptune barycenter',
    'pluto': 'pluto barycenter',
}