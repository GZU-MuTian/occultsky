from typing import Dict, Tuple

# Moon mean radius in kilometers (IAU standard)
MOON_RADIUS_KM: float = 1737.4

# Ancient Chinese capitals: (longitude, latitude) in degrees, elevation in meters.
# Coordinates from CHGIS/TGAZ seats of government; elevations sampled on Copernicus DEM GLO-90.
CAPITALS: Dict[str, Tuple[float, float, float]] = {
    'xian':      (108.91, 34.25, 411),   # Western Jin, Tang
    'luoyang':   (112.60, 34.73, 125),   # Cao Wei, Western Jin, Northern Wei (late), Later Tang
    'nanjing':   (118.77, 32.05, 19),    # Eastern Jin, Liu Song, Southern Qi, Ming (early)
    'datong':    (113.30, 40.09, 1053),  # Northern Wei (early, Pingcheng)
    'kaifeng':   (114.34, 34.79, 75),    # Later Jin, Northern Song, Jin (Nanjing, late)
    'hangzhou':  (120.17, 30.29, 10),    # Southern Song (Lin'an)
    'beijing':   (116.37, 39.93, 52),    # Jin (Zhongdu), Yuan, Ming (from 1420), Qing
    'linzhang':  (114.41, 36.27, 78),    # Northern Qi, Eastern Wei (Ye)
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