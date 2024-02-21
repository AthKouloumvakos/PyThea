from astropy.time import TimeDelta
from sunpy.coordinates.ephemeris import get_horizons_coord

__all__ = ['prep_lasco']


def prep_lasco(map_sequence):
    for map_ in map_sequence:
        map_.meta['DATE-AVG'] = (map_.date + TimeDelta(map_.meta['exptime']/2, format='sec')).value

    dates = [tmap.date_average for tmap in map_sequence]

    observer = get_horizons_coord('SOHO', dates)
    for map_, soho_ in zip(map_sequence, observer):
        map_.meta['HGLN_OBS'] = soho_.lon.to('deg').value
        map_.meta['HGLT_OBS'] = soho_.lat.to('deg').value
        map_.meta['DSUN_OBS'] = soho_.radius.to('m').value

    return map_sequence
