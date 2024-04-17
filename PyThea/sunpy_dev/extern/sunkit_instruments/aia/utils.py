from astropy.time import TimeDelta
from aiapy.calibrate import fix_observer_location, update_pointing

__all__ = ['prep_aia']


def prep_aia(map_sequence):
    map_sequence = [update_pointing(tmap) for tmap in map_sequence]
    map_sequence = [fix_observer_location(tmap) for tmap in map_sequence]

    for map_ in map_sequence:
        map_.meta['DATE-AVG'] = (map_.date + TimeDelta(map_.meta['exptime']/2, format='sec')).value

    return map_sequence
