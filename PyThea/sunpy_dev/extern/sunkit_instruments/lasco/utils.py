from sunpy.coordinates.ephemeris import get_horizons_coord

__all__ = ['prep_lasco']


def prep_lasco(map_sequence):
    dates = [tmap.date for tmap in map_sequence]
    observer = get_horizons_coord('SOHO', dates)
    for map_, soho_  in zip(map_sequence, observer):
        map_.meta['HGLN_OBS'] = soho_.lon.to('deg').value
        map_.meta['HGLT_OBS'] = soho_.lat.to('deg').value
        map_.meta['DSUN_OBS'] = soho_.radius.to('m').value

    return map_sequence
