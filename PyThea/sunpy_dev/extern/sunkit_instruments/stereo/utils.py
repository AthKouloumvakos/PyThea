import numpy as np
import sunpy.map

from sunpy_dev.map import maputils

__all__ = ['cor_polariz']


def cor_polariz(map_sequence):
    """
    Returns the total brightness image from three polarized images.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A SunPy map.

    Returns
    -------
    `~sunpy.map.GenericMap`
        A SunPy map.

    Note
    -------
    This script is partially based on cor_polariz.pro
    """
    sequence_polar_0 = [tmap for tmap in map_sequence if tmap.meta['polar'] == 0]
    sequence_polar_120 = [tmap for tmap in map_sequence if tmap.meta['polar'] == 120]
    sequence_polar_240 = [tmap for tmap in map_sequence if tmap.meta['polar'] == 240]

    # Calculate Mueller matrix
    s, d = 0.5, 0.5
    ang_0, ang_120, ang_240 = 0, 120, 240

    m = [[s, d*np.cos(2*np.deg2rad(ang_0)), d*np.sin(2*np.deg2rad(ang_0))],
         [s, d*np.cos(2*np.deg2rad(ang_120)), d*np.sin(2*np.deg2rad(ang_120))],
         [s, d*np.cos(2*np.deg2rad(ang_240)), d*np.sin(2*np.deg2rad(ang_240))]]
    m = np.linalg.inv(m)

    # Calculate the Stokes parameters for incident intensity
    I = m[0][0] + m[0][1] + m[0][2]
    # Q = m[1][0] + m[1][1] + m[1][2]
    # U = m[2][0] + m[2][1] + m[2][2]

    smap = []
    for map_0 in sequence_polar_0:
        map_120 = maputils.get_closest(sequence_polar_120, map_0.date)
        map_240 = maputils.get_closest(sequence_polar_240, map_0.date)

        map_0 = maputils.normalize_exposure(map_0)
        map_120 = maputils.normalize_exposure(map_120)
        map_240 = maputils.normalize_exposure(map_240)

        # Calculate the Stokes parameters for incident intensity
        I = map_0.data * m[0][0] + map_120.data * m[0][1] + map_240.data * m[0][2]
        # Q = map_0.data * m[1][0] + map_120.data * m[1][1] + map_240.data * m[1][2]
        # U = map_0.data * m[2][0] + map_120.data * m[2][1] + map_240.data * m[2][2]

        t = [map_0.date, map_120.date, map_240.date]
        obs_time = min(t) + (max(t)-min(t))/2
        crota = np.mean([map_0.meta['crota'], map_120.meta['crota'], map_240.meta['crota']])

        new_map_meta = map_120.meta
        # TODO: Check if any other keys need to change
        new_map_meta['date-obs'] = obs_time.strftime('%Y-%m-%dT%H:%M:%S.%f')
        new_map_meta['polar'] = '1001'
        new_map_meta['crval1'] = np.mean([map_0.meta['crval1'], map_120.meta['crval1'], map_240.meta['crval1']])
        new_map_meta['crval2'] = np.mean([map_0.meta['crval2'], map_120.meta['crval2'], map_240.meta['crval2']])
        new_map_meta['crpix1'] = np.mean([map_0.meta['crpix1'], map_120.meta['crpix1'], map_240.meta['crpix1']])
        new_map_meta['crpix2'] = np.mean([map_0.meta['crpix2'], map_120.meta['crpix2'], map_240.meta['crpix2']])
        new_map_meta['crota'] = crota
        new_map_meta['PC1_1'] = np.cos(np.deg2rad(crota))
        new_map_meta['PC1_2'] = -np.sin(np.deg2rad(crota))
        new_map_meta['PC2_1'] = np.sin(np.deg2rad(crota))
        new_map_meta['PC2_2'] = np.cos(np.deg2rad(crota))

        map_ = sunpy.map.Map(I, new_map_meta)
        smap.append(map_)

    return sunpy.map.Map(smap, sequence=True)
