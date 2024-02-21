import copy
import numpy as np
import sunpy.map
import scipy

from PyThea.sunpy_dev.map import maputils

__all__ = ['cor_polariz', 'euvi_prep']


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
        t = [map_0.date_average, map_120.date_average, map_240.date_average]
        time_avg = min(t) + (max(t)-min(t))/2

        crota = np.mean([map_0.meta['crota'], map_120.meta['crota'], map_240.meta['crota']])

        new_map_meta = map_120.meta
        # TODO: Check if any other keys need to change
        new_map_meta['date-obs'] = obs_time.strftime('%Y-%m-%dT%H:%M:%S.%f')
        new_map_meta['date-avg'] = time_avg.strftime('%Y-%m-%dT%H:%M:%S.%f')
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


def euvi_prep(map_sequence):
    smap = []
    for map_ in map_sequence:

        image = map_.data
        new_map_meta = copy.deepcopy(map_.meta)

        # Bias Subtraction
        if map_.meta['offsetcr'] == 0:
            bias = map_.meta['BIASMEAN']
            if map_.meta['ipsum'] > 1:
                bias *= ((2.0**(map_.meta['ipsum'] - 1))**2.0)
            if map_.meta['IP_PROG3'] == 95:  # Added manually
                bias *= 1.995  # Added manually
            image = image - bias

        # Normalize to Open filter position
        if map_.meta['FILTER'] == 'OPEN':
            filter_factor = 1  # no filter
        elif map_.meta['FILTER'] == 'S1':
            filter_factor = 0.49  # 1500A filter
        elif map_.meta['FILTER'] == 'S2':
            filter_factor = 0.49  # 1500A filter
        elif map_.meta['FILTER'] == 'DBL':
            filter_factor = 0.41  # 3000A filter
        image = image / filter_factor

        # Image Dejitter
        offset = [0., 0.]
        nsum = 2**(map_.meta['summed']-1)

        scale = map_.dimensions[0].value/2048

        offset[0] = scale * round(map_.meta['fpsoffz'] / 38.) / nsum
        offset[1] = scale * round(map_.meta['fpsoffy'] / 38.) / nsum

        if map_.meta['obsrvtry'] == 'STEREO_B':
            offset = [-offset[0], -offset[1]]

        # TODO: if hdr.date_obs lt '2015-05-19' then offset = -offset

        new_map_meta['crpix1'] += offset[0]
        new_map_meta['crpix2'] += offset[1]
        new_map_meta['crpix1a'] += offset[0]
        new_map_meta['crpix2a'] += offset[1]
        new_map_meta['xcen'] -= map_.meta['cdelt1'] * offset[0]
        new_map_meta['ycen'] -= map_.meta['cdelt2'] * offset[1]

        image = scipy.ndimage.interpolation.affine_transform(
            np.nan_to_num(image).T, [[1, 0], [0, 1]], offset=[-offset[0], -offset[1]], order=3,
            mode='constant', cval=0.0).T

        smap.append(sunpy.map.Map(image, new_map_meta))

    return sunpy.map.Map(smap, sequence=True)
