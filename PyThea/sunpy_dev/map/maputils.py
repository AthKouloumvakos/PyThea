"""
This submodule provides utility functions to act on `sunpy.map.GenericMap` instances.
"""

import copy
import warnings

import astropy.units as u
import numpy as np
import sunpy.map
from sunpy.map.mapsequence import MapSequence

import PyThea.sunpy_dev.extern.sunkit_instruments.aia.utils  # noqa
import PyThea.sunpy_dev.extern.sunkit_instruments.lasco.utils  # noqa
import PyThea.sunpy_dev.extern.sunkit_instruments.stereo.utils  # noqa

__all__ = ['maps_sequence_processing', 'get_closest', 'normalize_exposure',
           'prepare_maps', 'difference_maps', 'mask_occulter']


def check_maps_sequence(map_sequence):
    if map_sequence:
        if isinstance(map_sequence, MapSequence):
            return map_sequence
        else:
            return sunpy.map.Map(map_sequence, sequence=True)
    else:
        return []


def maps_sequence_processing(map_sequence, **kwargs):
    """
    Returns a sequence of maps which is processed as plain images, running or base difference images.

    Parameters
    ----------
    map_sequence : `~sunpy.map.GenericMap`
        A list SunPy maps.

    image_mode : '~str'
        The type of sequence processing: 'Running Diff.', 'Base Diff.', 'Plain'

    Returns
    -------
    `~sunpy.map.GenericMap`
        A SunPy map.
    """

    map_sequence = check_maps_sequence(map_sequence)
    if len(map_sequence) == 0:
        return []

    seq_type = kwargs.get('image_mode', 'Plain')
    diff_num = kwargs.get('diff_num', 1 if seq_type == 'Running Diff.' else 0 if seq_type == 'Base Diff.' else None)

    normalized = True if True in [tmap.exposure_time == 1.0*u.second for tmap in map_sequence] else False
    if not normalized:
        warnings.warn('Warning [maps_sequence_processing]: The exposure time of the maps is not normalized.')

    smap = []
    if seq_type == 'Running Diff.':
        for i in range(diff_num, len(map_sequence)):
            smap_diff = difference_maps(map_sequence[i], map_sequence[i-diff_num])
            smap.append(smap_diff)
    if seq_type == 'Base Diff.':
        for i in range(1, len(map_sequence)):
            smap_diff = difference_maps(map_sequence[i], map_sequence[diff_num])
            smap.append(smap_diff)
    if seq_type == 'Plain':
        for i in range(0, len(map_sequence)):
            smap.append(map_sequence[i])

    if len(smap) != 0:
        sequence_final = sunpy.map.Map(smap, sequence=True)
    else:
        sequence_final = []

    return sequence_final


def get_closest(smap, date):
    """
    Returns the map closest to the date provided from a sequence
    or a list of smaps.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A SunPy map.

    date : 'astropy.time.Time'
        A time object.

    Returns
    -------
    `~sunpy.map.GenericMap`
        A SunPy map.
    """
    delta_time = [abs(map_.date_average-date) for map_ in smap]
    map_closest = smap[delta_time.index(min(delta_time))]
    return map_closest


def normalize_exposure(smap):
    """
    Normalize the exposure time of a map.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A SunPy map.

    Returns
    -------
    `~sunpy.map.GenericMap`
        A SunPy map.
    """
    if smap.exposure_time == 1*u.second:
        return smap
    new_meta = copy.deepcopy(smap.meta)
    if 'exptime' in new_meta:
        new_meta['exptime'] = 1.0
    elif 'xposure' in new_meta:
        new_meta['xposure'] = 1.0
    smap_new = sunpy.map.Map(smap.data/smap.exposure_time.to_value(u.second), new_meta)
    if smap_new.exposure_time != 1*u.second:
        raise Exception('Could not normalize Maps exposure time.')
    return smap_new


def filter_maps(map_sequence, **kwargs):
    '''
    Returns filtered maps.

    For a map (or a sequence of maps) it performs a filtering, namely:
    removes map dublicates, maps with non-default exposure or dimensions, or polarized maps,
    combine the polarization images to total brightness,

    Parameters
    ----------
    map_sequence : `~sunpy.map.GenericMap`
        A SunPy map.

    kwargs : A list with with the preparation arguments

    Returns
    -------
    `~sunpy.map.GenericMap`
        A SunPy map.

    '''

    map_sequence = check_maps_sequence(map_sequence)
    if len(map_sequence) == 0:
        return []

    if 'processing_level' in kwargs:
        map_sequence = [tmap for tmap in map_sequence if tmap.processing_level == kwargs['processing_level']]

    if 'exposure' in kwargs.keys():
        map_sequence = [tmap for tmap in map_sequence if tmap.exposure_time > kwargs['exposure']*u.second]

    if 'dimensions' in kwargs.keys():
        map_sequence = [tmap for tmap in map_sequence if (tmap.dimensions[0], tmap.dimensions[1]) == kwargs['dimensions']]

    if 'polar' in kwargs.keys():
        map_sequence = [tmap for tmap in map_sequence if tmap.meta['polar'] == kwargs['polar']]

    if len(map_sequence) != 0:
        sequence_final = check_maps_sequence(map_sequence)
    else:
        sequence_final = []

    return sequence_final


def prepare_maps(map_sequence, **kwargs):
    '''
    Returns prepared maps.

    For a map (or a sequence of maps) it performs basic preparations, namely:
    exposure time normalization, mask the occulters, downsample with superpixel,
    combine the polarization images to total brightness,

    Parameters
    ----------
    map_sequence : `~sunpy.map.GenericMap`
        A SunPy map.

    kwargs : A list with with the preparation arguments

    Returns
    -------
    `~sunpy.map.GenericMap`
        A SunPy map.

    '''

    map_sequence = check_maps_sequence(map_sequence)
    if len(map_sequence) == 0:
        return []

    detector = map_sequence[0].detector
    print(f'Preparing image sequence for {detector}. This could take a while...')

    # Prepare the maps before anything else
    if detector == 'AIA':
        map_sequence = PyThea.sunpy_dev.extern.sunkit_instruments.aia.utils.prep_aia(map_sequence)
    elif detector == 'COR1':
        if 'polar' not in kwargs.keys():
            map_sequence = PyThea.sunpy_dev.extern.sunkit_instruments.stereo.utils.cor_polariz(map_sequence)
    elif detector == 'EUVI':
        map_sequence = PyThea.sunpy_dev.extern.sunkit_instruments.stereo.utils.euvi_prep(map_sequence)
    elif detector in ['C2', 'C3']:
        map_sequence = PyThea.sunpy_dev.extern.sunkit_instruments.lasco.utils.prep_lasco(map_sequence)

    map_sequence = [normalize_exposure(tmap) for tmap in map_sequence]

    map_sequence = [mask_occulter(tmap) for tmap in map_sequence]

    if 'superpixel' in kwargs.keys():
        nsuper = kwargs['superpixel']
        super_dim = u.Quantity([nsuper, nsuper] * u.pixel)
        map_sequence = [tmap.superpixel(super_dim) for tmap in map_sequence]

    # map_sequence = [tmap.rotate(recenter=True) for tmap in map_sequence]

    return map_sequence


def difference_maps(smapi, smapm):
    """
    Returns a running or base difference map.

    Parameters
    ----------
    smapi : `~sunpy.map.GenericMap`
        A SunPy map.

    smapm : `~sunpy.map.GenericMap`
        A SunPy map.

    Returns
    -------
    `~sunpy.map.GenericMap`
        A SunPy map.

    Note
    -------
    The two maps have to be coalligned.
    This is something that it is not checked here.
    """
    if smapi.exposure_time != 1*u.second:
        smapi = normalize_exposure(smapi)
    if smapm.exposure_time != 1*u.second:
        smapm = normalize_exposure(smapm)
    smap_difference = smapi.data - smapm.data

    return sunpy.map.Map(smap_difference, smapi.meta)


def mask_occulter(smap, apply_mask=True, mask_value=0):
    """
    For a given map it produces a mask for the coronagraphic occulters
    and returns the map with the mash applied or the mask. If the map
    is not from coronagraphs or it is not registered returns the initial map

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A SunPy map.

    Returns
    -------
    `~sunpy.map.GenericMap`
        A SunPy map.
    """
    if smap.detector == 'COR1':
        xcen, ycen = smap.meta['crpix1'], smap.meta['crpix2']
        in_fac, out_fac = (100/512), (280/512)
    elif smap.detector == 'COR2':
        xcen, ycen = smap.meta['crpix1'], smap.meta['crpix2']
        in_fac, out_fac = (150/2048), (1030/2048)
    elif smap.detector == 'C2':
        xcen = (smap.meta['crval1'] / smap.meta['cdelt1']) +  smap.meta['crpix1']
        ycen = (smap.meta['crval2'] / smap.meta['cdelt2']) + smap.meta['crpix2']
        in_fac, out_fac = (170/1024), (670/1024)
    elif smap.detector == 'C3':
        xcen = (smap.meta['crval1'] / smap.meta['cdelt1']) + smap.meta['crpix1']
        ycen = (smap.meta['crval2'] / smap.meta['cdelt2']) + smap.meta['crpix2']
        in_fac, out_fac = (60/1024), (560/1024)
    else:
        if apply_mask:
            return smap
        return []

    inner_distance = in_fac * smap.dimensions.x
    outer_distance = out_fac * smap.dimensions.x

    x, y = np.meshgrid(*[np.arange(v.value) for v in smap.dimensions]) * u.pixel
    xprime_sq = (x-xcen * u.pixel) ** 2
    yprime_sq = (y-ycen * u.pixel) ** 2
    mask_inner = np.sqrt(xprime_sq + yprime_sq) < inner_distance
    mask_outer = np.sqrt(xprime_sq + yprime_sq) > outer_distance

    if apply_mask:
        smap.data[mask_inner+mask_outer] = mask_value
        return sunpy.map.Map(smap.data, smap.meta)
    return mask_inner+mask_outer


def maps_clims(images):
    i = 0 if len(images) < 2 else 1

    if images[i].instrument == 'WISPR':
        if images[i].detector == 'Outer':
            return [14., 14.]
        elif images[i].detector == 'Inner':
            return [13., 13.]
    elif images[i].instrument == 'Metis':
        return [12.80, 13.00]
    elif images[i].instrument.startswith('SoloHI'):
        return [13.95, 14.05]

    return [np.nanquantile(images[i].data, 0.20)-10, np.nanquantile(images[i].data, 0.80)+10]
