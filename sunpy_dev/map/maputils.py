"""
This submodule provides utility functions to act on `sunpy.map.GenericMap` instances.
"""

import copy
import sunpy.map
import astropy.units as u
import numpy as np


__all__ = ['get_closest', 'normalize_exposure', 'difference_maps', 'mask_occulter']


def get_closest(smap, date):
    """
    Returns the map closest to the date provided from a sequence
    or a list of smaps.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A SunPy map.
    
    date : '~'

    Returns
    -------
    `~sunpy.map.GenericMap`
        A SunPy map.
    """
    delta_time = [abs(map_.date-date) for map_ in smap]
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
    exp_time = smap.exposure_time
    new_meta = copy.deepcopy(smap.meta)
    if 'exptime' in new_meta:
        new_meta['exptime'] = 1.0
    elif 'xposure' in new_meta:
        new_meta['xposure'] = 1.0
    smap_new = sunpy.map.Map(smap.data/smap.exposure_time.to_value(u.second), new_meta)
    if smap_new.exposure_time != 1*u.second:
        raise Exception('Could not normalize Maps exposure time.')
    return smap_new

def difference_maps(smapi, smapm):
    """
    Make a running or base difference map.

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
    smapi = normalize_exposure(smapi)
    smapm = normalize_exposure(smapm)
    smap_difference = smapi.data - smapm.data
    return sunpy.map.Map(smap_difference, smapi.meta)

def mask_occulter(smap, apply_mask = True, mask_value = 0):
    """
    Returns a mask of the coronagraphic occulters
    """
    if smap.detector == 'COR1':
        xcen, ycen = smap.meta['crpix1'], smap.meta['crpix2']
        in_fac, out_fac = (100/512), (280/512)
    elif smap.detector == 'COR2':
        xcen, ycen = smap.meta['crpix1'], smap.meta['crpix2']
        in_fac, out_fac = (150/2048), (1030/2048)
    elif smap.detector == 'C2':
        xcen = smap.meta['crval1'] / smap.meta['cdelt1'] + ( smap.meta['crpix1'])
        ycen = smap.meta['crval2'] / smap.meta['cdelt2'] + ( smap.meta['crpix2'])
        in_fac, out_fac = (170/1024), (670/1024)
    elif smap.detector == 'C3':
        xcen = smap.meta['crval1'] / smap.meta['cdelt1'] + ( smap.meta['crpix1'])
        ycen = smap.meta['crval2'] / smap.meta['cdelt2'] + ( smap.meta['crpix2'])
        in_fac, out_fac = (60/1024), (560/1024)
    else:
        if apply_mask:
            return smap
        else:
            return []
    
    inner_distance = in_fac * smap.dimensions.x
    outer_distance = out_fac * smap.dimensions.x
    
    x, y = np.meshgrid(*[np.arange(v.value) for v in smap.dimensions]) * u.pixel
    xprime_sq = (x-xcen* u.pixel) ** 2
    yprime_sq = (y-ycen* u.pixel) ** 2
    mask_inner = np.sqrt( xprime_sq + yprime_sq ) < inner_distance
    mask_outer = np.sqrt( xprime_sq + yprime_sq ) > outer_distance
    
    if apply_mask:
        smap.data[mask_inner+mask_outer] = mask_value
        return sunpy.map.Map(smap.data, smap.meta)
    else:
        return mask_inner+mask_outer
