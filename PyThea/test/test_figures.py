import astropy.units as u
import matplotlib.pyplot as plt
import pytest
import sunpy.map
from astropy.coordinates import SkyCoord, solar_system_ephemeris
from sunpy.coordinates import get_body_heliographic_stonyhurst

from PyThea.data.sample_data import aia_sample_data
from PyThea.geometrical_models import ellipsoid
from PyThea.sunpy_dev.map.maputils import prepare_maps

solar_system_ephemeris.set('de432s')


@pytest.mark.mpl_image_compare()
def test_ellipsoid_on_AIA_venus():
    """
    Plot the position of Venus as it transited in front
    of the Sun as observed by SDO/AIA.
    """
    venus_radius = 6051.8 * u.km  # From here: https://nssdc.gsfc.nasa.gov/planetary/factsheet/venusfact.html
    aia_fits_venus = aia_sample_data.fetch('aia_lev1_1600a_2012_06_06t04_07_29_12z_image_lev1.fits')

    aiamap = sunpy.map.Map(aia_fits_venus, sequence=True)
    aiamap = prepare_maps(aiamap)[0]

    venus = get_body_heliographic_stonyhurst('venus', aiamap.date_average, observer=aiamap.observer_coordinate)

    center = SkyCoord(venus)
    model_shock = ellipsoid(center, venus_radius, venus_radius, venus_radius, 0 * u.degree)

    venus_hpc = venus.transform_to(aiamap.coordinate_frame)
    fov = 200 * u.arcsec
    bottom_left = SkyCoord(venus_hpc.Tx - fov/2, venus_hpc.Ty - fov/2, frame=aiamap.coordinate_frame)
    smap = aiamap.submap(bottom_left, width=fov, height=fov)

    fig = plt.figure(dpi=200)
    ax = fig.add_subplot(projection=smap)
    smap.plot(axes=ax)
    smap.draw_limb(axes=ax)
    ax.grid(False)
    model_shock.plot(ax, mode='Full')
    ax.set_xlim([0, smap.data.shape[0]])
    ax.set_ylim([0, smap.data.shape[1]])

    return fig


@pytest.mark.mpl_image_compare()
def test_ellipsoid_on_AIA_mercury():
    """
    Plot the position of Venus as it transited in front
    of the Sun as observed by SDO/AIA.
    """
    mercury_radius = 2440.5 * u.km  # From here: https://nssdc.gsfc.nasa.gov/planetary/factsheet/venusfact.html
    aia_fits_mercury = aia_sample_data.fetch('aia_lev1_1600a_2016_05_09t11_35_51_12z_image_lev1.fits')

    aiamap = sunpy.map.Map(aia_fits_mercury, sequence=True)
    aiamap = prepare_maps(aiamap)[0]

    mercury = get_body_heliographic_stonyhurst('mercury', aiamap.date_average, observer=aiamap.observer_coordinate)

    center = SkyCoord(mercury)
    model_shock = ellipsoid(center, mercury_radius, mercury_radius, mercury_radius, 0 * u.degree)

    mercury_hpc = mercury.transform_to(aiamap.coordinate_frame)
    fov = 100 * u.arcsec
    bottom_left = SkyCoord(mercury_hpc.Tx - fov/2, mercury_hpc.Ty - fov/2, frame=aiamap.coordinate_frame)
    smap = aiamap.submap(bottom_left, width=fov, height=fov)

    fig = plt.figure(dpi=200)
    ax = fig.add_subplot(projection=smap)
    smap.plot(axes=ax)
    smap.draw_limb(axes=ax)
    ax.grid(False)
    model_shock.plot(ax, mode='Full')
    ax.set_xlim([0, smap.data.shape[0]])
    ax.set_ylim([0, smap.data.shape[1]])

    return fig
