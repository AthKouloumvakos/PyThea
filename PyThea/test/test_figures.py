import astropy.units as u
import matplotlib.pyplot as plt
import pytest
import sunpy.map
from astropy.coordinates import SkyCoord, solar_system_ephemeris
from sunpy.coordinates import get_body_heliographic_stonyhurst

from PyThea.data.sample_data import (aia_sample_data,
                                     json_fitting_file_sample_data,
                                     stereo_sample_data, wispr_sample_data)
from PyThea.geometrical_models import ellipsoid
from PyThea.sunpy_dev.map.maputils import prepare_maps
from PyThea.utils import model_fittings, plot_fitting_model

solar_system_ephemeris.set('de432s')

bodies_color_dict = {'mercury': 'darkturquoise',
                     'venus': 'darkorchid',
                     'earth': 'green',
                     'mars': 'maroon',
                     'jupiter': 'navy',
                     'saturn': 'dodgerblue'}

body_radious = {'mercury': 2440.5 * u.km,
                'venus': 6051.8 * u.km,
                'mars': 3396.2 * u.km,
                'jupiter': 71492 * u.km,
                'saturn': 60268 * u.km}  # From here: https://nssdc.gsfc.nasa.gov/planetary/factsheet/venusfact.html


def model_from_body(body, date, observer_coordinate):
    body_coord = get_body_heliographic_stonyhurst(body, date, observer=observer_coordinate)
    radius =  body_radious[body]
    center = SkyCoord(body_coord)
    return ellipsoid(center, radius, radius, radius, 0 * u.degree)


def make_submap(map_, center, fov):
    center_hpc = center.transform_to(map_.coordinate_frame)
    print(center_hpc)
    bottom_left = SkyCoord(center_hpc.Tx - fov/2, center_hpc.Ty - fov/2, frame=map_.coordinate_frame)
    return map_.submap(bottom_left, width=fov, height=fov)


@pytest.mark.mpl_image_compare()
def test_ellipsoid_on_AIA_venus():
    """
    Tests that the position of Venus is ploted correct at SDO/AIA images when it transited in front of the Sun.
    """
    fits_file = aia_sample_data.fetch('aia_lev1_1600a_2012_06_06t04_07_29_12z_image_lev1.fits')

    map_ = sunpy.map.Map(fits_file, sequence=True)
    map_ = prepare_maps(map_)[0]

    model_shock = model_from_body('venus', map_.date_average, map_.observer_coordinate)
    smap = make_submap(map_, model_shock.center, fov=200*u.arcsec)

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
    Tests that the position of Mercury is ploted correct at SDO/AIA images when it transited in front of the Sun.
    """
    fits_file = aia_sample_data.fetch('aia_lev1_1600a_2016_05_09t11_35_51_12z_image_lev1.fits')

    map_ = sunpy.map.Map(fits_file, sequence=True)
    map_ = prepare_maps(map_)[0]

    model_shock = model_from_body('mercury', map_.date_average, map_.observer_coordinate)
    smap = make_submap(map_, model_shock.center, fov=100*u.arcsec)

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
def test_ellipsoid_on_STEREO_COR1_venus():
    """
    Tests that the position of Venus is correct at STEREO COR1 images. (https://stereo.gsfc.nasa.gov/browse/2019/03/08/planets/)
    """
    fits_file = [stereo_sample_data.fetch('20070503_102018_s4c1b.fts'),
                 stereo_sample_data.fetch('20070503_102018_s4c1a.fts')]

    maps = sunpy.map.Map(fits_file, sequence=True)
    maps = prepare_maps(maps, polar=True)

    fig = plt.figure(figsize=(10, 5), dpi=200)

    for i, map_ in enumerate(maps):
        axis = fig.add_subplot(int(120+i+1), projection=map_)

        map_.plot(axes=axis)
        map_.draw_limb(axes=axis)
        axis.grid(False)

        model_shock = model_from_body('mercury', map_.date_average, map_.observer_coordinate)
        axis.plot_coord(model_shock.center, 'o', color='tab:red', fillstyle='none', markersize=8)
        x = map_.world_to_pixel(model_shock.center)[0].value
        y = map_.world_to_pixel(model_shock.center)[1].value+40
        axis.annotate('Mercury', (x, y), xytext=(x, y),
                      backgroundcolor='none', color='Black', fontsize=8,
                      horizontalalignment='center', verticalalignment='bottom')

        lon, lat = axis.coords
        lon .set_major_formatter('d.dd')
        lat.set_major_formatter('d.dd')

        axis.set_xlim([0, map_.data.shape[0]])
        axis.set_ylim([0, map_.data.shape[1]])

    fig.subplots_adjust(hspace=0.5, wspace=0.35)

    return fig


@pytest.mark.mpl_image_compare()
def test_ellipsoid_on_STEREO_COR2_three_planets():
    """
    Tests that the position of three planets is correct at STEREO COR2 images. (https://stereo.gsfc.nasa.gov/browse/2021/06/22/planets/)
    """
    fits_file = [stereo_sample_data.fetch('20120622_235400_d4c2a.fts'),
                 stereo_sample_data.fetch('20140221_235400_d4c2b.fts')]

    maps = sunpy.map.Map(fits_file, sequence=True)
    maps = prepare_maps(maps)

    fig = plt.figure(figsize=(10, 5), dpi=200)

    for i, map_ in enumerate(maps):
        axis = fig.add_subplot(int(120+i+1), projection=map_)

        map_.plot(axes=axis, vmin=300, vmax=3000)
        map_.draw_limb(axes=axis)
        axis.grid(False)

        for body in ['mercury', 'venus', 'mars', 'saturn']:
            model_shock = model_from_body(body, map_.date_average, map_.observer_coordinate)
            axis.plot_coord(model_shock.center, 'o', color=bodies_color_dict[body], fillstyle='none', markersize=8)
            x = map_.world_to_pixel(model_shock.center)[0].value
            y = map_.world_to_pixel(model_shock.center)[1].value+40
            axis.annotate(body.capitalize(), (x, y), xytext=(x, y),
                          backgroundcolor='none', color='Black', fontsize=8,
                          horizontalalignment='center', verticalalignment='bottom')
        lon, lat = axis.coords
        lon .set_major_formatter('d.dd')
        lat.set_major_formatter('d.dd')

        axis.set_xlim([0, map_.data.shape[0]])
        axis.set_ylim([0, map_.data.shape[1]])

    fig.subplots_adjust(hspace=0.5, wspace=0.35)

    return fig


@pytest.mark.mpl_image_compare()
def test_kinematics_figure():
    json_fitting_file = json_fitting_file_sample_data.fetch('FLX1p0D20211028T153500MEllipsoid.json')
    model_fittings_class = model_fittings.load_from_json(json_fitting_file)

    fit_method = {'type': 'polynomial', 'order': 2}

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5), tight_layout=True)

    fig, ax1 = plot_fitting_model(model_fittings_class,
                                  fit_args=fit_method,
                                  plt_type='HeightT',
                                  fig=fig, axis=ax1)
    fig, ax2 = plot_fitting_model(model_fittings_class,
                                  fit_args=fit_method,
                                  plt_type='SpeedT',
                                  fig=fig, axis=ax2)
    fig, ax3 = plot_fitting_model(model_fittings_class,
                                  fit_args=fit_method,
                                  plt_type='AccelerationT',
                                  fig=fig, axis=ax3)
    return fig


@pytest.mark.mpl_image_compare()
def test_ellipsoid_on_WISPR_Inner_Venus():
    """
    Tests that the position of Venus is ploted correct at WISPR-Inner telescope images.
    """
    fits_file = wispr_sample_data.fetch('psp_l3_wispr_20210808t103707_v1_1211.fits')

    map_ = sunpy.map.Map(fits_file, sequence=True)
    map_ = prepare_maps(map_)[0]

    model_shock = model_from_body('venus', map_.date_average, map_.observer_coordinate)

    fig = plt.figure(dpi=200)
    ax = fig.add_subplot(projection=map_)
    map_.plot(axes=ax)
    map_.draw_limb(axes=ax)
    ax.grid(False)

    model_shock = model_from_body('venus', map_.date_average, map_.observer_coordinate)
    ax.plot_coord(model_shock.center, 'o', color=bodies_color_dict['venus'], fillstyle='none', markersize=6, markeredgewidth=2)
    x = map_.world_to_pixel(model_shock.center)[0].value
    y = map_.world_to_pixel(model_shock.center)[1].value+20
    ax.annotate('Venus', (x, y), xytext=(x, y),
                backgroundcolor='none', color='white', fontsize=8,
                horizontalalignment='center', verticalalignment='bottom')
    lon, lat = ax.coords
    lon .set_major_formatter('d.dd')
    lat.set_major_formatter('d.dd')

    ax.set_xlim([0, map_.data.shape[0]-80])
    ax.set_ylim([0, map_.data.shape[1]])

    return fig


@pytest.mark.mpl_image_compare()
def test_ellipsoid_on_WISPR_Outer_Mercury():
    """
    Tests that the position of Mercury is ploted correct at WISPR-Outer telescope images.
    """
    fits_file = wispr_sample_data.fetch('psp_l3_wispr_20210808t104010_v1_2222.fits')

    map_ = sunpy.map.Map(fits_file, sequence=True)
    map_ = prepare_maps(map_)[0]

    model_shock = model_from_body('mercury', map_.date_average, map_.observer_coordinate)

    fig = plt.figure(dpi=200)
    ax = fig.add_subplot(projection=map_)
    map_.plot(axes=ax)
    map_.draw_limb(axes=ax)
    ax.grid(False)

    model_shock = model_from_body('mercury', map_.date_average, map_.observer_coordinate)
    ax.plot_coord(model_shock.center, 'o', color=bodies_color_dict['mercury'], fillstyle='none', markersize=6, markeredgewidth=2)
    x = map_.world_to_pixel(model_shock.center)[0].value
    y = map_.world_to_pixel(model_shock.center)[1].value+20
    ax.annotate('Mercury', (x, y), xytext=(x, y),
                backgroundcolor='none', color='white', fontsize=8,
                horizontalalignment='center', verticalalignment='bottom')
    lon, lat = ax.coords
    lon .set_major_formatter('d.dd')
    lat.set_major_formatter('d.dd')

    ax.set_xlim([0, map_.data.shape[0]-80])
    ax.set_ylim([0, map_.data.shape[1]])

    return fig
