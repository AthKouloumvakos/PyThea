"""
Test the extensions utilities
"""


import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose
from sunpy.coordinates import frames, get_horizons_coord

from PyThea.config.selected_bodies import bodies_dict
from PyThea.extensions import Parker_spirals


def test_parker_spiral():
    """
    Tests that the Parker spiral function returns consistent values.
    """
    time = '2022-01-01T00:00:00'
    body = bodies_dict['Earth'][0]

    pos = get_horizons_coord(body, time, 'id')
    pos = pos.transform_to(frames.HeliographicCarrington(observer='Earth', obstime=time))

    spiral_coord = Parker_spirals.spiral(pos, 350 * (u.km/u.second), time)

    test_coordinates = SkyCoord(211.21569163 * u.degree, -2.99926774 * u.degree, 1*u.Rsun,
                                frame=frames.HeliographicCarrington(observer='Earth', obstime=time))

    assert_quantity_allclose(spiral_coord[0].separation(test_coordinates), 0 * u.degree, rtol=1e-5, atol=1e-4 * u.degree)
    assert_quantity_allclose(spiral_coord[-1].separation(pos), 0 * u.degree, rtol=1e-5, atol=1e-4 * u.degree)

    spiral_coord_ = Parker_spirals.spiral(pos, 450*(u.km/u.second), time)

    assert_quantity_allclose(spiral_coord[0].separation(spiral_coord_[0]), 15.24035306 * u.degree, rtol=1e-5, atol=1e-4 * u.degree)
    assert_quantity_allclose(spiral_coord[-1].separation(spiral_coord_[-1]), 0 * u.degree, rtol=1e-5, atol=1e-4 * u.degree)


def test_footpoint():
    """
    Tests that the Parker spiral footpoint function returns consistent values.
    """
    time = '2022-01-01T00:00:00'
    body = bodies_dict['Earth'][0]

    pos = get_horizons_coord(body, time, 'id')
    pos = pos.transform_to(frames.HeliographicCarrington(observer='Earth', obstime=time))

    footpoint_coord = Parker_spirals.footpoint(pos, 350 * (u.km/u.second), time)
    spiral_coord = Parker_spirals.spiral(pos, 350 * (u.km/u.second), time)

    assert_quantity_allclose(footpoint_coord.separation(spiral_coord[0]), 0 * u.degree, rtol=1e-5, atol=1e-4 * u.degree)

    footpoint_coord_ = Parker_spirals.footpoint(pos, 450 * (u.km/u.second), time)

    assert_quantity_allclose(footpoint_coord.separation(footpoint_coord_), 15.24035306 * u.degree, rtol=1e-5, atol=1e-4 * u.degree)
