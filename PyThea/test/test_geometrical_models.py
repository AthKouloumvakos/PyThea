"""
Test the geometrical models
"""

import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose

from PyThea.geometrical_models import ellipsoid, sphere, spheroid


@pytest.fixture
def center():
    return SkyCoord(0*u.deg, 0*u.deg, 0*u.R_sun,
                    frame='heliographic_stonyhurst',
                    observer='earth',
                    obstime='2020/01/01T00:00:30')


def test_sphere_shape():
    """
    Tests if sphere produces the right matrix dimensions
    """
    x, y, z = sphere(n=20)
    assert (x.shape, y.shape, z.shape) == ((21, 21), (21, 21), (21, 21))


def test_sphere_isunit():
    """
    Tests that sphere is unit shpere
    """
    x, y, z = sphere(n=20)
    assert (x**2 + y**2 + z**2).all() == 1


def test_spheroid_pointlike(center):
    """
    Tests that spheroid is point-like when both semi-axes length is zero
    """
    sph = spheroid(center, 0.*u.R_sun, 0.*u.R_sun)
    assert (sph.coordinates == center).all()


def test_spheroid_isunit(center):
    """
    Tests that spheroid is unit sphere when both semi-axes are equal to one
    """
    sph = spheroid(center, 1.*u.R_sun, 1.*u.R_sun)
    assert_quantity_allclose(sph.coordinates.radius, 1.*u.R_sun, atol=5e-7*u.R_sun)


# TODO def test_spheroid_isunit_afrer_rotation(center):


def test_spheroid_apex_location():
    """
    Tests that the spheroid apex location is correct
    """
    center_ = SkyCoord(0*u.deg, 0*u.deg, 1*u.R_sun,
                       frame='heliographic_stonyhurst',
                       observer='earth',
                       obstime='2020/01/01T00:00:30')
    apex_ = SkyCoord(0*u.deg, 0*u.deg, 2*u.R_sun,
                       frame='heliographic_stonyhurst',
                       observer='earth',
                       obstime='2020/01/01T00:00:30')
    sph = spheroid(center_, 1.*u.R_sun, 1.*u.R_sun)
    assert_quantity_allclose(sph.apex.radius, apex_.radius, atol=5e-7*u.R_sun)


# TODO @pytest.mark.parametrize('num, expected', [(-1, 1), (0, 2), (2, 4)])
def test_spheroid_parameters_conversions():
    """
    Tests that the spheroid parameters conversion is correct
    """
    he, ep, ka = spheroid.hek_from_rab(1*u.R_sun, 2*u.R_sun, 3*u.R_sun)
    r, a, b = spheroid.rab_from_hek(he, ep, ka)
    he_, ep_, ka_ = spheroid.hek_from_rab(r, a, b)
    assert_quantity_allclose(he, he_, atol=5e-7*u.R_sun)
    assert_quantity_allclose(ep, ep_, atol=5e-7)
    assert_quantity_allclose(ka, ka_, atol=5e-7)


def test_ellipsoid_semiaxis_alpha(center):
    """
    Tests that the ellipsoid alpha parameter is calculated correct
    """
    ell = ellipsoid(center, 1.*u.R_sun, 1.*u.R_sun, 1.*u.R_sun, 0.*u.degree)
    assert_quantity_allclose(ell.alpha, 1, atol=5e-7)
    ell = ellipsoid(center, 1.*u.R_sun, 1.*u.R_sun, 2.*u.R_sun, 0.*u.degree)
    assert_quantity_allclose(ell.alpha, 0.5, atol=5e-7)
    ell = ellipsoid(center, 1.*u.R_sun, 2.*u.R_sun, 1.*u.R_sun, 0.*u.degree)
    assert_quantity_allclose(ell.alpha, 2, atol=5e-7)


def test_ellipsoid_parameters_conversions():
    """
    Tests that the ellipsoid parameters conversion is correct
    """
    he, ep, ka, al = ellipsoid.heka_from_rabc(1*u.R_sun, 2*u.R_sun, 3*u.R_sun, 4*u.R_sun)
    r, a, b, c = ellipsoid.rabc_from_heka(he, ep, ka, al)
    he_, ep_, ka_, al_ = ellipsoid.heka_from_rabc(r, a, b, c)
    assert_quantity_allclose(he, he_, atol=5e-7*u.R_sun)
    assert_quantity_allclose(ep, ep_, atol=5e-7)
    assert_quantity_allclose(ka, ka_, atol=5e-7)
    assert_quantity_allclose(al, al_, atol=5e-7)


def test_ellipsoid_unitsphere_equivalency(center):
    """
    Tests that an ellipsoid with all the semi-axis with length one is unit sphere
    """
    ell = ellipsoid(center, 1.*u.R_sun, 1.*u.R_sun, 1.*u.R_sun, 0.*u.degree)
    assert_quantity_allclose(ell.coordinates.radius, 1.*u.R_sun, atol=5e-7*u.R_sun)


def test_spheroid_ellipsoid_equivalency(center):
    """
    Tests that an ellipsoid and a spheroid are equivalent
    """
    sph = spheroid(center, 1.5*u.R_sun, 1.5*u.R_sun)
    ell = ellipsoid(center, 1.5*u.R_sun, 1.5*u.R_sun, 1.5*u.R_sun, 0.*u.degree)
    assert_quantity_allclose(sph.coordinates.radius, ell.coordinates.radius, atol=5e-7*u.R_sun)
    sph = spheroid(center, 2*u.R_sun, 1.5*u.R_sun)
    ell = ellipsoid(center, 2*u.R_sun, 1.5*u.R_sun, 1.5*u.R_sun, 0.*u.degree)
    assert_quantity_allclose(sph.coordinates.radius, ell.coordinates.radius, atol=5e-7*u.R_sun)

# TODO Compare the stony and carr center values
