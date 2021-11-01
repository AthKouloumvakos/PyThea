import pytest

from sunpy.net import attr, attrs, hek
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose
import astropy.units as u
import sunpy.coordinates

from PyThea.geometrical_models import sphere, spheroid, ellipsoid
from PyThea.utils import get_hek_flare

'''
Before you start the test the PyThea pkg should be in the python path
export PYTHONPATH="${PYTHONPATH}:{top_level_dir_that_pythea_lives}/PyThea"
'''

@pytest.mark.remote_data
def test_hek_client():
    startTime = '2011/08/09 07:23:56'
    endTime = '2011/08/09 12:40:29'
    eventType = 'FL'

    hekTime = attrs.Time(startTime, endTime)
    hekEvent = attrs.hek.EventType(eventType)

    h = hek.HEKClient()
    hek_query = h.search(hekTime, hekEvent)
    assert type(hek_query) == hek.hek.HEKTable
    
# ------------------
# Test Geometrical Models

def test_sphere_shape():
    x, y, z = sphere(n=20)
    assert (x.shape, y.shape, z.shape) == ((21, 21), (21, 21), (21, 21))

def test_sphere_isunit():
    x, y, z = sphere(n=20)
    assert (x**2 + y**2 + z**2).all() == 1

@pytest.fixture
def center():
    return SkyCoord(0*u.deg, 0*u.deg, 0*u.R_sun,
                    frame="heliographic_stonyhurst",
                    observer="earth",
                    obstime="2020/01/01T00:00:30")

def test_spheroid_pointlike(center):
    sph = spheroid(center, 0.*u.R_sun, 0.*u.R_sun)
    assert (sph.coordinates == center).all()

def test_spheroid_isunit(center):
    sph = spheroid(center, 1.*u.R_sun, 1.*u.R_sun)
    assert_quantity_allclose(sph.coordinates.radius, 1.*u.R_sun, atol=5e-7*u.R_sun)
    
# TODO def test_spheroid_isunit_afrer_rotation(center):

def test_spheroid_apex_location():
    center_ = SkyCoord(0*u.deg, 0*u.deg, 1*u.R_sun,
                       frame="heliographic_stonyhurst",
                       observer="earth",
                       obstime="2020/01/01T00:00:30")
    apex_ = SkyCoord(0*u.deg, 0*u.deg, 2*u.R_sun,
                       frame="heliographic_stonyhurst",
                       observer="earth",
                       obstime="2020/01/01T00:00:30")
    sph = spheroid(center_, 1.*u.R_sun, 1.*u.R_sun)
    assert_quantity_allclose(sph.apex.radius, apex_.radius, atol=5e-7*u.R_sun)

def test_spheroid_apex_location():
    center_ = SkyCoord(0*u.deg, 0*u.deg, 1*u.R_sun,
                       frame="heliographic_stonyhurst",
                       observer="earth",
                       obstime="2020/01/01T00:00:30")
    apex_ = SkyCoord(0*u.deg, 0*u.deg, 2*u.R_sun,
                       frame="heliographic_stonyhurst",
                       observer="earth",
                       obstime="2020/01/01T00:00:30")
    sph = spheroid(center_, 1.*u.R_sun, 1.*u.R_sun)
    assert_quantity_allclose(sph.apex.radius, apex_.radius, atol=5e-7*u.R_sun)

# TODO @pytest.mark.parametrize('num, expected', [(-1, 1), (0, 2), (2, 4)])
def test_spheroid_parameters_conversions():
    he, ep, ka = spheroid.hek_from_rab(1*u.R_sun, 2*u.R_sun, 3*u.R_sun)
    r, a, b = spheroid.rab_from_hek(he, ep, ka)
    he_, ep_, ka_ = spheroid.hek_from_rab(r, a, b)
    assert_quantity_allclose(he, he_, atol=5e-7*u.R_sun)
    assert_quantity_allclose(ep, ep_, atol=5e-7)
    assert_quantity_allclose(ka, ka_, atol=5e-7)

def test_ellipsoid_semiaxis_alpha(center):
    ell = ellipsoid(center, 1.*u.R_sun, 1.*u.R_sun, 1.*u.R_sun, 0.*u.degree)
    assert_quantity_allclose(ell.alpha, 1, atol=5e-7)
    ell = ellipsoid(center, 1.*u.R_sun, 1.*u.R_sun, 2.*u.R_sun, 0.*u.degree)
    assert_quantity_allclose(ell.alpha, 0.5, atol=5e-7)
    ell = ellipsoid(center, 1.*u.R_sun, 2.*u.R_sun, 1.*u.R_sun, 0.*u.degree)
    assert_quantity_allclose(ell.alpha, 2, atol=5e-7)

def test_ellipsoid_parameters_conversions():
    he, ep, ka, al = ellipsoid.heka_from_rabc(1*u.R_sun, 2*u.R_sun, 3*u.R_sun, 4*u.R_sun)
    r, a, b, c = ellipsoid.rabc_from_heka(he, ep, ka, al)
    he_, ep_, ka_, al_ = ellipsoid.heka_from_rabc(r, a, b, c)
    assert_quantity_allclose(he, he_, atol=5e-7*u.R_sun)
    assert_quantity_allclose(ep, ep_, atol=5e-7)
    assert_quantity_allclose(ka, ka_, atol=5e-7)
    assert_quantity_allclose(al, al_, atol=5e-7)

def test_spheroid_ellipsoid_equivalency(center):
    sph = spheroid(center, 1.5*u.R_sun, 1.5*u.R_sun)
    ell = ellipsoid(center, 1.5*u.R_sun, 1.5*u.R_sun, 1.5*u.R_sun, 0.*u.degree)
    assert_quantity_allclose(sph.coordinates.radius, ell.coordinates.radius, atol=5e-7*u.R_sun)

# ------------------
# Test utilities
def test_get_hek_flare_HEKTable(center):
    from datetime import datetime
    day = datetime(2017, 9, 10)
    _, flare_list_ = get_hek_flare(day)

    assert type(flare_list_) == hek.hek.HEKTable


    
# TODO Compare the stony and carr center values
