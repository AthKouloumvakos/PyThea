"""
Test the utilities
"""

from datetime import datetime

import matplotlib.dates as mdates
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose
from sunpy.coordinates import get_horizons_coord
from sunpy.coordinates.frames import HeliographicStonyhurst
from sunpy.net import hek

from PyThea.config.selected_bodies import bodies_dict
from PyThea.utils import get_hek_flare, parameter_fit


def test_get_hek_flare():
    """
    Test that the get_hek_flare returns a hek.hek.HEKTable and the results are formated correct.
    """
    # This test requests the flare info for a date that returns a list of flares and compares against previous results
    selectbox_list, flare_list = get_hek_flare(datetime(2017, 9, 10))
    assert type(flare_list) == hek.hek.HEKTable
    assert selectbox_list == ['FLM1.1|2017-09-09T23:53:00', 'FLC9.0|2017-09-10T03:09:00', 'FLC2.9|2017-09-10T09:20:00',
                              'FLC1.6|2017-09-10T14:23:00', 'FLC1.0|2017-09-10T15:26:00', 'FLX8.2|2017-09-10T16:06:00']

    # This test requests the flare info for a date that returns an empty list
    selectbox_list, flare_list = get_hek_flare(datetime(2117, 9, 10))
    assert flare_list == []
    assert selectbox_list == ['No events returned']


def test_get_horizons_coord():
    """
    Test that get_horizons_coord returns a Skycoord instance, in HeliographicStonyhurst frame, and returns the same coordinates every time.
    """
    bodies_coord = []
    for body in bodies_dict:
        bodies_coord.append(get_horizons_coord(bodies_dict[body][0], '2022-01-01T00:00:00'))

    coords = [[-106.56704928, 1.45508574, 0.3631443],
              [-4.65750614, -1.34575135, 0.71937524],
              [2.55987423e-06, -2.99926774, 0.98335562],
              [135.77581107, -2.68402546, 1.53663218],
              [-121.16524996, 6.08953762, 4.99274907],
              [32.43818351, -5.96794509, 1.0004962],
              [-34.99771814, 1.28725362, 0.96437702],
              [-11.53530598, -0.39313636, 0.99738376],
              [-138.81147107, 3.52846474, 0.74828612],
              [144.12380481, -2.64432759, 0.67510736]
              ]

    for i, body_coord in enumerate(bodies_coord):
        assert isinstance(body_coord, SkyCoord)
        assert isinstance(body_coord.frame, HeliographicStonyhurst)
        assert_quantity_allclose([body_coord.lon.value, body_coord.lat.value, body_coord.radius.value], coords[i], rtol=1e-5, atol=1e-4)


def test_parameter_fit_polynomial():
    x = [datetime(2020, 1, 1, 0, 0, 0), datetime(2020, 1, 1, 0, 10, 0), datetime(2020, 1, 1, 0, 20, 0), datetime(2020, 1, 1, 0, 30, 0)]
    xx = (mdates.date2num(x) - mdates.date2num(x[0]))

    y = [0, 1, 2, 3]

    fit = parameter_fit(x, y, {'type': 'polynomial', 'order': 1})
    assert_quantity_allclose([fit['popt'][0]*xx[1], fit['popt'][1]], [1, 0], atol=1e-7)

    fit = parameter_fit(x, y, {'type': 'polynomial', 'order': 2})
    assert_quantity_allclose([fit['popt'][0], fit['popt'][1]*xx[1], fit['popt'][2]], [0, 1, 0], atol=1e-5)

    y = [0, 2, 4, 6]

    fit = parameter_fit(x, y, {'type': 'polynomial', 'order': 1})
    assert_quantity_allclose([fit['popt'][0]*xx[1], fit['popt'][1]], [2, 0], atol=1e-7)

    fit = parameter_fit(x, y, {'type': 'polynomial', 'order': 2})
    assert_quantity_allclose([fit['popt'][0], fit['popt'][1]*xx[1], fit['popt'][2]], [0, 2, 0], atol=1e-5)
