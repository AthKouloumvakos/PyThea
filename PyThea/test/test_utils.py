"""
Test the utilities
"""

import json
from datetime import datetime, timedelta

import astropy.units as u
import matplotlib.dates as mdates
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose
from sunpy.coordinates import get_horizons_coord
from sunpy.coordinates.frames import HeliographicStonyhurst
from sunpy.net import attrs as a
from sunpy.net import hek

from PyThea.config.selected_bodies import bodies_dict
from PyThea.config.selected_imagers import imager_dict
from PyThea.data.sample_data import (json_fitting_file_sample_data,
                                     stereo_sample_data)
from PyThea.utils import (get_hek_flare, load_fits, maps_process,
                          model_fittings, parameter_fit)


def test_get_hek_flare():
    """
    Test that the get_hek_flare returns a hek.hek.HEKTable and the results are formated correct.
    """
    def make_selectbox_list(fl_list):
        if len(fl_list) == 0:
            selectbox_list = ['No events returned', ]
        else:
            selectbox_list = []
            for flares in fl_list:
                fl_ = flares['fl_goescls']
                t_ = flares['event_peaktime'].strftime('%Y-%m-%dT%H:%M:%S')
                selectbox_list.append((f'FL{fl_}|{t_}'))
        return selectbox_list

    # This test requests the flare info for a date that returns a list of flares and compares against previous results
    day = datetime(2017, 9, 10, 0, 0, 0)
    timerange = a.Time(day, day + timedelta(days=1))
    flare_list = get_hek_flare(timerange)
    assert type(flare_list) == hek.hek.HEKTable
    selectbox_list = make_selectbox_list(flare_list)
    assert selectbox_list == ['FLM1.1|2017-09-09T23:53:00', 'FLC9.0|2017-09-10T03:09:00', 'FLC2.9|2017-09-10T09:20:00',
                              'FLC1.6|2017-09-10T14:23:00', 'FLC1.0|2017-09-10T15:26:00', 'FLX8.2|2017-09-10T16:06:00']

    # This test requests the flare info for a date that returns an empty list
    day = datetime(2030, 9, 10, 0, 0, 0)
    timerange = a.Time(day, day + timedelta(days=1))
    flare_list = get_hek_flare(timerange)
    selectbox_list = make_selectbox_list(flare_list)
    assert flare_list == []
    assert selectbox_list == ['No events returned']


def test_load_fitting_json_file():
    """
    Tests that the load_fitting_json_file returns the model_fittings class and storing it to the model_fittings
    will return the same result as the model_fittings.to_json().
    """
    def compare_json_objects(obj1, obj2, path='', exclude=['date_created', 'version']):
        if exclude:
            for key in exclude:
                obj1.pop(key, None)
                obj2.pop(key, None)

        for key in set(obj1.keys()) | set(obj2.keys()):
            new_path = f'{path}.{key}' if path else key

            if key not in obj1 or key not in obj2:
                print(f"Key '{new_path}' is not present in one of the objects.")
                return False
            elif obj1[key] != obj2[key]:
                if isinstance(obj1[key], dict) and isinstance(obj2[key], dict):
                    check = compare_json_objects(obj1[key], obj2[key], path=new_path)
                    if not check:
                        return False
                else:
                    print(f"Value mismatch at key '{new_path}': {obj1[key]} != {obj2[key]}")
                    return False

        return True

    json_fitting_file = json_fitting_file_sample_data.fetch('FLX1p0D20211028T153500MEllipsoid.json')
    with open(json_fitting_file, 'r') as file:
        json_content = file.read()
        fitting_from_json = json.loads(json_content)

    model_fittings_class = model_fittings.load_from_json(json_fitting_file)
    model_fittings_class_dict = json.loads(model_fittings_class.to_json())
    assert isinstance(model_fittings_class, model_fittings)
    assert model_fittings_class.event_selected == 'FLX1.0|2021-10-28T15:35:00'
    assert model_fittings_class.date_process == '2021-10-28T15:35:00.000000'
    assert model_fittings_class.geometrical_model == 'Ellipsoid'
    assert compare_json_objects(model_fittings_class_dict, fitting_from_json)


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


def test_maps_process():
    sta_cor2 = load_fits([stereo_sample_data.fetch('20120622_235400_d4c2a.fts')])
    stb_cor2 = load_fits([stereo_sample_data.fetch('20140221_235400_d4c2b.fts')])
    maps, _ = maps_process({'COR2A': [sta_cor2], 'COR2B': [stb_cor2]}, ['COR2A', 'COR2B'], 'Plain')
    for key in maps.keys():
        assert_quantity_allclose(maps[key][0].exposure_time, 1*u.second)
        assert_quantity_allclose([maps[key][0].dimensions[0], maps[key][0].dimensions[1]],
                                 [d/imager_dict[key]['process']['superpixel'] for d in imager_dict[key]['process']['dimensions']])
