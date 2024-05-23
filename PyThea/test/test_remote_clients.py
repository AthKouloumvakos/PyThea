"""
Test the remote clients
"""

import numpy as np
import pytest
from astropy.time import Time
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net import hek

from PyThea.config.selected_imagers import imager_dict


def time_query_hek():
    for query in [a.Time('2011/08/09 07:23:56', '2011/08/09 12:40:29'),
                  a.Time('2012/09/10 08:24:57', '2012/09/10 13:41:30'),
                  a.Time('2013/10/11 09:25:58', '2013/10/11 14:42:31'), ]:
        yield query


@pytest.mark.remote_data
@pytest.mark.parametrize('time_query', time_query_hek())
def test_hek_client(time_query):
    """
    Tests that the HEK client returns a HEKTable and the hek_query results for EventType('FL') have consistent format
    """
    h = hek.HEKClient()
    hek_query = h.search(time_query, a.hek.EventType('FL'))
    assert type(hek_query) == hek.hek.HEKTable

    assert isinstance(hek_query['event_starttime'], Time)
    assert isinstance(hek_query['event_peaktime'], Time)
    assert isinstance(hek_query['event_endtime'], Time)
    assert all(isinstance(h, np.str_) for h in hek_query['fl_goescls'])


@pytest.mark.remote_data
def test_hek_client_flare_data():
    """
    Tests that the HEK client returns always the same results (flare class, ...)
    """

    startTime = '2011/08/09 07:23:56'
    endTime = '2011/08/09 12:40:29'
    eventType = 'FL'

    hekTime = a.Time(startTime, endTime)
    hekEvent = a.hek.EventType(eventType)

    h = hek.HEKClient()

    hek_query = h.search(hekTime, hekEvent, a.hek.FL.GOESCls > 'B1.0', a.hek.OBS.Observatory == 'GOES')
    assert all(hek_query['fl_goescls'] == ['C1.4', 'X6.9'])
    assert all(hek_query['event_starttime'] == Time(['2011-08-09 07:19:00.000', '2011-08-09 07:48:00.000']))
    assert all(hek_query['event_peaktime'] == Time(['2011-08-09 07:23:00.000', '2011-08-09 08:05:00.000']))
    assert all(hek_query['event_endtime'] == Time(['2011-08-09 07:27:00.000', '2011-08-09 08:08:00.000']))


def imager_query():
    for imager in imager_dict.keys():
        yield imager


@pytest.mark.remote_data
@pytest.mark.parametrize('imager', imager_query())
def test_vso_search(imager):
    """
    Tests that the Fido.search returns always the same number of files from VSO.
    """
    file_num = {'AIA-193': [60, 61], 'AIA-211': [60, 61],
                'LC2': [5], 'LC3': [5],
                'COR2A': [7], 'COR2B': [7],
                'EUVIA': [11], 'EUVIB': [11],
                'COR1A': [18, 34], 'COR1B': [17, 66],
                'HI1A': [1, 2, 3], 'HI1B': [1, 2],
                'HI2A': [1, 2], 'HI2B': [1], }

    for query in [a.Time('2011/01/01 00:00:00', '2011/01/01 01:00:00'),
                  a.Time('2012/01/01 00:00:00', '2012/01/01 01:00:00'),
                  a.Time('2013/01/01 00:00:00', '2013/01/01 01:00:00'), ]:
        args = imager_dict[imager]['fido']
        result = Fido.search(query, *args)

        assert len(result) == 1
        assert 'vso' in result.keys()

        assert result.file_num in file_num[imager]
