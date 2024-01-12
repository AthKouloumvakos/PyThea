"""
Test the utilities
"""

from datetime import datetime

from sunpy.net import hek

from PyThea.utils import get_hek_flare


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
