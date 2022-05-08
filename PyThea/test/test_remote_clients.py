"""
Test the remote clients
"""

import pytest
from sunpy.net import attrs, hek


@pytest.mark.remote_data
def test_hek_client():
    """
    Tests that the HEK client returns a HEKTable
    """
    startTime = '2011/08/09 07:23:56'
    endTime = '2011/08/09 12:40:29'
    eventType = 'FL'

    hekTime = attrs.Time(startTime, endTime)
    hekEvent = attrs.hek.EventType(eventType)

    h = hek.HEKClient()
    hek_query = h.search(hekTime, hekEvent)
    assert type(hek_query) == hek.hek.HEKTable

    # Test that the hek_query['event_peaktime'] is <class 'astropy.time.core.Time'>
