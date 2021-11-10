"""
Test the utilities
"""

from sunpy.net import hek

from PyThea.utils import get_hek_flare


# ------------------
# Test utilities
def test_get_hek_flare_HEKTable():
    """
    Test that the get_hek_flare returns a hek.hek.HEKTable
    """
    from datetime import datetime
    day = datetime(2017, 9, 10)
    _, flare_list_ = get_hek_flare(day)

    assert type(flare_list_) == hek.hek.HEKTable
