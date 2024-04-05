"""
A place for the more generic tests that do not fit in the other categories

Note
----
Before you start the test the PyThea pkg should be in the python path
export PYTHONPATH="${PYTHONPATH}:{top_level_dir_that_pythea_lives}/PyThea"
"""

import os
from pathlib import Path


def test_database_dir_exists():
    """
    Tests that Pythea's database directory exists.
    """
    database_dir = os.path.join(Path.home(), 'PyThea')
    assert os.path.exists(database_dir), f"PyThea's database directory '{database_dir}' does not exist."
