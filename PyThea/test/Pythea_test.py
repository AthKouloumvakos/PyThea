"""
A place for the more generic tests that do not fit in the other categories

Note
----
Before you start the test the PyThea pkg should be in the python path
export PYTHONPATH="${PYTHONPATH}:{top_level_dir_that_pythea_lives}/PyThea"
"""

import os

from PyThea.config import database_dir_default


def test_database_dir_exists():
    """
    Tests that Pythea's database directory exists.
    """
    # Skip the test if running in GitHub Actions
    if os.environ.get('CI') == 'true':
        print("Skipping test as it's running in GitHub Actions.")
        return

    assert os.path.exists(database_dir_default), f"PyThea's database directory '{database_dir_default}' does not exist."
