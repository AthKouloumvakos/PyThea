import os
from pathlib import Path

from setuptools import setup

# create home directory
if not os.path.isdir(os.path.join(Path.home(), 'PyThea')):
    os.mkdir(os.path.join(Path.home(), 'PyThea'))
    os.mkdir(os.path.join(Path.home(), 'PyThea', 'data'))
    os.mkdir(os.path.join(Path.home(), 'PyThea', 'data_manager'))

setup()
