import os
from pathlib import Path

import pooch

database_dir = os.path.join(Path.home(), 'PyThea')
github_main_url = 'https://github.com/AthKouloumvakos/PyThea-sample-data'

aia_sample_data = pooch.create(
    path=os.path.join(database_dir, 'sample_data'),  # The cache folder
    base_url=f'{github_main_url}/raw/main/data/',  # The remote data url on Github
    registry={
        'aia_lev1_1600a_2012_06_06t04_07_29_12z_image_lev1.fits': 'sha256:7e88d8a277ad3dd111c1bcff3c3936d6d5ef5186e05b4bfa7749ee43c05577d5',
        'aia_lev1_1600a_2016_05_09t11_35_51_12z_image_lev1.fits': 'sha256:1c3bb938e9bfb3178cecead1319a92fa661eb70469715a441a849a98039f0b61'
    },
)

json_fitting_file_sample_data = pooch.create(
    path=os.path.join(database_dir, 'sample_data'),  # The cache folder
    base_url=f'{github_main_url}/raw/main/data/',  # The remote data url on Github
    registry={
        'FLX1p0D20211028T153500MEllipsoid.json': 'sha256:700687ed982317400ec409390571747694237de9688a2816baa7b1770b75c9b1',
    },
)
