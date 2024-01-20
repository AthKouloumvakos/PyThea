import os
from pathlib import Path

import pooch

database_dir = os.path.join(Path.home(), 'PyThea')
github_main_url = 'https://github.com/AthKouloumvakos/PyThea-sample-data'

aia_sample_data = pooch.create(
    path=os.path.join(database_dir, 'sample_data'),  # The cache folder
    base_url=f'{github_main_url}/raw/main/data/',  # The remote data url on Github
    registry={
        'aia_lev1_1600a_2012_06_06t04_07_29_12z_image_lev1_lowres.fits': 'sha256:396b9ef5cbc1cbe2f2af806758f449cfd36f352826b61295e6330859ac7f7652',
    },
)
