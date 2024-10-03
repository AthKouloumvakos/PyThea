import os

import pooch

from PyThea.config import database_dir_default

github_main_url = 'https://github.com/AthKouloumvakos/PyThea-sample-data'

aia_sample_data = pooch.create(
    path=os.path.join(database_dir_default, 'sample_data'),  # The cache folder
    base_url=f'{github_main_url}/raw/main/data/',  # The remote data url on Github
    registry={
        'aia_lev1_1600a_2012_06_06t04_07_29_12z_image_lev1.fits': 'sha256:7e88d8a277ad3dd111c1bcff3c3936d6d5ef5186e05b4bfa7749ee43c05577d5',
        'aia_lev1_1600a_2016_05_09t11_35_51_12z_image_lev1.fits': 'sha256:1c3bb938e9bfb3178cecead1319a92fa661eb70469715a441a849a98039f0b61'
    },
)

stereo_sample_data = pooch.create(
    path=os.path.join(database_dir_default, 'sample_data'),  # The cache folder
    base_url=f'{github_main_url}/raw/main/data/',  # The remote data url on Github
    registry={
        '20120622_235400_d4c2a.fts': 'sha256:940680fe8bea754c9859170585d8299b9b049b631155435a31ecacf5b702d9cf',
        '20140221_235400_d4c2b.fts': 'sha256:13638dedeb2e510c9e2ef11d46b372acb16f9b9645ca1cabf4d0ef8444353adc',
        '20070503_102018_s4c1a.fts': 'sha256:55cf6e3741833a82e41b46e2345bf85a0e5580a2832900eccc41f235d36fc31e',
        '20070503_102018_s4c1b.fts': 'sha256:ca8f4aabaeb7e4bc6f4e5d4819057ddfe3f9f621716224d11acc80d701e3beb8'
    },
)

wispr_sample_data = pooch.create(
    path=os.path.join(database_dir_default, 'sample_data'),  # The cache folder
    base_url=f'{github_main_url}/raw/main/data/',  # The remote data url on Github
    registry={
        'psp_l3_wispr_20210808t103707_v1_1211.fits': 'sha256:85035693c4677e8fa8a3e61b8051d6f0a46e96806b5b0cbdb3efa77652308bda',
        'psp_l3_wispr_20210808t104010_v1_2222.fits': 'sha256:5e55e17f1aed7dc96f88b79f092947c58e6ee0def1b72822dbf1a696d8df7e6f'
    },
)

json_fitting_file_sample_data = pooch.create(
    path=os.path.join(database_dir_default, 'sample_data'),  # The cache folder
    base_url=f'{github_main_url}/raw/main/data/',  # The remote data url on Github
    registry={
        'FLX1p0D20211028T153500MEllipsoid.json': 'sha256:34fced8530875110117ba27a73541d3acd0c90bc8fca99367c1ec4cdfc05ce86',
    },
)
