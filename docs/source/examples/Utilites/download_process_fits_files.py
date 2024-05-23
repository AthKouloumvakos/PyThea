r"""
Download and process fits files
-------------------------------

In this example, using the utility `download_fits` we download data for a selected imager and then process the fits files.
"""
# %%
# Import Required Modules
from datetime import datetime, timedelta

from sunpy.net import attrs as a

from PyThea.config import selected_imagers
from PyThea.utils import download_fits, single_imager_maps_process

# %%
# Retrieve the list of available imagers and their corresponding keys using the following,
for key in selected_imagers.imager_dict.keys():
    imager = selected_imagers.imager_dict[key]
    detector_or_wavelenght = imager['detector'] if 'detector' in imager else imager['wavelength']
    print(f'{imager["source"]}/{imager["instrument"]}-{detector_or_wavelenght}:  {key}')


# %%
# First, select the imager and specify the time range of the query. For this example we will download data from SOHO LASCO coronagraph.

imager = 'LC2'
date_process = datetime.strptime('2021-10-28T16:30:00', '%Y-%m-%dT%H:%M:%S')
time_range = [-1, 1]
timerange = a.Time(date_process + timedelta(hours=time_range[0]),
                   date_process + timedelta(hours=time_range[1]))

# %%
# Download the fits files from VSO using the ``download_fits`` utility.

maps = download_fits(timerange, imager)

# %%
# Process the downloaded fits files using the ``single_imager_maps_process`` utility. This function filters, prepares, and process the maps.
#
# The default options for the maps processing, for each imager, can be found in,

print(selected_imagers.imager_dict[imager]['process'])

# For LASCO-C2 select only the images with dimensions 1024x1024 and only total brightness images. Then the filtered images are prepared. Depending on the instrument, this includes pointing corrections, calibrations, observer location corrections, exposure time normalization, and others.

# %%
# At the last step the images are resampled using SynPy's ``superpixel`` method and the final maps are processed  into running/base difference or plain image sequence maps.

processed_images = single_imager_maps_process(maps,
                                              **selected_imagers.imager_dict[imager]['process'],
                                              image_mode='Running Diff.',
                                              diff_num=1)

# These images can now be used in the fitting process or just plot them (see example).
