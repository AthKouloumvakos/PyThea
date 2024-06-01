r"""
Plots fits files for a selected imager
--------------------------------------

This example shows how to download data for a selected imager, process the fits files, and plot an image.
"""

# %%
# Import Required Modules
from datetime import datetime, timedelta

import matplotlib.pyplot as plt
from sunpy.net import attrs as a

from PyThea.config import selected_imagers
from PyThea.utils import (download_fits, load_fits, make_figure,
                          single_imager_maps_process)

# %%
# Select the imager and specify the time range of the query.

imager = 'LC2'
date_process = datetime.strptime('2021-10-28T16:30:00', '%Y-%m-%dT%H:%M:%S')
time_range = [-1, 1]
timerange = a.Time(date_process + timedelta(hours=time_range[0]),
                   date_process + timedelta(hours=time_range[1]))

# %%
# Download the fits files from VSO using the ``download_fits`` utility.

files = download_fits(timerange, imager)
maps = load_fits(files)

# %%
# Process the downloaded fits files using the ``single_imager_maps_process`` utility. We select to process the maps into running differences images.

processed_images = single_imager_maps_process(maps,
                                              **selected_imagers.imager_dict[imager]['process'],
                                              image_mode='Running Diff.',
                                              diff_num=1)

# %%
# Plot the first image using the ``make_figure`` utility.

fig, ax = make_figure(processed_images[0], 'Running Diff.', clim=[-20, 20], clip_model=True)

plt.show()
