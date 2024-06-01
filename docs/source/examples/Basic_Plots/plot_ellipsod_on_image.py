r"""
Plot geometrical model to an image
----------------------------------

In this example, we use PyThea's utilities to plot the ellipsoid model to a coronagraphic image from LASCO-C2. This example is a combination of two examples that show how to load, process, and plot an image and how to load a fitting file and construct the geometrical model.
"""

# %%
# Import Required Modules
from datetime import datetime, timedelta

import astropy.units as u
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
from sunpy.net import attrs as a

from PyThea.config import selected_imagers
from PyThea.data.sample_data import json_fitting_file_sample_data
from PyThea.geometrical_models import ellipsoid
from PyThea.utils import (download_fits, load_fits, make_figure,
                          model_fittings, single_imager_maps_process)

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
# Import the sample JSON fitting file using ``json_fitting_file_sample_data.fetch()`` method. This sample data contains a series of fitted ellipsoids for a selected event.

json_fitting_file = json_fitting_file_sample_data.fetch('FLX1p0D20211028T153500MEllipsoid.json')

# %%
# The ``model_fittings.load_from_json(json_fitting_file)`` loads the model fittings from the JSON file. The result is stored in the variable model_fittings_class.

model_fittings_class = model_fittings.load_from_json(json_fitting_file)

# %%
# Select the second fitting that corresponds to the first image on the list of maps sequence.

model_parameters = model_fittings_class.parameters.iloc[1]

# %%
# Construct the ellipsoid model. Define first the center of the ellipsoid and then use the ``ellipsoid`` class to make the model
obstime = model_parameters.name

center = SkyCoord(model_parameters['hgln']*u.degree,
                  model_parameters['hglt']*u.degree,
                  model_parameters['rcenter']*u.R_sun,
                  obstime=obstime, observer='earth',
                  frame=frames.HeliographicStonyhurst)

model_shock = ellipsoid(center,
                        model_parameters['radaxis']*u.R_sun,
                        model_parameters['orthoaxis1']*u.R_sun,
                        model_parameters['orthoaxis2']*u.R_sun,
                        model_parameters['tilt']*u.degree)

# %%
# Plot the first image using the ``make_figure`` utility and over-plot at the same axis the ellipsoid model.

fig, ax = make_figure(processed_images[0], 'Running Diff.', clim=[-20, 20], clip_model=True)
model_shock.plot(ax, mode='Full')

plt.show()
