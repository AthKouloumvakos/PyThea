r"""
Plot geometrical models on empty Map
------------------------------------

This example shows how to plot an ellipsoid on a blank map.
"""
import astropy.units as u
# %%
# Import Required Modules
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
from astropy.coordinates import SkyCoord
from IPython.display import display
from sunpy.coordinates import frames

from PyThea.data.sample_data import json_fitting_file_sample_data
from PyThea.geometrical_models import ellipsoid
from PyThea.utils import model_fittings

# %%
# Create a blank map using with an array of zeros and construct a header to pass to Map.
# This is based on the `SunPy's example <https://docs.sunpy.org/en/stable/generated/gallery/plotting/plotting_blank_map.html>`_


def create_blank_map(obstime):
    data = np.full((10, 10), np.nan)

    # Define a reference coordinate and create a header using sunpy.map.make_fitswcs_header
    skycoord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime=obstime,
                        observer='earth', frame=frames.Helioprojective)

    # Scale set to the following for solar limb to be in the field of view
    header = sunpy.map.make_fitswcs_header(data, skycoord, scale=[220, 220]*u.arcsec/u.pixel)

    # Use sunpy.map.Map to create the blank map
    blank_map = sunpy.map.Map(data, header)

    return blank_map


# %%
# User defined ellipsoid
# ^^^^^^^^^^^^^^^^^^^^^^
#
# In this part of the example we construct an ellipsoid model using ``ellipsoid`` class from PyThea's geometrical_models.
#

# %%
# First create an ellipsoid providing an observation time, the ellipsoid center coordinates, and the geomertical parameters.

obstime = '2021-10-28'
center = SkyCoord(90*u.degree, 0*u.degree, 1*u.R_sun, obstime=obstime,
                  observer='earth', frame=frames.HeliographicStonyhurst)
model_shock = ellipsoid(center, 1*u.R_sun, 1*u.R_sun, 1*u.R_sun, 0 * u.degree)

# %%
# Create a blank map with WCS defined by a helioprojective frame as observed from Earth at the observation time.
# We use this map to overplot the ellipsoid.

blank_map = create_blank_map(obstime)

# %%
# Create a figure and plot the map and the ellipsoid.

fig = plt.figure()
ax = fig.add_subplot(projection=blank_map)
blank_map.plot(axes=ax)
blank_map.draw_limb(axes=ax, color='k')
blank_map.draw_grid(axes=ax, color='k')
model_shock.plot(ax, mode='Full')

ax.set_title('Plotting ellipsoid on a map')

plt.show()

# %%
# Ellipsoid from fitting file
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# In this part we will construct an ellipsoid model using the parameters from a fitting file.

# %%
# Import a sample JSON fitting file using ``json_fitting_file_sample_data.fetch()`` method. This sample data contains a series of fitted ellipsoids for a selected event.
# Then use the ``model_fittings.load_from_json(json_fitting_file)`` to load the model parameters.

json_fitting_file = json_fitting_file_sample_data.fetch('FLX1p0D20211028T153500MEllipsoid.json')
model_fittings_class = model_fittings.load_from_json(json_fitting_file)

# %%
# Select one of the fittings and display the parameters of the geometrical model.

model_parameters = model_fittings_class.parameters.iloc[0]
display(model_parameters)

# %%
# Create the ellipsoid using the observation time, the ellipsoid center coordinates, and the geomertical parameters from the fitting.

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
# Create a figure and plot the map and the ellipsoid.

fig = plt.figure()
ax = fig.add_subplot(projection=blank_map)
blank_map.plot(axes=ax)
blank_map.draw_limb(axes=ax, color='k')
blank_map.draw_grid(axes=ax, color='k')
model_shock.plot(ax, mode='Full')

ax.set_title('Plotting ellipsoid on a map')

plt.show()
