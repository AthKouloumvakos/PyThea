r"""
Construct a geometrical model
-----------------------------

This example shows how to construct a geometrical model using either user defined parameters or the results from a fitting file.
"""

# %%
# Import Required Modules
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames

from PyThea.data.sample_data import json_fitting_file_sample_data
from PyThea.geometrical_models import ellipsoid
from PyThea.utils import model_fittings

# %%
# User defined ellipsoid
# ^^^^^^^^^^^^^^^^^^^^^^
#
# You can create the a geometrical model manualy by providing an observation time, the ellipsoid center coordinates, and the geomertical parameters and use the ``ellipsoid`` class from PyThea's geometrical_models.

obstime = '2021-10-28'

center = SkyCoord(90*u.degree, 0*u.degree, 1*u.R_sun, obstime=obstime,
                  observer='earth', frame=frames.HeliographicStonyhurst)

model_shock = ellipsoid(center, 1*u.R_sun, 1*u.R_sun, 1*u.R_sun, 0 * u.degree)

# %%
# Ellipsoid from fitting file
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Instead of manualy constructing the geometrical model you can also load the model parameters from a fitting file and then construct the model using the ``ellipsoid`` class like previously shown.

# %%
# Import a sample JSON fitting file using ``json_fitting_file_sample_data.fetch()`` method. This sample data contains a series of fitted ellipsoids for a selected event.

json_fitting_file = json_fitting_file_sample_data.fetch('FLX1p0D20211028T153500MEllipsoid.json')

# %%
# Then use the ``model_fittings.load_from_json(json_fitting_file)`` to load the model parameters.

model_fittings_class = model_fittings.load_from_json(json_fitting_file)

# %%
# Select one of the fittings and create the ellipsoid using the observation time, the ellipsoid center coordinates, and the geomertical parameters from the fitting.

model_parameters = model_fittings_class.parameters.iloc[0]

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
