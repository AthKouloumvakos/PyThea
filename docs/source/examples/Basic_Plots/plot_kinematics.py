r"""
Plot model kinematics
---------------------

In this example, we use PyThea's utilities to plot the kinematics of the geometrical model from a fitting file.
"""

# %%
# Import Required Modules
from datetime import datetime

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from IPython.display import display

from PyThea.data.sample_data import json_fitting_file_sample_data
from PyThea.utils import model_fittings, plot_fitting_model

# %%
# Import a sample JSON fitting file using ``json_fitting_file_sample_data.fetch()`` method. This sample data contains a series of fitted ellipsoids for a selected event.
# Then use the ``model_fittings.load_from_json(json_fitting_file)`` to load the model parameters.

json_fitting_file = json_fitting_file_sample_data.fetch('FLX1p0D20211028T153500MEllipsoid.json')
model_fittings_class = model_fittings.load_from_json(json_fitting_file)

# %%
# Print the information about the selected event and geometrical model using ``print()`` statements and the parameters of the model fittings using ``display(model_fittings_class.parameters)``

print(f'Selected Event: {model_fittings_class.event_selected}')
print(f'Geometrical Model: {model_fittings_class.geometrical_model}')
display(model_fittings_class.parameters)


# %%
# Use the stored curve fitting parameters
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# This part uses the curve fitting method parameters that are stored in the json fitting file. The curve fitting parameters are selected by the user during the fitting process in the application.
# View the curve fitting method parameters. These consist of the method (``type``) of the curve fitting ('polynomial' in this example) and the ``order`` (2nd order in this example).

print(f'Fitting Parameters: {model_fittings_class.kinematics["fit_method"]}')

# %%
# Use the utility function ``plot_fitting_model`` to plot the kinematic figures. Select for ``plt_type`` the ``HeightT`` option to plot the height-time profile of the geometrical fitting and the ``SpeedT`` to plot the speed-time profile.

fig, axis = plot_fitting_model(model_fittings_class,
                               fit_args=model_fittings_class.kinematics['fit_method'],
                               plt_type='HeightT')
axis.set_xlim([datetime(2021, 10, 28, 15, 40, 0), datetime(2021, 10, 28, 16, 40, 0)])
axis.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=range(60), interval=1))
plt.show()

fig, axis = plot_fitting_model(model_fittings_class,
                               fit_args=model_fittings_class.kinematics['fit_method'],
                               plt_type='SpeedT')
axis.set_xlim([datetime(2021, 10, 28, 15, 40, 0), datetime(2021, 10, 28, 16, 40, 0)])
axis.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=range(60), interval=1))
axis.set_ylim([1000, 2200])
plt.show()

# %%
# Manual fitting parameters
# ^^^^^^^^^^^^^^^^^^^^^^^^^
#
# To use other curve fitting parameters parse different 'type' and 'order'. In the example bellow we plot the height time profile using a linear fit.

fit_method = {'type': 'polynomial', 'order': 1}

fig, axis = plot_fitting_model(model_fittings_class,
                               fit_args=fit_method,
                               plt_type='HeightT')
axis.set_xlim([datetime(2021, 10, 28, 15, 40, 0), datetime(2021, 10, 28, 16, 40, 0)])
axis.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=range(60), interval=1))
plt.show()
