r"""
Load a fiting file
------------------

This script demonstrates how to use functionalities from the PyThea library to work with fitting data stored in JSON format.
"""

# %%
# Import Required Modules
from IPython.display import display

from PyThea.data.sample_data import json_fitting_file_sample_data
from PyThea.utils import model_fittings

# %%
# The script imports a sample JSON fitting file using ``json_fitting_file_sample_data.fetch()`` method. This sample data contains a series of fitted ellipsoids for a selected event.

json_fitting_file = json_fitting_file_sample_data.fetch('FLX1p0D20211028T153500MEllipsoid.json')

# %%
# The ``model_fittings.load_from_json(json_fitting_file)`` loads the model fittings from the JSON file. The result is stored in the variable model_fittings_class.

model_fittings_class = model_fittings.load_from_json(json_fitting_file)

# %%
# Print the information about the selected event and geometrical model using ``print()`` statements and the parameters of the model fittings using ``display(model_fittings_class.parameters)``

print(f'Selected Event: {model_fittings_class.event_selected}')
print(f'Geometrical Model: {model_fittings_class.geometrical_model}')
display(model_fittings_class.parameters)


# %%
