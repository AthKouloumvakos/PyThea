# v0.9.0 (8-Feb-2024)

## Features
- Adds sample data methods for testing and document builds
- Adds tests: test_get_horizons_coord, test_vso_search, test_hek_client, and test_parameter_fit_polynomial
- Adds calibration for AIA, LASCO, and STEREO EUVI

## Major Changes
- Changes the fitting time input from map.date to map.date_average when exist (see #24)

## Minor Changes
- Improves the following utilities: make_figure, maps_process, plot_fitting_model
- Improves the test_get_hek_flare and test_hek_client tests
- Changes how best_fit_x is defined in the parameters_fit utility
- Changes the pipeline of imaging data download and processing in the app

# v0.8.1 (13-Jan-2024)

## Features
- Implements PyThea test to cli

## Minor Changes
- Includes Python version 3.10

## Bug Fixes
- Fixes a runtime bug with streamlit 1.30.0

# v0.8.0 (11-Jan-2024)

## Major Changes
- Switches download of fits data to PyThea directory

## Minor Changes
- Adds in the README how to update PyThea and a note for the further development of PyThea
- Adds a CHANGELOG file
- Removes border lines from forms to improve the layout of components
- Replaces deprecated st.experimental_rerun with st.rerun
- Improves the datetick format of the kinematic plots

# v0.7.4 (09-Jul-2023)

## Minor Changes
- Changes the intersection calculation from Vedo to PyVista to reduce redering time
- Lowers the resolution of draw_limb to reduce redering time
- Removes unused code for duplicate fits file filtering

## Bug Fixes
- Fixes the css code that removes padding

# v0.7.3 (15-Feb-2023)

## Minor Changes
- Changes the relative imports to absolute
- Don't remove imports from __init__.py with autoflake

## Bug Fixes
- Fixes a missing package in custom mode in kinematics
- Fixes the github workflows when push
- Fixes a deprecation in distutils (use the packaging.version instead)

# v0.7.2 (28-Nov-2022)

## Minor Changes
- Improves sliders behavior
- Improves the main app. script readability and functionality

## Bug Fixes
- pip installation fails on Windows because of non-existing HOME (closed #17)

# v0.7.1 (18-Oct-2022)

## Minor Changes
- Changes the python version in GitHub workflow.
- A minor layout change on the main page.
- Improves the sliders' behaviour.
- Updates the pre-commit package versions.
- Removes the files related to Herocu (since it is deprecated).
- Reverse the coronagraphic images when crota was far from zero (closed #3).
- Includes the license file of the gcs_python package and improves the docstrings on gcs model (closed #16).

## Bug Fixes
- Fixes a bug with windows install (closed #15)

# v0.7.0 (01-Sep-2022)

## Minor Changes
- Changes the default minimum limit for the images download time range.
- Removes stqdm extension since this is now available from pip.
- General changes in the Readme page (Improves the installation part and also includes the information for the paper in frontiers)

## Bug Fixes
- Fixes a bug with streamlit>1.12 when the slider values is integers. Changed to float.
- Pin down the dependencies to parfive==1.5.1 (see #13 where data downloads failing)
- Fixes PyThea's runtime script since streamlit dropped the support for streamlit.cli
- Fixes an issue when no images found from VSO and resulted to error with sunpy>4.0.0
