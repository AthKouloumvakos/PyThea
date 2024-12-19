# v1.0.0 (19-Nov-2024)

## Features
 - Adds imaging data from Solar Orbiter's SOLOHi.
 - Implements offline fits file loading from local database.

## Minor Changes
 - Updates Test figure hashes.
 - Final changes in README before the major release.

 ## Bug Fixes
 - Fixes a temporary bug with JSOC and AIA prep
 - Fixes potential bug with maps loading into map sequence.
 - Fixes a bug when fits database is selected but no files downloaded.
 - Fixes a bug with the maps clims.

# v0.14.0 (03-Oct-2024)

## Features
 - Adds imaging data from Solar Orbiter EUI and METIS.
 - Implements offline fits file loading from local database.

## Minor Changes
 - Simplifies the configuration of the database directory.
 - Improves Parker Spiral visualization.

# v0.13.0 (14-Jul-2024)

## Features
 - Adds imaging data from PSP/WISPR.
 - Implements new feature of Parker spiral and HEK events visualization on the images.
 - Adds a method to get directly the geometrical model from model_fittings.
 - Adds tests for extensions, WISPR imaging.
 - Adds and improves docstrings in utilities, modules, and geometrical models.
 - Adds GCS and kinematic plots in documentation
 - Adds function reference in the documentation

## Minor Changes
 - Changes some optional inputs for get_hek_flare, make_figure, plot_bodies, and plot_solar_reference_lines to be better utilised.

# v0.12.0 (09-Jun-2024)

## Features
 - Adds more imaging data from SDO/AIA.
 - Adds acceleration calculation in the kinematic plots.
 - Implements acceleration plotting in the app.
 - Includes fits file name in fitting files and dataframes.
 - Highlights a row in the fitting table when fitting exist in table.

## Major Changes
 - Decouples the loading of fits files from download_fits function.

## Minor Changes
 - Improves configuration dictionary for selected imagers.
 - Improves json warning message in the main app page.
 - Adds documentation page link to the main app page.
 - Updates README files.
 - Improves plot_fitting_model and adds figure test.
 - Improves application layout.

## Bug Fixes
 - Fixes a potential bug in filter_maps when no filtering is applied.
 - Fixes a bug with progress bar with stqdm.
 - Fixes a bug with image processing mode selection.

# v0.11.0 (19-May-2024)

## Features
 - Adds tests for ellipsoid location on AIA and STEREO COR images.
 - Adds test_load_fitting_json_file to check the fitting file load and save.
 - Improve PyThea tests on cli.

## Major Changes
- Change the inputs for download_fits from string dates to timerange.

# v0.10.0 (17-April-2024)

## Features
 - Option added to choose the step for Running Difference or the background image for the Base Difference image processing.
 - Fetures added in the app to plot limb or meridians from different observers and grid.
 - Feture added in the app for manual change of the images colorbar limits.
 - Feture added in the app to use median filter images processing.
 - Tests added: A test to verify the existence of Pythea's database directory

## Major Changes
- Updates sunpy dependencies to sunpy 5.1.2

## Minor Changes
- Improvements on the pipeline of image download and processing in the app.
- Improvements and simplifications on the figures production in the app.
- Values of the colormap limits in the app are stored and are reused when changing the imager.

## Bug Fixes
- Fixes a bug with test_database_dir_exists.
- Fixes a bug when corrupted fits files loaded and imager load skipped.
- Fixes a bug with AIA visuallization from missing date_average.
- Fixes a bug with supplementary imaging selection.

# v0.9.1 (5-March-2024)

## Features
 - A utility added to load the json fitting files

## Minor Changes
- Simplifies the code structure for loading the fitting files in the app

## Bug Fixes
- Fixes a bug when passing no maps in maputils
- Fixes a bug when corrupted fits files loaded and imager load skipped

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
