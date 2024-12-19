"""
    PyThea: A software package to reconstruct the 3D structure of CMEs and
    shock waves using multi-viewpoint remote-sensing observations.
    Copyright (C) 2021  Athanasios Kouloumvakos

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


import datetime
import io
import json
import os
import re
import warnings
from copy import copy

import astropy.units as u
import matplotlib.colors as colors
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numexpr
import numpy as np
import pandas as pd
import seaborn as sns
import sunpy.map
from astropy.coordinates import Distance, SkyCoord, SphericalRepresentation
from astropy.visualization.wcsaxes.wcsapi import wcsapi_to_celestial_frame
from scipy.interpolate import UnivariateSpline
from scipy.ndimage import median_filter
from scipy.optimize import curve_fit
from sunpy.coordinates import frames, get_horizons_coord
from sunpy.map.maputils import contains_coordinate
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import parse_time
from sunpy.visualization import drawing

from PyThea.config import database_dir_default
from PyThea.config.selected_bodies import bodies_dict as bodies_dict_default
from PyThea.config.selected_imagers import imager_dict as imager_dict_default
from PyThea.geometrical_models import ellipsoid, gcs, spheroid
from PyThea.sunpy_dev.map.maputils import (filter_maps,
                                           maps_sequence_processing,
                                           prepare_maps)
from PyThea.version import version


def get_hek_flare(timerange, thresshold='B1.0'):
    """
    Returns the flare list for a given day from HEK.

    This function will download the flare data from HEK and return a list of flares.

    Parameters
    ----------
    timerange: sunpy.net.attrs.Time
        The time range of the data search.
    thresshold: str
        The flare class above which to do the search (default 'B1.0').

    thresshold: str, optional
        The flare class above which the search will return records (default: 'B1.0').

    Returns
    -------
    parfive.Results
    """
    flare_list = Fido.search(timerange,
                             a.hek.EventType('FL'),
                             a.hek.FL.GOESCls > thresshold,
                             a.hek.OBS.Observatory == 'GOES')
    if len(flare_list['hek']) == 0:
        flare_list_ = []
    else:
        flare_list_ = flare_list['hek']['event_starttime', 'event_peaktime',
                                        'event_endtime', 'fl_goescls',
                                        'fl_peakflux', 'hgs_x', 'hgs_y', 'ar_noaanum']

    return flare_list_


def make_figure(map, cmap='Greys_r', clim=[-20, 20], clip_model=True, **kwargs):
    """
    Creates the main imager figure and returns the figure and axis handles.

    Parameters
    ----------
    map : sunpy.map.Map
        The SunPy map to be plotted.
    cmap : str, optional
        The colormap to be used for the plot (default is 'Greys_r').
    clim : list, optional
        The color limits for the plot as [min, max] (default is [-20, 20]).
    clip_model : bool, optional
        Whether to clip the axis limits to the map's data shape (default is True).
    **kwargs : dict
        Additional keyword arguments for customization, including:
            - fig: The figure object (default creates a new figure).
            - axis: The axis object (default creates a new axis with the map's projection).
            - median_filter: The size of the median filter to apply (default is 1, which means no filtering).

    Returns
    -------
    tuple
        The figure and axis handles.

    Notes
    -----
    This function customizes the plot by optionally applying a median filter,
    setting color limits, drawing the solar limb, and adjusting axis properties.
    """
    fig = kwargs.get('fig', plt.figure())
    axis = kwargs.get('axis', plt.subplot(projection=map))

    median_filter_value = kwargs.get('median_filter', 1)
    if median_filter_value != 1:
        map = sunpy.map.Map(median_filter(map.data, size=int(median_filter_value)), map.meta)

    if map.instrument in ['WISPR', 'Metis'] or map.instrument.startswith('SoloHI'):
        clim = [-10**-clim[0], 10**-clim[1]]

    if cmap == 'default':
        # TODO: For plain images or when EUVIA-B are used, this does not work very well.
        map.plot(norm=colors.Normalize(vmin=clim[0], vmax=clim[1]))
    else:
        map.plot(cmap=cmap, norm=colors.Normalize(vmin=clim[0], vmax=clim[1]))

    map.draw_limb(resolution=90)

    yax = axis.coords[1]
    yax.set_ticklabel(rotation=90)

    if clip_model:
        axis.set_xlim([0, map.data.shape[0]])
        axis.set_ylim([0, map.data.shape[1]])

    cref = map.pixel_to_world(0*u.pix, 0*u.pix)
    if cref.Tx > 0 and (map.instrument != 'WISPR'):
        axis.invert_xaxis()
    if cref.Ty > 0:
        axis.invert_yaxis()

    if map.instrument == 'SoloHI':
        title = 'SoloHI' + f' Tile-{map.detector}' ' $T_{AGV}:$' + parse_time(map.date_average).strftime('%Y-%m-%d %H:%M:%S')
    elif map.instrument == 'Metis':
        title = re.sub(r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}',
                       ' $T_{AGV}:$' + parse_time(map.date_average).strftime('%Y-%m-%d %H:%M:%S'),
                       map.latex_name.replace('VLD', 'METIS-VDL'))
    else:
        title = re.sub(r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}',
                       ' $T_{AGV}:$' + parse_time(map.date_average).strftime('%Y-%m-%d %H:%M:%S'),
                       map.latex_name)
    print(map.instrument)
    print(map.latex_name)
    axis.set_title(title,
                   fontsize=10, pad=8)

    return fig, axis


def plot_bodies(axis, bodies_list, smap, bodies_dict=None):
    """
    Plots the positions of pre-configured bodies (e.g., Earth, STA, Venus) on the given axis.

    Parameters
    ----------
    axis : matplotlib.axes._subplots.AxesSubplot
        The matplotlib axis to plot on.
    bodies_list : list
        A list of bodies to be plotted.
    smap : sunpy.map.Map
        The SunPy map object containing the data.
    bodies_dict : dict
        A dictionary where keys are body names and values are tuples containing
        the body ID for Horizons coordinates and the color for plotting.
        The default bodies_dict is from PyThea.config.selected_bodies import bodies_dict

    Returns
    -------
    None
    """
    if bodies_dict is None:
        bodies_dict = bodies_dict_default

    for body in bodies_list:
        if body not in bodies_dict:
            RuntimeWarning(f'Body information was missing in bodies_dict for {body}')
            continue

        body_coo = get_horizons_coord(bodies_dict[body][0], smap.date_average)
        if contains_coordinate(smap, body_coo):
            axis.plot_coord(body_coo, 'o', color=bodies_dict[body][1],
                            fillstyle='none', markersize=6, label=body)


def plot_solar_reference_lines(axis, bodies_list, smap, mode='Limb from Obs.', bodies_dict=None):
    """
    Plots solar reference lines (e.g., solar limb, central meridians, equator) on a WCS axis.

    Parameters
    ----------
    axis : matplotlib.axes.Axes
        The axis to plot on.
    bodies_list : list
        List of bodies to plot reference lines for (e.g., 'Earth', 'Venus').
    smap : sunpy.map.Map
        The SunPy map for which the reference lines are plotted.
    mode : str, optional
        The mode of reference lines to plot. Options are 'Limb from Obs.', 'Central Meridian from Obs.',
        'Carr. Prime Meridian+Solar Equator', 'Stonyhurst Grid', 'Carrington Grid' (default is 'Limb from Obs.').
    bodies_dict : dict, optional
        Dictionary containing body names as keys and tuples of (body ID, color) as values.
        Default is `bodies_dict_default`.

    Notes
    -----
    The function uses the Horizons system to get body coordinates and supports various modes for plotting
    reference lines including limb, central meridian, solar equator, and grids in different coordinate systems.
    """
    def draw_central_meridian(axis, body_coo, date, color):
        body_coo = body_coo.transform_to(frames.HeliographicStonyhurst)
        constant_lon = SkyCoord(body_coo.lon, np.linspace(-90, 90, 90) * u.deg,
                                frame=frames.HeliographicStonyhurst(obstime=date))
        v, _ = drawing._plot_vertices(constant_lon, axis, wcsapi_to_celestial_frame(axis.wcs), None,
                                      color=color, close_path=False, linewidth=1)

    if bodies_dict is None:
        bodies_dict = bodies_dict_default

    for body in bodies_list:
        if body not in bodies_dict:
            RuntimeWarning(f'Body information was missing in bodies_dict for {body}')
            continue

        try:
            body_coo = get_horizons_coord(bodies_dict[body][0], smap.date_average)
        except ValueError as ve:
            print(f'Error processing {body}: {ve}')
            continue

        if mode in ['Limb from Obs.', 'limb']:
            drawing.limb(axis, body_coo, color=bodies_dict[body][1],
                         linewidth=1, resolution=90, label=f'Limb seen from {body}')
        if mode in ['Central Meridian from Obs.', 'meridian']:
            draw_central_meridian(axis, body_coo, smap.date_average, bodies_dict[body][1])

    if mode in ['Carr. Prime Meridian+Solar Equator', 'carr-me']:
        drawing.equator(axis, linewidth=1, resolution=90)
        drawing.prime_meridian(axis, linewidth=1, resolution=90)
    elif mode in ['Stonyhurst Grid', 'stony-grid']:
        earth = get_horizons_coord(3, smap.date_average)
        draw_central_meridian(axis, earth, smap.date_average, 'white')
        smap.draw_grid(linewidth=1, color='red', system='stonyhurst', alpha=0.8)
    elif mode in ['Carrington Grid', 'carr-grid']:
        drawing.prime_meridian(axis, linewidth=1, resolution=90)
        smap.draw_grid(linewidth=1, color='red', system='carrington', alpha=0.8)


def download_fits(timerange, imager, imager_dict=None, database_dir=None):
    """
    Downloads the imaging data (FITS files) from the Virtual Solar Observatory (VSO).

    Parameters
    ----------
    timerange : sunpy.net.attrs.Time
        The time range of the data search.
    imager : str
        The name of the imager (e.g., 'AIA', 'HMI').
    imager_dict : dict, optional
        A dictionary containing imager properties where keys are imager names and
        values are dictionaries with properties including 'fido', 'source', 'instrument',
        'detector' (optional), and 'wavelength' (optional). The default imager_dict is from
        PyThea.config.selected_imagers.import imager_dict.
    database_dir : str, optional
        The base directory where the data will be saved. The default is database_dir_default from PyThea.config.

    Returns
    -------
    list
        A list of downloaded FITS file paths.
    """

    if imager_dict is None:
        imager_dict = imager_dict_default[imager]
    else:
        if imager in imager_dict:
            imager_dict = imager_dict[imager]
        else:
            if 'fido' not in imager_dict:
                raise ValueError('Provided imager_dict does not contain fido information')

    if database_dir is None:
        database_dir = database_dir_default
    else:
        if not os.path.isdir(database_dir):
            raise ValueError(f"The provided database_dir '{database_dir}' is not a valid directory")

    result = Fido.search(timerange, *imager_dict['fido'])
    print(result)

    if result:
        if 'detector' in imager_dict:
            sub_path = imager_dict['detector']
        elif 'wavelength' in imager_dict:
            sub_path = imager_dict['wavelength']
        else:
            sub_path = ''

        if sub_path:
            path_str = f'{database_dir}/data/{imager_dict["source"]}/{imager_dict["instrument"]}/{sub_path}' + '/{file}'
        else:
            path_str = f'{database_dir}/data/{imager_dict["source"]}/{imager_dict["instrument"]}' + '/{file}'
        downloaded_files = Fido.fetch(result, path=path_str)
    else:
        downloaded_files = []

    return downloaded_files


def load_fits(files):
    """
    Loads imaging data (FITS files) from a list of files and returns them as SunPy map objects.

    Parameters
    ----------
    files : list
        A list of file paths to the FITS files.

    Returns
    -------
    sunpy.map.Map or list of sunpy.map.Map
        A SunPy map sequence if multiple maps are loaded, or a single SunPy map if only one map is loaded.
    """
    maps_ = []

    if files:
        for file_path in files:
            try:
                map_ = sunpy.map.Map(file_path)
                if isinstance(map_, list):  # Added because of METIS
                    map_ = map_[0]
                map_.meta['fits_file'] = os.path.basename(file_path)
                maps_.append(map_)
            except RuntimeError as err:
                print('Handling RuntimeError error:', err)
            except OSError as err:
                print('Handling OSError error:', err)
                os.remove(file_path)
                print(f"File '{file_path}' has been removed.")

        if maps_:
            maps_ = sunpy.map.Map(maps_, sequence=True) if len(maps_) > 1 else maps_[0]

    return maps_


def maps_process(maps_dict_in, imagers_list_in, image_mode, **kwargs):
    """
    Parameters
    ----------
    maps_dict_in : dict
        Dictionary of input maps, where keys are imager names and values are lists of maps.
    imagers_list_in : list
        List of imagers to process.
    image_mode : str
        The mode in which to process the images.
    **kwargs : dict, optional
        Additional processing settings for specific imagers.

    Returns
    -------
    dict
        Dictionary of processed maps, where keys are imager names and values are lists of processed maps.
    list
        List of imagers that were successfully processed and have output maps.

    Notes
    -----
    The `maps_dict_in` is the `session_state.map_` when used from the application.
    """
    maps_dict_out = {}
    imagers_list_out = []

    for imager in imagers_list_in:
        if imager in maps_dict_in and maps_dict_in[imager]:
            if not kwargs:
                extras = imager_dict_default[imager]['process']
            else:
                if imager in kwargs:
                    extras = kwargs[imager]
                else:
                    print(f'Warning [maps_process]: No extras provided for {imager}.')
            maps_dict_out[imager] = single_imager_maps_process(maps_dict_in[imager],
                                                               image_mode=image_mode,
                                                               **extras)
            if maps_dict_out[imager]:
                imagers_list_out.append(imager)

    return maps_dict_out, imagers_list_out


def single_imager_maps_process(map_list, skip=None, **kwargs):
    """
    Processes the images for a single imager and returns the final maps after applying specified operations.

    Parameters
    ----------
    map_list : list
        List of maps to be processed.
    skip : str or None, optional
        Operations to skip during processing (default is None).
    **kwargs : dict
        Additional keyword arguments for processing operations.

    Returns
    -------
    list
        List of processed maps after applying the specified operations.

    Notes
    -----
    Operations can be skipped by specifying them in `skip`, e.g., 'filter', 'prepare', 'sequence_processing'.
    """
    if 'filter' not in str(skip or ''):
        map_list = filter_maps(map_list, **kwargs)

    if 'prepare' not in str(skip or ''):
        map_list = prepare_maps(map_list, **kwargs)

    if 'sequence_processing' not in str(skip or ''):
        map_list = maps_sequence_processing(map_list, **kwargs)

    return map_list


# TODO: Implement units here
class model_fittings:
    """
    A class to store and manage the fittings of the geometrical model.

    Parameters
    ----------
    event_selected : str
        The selected event for the fitting.
    date_process : str
        The date of the fitting process.
    geometrical_model : str
        The type of the geometrical model used.
    model_parameters : pandas.DataFrame
        The parameters of the geometrical model.
    kinematics : dict, optional
        The kinematic parameters of the model (default is {'fit_method': None}).

    Methods
    -------
    load_from_json(json_file)
        Loads model fittings from a JSON file.
    model_id()
        Generates a unique ID for the model.
    get_geomertical_model(index)
        Retrieves the geometrical model at a given index.
    to_dict()
        Converts the model fittings to a dictionary.
    to_json(buffer=False)
        Converts the model fittings to JSON format.
    """

    def __init__(self, event_selected, date_process, geometrical_model, model_parameters, kinematics={'fit_method': None}):
        self.event_selected = event_selected
        self.date_process = date_process
        self.geometrical_model = geometrical_model
        self.parameters = model_parameters
        self.kinematics = kinematics

    @staticmethod
    def load_from_json(json_file):
        """
        Loads the model fittings from a JSON file.

        Parameters
        ----------
        json_file : str or file object
            The JSON file containing the model fittings.

        Returns
        -------
        model_fittings
            An instance of the model_fittings class.
        """
        if isinstance(json_file, str):
            with open(json_file, 'r') as file:
                json_content = file.read()
                fitting = json.loads(json_content)
        else:
            fitting = json.loads(json_file.read())

        table_indx = [datetime.datetime.strptime(t, '%Y-%m-%dT%H:%M:%S.%f') for t in fitting['geometrical_model']['parameters']['time']]
        parameters = pd.DataFrame(fitting['geometrical_model']['parameters'], index=table_indx)
        parameters = parameters.drop(['time'], axis=1)
        if 'kinematics' in fitting:
            kinematics = fitting['kinematics']
        else:
            # TODO: Remove this in version 1.0.0
            warnings.warn('The .json fitting file does not contain the "kinematics" information. \n \
                            This means that the file was prodused using Pythea with <V0.6.0. \n \
                            To resolve this: Just do a save of the loaded fitting now and replace \n \
                            the old file with the new one. This will not alter your fittings.')
            kinematics = {'fit_method': {'type': 'polynomial', 'order': 1}}

        if fitting['version'] <= '0.11.0':
            # In version >0.11.0 more AIA channels added so this will update the fitting files
            # prodused at earlier versions by replasing the imager values from "AIA" to "AIA-193".
            parameters['imager'] = parameters['imager'].replace('AIA', 'AIA-193')
            # In version >0.11.0 the fits file name is added to the single fits
            if 'fits_file' not in parameters.columns:
                parameters['fits_file'] = ''

        model_fittings_class = model_fittings(fitting['event_selected'],
                                              fitting['date_process'],
                                              fitting['geometrical_model']['type'],
                                              parameters,
                                              kinematics=kinematics)

        return model_fittings_class

    def model_id(self):
        """
        Generates a unique ID for the model based on event_selected and geometrical_model.

        Returns
        -------
        str
            The unique ID of the model.
        """
        str_id = self.event_selected.replace('-', '').replace(':', '').replace('|', 'D').replace('.', 'p') + 'M' + self.geometrical_model
        return str_id

    def get_geomertical_model(self, index):
        """
        Retrieves the geometrical model at the given index.

        Parameters
        ----------
        index : int or str
            The index of the model parameters.

        Returns
        -------
        object
            The geometrical model object.
        """
        if isinstance(index, int):
            model_parameters = self.parameters.iloc[index]
        elif isinstance(index, str):
            model_parameters = self.parameters.loc[index]

        obstime = model_parameters.name

        Spher_rep = SphericalRepresentation(model_parameters['hgln']*u.degree,
                                            model_parameters['hglt']*u.degree,
                                            Distance(np.abs(model_parameters['rcenter']*u.R_sun)))
        center = SkyCoord(np.sign(model_parameters['rcenter']*u.R_sun) * Spher_rep.to_cartesian(),
                          frame=frames.HeliographicStonyhurst,
                          observer='earth',
                          obstime=obstime)

        if self.geometrical_model == 'Spheroid':
            model_shock = spheroid(center,
                                   model_parameters['radaxis']*u.R_sun,
                                   model_parameters['orthoaxis1']*u.R_sun)
        elif self.geometrical_model == 'Ellipsoid':
            model_shock = ellipsoid(center,
                                    model_parameters['radaxis']*u.R_sun,
                                    model_parameters['orthoaxis1']*u.R_sun,
                                    model_parameters['orthoaxis2']*u.R_sun,
                                    model_parameters['tilt']*u.degree)
        elif self.geometrical_model == 'GCS':
            model_shock = gcs(center,
                              model_parameters['height']*u.R_sun,
                              model_parameters['alpha']*u.degree,
                              model_parameters['kappa'],
                              model_parameters['tilt']*u.degree)

        return model_shock

    def to_dict(self):
        """
        Converts the model fittings to a dictionary.

        Returns
        -------
        dict
            The dictionary representation of the model fittings.
        """
        parameters = copy(self.parameters)
        parameters['time'] = parameters.index.strftime('%Y-%m-%dT%H:%M:%S.%f')
        parameters = parameters.to_dict(orient='list')
        model_fittings_dict = {'event_selected': self.event_selected,
                               'date_process': self.date_process,
                               'geometrical_model': {'type': self.geometrical_model,
                                                     'parameters': parameters},
                               'kinematics': self.kinematics
                               }
        return model_fittings_dict

    def to_json(self, buffer=False):
        """
        Converts the model fittings to JSON format.

        Parameters
        ----------
        buffer : bool, optional
            If True, returns a BytesIO buffer containing the JSON data. If False, returns a JSON string (default is False).

        Returns
        -------
        str or io.BytesIO
            The JSON representation of the model fittings. Returns a JSON string if buffer is False,
            otherwise returns a BytesIO buffer containing the JSON data.
        """
        model_dict = self.to_dict()
        model_dict.update({'date_created': (datetime.datetime.utcnow()).strftime('%Y-%m-%dT%H:%M:%S')})
        model_dict.update({'version': version})
        if buffer:
            json_buffer = io.BytesIO()
            return json_buffer.write(json.dumps(model_dict, indent=' ').encode())
        else:
            return json.dumps(model_dict, indent=' ')


def plot_fitting_model(model, fit_args, plt_type='HeightT', fig=None, axis=None):
    """
    Plot the height (or speed/acceleration) evolution of the fitting parameters over time.

    Parameters
    ----------
    model : model_fittings
        The model_fittings object containing the geometrical model and its parameters.
    fit_args : dict
        A dictionary containing the fitting parameters:
            - 'type': str, type of fit ('polynomial', 'spline', or 'custom')
            - 'order': int, order of the polynomial or spline
            - 'smooth': float, smoothing factor for spline fitting (only for 'spline')
            - 'expression': str, custom expression for fitting (only for 'custom')
            - 'bounds': tuple, bounds for the custom fitting parameters (only for 'custom')
    plt_type : str, optional
        The type of plot to generate ('HeightT', 'SpeedT', 'AccelerationT', 'LongT', 'LatT').
        Default is 'HeightT'.
    fig : matplotlib.figure.Figure, optional
        The figure object to plot on. If None, a new figure is created.
    axis : matplotlib.axes.Axes, optional
        The axis object to plot on. If None, a new axis is created.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object with the plot.
    axis : matplotlib.axes.Axes
        The axis object with the plot.
    """
    palete = sns.color_palette('deep')

    Rs2km = (1 * u.R_sun).to_value(u.km)
    sec = 24*60*60

    parameters = {
        'Spheroid': {
            'height': ['+', '', palete[3], 'h-apex'],
            'orthoaxis1': ['x', '', palete[0], 'r-axis1'],
        },
        'Ellipsoid': {
            'height': ['+', '', palete[3], 'h-apex'],
            'orthoaxis1': ['x', '', palete[0], 'r-axis1'],
            'orthoaxis2': ['x', '', palete[2], 'r-axis2'],
        },
        'GCS': {
            'height': ['+', '', palete[3], 'h-apex'],
            'rapex': ['x', '', palete[0], 'r-apex'],
        },
    }

    parameters = parameters[model.geometrical_model]

    if fig is None and axis is None:
        fig, axis = plt.subplots(figsize=(5.5, 5.5), tight_layout=True)
    else:
        fig, axis = fig, axis

    if plt_type == 'HeightT':
        for p in parameters.keys():
            axis.plot(
                model.parameters.index,
                model.parameters[p],
                marker=parameters[p][0],
                linestyle=parameters[p][1],
                color=parameters[p][2],
                label=parameters[p][3]
            )
            if len(model.parameters[p]) - 1 > fit_args['order']:
                fit = parameter_fit(model.parameters.index, model.parameters[p], fit_args)
                axis.plot(fit['best_fit_x'], fit['best_fit_y'], '-', color=parameters[p][2])
                axis.fill_between(
                    fit['best_fit_x'],
                    fit['sigma_bounds']['up'],
                    fit['sigma_bounds']['low'],
                    color=parameters[p][2],
                    alpha=0.20
                )
                if fit_args['type'] == 'spline':
                    axis.fill_between(
                        fit['best_fit_x'],
                        fit['sigv_bounds']['up'],
                        fit['sigv_bounds']['low'],
                        color=parameters[p][2],
                        alpha=0.05
                    )
            else:
                axis.plot(model.parameters.index, model.parameters[p], '--', color=parameters[p][2])
        ylabel = 'Height of apex and length of flanks [Rsun]'
        axis.set_ylim(bottom=0)
    elif plt_type in ('SpeedT', 'AccelerationT'):
        gradient_power = {'SpeedT': 1, 'AccelerationT': 2}
        const = {'SpeedT': Rs2km / sec, 'AccelerationT': Rs2km / sec**2}
        for p in parameters.keys():
            if len(model.parameters[p]) - gradient_power[plt_type] > fit_args['order']:
                fit = parameter_fit(model.parameters.index, model.parameters[p], fit_args)
                if plt_type == 'SpeedT':
                    best_fit = const[plt_type] * np.gradient(fit['best_fit_y'], fit['best_fit_x_num'])
                    upper_bound = const[plt_type] * np.gradient(fit['sigma_bounds']['up'], fit['best_fit_x_num'])
                    lower_bound = const[plt_type] * np.gradient(fit['sigma_bounds']['low'], fit['best_fit_x_num'])
                    if fit_args['type'] == 'spline':
                        axis.fill_between(
                            fit['best_fit_x'],
                            const[plt_type] * fit['sigv_bounds']['dlow'],
                            const[plt_type] * fit['sigv_bounds']['dup'],
                            color=parameters[p][2],
                            alpha=0.05
                        )
                    ylabel = 'Speed [km/s]'
                elif plt_type == 'AccelerationT':
                    best_fit = const[plt_type] * np.gradient(
                        np.gradient(fit['best_fit_y'], fit['best_fit_x_num'], edge_order=2),
                        fit['best_fit_x_num'],
                        edge_order=2
                    )
                    upper_bound = const[plt_type] * np.gradient(
                        np.gradient(fit['sigma_bounds']['up'], fit['best_fit_x_num'], edge_order=2),
                        fit['best_fit_x_num'],
                        edge_order=2
                    )
                    lower_bound = const[plt_type] * np.gradient(
                        np.gradient(fit['sigma_bounds']['low'], fit['best_fit_x_num'], edge_order=2),
                        fit['best_fit_x_num'],
                        edge_order=2
                    )
                    if fit_args['type'] == 'spline':
                        axis.fill_between(
                            fit['best_fit_x'],
                            const[plt_type] * fit['sigv_bounds']['ddlow'],
                            const[plt_type] * fit['sigv_bounds']['ddup'],
                            color=parameters[p][2],
                            alpha=0.05
                        )
                    ylabel = 'Acceleration [km/s$^2$]'

                axis.plot(fit['best_fit_x'], best_fit, '-', color=parameters[p][2], label=parameters[p][3])
                axis.fill_between(fit['best_fit_x'], lower_bound, upper_bound, color=parameters[p][2], alpha=0.20)
            else:
                ylabel = ' '
                axis.text(
                    0.5, 0.5,
                    f'Not enough points for \n fitting with order {fit_args["order"]}.',
                    transform=axis.transAxes,
                    fontsize=20, color='gray', alpha=0.5,
                    ha='center', va='center', rotation=30
                )
                break

    if plt_type != 'AccelerationT':
        axis.set_ylim(bottom=0)

    if plt_type == 'LongT' or plt_type == 'LatT':
        parameters = {
            'LongT': ['hgln', '+', palete[3], 'Longitude'],
            'LatT': ['hglt', '+', palete[2], 'Latitude']
        }
        axis.plot(
            model.parameters.index,
            model.parameters[parameters[plt_type][0]],
            marker=parameters[plt_type][1],
            linestyle='',
            color=parameters[plt_type][2],
            label=parameters[plt_type][3]
        )
        if len(model.parameters[parameters[plt_type][0]]) - 1 > 3:
            fit = parameter_fit(model.parameters.index, model.parameters[parameters[plt_type][0]], fit_args)
            axis.plot(fit['best_fit_x'], fit['best_fit_y'], '-', color=parameters[plt_type][2])
            axis.fill_between(
                fit['best_fit_x'],
                fit['sigma_bounds']['up'],
                fit['sigma_bounds']['low'],
                color=parameters[plt_type][2],
                alpha=0.05
            )
        ylabel = parameters[plt_type][3] + ' [degrees]'

    axis.set_xlabel('Time [UTC]')
    axis.set_ylabel(ylabel)
    if fit_args['type'] == 'polynomial':
        title = 'Event: ' + model.event_selected + ' | ' + fit_args['type'] + str(fit_args['order'])
    elif fit_args['type'] == 'spline':
        title = 'Event: ' + model.event_selected + ' | ' + fit_args['type'] + str(fit_args['order']) + ' (' + str(fit_args['smooth']) + ')'
    elif fit_args['type'] == 'custom':
        title = 'Event: ' + model.event_selected + ' | Funct: ' + fit_args['expression']

    axis.set_title(title)
    locator = mdates.AutoDateLocator(minticks=4, maxticks=8)
    formatter = mdates.ConciseDateFormatter(locator)
    axis.xaxis.set_major_locator(locator)
    axis.xaxis.set_major_formatter(formatter)
    axis.minorticks_on()
    xlim = axis.get_xlim()

    hour_threshold = 24 * (xlim[1] - xlim[0])

    if hour_threshold <= 2:
        axis.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=np.arange(0, 61, 1)))
    elif hour_threshold <= 6:
        axis.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=np.arange(0, 61, 5)))
    else:
        axis.xaxis.set_minor_locator(mdates.HourLocator(byhour=np.arange(0, 25, 1)))

    fig.autofmt_xdate(bottom=0, rotation=0, ha='center')
    axis.legend(loc='lower right')

    return fig, axis


def parameter_fit(x, y, fit_args):
    """
    Performs a polynomial or spline fit to a set of x, y parameters.

    Parameters
    ----------
    x : array-like
        The independent variable, typically datetime objects.
    y : array-like
        The dependent variable, numerical data to be fitted.
    fit_args : dict
        A dictionary containing the fitting parameters:
            - 'type': str, type of fit ('polynomial', 'spline', or 'custom')
            - 'order': int, order of the polynomial or spline
            - 'smooth': float, smoothing factor for spline fitting (only for 'spline')
            - 'expression': str, custom expression for fitting (only for 'custom')
            - 'bounds': tuple, bounds for the custom fitting parameters (only for 'custom')

    Returns
    -------
    fitting : dict
        A dictionary containing the fitting results:
            - 'popt': optimal parameters for the fit (only for 'polynomial' and 'custom')
            - 'pcov': covariance of the optimal parameters (only for 'polynomial' and 'custom')
            - 'sigma': standard deviation of the residuals
            - 'best_fit_x_num': array, x-values for the best fit
            - 'best_fit_x': array, datetime index corresponding to best_fit_x_num
            - 'best_fit_y': array, y-values for the best fit
            - 'sigma_bounds': dict, upper and lower sigma bounds for the fit
            - 'sigv_bounds': dict, velocity and acceleration sigma bounds for the spline fit (only for 'spline')
    """
    minx = x[0].replace(second=0, microsecond=0)
    xx = (mdates.date2num(x) - mdates.date2num(minx))  # fractional days from minx
    step = 1/(60*24)  # one minute timestep in units of xx (fractional days)
    xxx = np.arange(0, xx.max() + step, step)  # Makes an array with a timestep of one minute in fractional days
    dd = pd.DatetimeIndex(mdates.num2date(xxx + mdates.date2num(minx)))

    if fit_args['type'] == 'polynomial':
        popt, pcov = np.polyfit(xx, y, fit_args['order'], full=False, cov=True)
        sigma = np.sqrt(np.diagonal(pcov))  # calculate sigma from covariance matrix
        best_fit = np.polyval(popt, xxx)
        sigma_bound_up = np.polyval((popt + sigma), xxx)
        sigma_bound_low = np.polyval((popt - sigma), xxx)
        fitting = {
            'popt': popt,
            'pcov': pcov,
            'sigma': sigma,
            'best_fit_x_num': xxx,
            'best_fit_x': dd,
            'best_fit_y': best_fit,
            'sigma_bounds': {
                'up': sigma_bound_up,
                'low': sigma_bound_low
            },
        }
    elif fit_args['type'] == 'spline':
        spl = UnivariateSpline(xx, y, s=fit_args['smooth'], k=fit_args['order'])
        resid = y - spl(xx)  # true - prediction
        sigma = np.nanstd(resid, axis=0)  # sigma = std(res)
        best_fit = spl(xxx)
        sigma_bound_up = best_fit + sigma
        sigma_bound_low = best_fit - sigma

        sv_bound_up, sv_bound_low, sv_bound_dup, sv_bound_dlow, sv_bound_ddup, sv_bound_ddlow = \
            best_fit, best_fit, np.gradient(best_fit, xxx), np.gradient(best_fit, xxx), \
            np.gradient(np.gradient(best_fit, xxx)), np.gradient(np.gradient(best_fit, xxx))
        for i in range(2, 100):
            spl = UnivariateSpline(xx, y, s=i/100, k=fit_args['order'])
            sv_bound_up = np.maximum(sv_bound_up, spl(xxx))
            sv_bound_low = np.minimum(sv_bound_low, spl(xxx))
            sv_bound_dup = np.maximum(sv_bound_dup, np.gradient(spl(xxx), xxx))
            sv_bound_dlow = np.minimum(sv_bound_dlow, np.gradient(spl(xxx), xxx))
            sv_bound_ddup = np.maximum(sv_bound_ddup, np.gradient(np.gradient(spl(xxx), xxx)))
            sv_bound_ddlow = np.minimum(sv_bound_ddlow, np.gradient(np.gradient(spl(xxx), xxx)))

        fitting = {
            'spl': spl,
            'sigma': sigma,
            'best_fit_x_num': xxx,
            'best_fit_x': dd,
            'best_fit_y': best_fit,
            'sigma_bounds': {
                'up': sigma_bound_up,
                'low': sigma_bound_low
            },
            'sigv_bounds': {
                'up': sv_bound_up,
                'low': sv_bound_low,
                'dup': sv_bound_dup,
                'dlow': sv_bound_dlow,
                'ddup': sv_bound_ddup,
                'ddlow': sv_bound_ddlow
            },
        }
    elif fit_args['type'] == 'custom':
        expression = fit_args['expression']  # 'a * exp(-b * x) + c'

        def func(x, a, b, c):
            return numexpr.evaluate(expression)

        popt, pcov = curve_fit(func, xx, y, bounds=fit_args['bounds'], maxfev=2500)
        sigma = np.sqrt(np.diagonal(pcov))  # calculate sigma from covariance matrix
        best_fit = func(xxx, *popt)
        sigma_bound_up = func(xxx, *(popt + sigma))
        sigma_bound_low = func(xxx, *(popt - sigma))
        fitting = {
            'popt': popt,
            'pcov': pcov,
            'sigma': sigma,
            'best_fit_x_num': xxx,
            'best_fit_x': dd,
            'best_fit_y': best_fit,
            'sigma_bounds': {
                'up': sigma_bound_up,
                'low': sigma_bound_low
            },
        }
    return fitting
