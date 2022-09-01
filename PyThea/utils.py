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
from copy import copy

import astropy.units as u
import matplotlib.colors as colors
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sunpy.map
from config.selected_bodies import bodies_dict
from config.selected_imagers import imager_dict
from scipy.interpolate import UnivariateSpline
from sunpy.coordinates import get_horizons_coord
from sunpy.map.maputils import contains_coordinate
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy_dev.map.maputils import (filter_maps, maps_sequence_processing,
                                    prepare_maps)
from version import version


def get_hek_flare(day):
    '''
    Returns from HEK the flare list for a given day
    '''
    flare_list = Fido.search(a.Time(day, day + datetime.timedelta(days=1)),
                             a.hek.EventType('FL'),
                             a.hek.FL.GOESCls > 'B1.0',
                             a.hek.OBS.Observatory == 'GOES')
    if len(flare_list['hek']) == 0:
        selectbox_list = ['No events returned', ]
        flare_list_ = []
    else:
        flare_list_ = flare_list['hek']['event_starttime', 'event_peaktime',
                                        'event_endtime', 'fl_goescls',
                                        'fl_peakflux', 'hgs_x', 'hgs_y', 'ar_noaanum']
        selectbox_list = []
        for flares in flare_list_:
            fl_ = flares['fl_goescls']
            t_ = flares['event_peaktime'].strftime('%Y-%m-%dT%H:%M:%S')
            selectbox_list.append((f'FL{fl_}|{t_}'))

    return selectbox_list, flare_list_


def make_figure(map, image_mode, clim=[-20, 20], clip_model=True):
    '''
    Makes the main imager figure and returns the figure and axis handle.
    '''
    fig = plt.figure()
    axis = plt.subplot(projection=map)
    # TODO: For plain images or when EUVIA-B are used, this does not work very well.
    if image_mode == 'Plain':
        map.plot()
    else:
        map.plot(cmap='Greys_r',
                 norm=colors.Normalize(vmin=clim[0], vmax=clim[1]))
    map.draw_limb()
    # map.draw_grid(linewidth=2, color='red') # TODO: This takes too much computation time. Maybe for AIA or EUVI?
    yax = axis.coords[1]
    yax.set_ticklabel(rotation=90)
    if clip_model:
        axis.set_xlim([0, map.data.shape[0]])
        axis.set_ylim([0, map.data.shape[1]])

    return fig, axis


def plot_bodies(axis, bodies_list, smap):
    '''
    Plots in the images the possition of the pre-configured bodies (Earth, STA, Venus etc.)
    '''
    for body in bodies_list:
        body_coo = get_horizons_coord(bodies_dict[body][0], smap.date)
        if contains_coordinate(smap, body_coo):
            axis.plot_coord(body_coo, 'o', color=bodies_dict[body][1],
                            fillstyle='none', markersize=6, label=body)


def download_fits(date_process, imager, time_range=[-1, 1]):
    '''
    Downloads the imaging data (fits files) from VSO
    '''
    timerange = a.Time(date_process + datetime.timedelta(hours=time_range[0]),
                       date_process + datetime.timedelta(hours=time_range[1]))

    map_ = {}
    args = imager_dict[imager][0]
    result = Fido.search(timerange, *args)
    print(result)
    if result:
        downloaded_files = Fido.fetch(result)
        try:
            map_ = sunpy.map.Map(downloaded_files)
        except RuntimeError as err:
            print('Handling RuntimeError error:', err)
            map_ = []
        except OSError as err:
            print('Handling OSError error:', err)
            map_ = []
    else:
        map_ = []

    return map_


def maps_process(ninstr_map_in, imagers_list_in, image_mode):
    '''
    Process the images for the selected imagers and return the final maps.

    Note
    ----
    Here the ninstr_map_in is the session_state.map_ when used from the application.
    '''
    ninstr_map_out = {}
    imagers_list_out = []

    for imager in imagers_list_in:
        extra = imager_dict[imager][1]
        if imager in ninstr_map_in and ninstr_map_in[imager] != []:
            ninstr_map_out[imager] = filter_maps(ninstr_map_in[imager], extra)
            ninstr_map_out[imager] = prepare_maps(ninstr_map_out[imager], extra)
            ninstr_map_out[imager] = maps_sequence_processing(ninstr_map_out[imager],
                                                              seq_type=image_mode)
            if ninstr_map_out[imager] != []:
                imagers_list_out.append(imager)

    return ninstr_map_out, imagers_list_out


# TODO: Implement units here
class model_fittings:
    '''
    A class to store the fittings of the geometrical model.
    '''

    def __init__(self, event_selected, date_process, geometrical_model, model_parameters, kinematics={'fit_method': None}):
        self.event_selected = event_selected
        self.date_process = date_process
        self.geometrical_model = geometrical_model
        self.parameters = model_parameters
        self.kinematics = kinematics

    def model_id(self):
        str_id = self.event_selected.replace('-', '').replace(':', '').replace('|', 'D').replace('.', 'p') + 'M' + self.geometrical_model
        return str_id

    def to_dict(self):
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
        '''
        Returns the fittings of the geometrical model as a json format
        '''
        model_dict = self.to_dict()
        model_dict.update({'date_created': (datetime.datetime.utcnow()).strftime('%Y-%m-%dT%H:%M:%S')})
        model_dict.update({'version': version})
        if buffer:
            json_buffer = io.BytesIO()
            return json_buffer.write(json.dumps(model_dict, indent=' ').encode())
        else:
            return json.dumps(model_dict, indent=' ')


def plot_fitting_model(model, fit_args, plt_type='HeightT'):
    '''
    Plot the height(speed)--time evolution of the fitting parameters.
    '''
    palete = sns.color_palette('deep')
    parameters = {
        'Spheroid':
            {'height': ['+', '', palete[3], 'h-apex'],
             'orthoaxis1': ['x', '', palete[0], 'r-axis1'],
             },
        'Ellipsoid':
            {'height': ['+', '', palete[3], 'h-apex'],
             'orthoaxis1': ['x', '', palete[0], 'r-axis1'],
             'orthoaxis2': ['x', '', palete[2], 'r-axis2'],
             },
        'GCS':
            {'height': ['+', '', palete[3], 'h-apex'],
             'rapex': ['x', '', palete[0], 'r-apex'],
             },
    }
    parameters = parameters[model.geometrical_model]
    # Height vs time
    fig = plt.figure(figsize=(5.5, 5.5), tight_layout=True)
    axis = plt.subplot()
    if plt_type == 'HeightT':
        for p in parameters.keys():
            plt.plot(model.parameters.index,
                     model.parameters[p],
                     marker=parameters[p][0],
                     linestyle=parameters[p][1],
                     color=parameters[p][2],
                     label=parameters[p][3])  # label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt)
            if len(model.parameters[p])-1 > fit_args['order']:
                # How to get confidence intervals from curve_fit?
                fit = parameter_fit(model.parameters.index, model.parameters[p], fit_args)
                plt.plot(fit['best_fit_x'], fit['best_fit_y'], '-', color=parameters[p][2])
                plt.fill_between(fit['best_fit_x'], fit['sigma_bounds']['up'], fit['sigma_bounds']['low'],
                                 color=parameters[p][2], alpha=0.20)
                if fit_args['type'] == 'spline':
                    plt.fill_between(fit['best_fit_x'], fit['sigv_bounds']['up'], fit['sigv_bounds']['low'],
                                     color=parameters[p][2], alpha=0.05)
            else:
                plt.plot(model.parameters.index, model.parameters[p], '--', color=parameters[p][2])
        ylabel = 'Height or Length [Rsun]'
        plt.gca().set_ylim(bottom=0)
    elif plt_type == 'SpeedT':
        for p in parameters.keys():
            if len(model.parameters[p])-1 > fit_args['order']:
                # How to get confidence intervals from curve_fit?
                fit = parameter_fit(model.parameters.index, model.parameters[p], fit_args)
                Rs2km = (1 * u.R_sun).to_value(u.km)
                sec = 24*60*60
                speed_best_fit = (Rs2km/sec) * np.gradient(fit['best_fit_y'], fit['best_fit_x_num'])
                speed_bound_upper = (Rs2km/sec) * np.gradient(fit['sigma_bounds']['up'], fit['best_fit_x_num'])
                speed_bound_lower = (Rs2km/sec) * np.gradient(fit['sigma_bounds']['low'], fit['best_fit_x_num'])
                plt.plot(fit['best_fit_x'], speed_best_fit, '-', color=parameters[p][2], label=parameters[p][3])
                plt.fill_between(fit['best_fit_x'], speed_bound_lower, speed_bound_upper,
                                 color=parameters[p][2], alpha=0.20)
                if fit_args['type'] == 'spline':
                    plt.fill_between(fit['best_fit_x'], (Rs2km/sec) * fit['sigv_bounds']['dlow'], (Rs2km/sec) * fit['sigv_bounds']['dup'],
                                     color=parameters[p][2], alpha=0.05)
            else:
                pass
        ylabel = 'Speed [km/s]'
        plt.gca().set_ylim(bottom=0)
    if plt_type == 'LongT' or plt_type == 'LatT':
        parameters = {'LongT': ['hgln', '+', palete[3], 'Longitude'],
                      'LatT': ['hglt', '+', palete[2], 'Latitude']}
        plt.plot(model.parameters.index,
                 model.parameters[parameters[plt_type][0]],
                 marker=parameters[plt_type][1],
                 linestyle='',
                 color=parameters[plt_type][2],
                 label=parameters[plt_type][3])
        if len(model.parameters[parameters[plt_type][0]])-1 > 3:
            fit = parameter_fit(model.parameters.index, model.parameters[parameters[plt_type][0]], fit_args)
            plt.plot(fit['best_fit_x'], fit['best_fit_y'], '-', color=parameters[plt_type][2])
            plt.fill_between(fit['best_fit_x'], fit['sigma_bounds']['up'], fit['sigma_bounds']['low'],
                             color=parameters[plt_type][2], alpha=0.05)
        ylabel = parameters[plt_type][3] + ' [degrees]'

    plt.xlabel('Time [UT]')
    plt.ylabel(ylabel)
    if fit_args['type'] == 'polynomial':
        title = 'Event: ' + model.event_selected + ' | ' + fit_args['type'] + str(fit_args['order'])
    elif fit_args['type'] == 'spline':
        title = 'Event: ' + model.event_selected + ' | ' + fit_args['type'] + str(fit_args['order']) + ' (' + str(fit_args['smooth']) + ')'
    elif fit_args['type'] == 'custom':
        title = 'Event: ' + model.event_selected + ' | Funct: ' + fit_args['expression']

    plt.title(title)
    axis.xaxis.set_major_formatter(mdates.DateFormatter('%Y\n%b-%d\n%H:%M'))
    axis.minorticks_on()
    fig.autofmt_xdate(bottom=0, rotation=0, ha='center')
    plt.legend(loc='lower right')

    return fig


def parameter_fit(x, y, fit_args):
    '''
    Performs a polynomial or spline fit to a set of x, y parameters.
    '''
    xx = (mdates.date2num(x) - mdates.date2num(x[0]))
    xxx = np.linspace(xx.min(), xx.max(), 120)
    dd = mdates.num2date(xxx + mdates.date2num(x[0]))

    if fit_args['type'] == 'polynomial':
        # scipy.optimize.curve_fit and numpy.polyfit
        popt, pcov = np.polyfit(xx, y, fit_args['order'], full=False, cov=True)
        sigma = np.sqrt(np.diagonal(pcov))  # calculate sigma from covariance matrix
        best_fit = np.polyval(popt, xxx)
        sigma_bound_up = np.polyval((popt + sigma), xxx)
        sigma_bound_low = np.polyval((popt - sigma), xxx)
        fitting = {'popt': popt,
                   'pcov': pcov,
                   'sigma': sigma,
                   'best_fit_x_num': xxx,
                   'best_fit_x': dd,
                   'best_fit_y': best_fit,
                   'sigma_bounds': {'up': sigma_bound_up,
                                    'low': sigma_bound_low},
                   }
    elif fit_args['type'] == 'spline':
        spl = UnivariateSpline(xx, y, s=fit_args['smooth'], k=fit_args['order'])
        resid = y - spl(xx)  # true - prediction
        sigma = np.nanstd(resid, axis=0)  # sigma = std(res)
        best_fit = spl(xxx)
        sigma_bound_up = best_fit + sigma
        sigma_bound_low = best_fit - sigma

        sv_bound_up, sv_bound_low, sv_bound_dup, sv_bound_dlow = best_fit, best_fit, np.gradient(best_fit, xxx), np.gradient(best_fit, xxx)
        for i in range(2, 100):
            spl = UnivariateSpline(xx, y, s=i/100, k=fit_args['order'])
            sv_bound_up = np.maximum(sv_bound_up, spl(xxx))
            sv_bound_low = np.minimum(sv_bound_low, spl(xxx))
            sv_bound_dup = np.maximum(sv_bound_dup, np.gradient(spl(xxx), xxx))
            sv_bound_dlow = np.minimum(sv_bound_dlow, np.gradient(spl(xxx), xxx))

        fitting = {'spl': spl,
                   'sigma': sigma,
                   'best_fit_x_num': xxx,
                   'best_fit_x': dd,
                   'best_fit_y': best_fit,
                   'sigma_bounds': {'up': sigma_bound_up,
                                    'low': sigma_bound_low},
                   'sigv_bounds': {'up': sv_bound_up,
                                   'low': sv_bound_low,
                                   'dup': sv_bound_dup,
                                   'dlow': sv_bound_dlow},
                   }
    elif fit_args['type'] == 'custom':
        import numexpr
        from scipy.optimize import curve_fit

        expression = fit_args['expression']  # 'a * exp(-b * x) + c'

        def func(x, a, b, c):
            return numexpr.evaluate(expression)

        popt, pcov = curve_fit(func, xx, y, bounds=fit_args['bounds'], maxfev=2500)
        sigma = np.sqrt(np.diagonal(pcov))  # calculate sigma from covariance matrix
        best_fit = func(xxx, *popt)
        sigma_bound_up = func(xxx, *(popt + sigma))
        sigma_bound_low = func(xxx, *(popt - sigma))
        fitting = {'popt': popt,
                   'pcov': pcov,
                   'sigma': sigma,
                   'best_fit_x_num': xxx,
                   'best_fit_x': dd,
                   'best_fit_y': best_fit,
                   'sigma_bounds': {'up': sigma_bound_up,
                                    'low': sigma_bound_low},
                   }
    return fitting
