"""
    PyThea: A software package to perform forward modeling of CMEs and
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


import io
from operator import attrgetter
import numpy as np
import pandas as pd
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import sunpy.map
from sunpy.map.maputils import contains_coordinate
from sunpy.coordinates import frames
from sunpy.net import Fido, hek
from sunpy.net import attrs as a
from sunpy.coordinates import get_horizons_coord
import datetime
import astropy.units as u
from astropy.coordinates import (
    SkyCoord,
    Distance,
    SphericalRepresentation,
    CartesianRepresentation,
    concatenate
)
import json
from copy import copy
from scipy.optimize import curve_fit
import seaborn as sns

from config.selected_imagers import imager_dict
from config.selected_bodies import bodies_dict
from sunpy_dev.map.maputils import (difference_maps, normalize_exposure)

def get_hek_flare(day):
    flare_list = Fido.search(a.Time(day, day + datetime.timedelta(days=1)),
                            a.hek.EventType("FL"),
                            a.hek.FL.GOESCls > "B1.0",
                            a.hek.OBS.Observatory == "GOES")
    if len(flare_list["hek"]) == 0:
        selectbox_list = ['No events returned',]
        flare_list_ = []
    else:
        flare_list_ = flare_list["hek"]["event_starttime", "event_peaktime",
                                        "event_endtime", "fl_goescls",
                                        "fl_peakflux","hgs_x","hgs_y","ar_noaanum"]
        selectbox_list = []
        for flares in flare_list_:
            selectbox_list.append((f'FL{flares["fl_goescls"]}|{flares["event_peaktime"]}'))

    return selectbox_list, flare_list_

def make_figure(map, image_mode, clim=[-20,20]):
    fig = plt.figure()
    axis = plt.subplot(projection=map)
    #TODO: For plain images or when EUVIA-B are used this does not work very well.
    if image_mode=='Plain':
        map.plot()
    else:
        map.plot(cmap='Greys_r',
                 norm=colors.Normalize(vmin=clim[0], vmax=clim[1]))
    map.draw_limb()
    #map.draw_grid(linewidth=2, color='red') # TODO: This takes too much computation time. Maybe for AIA or EUVI?
    yax = axis.coords[1]
    yax.set_ticklabel(rotation=90)

    return fig, axis

def plot_bodies(axis, bodies_list, smap):
    for body in bodies_list:
        body_coo = get_horizons_coord(bodies_dict[body][0], smap.date)
        if contains_coordinate(smap, body_coo):
            axis.plot_coord(body_coo, 'o', color=bodies_dict[body][1], 
                fillstyle='none', markersize=6, label=body)


def download_fits(date_process, imager, time_range=[-1,1]):
    timerange = a.Time(date_process + datetime.timedelta(hours=time_range[0]),
                       date_process + datetime.timedelta(hours=time_range[1]))
    
    map_ = {}
    args = imager_dict[imager][0]
    extra = imager_dict[imager][1]
    result = Fido.search(timerange, *args)
    print(result)
    if len(result['vso'])!=0:
        nsuper = 2; # Superpixel by x2
        super_dim = u.Quantity([nsuper, nsuper] * u.pixel)
        downloaded_files = Fido.fetch(result)
        map_ = sunpy.map.Map(downloaded_files)
        map_ = [mask_occulter(smap) for smap in map_]
        map_ = [smap.superpixel(super_dim) for smap in map_]
        map_ = filter_maps(map_, extra)
        #map_ = [smap.rotate(recenter=True) for smap in map_] # Keep this afrer filter because images resize   
    else:
        map_ = []

    return map_

def maps_process(session_state, imagers_list, image_mode):
    session_state.map = {}
    session_state.imagers_list_ = []
    for imager in imagers_list:
        if imager not in session_state.map_ or session_state.map_[imager]==[]:
            pass
        else:
            session_state.map[imager] = map_diff(session_state.map_[imager], image_mode=image_mode)
            session_state.imagers_list_.append(imager)
    
    return session_state

def maps_clims(session_state, imagers_list):
    session_state.map_clim = {}
    for imager in imagers_list:
        if imager not in session_state.map or session_state.map[imager]==[]:
            pass
        else:
            map_ = session_state.map[imager][0]
            session_state.map_clim[imager] = [np.nanquantile(map_.data, 0.20), np.nanquantile(map_.data, 0.80)]
    
    return session_state

def filter_maps(map_sequence, extra):
    
    indices = []
    # TODO: This has been added because VSO search returns duplicates for WISPR
    #       so we manualy filter the dublicates until this problem is solved https://github.com/sunpy/sunpy/issues/5481
    if 'dublicates' in extra: 
        map_sequence.sort(key=attrgetter('date'))
        for i in range(0,len(map_sequence)-1):
            if map_sequence[i].date == map_sequence[i+1].date:
                indices.append(i)
        for i in sorted(indices, reverse=True):
            map_sequence.pop(i)

    if 'dimensions' in extra:
        map_sequence = [tmap for tmap in map_sequence if (2*tmap.dimensions[0],2*tmap.dimensions[1]) == extra['dimensions']]

    if 'polar' in extra:
        map_sequence = [tmap for tmap in map_sequence if tmap.meta['polar'] == extra['polar']]

    if len(map_sequence)!=0:
        sequence_final = sunpy.map.Map(map_sequence, sequence=True)
    else:
        sequence_final = []

    return sequence_final

def map_diff(map_sequence, image_mode='Plain'):
    smap = []
    if image_mode == 'Running Diff.':
        for i in range(1,len(map_sequence)):
            smap_diff = difference_maps(map_sequence[i], map_sequence[i-1])
            smap.append(smap_diff)
    if image_mode == 'Base Diff.':
        for i in range(1,len(map_sequence)):
            smap_diff = difference_maps(map_sequence[i], map_sequence[0])
            smap.append(smap_diff)
    if image_mode == 'Plain':
        for i in range(0,len(map_sequence)):
            smap.append(normalize_exposure(map_sequence[i]))

    if len(smap)!=0:
        sequence_final = sunpy.map.Map(smap, sequence=True)
    else:
        sequence_final = []

    return sequence_final

def mask_occulter(smap, apply_mask = True, mask_value = 0):
    """
    Returns a mask of the coronagraphic occulters
    """
    if smap.detector == 'COR1':
        xcen, ycen = smap.meta['crpix1'], smap.meta['crpix2']
        in_fac, out_fac = (100/512), (280/512)
    elif smap.detector == 'COR2':
        xcen, ycen = smap.meta['crpix1'], smap.meta['crpix2']
        in_fac, out_fac = (150/2048), (1030/2048)
    elif smap.detector == 'C2':
        xcen = smap.meta['crval1'] / smap.meta['cdelt1'] + ( smap.meta['crpix1'])
        ycen = smap.meta['crval2'] / smap.meta['cdelt2'] + ( smap.meta['crpix2'])
        in_fac, out_fac = (170/1024), (670/1024)
    elif smap.detector == 'C3':
        xcen = smap.meta['crval1'] / smap.meta['cdelt1'] + ( smap.meta['crpix1'])
        ycen = smap.meta['crval2'] / smap.meta['cdelt2'] + ( smap.meta['crpix2'])
        in_fac, out_fac = (60/1024), (560/1024)
    else:
        if apply_mask:
            return smap
        else:
            return []
    
    inner_distance = in_fac * smap.dimensions.x
    outer_distance = out_fac * smap.dimensions.x
    
    x, y = np.meshgrid(*[np.arange(v.value) for v in smap.dimensions]) * u.pixel
    xprime_sq = (x-xcen* u.pixel) ** 2
    yprime_sq = (y-ycen* u.pixel) ** 2
    mask_inner = np.sqrt( xprime_sq + yprime_sq ) < inner_distance
    mask_outer = np.sqrt( xprime_sq + yprime_sq ) > outer_distance
    
    if apply_mask:
        smap.data[mask_inner+mask_outer] = mask_value
        return sunpy.map.Map(smap.data, smap.meta)
    else:
        return mask_inner+mask_outer

# TODO: Implement units here
class model_fittings:
    def __init__(self, event_selected, date_process, geometrical_model, model_parameters):
        self.event_selected = event_selected
        self.date_process = date_process
        self.geometrical_model = geometrical_model
        self.parameters = model_parameters

    def model_id(self):
        str_id = self.event_selected.replace('-','').replace(':','').replace('|','D').replace('.','p') + 'M' + self.geometrical_model
        return str_id

    def to_jsonbuffer(self):
        parameters = copy(self.parameters)
        parameters['time'] = parameters.index.strftime("%Y-%m-%dT%H:%M:%S.%f")
        parameters = parameters.to_dict(orient='list')
        fitting_full_file = {'event_selected': self.event_selected,
                             'date_process': self.date_process,
                             'geometrical_model': self.geometrical_model,
                             'parameters': parameters
                             }
        json_buffer = io.BytesIO()
        json_buffer.write(json.dumps(fitting_full_file,indent=' ').encode())

        return json_buffer

def plot_fitting_model(model, order=2, plt_type='HeightT'):
    palete = sns.color_palette("deep")
    parameters = {'Spheroid': 
                             {'height': ['+', '', palete[3], 'h-apex'],
                              'orthoaxis1': ['x', '', palete[0], 'r-axis1'],
                              },
                  'GCS': {'height': ['+', '', palete[3], 'h-apex'],
                          'rappex': ['x', '', palete[0], 'r-apex'],
                                },
                      }
    parameters = parameters[model.geometrical_model]
    # Height vs time
    fig = plt.figure(figsize=(5.5,5.5), tight_layout=True)
    axis = plt.subplot()
    for p in parameters.keys():
        if plt_type == 'HeightT':
            plt.plot(model.parameters.index,
                     model.parameters[p],
                     marker=parameters[p][0],
                     linestyle=parameters[p][1],
                     color=parameters[p][2],
                     label=parameters[p][3]) # label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt)
            if len(model.parameters[p])-1>order:
                # How to get confidence intervals from curve_fit?
                fit = parameter_fit(model.parameters.index, model.parameters[p], order)
                plt.plot(fit['best_fit_x'], fit['best_fit_y'], '-', color=parameters[p][2])
                plt.fill_between(fit['best_fit_x'], fit['bound_lower'], fit['bound_upper'],
                                color = parameters[p][2], alpha = 0.15)
            else:
                plt.plot(model.parameters.index, model.parameters[p], '--', color=parameters[p][2])
            ylabel = 'Height [Rsun]'
        elif plt_type == 'SpeedT':
            if len(model.parameters[p])-1>order:
                # How to get confidence intervals from curve_fit?
                fit = parameter_fit(model.parameters.index, model.parameters[p], order)
                Rs2km = (1 * u.R_sun).to_value(u.km)
                sec = 24*60*60
                speed_best_fit = (Rs2km/sec) * np.gradient(fit['best_fit_y'], fit['best_fit_x_num'])
                speed_bound_upper = (Rs2km/sec) * np.gradient(fit['bound_upper'], fit['best_fit_x_num'])
                speed_bound_lower = (Rs2km/sec) * np.gradient(fit['bound_lower'], fit['best_fit_x_num'])
                plt.plot(fit['best_fit_x'], speed_best_fit, '-', color=parameters[p][2], label=parameters[p][3])
                plt.fill_between(fit['best_fit_x'], speed_bound_lower, speed_bound_upper,
                                 color = parameters[p][2], alpha = 0.15)
            else:
                pass
            ylabel = 'Speed [km/s]'
    plt.xlabel('Time [UT]')
    plt.ylabel(ylabel)
    plt.title('Event: '+model.event_selected)
    plt.gca().set_ylim(bottom=0)
    axis.xaxis.set_major_formatter(mdates.DateFormatter("%Y\n%b-%d\n%H:%M"))
    fig.autofmt_xdate(bottom=0, rotation=0, ha='center')
    plt.legend(loc='lower right')

    return fig

def parameter_fit(x, y, order=2):
    def fit_func(x, a, b, c):
        return a * x**2 + b * x + c
    xx = (mdates.date2num(x) - mdates.date2num(x[0]))
    ## scipy.optimize.curve_fit and numpy.polyfit
    popt, pcov = np.polyfit(xx, y, order, full=False, cov=True) #curve_fit(fit_func, xx, y)
    sigma = np.sqrt(np.diagonal(pcov))
    xxx = np.linspace(xx.min(), xx.max(), 120)
    dd = mdates.num2date(xxx + mdates.date2num(x[0]))
    best_fit = np.polyval(popt, xxx) #fit_func(xxx, *popt)
    bound_upper = np.polyval((popt + sigma), xxx) #fit_func(xxx, *(popt + sigma))
    bound_lower = np.polyval((popt - sigma), xxx) #fit_func(xxx, *(popt - sigma))
    fitting = {'popt': popt,
               'pcov': pcov,
               'sigma': sigma,
               'best_fit_x_num': xxx,
               'best_fit_x': dd,
               'best_fit_y': best_fit,
               'bound_upper': bound_upper,
               'bound_lower': bound_lower
               }
    return fitting
