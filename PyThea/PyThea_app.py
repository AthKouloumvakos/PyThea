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
from copy import copy

import astropy.units as u
import numpy as np
import streamlit as st
from astropy.coordinates import Distance, SkyCoord, SphericalRepresentation
from sunpy.coordinates import frames
from sunpy.net import attrs as a

from PyThea.callbacks import load_or_delete_fittings
from PyThea.config import (app_styles, config_sliders, selected_bodies,
                           selected_imagers)
from PyThea.geometrical_models import ellipsoid, gcs, spheroid
from PyThea.modules import (date_and_event_selection, figure_streamlit,
                            final_parameters_gmodel,
                            fitting_and_slider_options_container,
                            fitting_sliders)
from PyThea.sunpy_dev.map.maputils import get_closest, maps_clims
from PyThea.utils import (download_fits, load_fits, model_fittings,
                          plot_fitting_model, single_imager_maps_process)
from PyThea.utils_database import get_fits_filenames_from_database
from PyThea.version import version


def delete_from_state(vars):
    for var in vars:
        if var in st.session_state:
            del st.session_state[var]


def highlight_row(row, row_index):
    if row.name.strftime('%Y-%m-%dT%H:%M:%S.%f') == row_index.strftime('%Y-%m-%dT%H:%M:%S.%f').values[0]:
        return ['background-color: lightpink']*len(row)
    else:
        return ['']*len(row)


def footer_text():
    st.subheader('About this application:')
    st.markdown("""
                   _PyThea_  is an open-source software package that can be used to
                   reconstruct the 3D structure of Coronal Mass Ejections (CMEs) and
                   shock waves and determine their kinematics using remote-sensing observations.
                   The tool implements the Graduate Cylindrical Shell (GCS) model that can be used
                   to reconstruct CMEs and two geometrical models, namely a spheroid and ellipsoid model
                   to reconstruct shock waves. It also implements remote-sensing observations
                   from multiple viewpoints such as the SOlar and Heliospheric Observatory (SOHO)
                   and Solar Terrestrial Relations Observatory (STEREO).
                """)
    right, left = st.columns((1, 1))
    right.markdown("""
                   **Github**: Find the latest version here
                               [![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)](https://github.com/AthKouloumvakos/PyThea) \n
                   **Documentation Page**: https://www.pythea.org/
                   """ +
                   f"""
                   **Version**: {version} (latest release [![Version](https://img.shields.io/github/v/release/AthKouloumvakos/PyThea)](https://github.com/AthKouloumvakos/PyThea/releases))
                   """)
    left.image('https://github.com/AthKouloumvakos/PyThea/blob/master/docs/logo/pythea_logo.png?raw=true')
    st.markdown("""
                **Citation**: Please cite the following paper [![https://www.frontiersin.org/articles/10.3389/fspas.2022.974137/](https://img.shields.io/static/v1?label=Paper&message=Frontiers&color=red)](https://www.frontiersin.org/articles/10.3389/fspas.2022.974137/) and [![https://doi.org/10.5281/zenodo.5713659](https://zenodo.org/badge/DOI/10.5281/zenodo.5713659.svg)](https://doi.org/10.5281/zenodo.5713659)
                """)
    st.info('''
            More imaging data have been added:
            - SDO/AIA images from 211A channel.
            - PSP/WISPR inner and outer telescope images.
            - SolO/EUI, METIS, and HI(tiles 1&2) images.
            ''', icon='ℹ️')
    st.warning('''
           **NOTE: From PyThea >0.8.1 the JSON fitting files will be slightly different from the old ones.**
           Due to a change in the fitting time input, the new fitting files may have slightly different times for the same images (see further information [here](https://github.com/AthKouloumvakos/PyThea/discussions/24)).
           ''', icon='⚠️')
    st.markdown('---')


def run():
    ###################################################################
    # IMPORTANT: DON'T CHANGE THE POSSITION OF THE COMPONENTS. NEVER! #
    ###################################################################

    #############################################################
    # set page config
    st.set_page_config(page_title='PyThea', page_icon=':rocket:',
                       initial_sidebar_state='expanded',
                       menu_items={
                           'Get Help': 'https://www.pythea.org/',
                           'Report a Bug': 'https://github.com/AthKouloumvakos/PyThea'})

    #############################################################
    # HTML Styles
    app_styles.apply(st)

    #############################################################
    # Main page information text
    st.header('PyThea: Reconstruct CMEs & Shocks')

    #############################################################
    # Startup Variables
    if 'startup' not in st.session_state:
        st.session_state.startup = {'fitting': True}

    #############################################################
    # Date and Event selection

    if 'date_process' not in st.session_state:
        date_and_event_selection(st)
    else:
        st.sidebar.markdown('## Processing Event|Date:')
        st.sidebar.info(f'{st.session_state.event_selected}')

    if 'date_process' not in st.session_state:
        st.markdown('---')
        footer_text()
        st.stop()

    #############################################################
    # 3D Model Selection
    st.sidebar.markdown('## 3D Fitting and Reconstruction')

    if 'geometrical_model' not in st.session_state:
        st.markdown(""" #### ⏻ Select the geometrical model you want to fit. """)
        geometrical_model = st.sidebar.selectbox('Geometrical model to fit',
                                                 options=['Select a model', 'Spheroid', 'Ellipsoid', 'GCS'])
        if geometrical_model != 'Select a model':
            st.session_state.geometrical_model = geometrical_model
            st.rerun()
    else:
        st.sidebar.info(f'Geometrical Model: \
                        {st.session_state.geometrical_model}')

    if 'geometrical_model' not in st.session_state:
        st.markdown('---')
        footer_text()
        st.stop()

    #############################################################
    # 3D Fitting and Reconstruction
    st.markdown('---')
    fitting_and_slider_options_container(st)

    if st.session_state.coord_system == 'HGC':
        long_val = [0., 360.]
    elif st.session_state.coord_system == 'HGS':
        long_val = [-180., 180.]

    longit = st.sidebar.slider(f'{st.session_state.coord_system} \
                               Longitude [deg.]:',
                               min_value=long_val[0],
                               max_value=long_val[1],
                               value=0.,
                               step=0.01, key='longit') * u.degree

    latitu = st.sidebar.slider(f'{st.session_state.coord_system} \
                               Latitude [deg.]:',
                               min_value=-90.,
                               max_value=90.,
                               value=0.,
                               step=0.01, key='latitu') * u.degree

    fitting_sliders(st)

    col1, col3, col2 = st.sidebar.columns(3)
    store_fit_button_pressed = st.sidebar.button('Store Fit')

    #############################################################
    # 3D Fitting and Reconstruction
    st.sidebar.markdown('## Imaging menu')
    with st.sidebar.expander('Download Options'):
        select_imagers_form = st.form(key='select_imagers_form', border=False)
        imagers_list = select_imagers_form.multiselect('Select Imagers',
                                                       options=selected_imagers.imager_dict.keys(),
                                                       default=['LC2', 'LC3', 'COR2A'],
                                                       key='imagers_list')
        select_imagers_form.form_submit_button(label='Submit')
        select_timerange_form = st.form(key='select_timerange_form', border=False)
        imaging_time_range = select_timerange_form.slider('Time Range [hours]',
                                                          min_value=-3., max_value=6.,
                                                          value=[-1., 1.], step=0.5,
                                                          key='imaging_time_range')
        select_timerange_form.form_submit_button(label='Submit',
                                                 on_click=delete_from_state,
                                                 kwargs={'vars': ['map', 'map_', 'imagers_list_', 'hek_responses']})

    with st.sidebar.expander('Processing Options'):
        procoption_container = st.container()
        if 'imagers_list_' not in st.session_state:
            # imagers_list_ is used later when we download or filter the images
            st.session_state['imagers_list_'] = []
        image_mode = procoption_container.selectbox('Map Sequence Processing',
                                                    options=['Running Diff.', 'Base Diff.', 'Plain'], key='image_mode',
                                                    on_change=delete_from_state,
                                                    kwargs={'vars': ['map', 'imagers_list_', 'clim', 'clim_manual_low', 'clim_manual_high']})
        if image_mode in ['Running Diff.', 'Base Diff.']:
            min_diff = 1 if image_mode == 'Running Diff.' else 0
            diff_value = procoption_container.number_input('Select Step/Image for Running/Base Diff.', min_diff, 5, key='diff_value',
                                                           on_change=delete_from_state,
                                                           kwargs={'vars': ['map', 'imagers_list_', 'clim', 'clim_manual_low', 'clim_manual_high']})
        else:
            diff_value = None

        procoption_container.number_input('Median Filter', 1, 6, 1, 1, key='images_median_filter')

    with st.sidebar.expander('Plot/View Options'):
        plotviewopt_container = st.container()
        plotviewopt_container.checkbox('Clip plot on image limits', value=True, key='clip_model')
        plt_supp_imagers = plotviewopt_container.checkbox('Supplementary Imaging', value=False)

    #############################################################
    # Magnetic Connectivity
    st.sidebar.markdown('## Overlays')
    with st.sidebar.expander('Bodies or s/c Location'):
        bodies_container = st.container()
        bodies_container.checkbox('View Bodies or s/c', key='star_field')
        bodies_container.multiselect('Select Bodies', options=selected_bodies.bodies_dict.keys(),
                                     default=['Mercury', 'Venus', 'Jupiter'],
                                     key='bodies_list',
                                     disabled=not st.session_state.star_field)
    with st.sidebar.expander('Limb and Meridians'):
        limb_meridians_container = st.container()
        limb_meridians_container.checkbox('View Limb and Meridians', key='plot_solar_reference_lines_')
        markdown = '''
            **Limb**: plots the solar limb as observed
            for the selected observers.
            **Central Meridian**: plots the central meridian
            observed for the selected observers.
            **CR Meridan+Equator**: plots the primary meridian
            and equator of Carrington coordinate system.
            '''.strip()
        limb_meridians_container.selectbox('Select plot option',
                                           options=['Limb from Obs.', 'Central Meridian from Obs.', 'Carr. Prime Meridian+Solar Equator',
                                                    'Stonyhurst Grid', 'Carrington Grid'],
                                           help=markdown,
                                           key='plot_solar_reference_lines_mode',
                                           disabled=not st.session_state.plot_solar_reference_lines_)
        disabled = False if (st.session_state.plot_solar_reference_lines_) and \
            (st.session_state.plot_solar_reference_lines_mode in ['Limb from Obs.', 'Central Meridian from Obs.']) else True
        limb_meridians_container.multiselect('Select Bodies',
                                             label_visibility='collapsed',
                                             options=selected_bodies.bodies_dict.keys(),
                                             default=['Earth'], disabled=disabled,
                                             key='plot_solar_reference_lines_bodies_list')
    with st.sidebar.expander('Magnetic Connectivity'):
        connectivity_container = st.container()
        connectivity_container.checkbox('Plot Parker spirals', key='plot_parker_spirals')
        connectivity_container.multiselect('Select bodies/spacecraft', options=selected_bodies.bodies_dict.keys(),
                                           default=['Earth'], key='mag_bodies_list',
                                           disabled=not st.session_state.plot_parker_spirals)
        st.session_state.sw_speed_select = {}
        for body in st.session_state.mag_bodies_list:
            connectivity_container.number_input(f'{body} solar wind speed', min_value=200, max_value=800, value=350, step=50,
                                                key=f'sw_speed_select_{body}',
                                                disabled=not st.session_state.plot_parker_spirals)
            st.session_state.sw_speed_select[body] = st.session_state[f'sw_speed_select_{body}'] * (u.km/u.second)
    with st.sidebar.expander('HEK feature/events'):
        hek_container = st.container()
        hek_container.checkbox('Plot HEK feature/events', key='plot_hek', on_change=delete_from_state,
                               kwargs={'vars': ['hek_responses', ]})
        hek_container.multiselect('Select HEK feature/events', options=['Active Regions', 'Coronal Holes', 'Flares'],
                                  default=['Flares'], key='hek_list', disabled=not st.session_state.plot_hek)

    #############################################################
    # Download and Process the Images
    # This part runs only if the map_ doesn't exits or the session_state.map_ does not contain all the imagers requested
    imager_added = list(set(imagers_list) - set(st.session_state['imagers_list_']))
    imager_removed = list(set(st.session_state['imagers_list_']) - set(imagers_list))
    if 'map_colormap_limits' not in st.session_state:
        st.session_state['map_colormap_limits'] = {}

    if 'map_' not in st.session_state or imager_added != []:
        st.session_state.map_ = {} if 'map_' not in st.session_state else st.session_state.map_
        st.session_state.map = {} if 'map' not in st.session_state else st.session_state.map
        st.session_state['imagers_list_'] = imagers_list
        progress_bar = st.progress(0, text='Preparing to Download data')
        for i, imager in enumerate(imager_added):
            progress_bar.progress(i/len(imager_added), text=f'Download {imager} images from VSO')
            if imager not in st.session_state.map_:
                timerange = a.Time(st.session_state.date_process + datetime.timedelta(hours=imaging_time_range[0]),
                                   st.session_state.date_process + datetime.timedelta(hours=imaging_time_range[1]))

                if st.session_state.offline_mode is False:
                    downloaded_files = download_fits(timerange, imager)
                elif st.session_state.offline_mode is True:
                    progress_bar.desc = f'Load {imager} images from local database.'
                    event_id = st.session_state.event_selected.replace('-', '').replace(':', '').replace('|', 'D').replace('.', 'p') \
                        + 'M' + st.session_state.geometrical_model
                    downloaded_files = get_fits_filenames_from_database(event_id, timerange, imager)

                st.session_state.map_[imager] = load_fits(downloaded_files)
                st.session_state.map_[imager] = single_imager_maps_process(st.session_state.map_[imager],
                                                                           **selected_imagers.imager_dict[imager]['process'],
                                                                           skip='sequence_processing')
            processed_images = single_imager_maps_process(st.session_state.map_[imager],
                                                          skip=['filter', 'prepare'],
                                                          image_mode=image_mode,
                                                          diff_num=diff_value)
            if processed_images != []:
                st.session_state.map[imager] = processed_images
                st.session_state.map_colormap_limits[imager] = maps_clims(st.session_state.map[imager])
            else:
                st.session_state.imagers_list_.remove(imager)
                st.session_state.map_colormap_limits.pop(imager, None)
        progress_bar.empty()

    if imager_removed != []:
        st.session_state['imagers_list_'] = imagers_list
        for imager in imager_removed:
            st.session_state.map_.pop(imager, None)

    if not st.session_state.imagers_list_:
        st.error('No images have been downloaded or processed.')  # TODO: Explain better
        st.stop()

    #############################################################
    # Selection for the primary image to plot
    st.markdown('### Multi-viewpoint imaging')
    col1, col2 = st.columns((1, 3))
    imager_select = col1.selectbox('Select an imager',
                                   options=st.session_state.imagers_list_,
                                   on_change=delete_from_state, kwargs={'vars': ['clim', 'clim_manual_low', 'clim_manual_high', 'hek_responses']})

    maps_date = [getattr(maps, 'date_average', None) or getattr(maps, 'date', None) for maps in st.session_state.map[imager_select]]
    if len(maps_date) > 1:
        running_map_date = col2.select_slider('Slide to the image time',
                                              options=maps_date, value=maps_date[0],
                                              key='running_map_date',
                                              on_change=delete_from_state, kwargs={'vars': ['hek_responses']})
    else:
        running_map_date = maps_date[0]

    running_map = st.session_state.map[imager_select][maps_date.index(running_map_date)]

    manual_clim = plotviewopt_container.checkbox('Provide clims values', on_change=delete_from_state, kwargs={'vars': ['clim']})
    if manual_clim:
        col1, col2 = plotviewopt_container.columns(2)
        if 'clim_manual_low' not in st.session_state:
            st.session_state.clim_manual_low = st.session_state.map_colormap_limits[imager_select][0]
        qmin = col1.number_input('Low clim value', label_visibility='collapsed', key='clim_manual_low')
        if 'clim_manual_high' not in st.session_state:
            st.session_state.clim_manual_high = st.session_state.map_colormap_limits[imager_select][1]
        qmax = col2.number_input('High clim value', label_visibility='collapsed', key='clim_manual_high')
        st.session_state.map_colormap_limits[imager_select] = [qmin, qmax]
    else:
        qmin = st.session_state.map_colormap_limits[imager_select][0]
        qmax = st.session_state.map_colormap_limits[imager_select][1]
        cmmin, cpmax = config_sliders.slider_image_pmclims[image_mode]
        if 'clim' not in st.session_state:
            st.session_state.clim = [float(qmin), float(qmax)]
        st.session_state.map_colormap_limits[imager_select] = \
            plotviewopt_container.slider('Images colorbar limits:',
                                         min_value=float(cmmin), max_value=float(cpmax),
                                         key='clim')

    #############################################################
    # Finalize the geometrical model

    if st.session_state.geometrical_model == 'Spheroid':
        rcenter, radaxis, orthoaxis1 = final_parameters_gmodel(st)
    elif st.session_state.geometrical_model == 'Ellipsoid':
        rcenter, radaxis, orthoaxis1, orthoaxis2, tilt = final_parameters_gmodel(st)
    elif st.session_state.geometrical_model == 'GCS':
        rcenter, height, alpha, kappa, tilt = final_parameters_gmodel(st)

    Spher_rep = SphericalRepresentation(longit, latitu,
                                        Distance(np.abs(rcenter)))
    if st.session_state.coord_system == 'HGC':
        center = SkyCoord(np.sign(rcenter) * Spher_rep.to_cartesian(),
                          frame=frames.HeliographicCarrington,
                          observer='earth',
                          obstime=running_map_date)
    elif st.session_state.coord_system == 'HGS':
        center = SkyCoord(np.sign(rcenter) * Spher_rep.to_cartesian(),
                          frame=frames.HeliographicStonyhurst,
                          observer='earth',
                          obstime=running_map_date)
    st.session_state.center = center

    if st.session_state.geometrical_model == 'Spheroid':
        model = spheroid(center, radaxis, orthoaxis1)
    elif st.session_state.geometrical_model == 'Ellipsoid':
        model = ellipsoid(center, radaxis, orthoaxis1, orthoaxis2, tilt)
    elif st.session_state.geometrical_model == 'GCS':
        model = gcs(center, height, alpha, kappa, tilt)

    #############################################################
    # Plot main and supplement figure images
    fig, axis = figure_streamlit(st, running_map, image_mode, imager_select, model)
    st.pyplot(fig)

    if st.session_state.plot_hek and st.session_state.hek_responses['Flares']:
        st.markdown('**HEK Flare List:**')
        st.write(st.session_state.hek_responses['Flares'].to_pandas())
        st.markdown('*Flares with no location have been removed from the table.')
    if st.session_state.plot_hek and st.session_state.hek_responses['Active Regions']:
        st.markdown('**HEK Active Regions List:**')
        st.write(st.session_state.hek_responses['Active Regions'].to_pandas())
        st.markdown('*Active Regions without NOAA number have been removed from the table.')

    if plt_supp_imagers:
        if len(st.session_state.imagers_list_) < 3:
            other_element = [element for element in st.session_state.imagers_list_ if element != imager_select][0]
            fig, axis = figure_streamlit(st, get_closest(st.session_state.map[other_element], running_map_date), image_mode, other_element, model)
            st.pyplot(fig)
        else:
            supl_imagers_list = copy(st.session_state.imagers_list_)
            supl_imagers_list.remove(imager_select)
            supl_imagers = st.select_slider('Select supplement imagers',
                                            options=supl_imagers_list,
                                            value=(supl_imagers_list[0],
                                                   supl_imagers_list[-1]),
                                            key='supl_imagers')
            col1, col2 = st.columns(2)
            fig, axis = figure_streamlit(st, get_closest(st.session_state.map[supl_imagers[0]], running_map_date), image_mode, supl_imagers[0], model)
            col1.pyplot(fig)
            fig, axis = figure_streamlit(st, get_closest(st.session_state.map[supl_imagers[1]], running_map_date), image_mode, supl_imagers[1], model)
            model.plot(axis, mode='Skeleton')
            col2.pyplot(fig)

    #############################################################
    # Store the fitting
    single_fit = model.to_dataframe()
    single_fit['imager'] = imager_select
    single_fit['fits_file'] = running_map.meta['fits_file']

    if store_fit_button_pressed:
        if 'model_fittings' not in st.session_state:
            st.session_state.model_fittings = model_fittings(st.session_state.event_selected,
                                                             st.session_state.date_process.strftime('%Y-%m-%dT%H:%M:%S.%f'),
                                                             st.session_state.geometrical_model,
                                                             single_fit)
        else:
            st.session_state.model_fittings.parameters = single_fit.combine_first(st.session_state.model_fittings.parameters)

    #############################################################
    # View the fittings table
    if st.session_state.fitting_table:
        st.markdown('---')
        st.markdown('### Parameters Table')
        st.markdown('**Running Fitting Table:**')
        col_formats = {col: '{:,.2f}'.format for col in single_fit.select_dtypes(include='number').columns}
        st.dataframe(single_fit.style.format(col_formats))
        if 'model_fittings' in st.session_state:
            st.markdown('**Stored Fitting Table:**')
            st.dataframe(st.session_state.model_fittings.parameters.style.apply(highlight_row,
                                                                                row_index=single_fit.index, axis=1).format(col_formats)
                         )
            col1, col2 = st.columns((2, 2))
            col2.selectbox('Select a fitting',
                           options=st.session_state.model_fittings.parameters.index,
                           key='fitting_select')
            if 'fit_action' not in st.session_state:
                st.session_state.fit_action = 'Select'
            col1.selectbox('Action',
                           options=['Select', 'Load', 'Delete'],
                           on_change=load_or_delete_fittings, args=[st],
                           key='fit_action')

    #############################################################
    # Plot the kinematics
    if st.session_state.plt_kinematics and \
       ('model_fittings' in st.session_state):
        st.markdown('---')
        st.markdown('### Plots of kinematics ')
        col1, col2, col3 = st.columns(3)
        plt_kinematics_select = col1.selectbox('Select Plots',
                                               options=['All', 'HeightT', 'SpeedT', 'AccelerationT', 'Long-LatT'])

        fit_mode = col2.selectbox('Select Fitting Mode',
                                  options=['polynomial', 'spline', 'custom'],
                                  key='fit_mode')
        if fit_mode == 'polynomial':
            polyfit_order = col3.slider('Polynomial order', value=2, min_value=1, max_value=4, step=1, key='polyfit_order')
            fit_args_ = {'type': 'polynomial', 'order': polyfit_order}
        elif fit_mode == 'spline':
            splinefit_order = col3.slider('Spline order', value=3, min_value=1, max_value=5, step=1, key='splinefit_order')
            splinefit_smooth = st.slider('Spline smooth', value=0.5, min_value=0., max_value=1., step=0.01, key='splinefit_smooth')
            fit_args_ = {'type': 'spline', 'order': splinefit_order, 'smooth': splinefit_smooth}
        elif fit_mode == 'custom':
            fitcustexpres_select = col3.selectbox('Select a custom function',
                                                  options=['a*exp(-b*x)+c'],  # 'a*sqrt(b*x)+c',  '((x+a)/(x+b))+c'
                                                  key='plt_fitcustexpres_select')
            fit_args_ = {'type': 'custom', 'expression': fitcustexpres_select,
                         'bounds': ([-np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf]), 'order': 3}
            st.info(f'A custom function {fit_args_["expression"]} was fitted when processing parameters.')

        if plt_kinematics_select == 'All':
            col1, col2 = st.columns(2)
            fig_ht, axis_ht = plot_fitting_model(st.session_state.model_fittings,
                                                 fit_args=fit_args_,
                                                 plt_type='HeightT')
            col1.pyplot(fig_ht)
            fig_vt, axis_vt = plot_fitting_model(st.session_state.model_fittings,
                                                 fit_args=fit_args_,
                                                 plt_type='SpeedT')
            col2.pyplot(fig_vt)
        elif plt_kinematics_select == 'Long-LatT':
            fit_args_['bounds'] = ([-np.inf, 0, -np.inf], [np.inf, np.inf, np.inf])
            col1, col2 = st.columns(2)
            fig_loT, axis_loT = plot_fitting_model(st.session_state.model_fittings,
                                                   fit_args=fit_args_,
                                                   plt_type='LongT')
            col1.pyplot(fig_loT)
            fig_laT, axis_laT = plot_fitting_model(st.session_state.model_fittings,
                                                   fit_args=fit_args_,
                                                   plt_type='LatT')
            col2.pyplot(fig_laT)
        else:
            fig, axis = plot_fitting_model(st.session_state.model_fittings,
                                           fit_args=fit_args_,
                                           plt_type=plt_kinematics_select)
            st.pyplot(fig)
        st.session_state.model_fittings.kinematics['fit_method'] = fit_args_
    else:
        st.session_state.startup['fitting'] = True

    #############################################################
    # Download Fitting and Figures
    st.sidebar.markdown('## Finalize and save results')
    if 'model_fittings' in st.session_state:
        st.sidebar.download_button('Download Fitting as .json file',
                                   st.session_state.model_fittings.to_json(),
                                   st.session_state.model_fittings.model_id()+'.json')
    else:
        st.sidebar.info('Store a fit to enable this feature.')

    st.sidebar.markdown('---')
    st.markdown('---')
    # footer_text()


if __name__ == '__main__':
    run()
