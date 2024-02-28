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


import astropy.units as u
import numpy as np
import stqdm  # See also https://github.com/tqdm/tqdm
import streamlit as st
from astropy.coordinates import Distance, SkyCoord, SphericalRepresentation
from sunpy.coordinates import frames

from PyThea.callbacks import load_or_delete_fittings
from PyThea.config import (app_styles, config_sliders, selected_bodies,
                           selected_imagers)
from PyThea.geometrical_models import ellipsoid, gcs, spheroid
from PyThea.modules import (date_and_event_selection, final_parameters_gmodel,
                            fitting_and_slider_options_container,
                            fitting_sliders, maps_clims)
from PyThea.sunpy_dev.map.maputils import get_closest
from PyThea.utils import (download_fits, make_figure, model_fittings,
                          plot_bodies, plot_fitting_model,
                          single_imager_maps_process)
from PyThea.version import version


def delete_from_state(vars):
    for var in vars:
        if var in st.session_state:
            del st.session_state[var]


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
                   **Citation**: Cite PyThea as [![https://doi.org/10.5281/zenodo.5713659](https://zenodo.org/badge/DOI/10.5281/zenodo.5713659.svg)](https://doi.org/10.5281/zenodo.5713659) \n
                   """ +
                   f"""
                   **Version**: {version} (latest release [![Version](https://img.shields.io/github/v/release/AthKouloumvakos/PyThea)](https://github.com/AthKouloumvakos/PyThea/releases))
                   """)
    left.image('https://github.com/AthKouloumvakos/PyThea/blob/master/docs/logo/pythea_logo.png?raw=true')
    st.warning('''
           ⚠️ **NOTE: From PyThea >0.8.1 the JSON fitting files will be slightly different from the old ones.** ⚠️
           * Due to a change in the fitting time input, the new fitting files may have slightly different times for the same images (see further information [here](https://github.com/AthKouloumvakos/PyThea/discussions/24)).
           ''')
    st.markdown('---')


def run():
    ###################################################################
    # IMPORTANT: DON'T CHANGE THE POSSITION OF THE COMPONENTS. NEVER! #
    ###################################################################

    #############################################################
    # set page config
    st.set_page_config(page_title='PyThea', page_icon=':rocket:',
                       initial_sidebar_state='expanded')

    #############################################################
    # HTML Styles
    app_styles.apply(st)

    #############################################################
    # Main page information text
    st.title('PyThea: Reconstruct CMEs & Shocks')

    #############################################################
    # Startup Variables
    if 'startup' not in st.session_state:
        st.session_state.startup = {'fitting': True}

    #############################################################
    # Date and Event selection

    if 'date_process' not in st.session_state:
        date_and_event_selection(st)
        st.markdown(""" #### ⏻ Select a day & solar event. """)
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
                                                 kwargs={'vars': ['map', 'map_', 'imagers_list_']})

    with st.sidebar.expander('Processing Options'):
        procoption_container = st.container()
        if 'imagers_list_' not in st.session_state:
            # imagers_list_ is used later when we download or filter the images
            st.session_state['imagers_list_'] = []
        image_mode = procoption_container.selectbox('Map Sequence Processing',
                                                    options=['Running Diff.', 'Base Diff.', 'Plain'], key='image_mode',
                                                    on_change=delete_from_state, kwargs={'vars': ['map', 'imagers_list_']})
        if image_mode in ['Running Diff.', 'Base Diff.']:
            diff_value = procoption_container.number_input('Select Step/Image for Running/Base Diff.', 1, 5, key='diff_value',
                                                           on_change=delete_from_state, kwargs={'vars': ['map', 'imagers_list_']})
        else:
            diff_value = None

    with st.sidebar.expander('Plot/View Options'):
        plotviewopt_container = st.container()
        clip_model = plotviewopt_container.checkbox('Clip plot on image limits', value=True)
        plt_supp_imagers = plotviewopt_container.checkbox('Supplementary Imaging', value=False)
        star_field = plotviewopt_container.checkbox('View Bodies or s/c')
        if star_field:
            bodies_list = plotviewopt_container.multiselect('Select Bodies', options=selected_bodies.bodies_dict.keys(),
                                                            default=['Mercury', 'Venus', 'Jupiter'])

    #############################################################
    # Download and Process the Images
    # This part runs only if the map_ doesn't exits or the session_state.map_ does not contain all the imagers requested
    imager_added = list(set(imagers_list) - set(st.session_state['imagers_list_']))
    imager_removed = list(set(st.session_state['imagers_list_']) - set(imagers_list))
    if 'map_' not in st.session_state or imager_added != []:
        st.session_state.map_ = {} if 'map_' not in st.session_state else st.session_state.map_
        st.session_state.map = {} if 'map' not in st.session_state else st.session_state.map
        st.session_state['imagers_list_'] = imagers_list
        progress_bar = stqdm.stqdm(imager_added, desc='Preparing to Download data')
        for imager in progress_bar:
            progress_bar.desc = f'Downloaded {imager} images from VSO'
            if imager not in st.session_state.map_:
                st.session_state.map_[imager] = download_fits(st.session_state.date_process,
                                                              imager, time_range=imaging_time_range)
                st.session_state.map_[imager] = single_imager_maps_process(st.session_state.map_[imager],
                                                                           **selected_imagers.imager_dict[imager][1],
                                                                           skip='sequence_processing')
            processed_images = single_imager_maps_process(st.session_state.map_[imager],
                                                          skip=['filter', 'prepare'],
                                                          image_mode=image_mode,
                                                          diff_num=diff_value)
            if processed_images != []:
                st.session_state.map[imager] = processed_images
            else:
                st.session_state.imagers_list_.remove(imager)
    if imager_removed != []:
        st.session_state['imagers_list_'] = imagers_list
        for imager in imager_removed:
            st.session_state.map_.pop(imager, None)

    maps_clims(st, st.session_state.imagers_list_)

    if not st.session_state.imagers_list_:
        st.error('No images have been downloaded or processed.')  # TODO: Explain better
        st.stop()

    #############################################################
    # Selection for the primary image to plot
    st.markdown('### Multi-viewpoint imaging')
    col1, col2 = st.columns((1, 3))
    imager_select = col1.selectbox('Select an imager',
                                   options=st.session_state.imagers_list_)

    maps_date = [getattr(maps, 'date_average', None) or getattr(maps, 'date', None) for maps in st.session_state.map[imager_select]]
    if len(maps_date) > 1:
        running_map_date = col2.select_slider('Slide to the image time',
                                              options=maps_date, value=maps_date[0],
                                              key='running_map_date')
    else:
        running_map_date = maps_date[0]

    running_map = st.session_state.map[imager_select][maps_date.index(running_map_date)]

    qmin = st.session_state.map_clim[imager_select][0]  # np.nanquantile(running_map.data, 0.20)
    qmax = st.session_state.map_clim[imager_select][1]  # np.nanquantile(running_map.data, 0.80)
    col1, col2 = st.columns((1, 3))
    cmmin, cpmax = config_sliders.slider_image_pmclims[image_mode]
    clim = plotviewopt_container.slider('Images climits:',
                                        min_value=float(cmmin), max_value=float(cpmax),
                                        value=(float(qmin-10), float(qmax+10)),
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
    fig, axis = make_figure(running_map, image_mode, clim=clim, clip_model=clip_model)
    if st.session_state.plot_mesh_mode == 'Skeleton':
        model.plot(axis, mode='Skeleton')
    elif st.session_state.plot_mesh_mode == 'Full':
        model.plot(axis, mode='Skeleton')
        model.plot(axis, mode='Full')
    elif st.session_state.plot_mesh_mode == 'Surface':
        model.plot(axis, only_surface=True)

    if star_field:
        plot_bodies(axis, bodies_list, running_map)
        axis.legend()
    st.pyplot(fig)

    if plt_supp_imagers and len(st.session_state.imagers_list_) > 2:
        supl_imagers = st.select_slider('Select supplement imagers',
                                        options=st.session_state.imagers_list_,
                                        value=(st.session_state.imagers_list_[1],
                                               st.session_state.imagers_list_[-1]),
                                        key='supl_imagers')
        col1, col2 = st.columns(2)
        fig, axis = make_figure(get_closest(st.session_state.map[supl_imagers[0]],
                                running_map_date), image_mode)
        model.plot(axis, mode='Skeleton')
        col1.pyplot(fig)
        fig, axis = make_figure(get_closest(st.session_state.map[supl_imagers[1]],
                                            running_map_date), image_mode)
        model.plot(axis, mode='Skeleton')
        col2.pyplot(fig)

    #############################################################
    # Store the fitting
    single_fit = model.to_dataframe()
    single_fit['imager'] = imager_select

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
        st.dataframe(single_fit)
        if 'model_fittings' in st.session_state:
            st.markdown('**Stored Fitting Table:**')
            st.dataframe(st.session_state.model_fittings.parameters)
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
                                               options=['All', 'HeightT', 'SpeedT', 'Long-LatT'])

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
