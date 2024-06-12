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

import astropy.units as u

from PyThea.callbacks import change_fitting_sliders, change_long_lat_sliders
from PyThea.config.config_sliders import sliders_dict as sd
from PyThea.geometrical_models import ellipsoid, gcs, spheroid
from PyThea.utils import (get_hek_flare, make_figure, model_fittings,
                          plot_bodies, plot_solar_reference_lines)


def date_and_event_selection(st):
    st.sidebar.markdown('## Date and event selection')
    col1, col2 = st.sidebar.columns(2)
    day = col1.date_input('Select a day to process',
                          value=datetime.date.today() - datetime.timedelta(days=1),
                          min_value=datetime.datetime(2010, 1, 1, 0, 0, 0),
                          max_value=datetime.date.today())
    initialisation = col2.select_slider('How to initialise?',
                                        options=('Manual', 'Event', 'File'), value='Event')
    if initialisation == 'Event':
        selectbox_list, flare_list_ = get_hek_flare(day)
        selectbox_list.insert(0, 'Select event')
        event_selected = st.sidebar.selectbox('Solar Events', options=selectbox_list)
        if event_selected != 'Select event' and event_selected != 'No events returned':
            st.session_state.date_process = flare_list_['event_peaktime'][selectbox_list.index(event_selected)-1]
            st.session_state.event_selected = event_selected
            st.rerun()
        elif event_selected == 'No events returned':
            st.sidebar.warning('Initiate manually')
    elif initialisation == 'Manual':
        with st.sidebar.form('my_form'):
            col1, col2 = st.columns(2)
            ev_id = col1.selectbox('Event ID',
                                   options=['Select', 'FL', 'CME', 'SHW'])
            time = col2.time_input('Time', datetime.time(12, 0))
            if ev_id != 'Select':
                st.session_state.date_process = datetime.datetime.combine(day, time)
                t_ = st.session_state.date_process.strftime('%Y-%m-%dT%H:%M:%S')
                st.session_state.event_selected = f'{ev_id}|{t_}'
                st.rerun()
            st.form_submit_button()
    elif initialisation == 'File':
        uploaded_file = st.sidebar.file_uploader('Provide a fitting file', accept_multiple_files=False)
        if uploaded_file:
            if uploaded_file.type != 'application/json':
                st.sidebar.error('Unsupported file type')
                st.stop()

            st.session_state.model_fittings = model_fittings.load_from_json(uploaded_file)

            st.session_state.event_selected = st.session_state.model_fittings.event_selected
            st.session_state.date_process  = \
                datetime.datetime.strptime(st.session_state.model_fittings.date_process, '%Y-%m-%dT%H:%M:%S.%f')
            st.session_state.geometrical_model = st.session_state.model_fittings.geometrical_model

            st.rerun()


def fitting_and_slider_options_container(st):
    container = st.sidebar.container()

    with container.expander('Options'):
        col1, col2 = st.columns(2)
        col1.radio('Coordinate system', options=['HGS', 'HGC'],
                   on_change=change_long_lat_sliders, args=[st], key='coord_system')

        if st.session_state.geometrical_model == 'Spheroid':
            col2.radio('Axes representation', options=['h, e, k', 'r, a, b'], key='sliders_repr_mode')
        elif st.session_state.geometrical_model == 'Ellipsoid':
            col2.radio('Axes representation', options=['h, e, k, a', 'r, a, b, c'], key='sliders_repr_mode')
        elif st.session_state.geometrical_model == 'GCS':
            col2.radio('Axes representation', options=['h, a, k, t'], key='sliders_repr_mode')
        else:
            st.error('Unrecognized geometrical model')

        st.selectbox('Adjustments mode',
                     options=['Standard', '<10Rsun', '>10Rsun', '>30Rsun'],
                     key='sliders_adjustment_mode')

        st.select_slider('Plot option for fitting mesh',
                         options=['No plot', 'Skeleton', 'Full', 'Surface'],
                         value='Skeleton', key='plot_mesh_mode')

        st.checkbox('View Fitting Table', value=False, key='fitting_table')

        st.checkbox('View Kinematics Plot', value=False, key='plt_kinematics')

        if 'fit_args_prime' not in st.session_state and\
           'model_fittings' in st.session_state:
            if st.session_state.model_fittings.kinematics['fit_method'] is None:
                # TODO: When initialization is 'Event' or 'Manual' the st.session_state.model_fittings.kinematics
                # does not take a default value. Explore if I can move this earlier.
                st.session_state.model_fittings.kinematics = {'fit_method': {'type': 'polynomial', 'order': 1}}
            st.session_state.fit_args_prime = st.session_state.model_fittings.kinematics['fit_method']

        if (('model_fittings' in st.session_state) and
           (st.session_state.plt_kinematics) and
           st.button('Load Fitting Parameters', key='load_kinematics_param')) or\
           (st.session_state.startup['fitting'] and
           'fit_args_prime' in st.session_state):
            change_fitting_sliders(st)


def fitting_sliders(st):
    options = {
        'Spheroid': {'h, e, k': ['height', 'kappa', 'epsilon'],
                     'r, a, b': ['rcenter', 'radaxis', 'orthoaxis1']
                     },
        'Ellipsoid': {'h, e, k, a': ['height', 'kappa', 'epsilon', 'alpha', 'tilt'],
                      'r, a, b, c': ['rcenter', 'radaxis', 'orthoaxis1', 'orthoaxis2', 'tilt']
                      },
        'GCS': {'h, a, k, t': ['height', 'alpha', 'kappa', 'tilt'],
                },
    }

    adjustments = st.session_state.sliders_adjustment_mode

    if st.session_state.geometrical_model == 'Spheroid' or \
       st.session_state.geometrical_model == 'Ellipsoid' or \
       st.session_state.geometrical_model == 'GCS':
        gmodel = st.session_state.geometrical_model
        rmode = st.session_state.sliders_repr_mode
        sliders = options[gmodel][rmode]
        for slider in sliders:
            st.sidebar.slider(f'{slider} {sd[gmodel][slider]["unit"]}:',  # Add units here
                              min_value=sd[gmodel][slider][adjustments]['min'],
                              max_value=sd[gmodel][slider][adjustments]['max'],
                              value=sd[gmodel][slider][adjustments]['def'],
                              step=sd[gmodel][slider][adjustments]['step'], key=slider)  # * u.R_sun


def final_parameters_gmodel(st):
    if st.session_state.geometrical_model == 'Spheroid':
        if st.session_state.sliders_repr_mode == 'h, e, k':
            rcenter, radaxis, orthoaxis1 = spheroid.rab_from_hek(st.session_state.height * u.R_sun,
                                                                 st.session_state.epsilon,
                                                                 st.session_state.kappa)
        elif st.session_state.sliders_repr_mode == 'r, a, b':
            rcenter, radaxis, orthoaxis1 = st.session_state.rcenter * u.R_sun, \
                st.session_state.radaxis * u.R_sun, \
                st.session_state.orthoaxis1 * u.R_sun
        return rcenter, radaxis, orthoaxis1
    elif st.session_state.geometrical_model == 'Ellipsoid':
        if st.session_state.sliders_repr_mode == 'h, e, k, a':
            rcenter, radaxis, orthoaxis1, orthoaxis2 = ellipsoid.rabc_from_heka(st.session_state.height * u.R_sun,
                                                                                st.session_state.epsilon,
                                                                                st.session_state.kappa,
                                                                                st.session_state.alpha)
        elif st.session_state.sliders_repr_mode == 'r, a, b, c':
            rcenter, radaxis, orthoaxis1, orthoaxis2 = st.session_state.rcenter * u.R_sun, \
                st.session_state.radaxis * u.R_sun, \
                st.session_state.orthoaxis1 * u.R_sun, \
                st.session_state.orthoaxis2 * u.R_sun
        tilt = st.session_state.tilt * u.degree
        return rcenter, radaxis, orthoaxis1, orthoaxis2, tilt
    elif st.session_state.geometrical_model == 'GCS':
        if st.session_state.sliders_repr_mode == 'h, a, k, t':
            height = st.session_state.height * u.R_sun
            alpha = st.session_state.alpha * u.degree
            kappa = st.session_state.kappa
            rcenter = gcs.rcenter_(height, alpha, kappa)
        tilt = st.session_state.tilt * u.degree
        return rcenter, height, alpha, kappa, tilt


def figure_streamlit(st, running_map, image_mode, imager, model):
    if image_mode == 'Plain':
        cmap = 'default'
    fig, axis = make_figure(running_map,
                            cmap=cmap,
                            clim=st.session_state.map_colormap_limits[imager],
                            clip_model=st.session_state.clip_model,
                            median_filter=st.session_state.images_median_filter)

    if st.session_state.plot_mesh_mode == 'Skeleton':
        model.plot(axis, mode='Skeleton')
    elif st.session_state.plot_mesh_mode == 'Full':
        model.plot(axis, mode='Skeleton')
        model.plot(axis, mode='Full')
    elif st.session_state.plot_mesh_mode == 'Surface':
        model.plot(axis, only_surface=True)

    if st.session_state.star_field:
        plot_bodies(axis, st.session_state.bodies_list, running_map)
        axis.legend()

    if st.session_state.plot_solar_reference_lines_:
        plot_solar_reference_lines(axis, st.session_state.plot_solar_reference_lines_bodies_list,
                                   running_map, mode=st.session_state.plot_solar_reference_lines_mode)

    return fig, axis
