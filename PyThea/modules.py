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
from sunpy.net import attrs as a

from PyThea.callbacks import change_fitting_sliders, change_long_lat_sliders
from PyThea.config import database_dir_default
from PyThea.config.config_sliders import sliders_dict as sd
from PyThea.extensions.hek.utils import plot_hek
from PyThea.extensions.Parker_spirals.utils import plot_parker_spiral
from PyThea.geometrical_models import ellipsoid, gcs, spheroid
from PyThea.utils import (get_hek_flare, make_figure, model_fittings,
                          plot_bodies, plot_solar_reference_lines)


def date_and_event_selection(st):
    """
    Handles date and event selection through Streamlit sidebar components.

    Parameters
    ----------
    st : Streamlit session state object
        Streamlit session state object for managing UI components.

    Returns:
    None
    """
    st.sidebar.markdown('## Date and event selection')
    col1, col2 = st.sidebar.columns(2)

    # Date input for selecting a day to process
    day = col1.date_input('Select a day to process',
                          value=datetime.date.today() - datetime.timedelta(days=1),
                          min_value=datetime.datetime(2010, 1, 1, 0, 0, 0),
                          max_value=datetime.date.today())

    # Selection slider for initialization method: Manual, Event, or File
    initialisation = col2.select_slider('How to initialise?',
                                        options=('Manual', 'Event', 'File'), value='Event')

    if initialisation == 'Event':
        # Retrieve solar events for the selected day
        timerange = a.Time(day, day + datetime.timedelta(days=1))
        flare_list_ = get_hek_flare(timerange)

        # Prepare options for the selectbox based on returned events
        if len(flare_list_) == 0:
            selectbox_list = ['No events returned', ]
        else:
            selectbox_list = []
            for flares in flare_list_:
                fl_ = flares['fl_goescls']
                t_ = flares['event_peaktime'].strftime('%Y-%m-%dT%H:%M:%S')
                selectbox_list.append((f'FL{fl_}|{t_}'))

        selectbox_list.insert(0, 'Select event')

        # Selectbox to choose a specific solar event
        event_selected = st.sidebar.selectbox('Solar Events', options=selectbox_list)

        # If a valid event is selected, update session state and rerun
        if event_selected != 'Select event' and event_selected != 'No events returned':
            st.session_state.date_process = flare_list_['event_peaktime'][selectbox_list.index(event_selected)-1]
            st.session_state.event_selected = event_selected
            st.rerun()
        elif event_selected == 'No events returned':
            st.sidebar.warning('Initiate manually')

    elif initialisation == 'Manual':
        # Manual form for event ID selection and time input
        with st.sidebar.form('my_form'):
            col1, col2 = st.columns(2)
            ev_id = col1.selectbox('Event ID',
                                   options=['Select', 'FL', 'CME', 'SHW'])
            time = col2.time_input('Time', datetime.time(12, 0))

            # If event ID is selected, update session state and rerun
            if ev_id != 'Select':
                st.session_state.date_process = datetime.datetime.combine(day, time)
                t_ = st.session_state.date_process.strftime('%Y-%m-%dT%H:%M:%S')
                st.session_state.event_selected = f'{ev_id}|{t_}'
                st.rerun()
            st.form_submit_button()

    elif initialisation == 'File':
        # File upload for loading fitting file (JSON)
        uploaded_file = st.sidebar.file_uploader('Provide a fitting file', accept_multiple_files=False)

        # Handle uploaded file
        if uploaded_file:
            if uploaded_file.type != 'application/json':
                st.sidebar.error('Unsupported file type')
                st.stop()

            # Load model fittings from JSON and update session state
            st.session_state.model_fittings = model_fittings.load_from_json(uploaded_file)

            st.session_state.event_selected = st.session_state.model_fittings.event_selected
            st.session_state.date_process  = \
                datetime.datetime.strptime(st.session_state.model_fittings.date_process, '%Y-%m-%dT%H:%M:%S.%f')
            st.session_state.geometrical_model = st.session_state.model_fittings.geometrical_model

            st.rerun()

    st.session_state.offline_mode = st.sidebar.checkbox('Work with offline mode active?', value=False)
    if st.session_state.offline_mode:
        st.sidebar.info(f'Fits files will be loaded from PyThea\'s database directory: {database_dir_default}', icon='ℹ️')


def fitting_and_slider_options_container(st):
    """
    Constructs a sidebar container with options for fitting and slider adjustments.

    Parameters
    ----------
    st : Streamlit session state object
        Streamlit session state object for managing UI components.

    Returns:
    None
    """
    container = st.sidebar.container()

    with container.expander('Options'):
        col1, col2 = st.columns(2)

        # Radio button for selecting coordinate system (HGS or HGC)
        col1.radio('Coordinate system', options=['HGS', 'HGC'],
                   on_change=change_long_lat_sliders, args=[st], key='coord_system')

        # Radio button for selecting axes representation based on geometrical model
        if st.session_state.geometrical_model == 'Spheroid':
            col2.radio('Axes representation', options=['h, e, k', 'r, a, b'], key='sliders_repr_mode')
        elif st.session_state.geometrical_model == 'Ellipsoid':
            col2.radio('Axes representation', options=['h, e, k, a', 'r, a, b, c'], key='sliders_repr_mode')
        elif st.session_state.geometrical_model == 'GCS':
            col2.radio('Axes representation', options=['h, a, k, t'], key='sliders_repr_mode')
        else:
            st.error('Unrecognized geometrical model')

        # Selectbox for adjusting mode
        st.selectbox('Adjustments mode',
                     options=['Standard', '<10Rsun', '>10Rsun', '>30Rsun'],
                     key='sliders_adjustment_mode')

        # Select slider for plotting option for fitting mesh
        st.select_slider('Plot option for fitting mesh',
                         options=['No plot', 'Skeleton', 'Full', 'Surface'],
                         value='Skeleton', key='plot_mesh_mode')

        # Checkbox for viewing fitting table
        st.checkbox('View Fitting Table', value=False, key='fitting_table')

        # Checkbox for viewing kinematics plot
        st.checkbox('View Kinematics Plot', value=False, key='plt_kinematics')

        # Load fitting parameters if conditions are met
        if ('fit_args_prime' not in st.session_state and
           'model_fittings' in st.session_state):
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
    """
    Constructs sliders in the Streamlit sidebar for adjusting fitting parameters
    based on the selected geometrical model and representation mode.

    Parameters
    ----------
    st : Streamlit session state object
        Streamlit session state object for managing UI components.

    Returns
    -------
    None
    """
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

    if st.session_state.geometrical_model in options:
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
    """
    Calculate final parameters based on the selected geometrical model and representation mode.

    Parameters
    ----------
    st : Streamlit session state object
        Streamlit session state object for managing UI components.

    Returns
    -------
    tuple
        Tuple of final parameters depending on the geometrical model and representation mode.
    """
    gmodel = st.session_state.geometrical_model
    rmode = st.session_state.sliders_repr_mode

    if gmodel == 'Spheroid':
        if rmode == 'h, e, k':
            rcenter, radaxis, orthoaxis1 = spheroid.rab_from_hek(st.session_state.height * u.R_sun,
                                                                 st.session_state.epsilon,
                                                                 st.session_state.kappa)
        elif rmode == 'r, a, b':
            rcenter, radaxis, orthoaxis1 = (st.session_state.rcenter * u.R_sun,
                                            st.session_state.radaxis * u.R_sun,
                                            st.session_state.orthoaxis1 * u.R_sun)
        return rcenter, radaxis, orthoaxis1

    elif gmodel == 'Ellipsoid':
        if rmode == 'h, e, k, a':
            rcenter, radaxis, orthoaxis1, orthoaxis2 = ellipsoid.rabc_from_heka(st.session_state.height * u.R_sun,
                                                                                st.session_state.epsilon,
                                                                                st.session_state.kappa,
                                                                                st.session_state.alpha)
        elif rmode == 'r, a, b, c':
            rcenter, radaxis, orthoaxis1, orthoaxis2 = (st.session_state.rcenter * u.R_sun,
                                                        st.session_state.radaxis * u.R_sun,
                                                        st.session_state.orthoaxis1 * u.R_sun,
                                                        st.session_state.orthoaxis2 * u.R_sun)
        tilt = st.session_state.tilt * u.degree
        return rcenter, radaxis, orthoaxis1, orthoaxis2, tilt

    elif gmodel == 'GCS':
        if rmode == 'h, a, k, t':
            height = st.session_state.height * u.R_sun
            alpha = st.session_state.alpha * u.degree
            kappa = st.session_state.kappa
            rcenter = gcs.rcenter_(height, alpha, kappa)
        tilt = st.session_state.tilt * u.degree
        return rcenter, height, alpha, kappa, tilt


def figure_streamlit(st, running_map, image_mode, imager, model):
    """
    Generates and plots a figure in Streamlit based on the given parameters and session state settings.

    This function will create a figure with various plotting options, including mesh modes, star fields,
    solar reference lines, Parker spirals, and HEK responses.

    Parameters
    ----------
    st : Streamlit session state object
        The session state object for managing UI components.

    running_map : map object
        The map object to be plotted.

    image_mode : str
        Mode of the image ('Plain' or 'Greys_r').

    imager : str
        The imager used.

    model : model object
        The model to be plotted on the figure.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The generated figure object.

    axis : matplotlib.axes._axes.Axes
        The axis of the generated figure.
    """
    cmap = 'default' if image_mode == 'Plain' else 'Greys_r'

    # Create figure with the specified parameters
    fig, axis = make_figure(running_map,
                            cmap=cmap,
                            clim=st.session_state.map_colormap_limits[imager],
                            clip_model=st.session_state.clip_model,
                            median_filter=st.session_state.images_median_filter)

    # Plot the model based on the selected mesh mode
    plot_mode = st.session_state.plot_mesh_mode
    if plot_mode == 'Skeleton':
        model.plot(axis, mode='Skeleton')
    elif plot_mode == 'Full':
        model.plot(axis, mode='Skeleton')
        model.plot(axis, mode='Full')
    elif plot_mode == 'Surface':
        model.plot(axis, only_surface=True)

    # Optionally plot star field, solar reference lines, and Parker spirals
    if st.session_state.star_field:
        plot_bodies(axis, st.session_state.bodies_list, running_map)
        axis.legend()

    if st.session_state.plot_solar_reference_lines_:
        plot_solar_reference_lines(axis, st.session_state.plot_solar_reference_lines_bodies_list,
                                   running_map, mode=st.session_state.plot_solar_reference_lines_mode)

    if st.session_state.plot_parker_spirals:
        plot_parker_spiral(axis, running_map,
                           st.session_state.mag_bodies_list,
                           sw_speed=st.session_state.sw_speed_select)

    # Optionally plot HEK (Heliophysics Event Knowledgebase) responses
    if st.session_state.plot_hek:
        if 'hek_responses' not in st.session_state:
            st.session_state.hek_responses = {'Active Regions': [], 'Coronal Holes': [], 'Flares': []}
        for mode in st.session_state.hek_list:
            st.session_state.hek_responses[mode] = plot_hek(axis, running_map, mode,
                                                            st.session_state.imaging_time_range,
                                                            hek_responses=st.session_state.hek_responses)

    return fig, axis
