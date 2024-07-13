import astropy.units as u
from sunpy.coordinates import frames


def load_or_delete_fittings(st):
    """
    Load or delete fitting parameters from the model fittings based on user selection.

    Parameters
    ----------
    st : Streamlit session state object
        Streamlit session state object for managing UI components.

    Returns
    -------
    None
    """
    selected_row = str(st.session_state.fitting_select)
    dataframe = st.session_state.model_fittings.parameters

    if st.session_state.fit_action == 'Select':
        pass

    elif st.session_state.fit_action == 'Load':
        if st.session_state.coord_system == 'HGS':
            st.session_state.longit = float(dataframe.loc[selected_row, 'hgln'])
            st.session_state.latitu = float(dataframe.loc[selected_row, 'hglt'])
        elif st.session_state.coord_system == 'HGC':
            st.session_state.longit = float(dataframe.loc[selected_row, 'crln'])
            st.session_state.latitu = float(dataframe.loc[selected_row, 'crlt'])

        if st.session_state.model_fittings.geometrical_model == 'Spheroid':
            keys = ['height', 'kappa', 'epsilon', 'rcenter', 'radaxis', 'orthoaxis1']
        elif st.session_state.model_fittings.geometrical_model == 'Ellipsoid':
            keys = ['height', 'kappa', 'epsilon', 'alpha', 'rcenter', 'radaxis', 'orthoaxis1', 'orthoaxis2', 'tilt']
        elif st.session_state.model_fittings.geometrical_model == 'GCS':
            keys = ['height', 'alpha', 'kappa', 'tilt']

        for key in keys:
            if key in st.session_state:
                # del st.session_state[key]
                st.session_state[key] = float(dataframe.loc[selected_row, key])
        del st.session_state.fit_action

    elif st.session_state.fit_action == 'Delete':
        st.session_state.model_fittings.parameters = dataframe.drop(st.session_state.fitting_select)
        del st.session_state.fitting_select
        if len(st.session_state.model_fittings.parameters) < 1:
            del st.session_state.model_fittings
        del st.session_state.fit_action


def change_long_lat_sliders(st):
    """
    Update longitude and latitude sliders based on the selected coordinate system.

    Parameters
    ----------
    st : Streamlit session state object
        Streamlit session state object for managing UI components.

    Returns
    -------
    None
    """
    if st.session_state.coord_system == 'HGS':
        center_ = st.session_state.center.transform_to(frames.HeliographicStonyhurst)
    elif st.session_state.coord_system == 'HGC':
        center_ = st.session_state.center.transform_to(frames.HeliographicCarrington(observer='Earth',
                                                                                     obstime=st.session_state.center.obstime))
    st.session_state.longit = float(center_.lon.to_value(u.degree))
    st.session_state.latitu = float(center_.lat.to_value(u.degree))


def change_fitting_sliders(st):
    """
    Update fitting sliders based on the selected fitting method.

    Parameters
    ----------
    st : Streamlit session state object
        Streamlit session state object for managing UI components.

    Returns
    -------
    None
    """
    st.session_state.startup['fitting'] = False
    fit_opt = st.session_state.fit_args_prime
    st.session_state.fit_mode = fit_opt['type']

    if fit_opt['type'] == 'polynomial':
        st.session_state.polyfit_order = fit_opt['order']
    elif fit_opt['type'] == 'spline':
        st.session_state.splinefit_order = fit_opt['order']
        st.session_state.splinefit_smooth = fit_opt['smooth']
    elif fit_opt['type'] == 'custom':
        st.session_state.fitcustexpres_select = fit_opt['expression']
        st.session_state.splinefit_order = fit_opt['order']
