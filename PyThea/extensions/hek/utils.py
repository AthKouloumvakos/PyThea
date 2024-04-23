import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
from sunpy.coordinates import frames
from sunpy.net import attrs as a
from sunpy.net import hek
from sunpy.physics.differential_rotation import solar_rotate_coordinate
from sunpy.time import parse_time


def plot_hek(axis, map, mode, time_range=[-2, 2], hek_responses=None):
    """
    Plots HEK (Heliophysics Event Knowledgebase) events on a given axis.

    This function fetches and plots HEK events, such as Active Regions, Coronal Holes, or Flares, on a provided
    matplotlib axis using data from a SunPy map object.

    Parameters
    ----------
    axis : matplotlib.axes._subplots.AxesSubplot
        The matplotlib axis to plot on.
    map : sunpy.map.Map
        The SunPy map object containing the data.
    mode : str
        The type of HEK events to plot ('Active Regions', 'Coronal Holes', or 'Flares').
    time_range : list, optional
        Time range in hours around the map date to search for events (default is [-2, 2]).
    hek_responses : dict, optional
        Pre-fetched HEK responses (default is None).

    Returns
    -------
    responses: list
        The list of HEK event responses.
    """
    hek_client = hek.HEKClient()
    start_time = map.date_average + TimeDelta(time_range[0] * u.hour)
    end_time = map.date_average + TimeDelta(time_range[1] * u.hour)

    if mode == 'Active Regions':
        if hek_responses is None or not hek_responses['Active Regions']:
            responses = hek_client.search(a.Time(start_time, end_time),
                                          a.hek.AR, a.hek.FRM.Name == 'HMI SHARP')

            responses.keep_columns(['ar_noaanum', 'ar_mcintoshcls', 'ar_mtwilsoncls', 'hgs_x', 'hgs_y', 'event_starttime'])
            indx = [i for i, x in enumerate(responses['ar_noaanum']) if x is None]
            responses.remove_rows(indx)
            responses = responses[['ar_noaanum', 'ar_mcintoshcls', 'ar_mtwilsoncls', 'hgs_x', 'hgs_y', 'event_starttime']]
            dates = Time(responses['event_starttime'])
            time_diff = np.abs(dates - map.date_average)
            time_diff_seconds = np.array([td.sec for td in time_diff])
            if len(np.unique(time_diff_seconds)) > 1:
                responses.remove_rows(np.where(time_diff_seconds == np.max(time_diff_seconds)))

            responses_ = hek_client.search(a.Time(start_time, end_time),
                                           a.hek.AR, a.hek.FRM.Name == 'NOAA SWPC Observer')
            responses_.keep_columns(['ar_noaanum', 'ar_mcintoshcls', 'ar_mtwilsoncls', 'hgs_x', 'hgs_y'])

            responses['ar_mcintoshcls'] = np.full(responses['ar_noaanum'].shape[0], None)
            responses['ar_mtwilsoncls'] = np.full(responses['ar_noaanum'].shape[0], None)

            for rn, clm, clw in zip(responses_['ar_noaanum'], responses_['ar_mcintoshcls'], responses_['ar_mtwilsoncls']):
                i = np.argwhere(responses['ar_noaanum'] == rn)
                responses['ar_mcintoshcls'][i] = clm
                responses['ar_mtwilsoncls'][i] = clw
            print(responses)
        else:
            responses = hek_responses['Active Regions']

        for i, response in enumerate(responses):
            axis.annotate(i, (float(response['hgs_x']), float(response['hgs_y'])),
                          xytext=(float(response['hgs_x']), float(response['hgs_y'])),
                          xycoords=axis.get_transform('heliographic_stonyhurst'),
                          backgroundcolor='none',
                          color='tab:blue',
                          fontsize=8,
                          horizontalalignment='center', verticalalignment='center')
    elif mode == 'Coronal Holes':
        if hek_responses is None or not hek_responses['Coronal Holes']:
            responses = hek_client.search(a.Time(start_time, end_time),
                                          a.hek.CH, a.hek.FRM.Name == 'SPoCA')
        else:
            responses = hek_responses['Coronal Holes']

        for response in responses:
            p1 = response['hpc_boundcc'][9:-2]
            p2 = p1.split(',')
            p3 = [v.split(' ') for v in p2]
            ch_date = parse_time(response['event_starttime'])
            ch_boundary = SkyCoord(
                [(float(v[0]), float(v[1])) * u.arcsec for v in p3],
                obstime=ch_date, observer='earth',
                frame=frames.Helioprojective)
            rotated_ch_boundary = solar_rotate_coordinate(ch_boundary, time=map.date_average)
            axis.plot_coord(rotated_ch_boundary, color='c')
    elif mode == 'Flares':
        if hek_responses is None or not hek_responses['Flares']:
            responses = hek_client.search(a.Time(start_time, end_time),
                                          a.hek.FL, a.hek.FRM.Name == 'SWPC')
            responses.keep_columns(['event_starttime', 'event_peaktime', 'fl_goescls', 'hgs_x', 'hgs_y', 'ar_noaanum'])
            responses.remove_rows(np.where((responses['hgs_x'] == 0) & (responses['hgs_y'] == 0)))
            responses = responses[['event_starttime', 'event_peaktime', 'fl_goescls', 'hgs_x', 'hgs_y', 'ar_noaanum']]
        else:
            responses = hek_responses['Flares']

        for i, response in enumerate(responses):
            axis.annotate(i, (float(response['hgs_x']), float(response['hgs_y'])),
                          xytext=(float(response['hgs_x']), float(response['hgs_y'])),
                          xycoords=axis.get_transform('heliographic_stonyhurst'),
                          backgroundcolor='none',
                          color='tab:red',
                          fontsize=8,
                          horizontalalignment='center', verticalalignment='center')
    else:
        raise ValueError(f"Unsupported mode: {mode}. Use 'Coronal Holes' or 'Flares'.")

    return responses
