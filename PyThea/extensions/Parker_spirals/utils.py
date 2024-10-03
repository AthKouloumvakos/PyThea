import astropy.units as u
from sunpy.coordinates import frames, get_horizons_coord

from PyThea.config.selected_bodies import bodies_dict
from PyThea.extensions import Parker_spirals


def plot_parker_spiral(axis, map, bodies, sw_speed=350 * (u.km / u.second)):
    """
    Plot Parker spiral for specified bodies on the given map axis.

    Parameters
    ----------
    axis : matplotlib.axes.Axes
        The axis on which to plot the Parker spiral.
    map : sunpy.map.Map
        The solar map on which the Parker spiral will be plotted.
    bodies : list of str
        List of bodies for which to plot the Parker spiral.
    sw_speed : Quantity, optional
        Solar wind speed for the Parker spiral (default: 350 km/s).

    Returns
    -------
    None
    """
    for body in bodies:
        pos = get_horizons_coord(bodies_dict[body][0], map.date_average, 'id')
        pos = pos.transform_to(frames.HeliographicCarrington(observer='Earth', obstime=map.date_average))

        spiral_coord = Parker_spirals.spiral(pos, sw_speed[body], map.date_average)

        hpc_frame = frames.Helioprojective(observer=map.observer_coordinate, obstime=map.date_average)

        spiral_coord_visible_index = spiral_coord.transform_to(hpc_frame).is_visible()

        axis.plot_coord(spiral_coord[spiral_coord_visible_index], markersize=0, linewidth=1, color=bodies_dict[body][1])
        if spiral_coord_visible_index[1]:
            axis.plot_coord(spiral_coord[0], markersize=16, linewidth=10, marker='+', color=bodies_dict[body][1])
        else:
            axis.plot_coord(spiral_coord[0], markersize=16, linewidth=10, marker='+', color='white')

        axis.plot_coord(spiral_coord[0], markersize=3, linewidth=4, marker='s', color=bodies_dict[body][1])

        axis.plot_coord(spiral_coord[-1], markersize=8, linewidth=6, marker='x', color=bodies_dict[body][1])
