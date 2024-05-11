r"""
Plot Model Kinematics
---------------------

In this example, we use PyThea's utilities to plot the kinematics of the geometrical model from a fitting file.
"""
from datetime import datetime

import astropy.units as u
# %%
# Import Required Modules
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import Longitude, SkyCoord
from matplotlib.ticker import MultipleLocator
from sunpy.coordinates import (HeliocentricEarthEcliptic, frames,
                               get_body_heliographic_stonyhurst,
                               get_horizons_coord)
from sunpy.time import parse_time

from PyThea.geometrical_models import ellipsoid

# %%
# Define utility functions what are used to plot the planet and spacecraft location, and the Parker spirals connected to each observer.

# This function returns the Parker spiral coordinates for an observer


def spiral(pos, sw_speed):
    omega = 360. * u.degree / (25.38 * 24 * 60 * 60 * u.second)  # rot-angle in rad/sec, sidereal period
    r = np.arange(1, pos.radius.to_value(u.R_sun), 0.1) * u.R_sun
    alpha = omega * ((pos.radius - r) / sw_speed)
    hg_coord = SkyCoord(pos.lon+alpha, pos.lat, r,
                        frame=pos.frame)
    return hg_coord

# A function to convert coordinates in HEE.


def coord_to_heexy(coord, obstime):
    hee_frame = HeliocentricEarthEcliptic(obstime=obstime)
    coord = coord.transform_to(hee_frame)
    coord.representation_type = 'cartesian'
    return coord.y.to_value('AU'), coord.x.to_value('AU')

# A function to make the planet and spacecraft orbit plot.
# This is a modified version of the example provided by SunPy.


def orbit_plot(time, fig, ax):
    # Define the time for the plot as the time when this script is run.
    obstime = parse_time(time)

    hee_frame = HeliocentricEarthEcliptic(obstime=obstime)

    # Define a convenience function to extract the first full orbit from a
    # trajectory, assuming that the trajectory moves in the direction of positive
    # ecliptic longitude.
    def get_first_orbit(coord):
        lon = coord.transform_to(hee_frame).spherical.lon
        shifted = Longitude(lon - lon[0])
        ends = np.flatnonzero(np.diff(shifted) < 0)
        if ends.size > 0:
            return coord[:ends[0]]
        return coord

    # Obtain the locations and trajectories of the various planets and spacecraft.
    # To ensure that each trajectory contains at least one full orbit, we request
    # 700 days for each planet and 1 year for each spacecraft.

    planets = ['Mercury', 'Venus', 'Earth']
    times = obstime + np.arange(700) * u.day
    planet_coords = {planet: get_first_orbit(get_body_heliographic_stonyhurst(planet, times))
                     for planet in planets}

    stereo_a = get_horizons_coord('STEREO-A', obstime)

    missions = ['Parker Solar Probe', 'Solar Orbiter', 'BepiColombo']  # , 'BepiColombo'
    mission_labels = {'Parker Solar Probe': 'PSP', 'Solar Orbiter': 'SolO', 'BepiColombo': 'Bepi'}
    mission_colors = {'Parker Solar Probe': 'purple', 'Solar Orbiter': 'dodgerblue', 'BepiColombo': 'orange'}
    mission_coords = {mission: get_first_orbit(get_horizons_coord(mission, {'start': obstime,
                                                                            'stop': obstime + 1 * u.yr,
                                                                            'step': '6h'}))
                      for mission in missions}

    # Set Matplotlib settings to the desired appearance and initialize the axes.
    ax.set_xlim(-1.15, 1.15)
    ax.set_xlabel('Y (HEE)')
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))

    ax.set_ylim(1.15, -1.15)
    ax.set_ylabel('X (HEE)')
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))

    ax.set_title(f'Orbit Plot: {obstime.strftime("%d-%b-%Y %H:%M UT")}')
    ax.set_aspect('equal')

    # sphinx_gallery_defer_figures

    # Draw the Sun-Earth line.

    ax.plot([0, 0], [0, 2], linestyle='dotted', linewidth=1, color='gray')

    # sphinx_gallery_defer_figures

    # Draw Mercury, Venus, Earth, and Mars, with Earth formatted differently.

    for planet, coord in planet_coords.items():
        ax.plot(*coord_to_heexy(coord, obstime), linestyle='dashed', linewidth=1, color='gray')
        if planet == 'Earth':
            color, markersize, offset = 'green', None, 0.05
        elif planet == 'Venus':
            color, markersize, offset = 'gray', None, -0.25
        else:
            color, markersize, offset = 'gray', None, 0.05

        x, y = coord_to_heexy(coord[0], obstime)
        ax.plot(x, y, 'o', markersize=markersize, color=color)
        ax.text(x + offset, y, planet, color=color, clip_on=True)

    # sphinx_gallery_defer_figures

    # Draw the STEREO spacecraft (without orbits), as well as Sun-STEREO lines.

    for stereo, label, color in [(stereo_a, 'STA', 'red'), ]:  # , (stereo_b, 'B', 'blue')
        x, y = coord_to_heexy(stereo, obstime)
        ax.plot([0, 5*x], [0, 5*y], linestyle='dotted', linewidth=1, color='gray')
        ax.plot(x, y, 'o', color=color)
        ax.text(x + 0.05, y, label, color=color, fontsize=12, clip_on=True)

    # sphinx_gallery_defer_figures

    # Draw the Sun, which is at the origin by definition.

    ax.plot(0, 0, 'o', markersize=5, color='orange', zorder=100)
    # ax.text(0.09, 0, 'Sun', color='orange')

    # sphinx_gallery_defer_figures

    # Finally, draw the various spacecraft, with Solar Orbiter colored differently.

    for mission, coord in mission_coords.items():
        ax.plot(*coord_to_heexy(coord, obstime), linestyle='dashed', color=mission_colors[mission])

        x, y = coord_to_heexy(coord[0], obstime)
        ax.plot(x, y, 'o', color=mission_colors[mission])

        spiral_coord = spiral(coord[0], 350*(u.km/u.second))
        ax.plot(*coord_to_heexy(spiral_coord, obstime), linestyle='-', linewidth=1.5, color=mission_colors[mission])

        ax.text(x + 0.05, y, mission_labels[mission], color=mission_colors[mission], clip_on=True)

    return fig, ax


# %%
# Create an ellipsoid providing an observation time, the ellipsoid center coordinates, and the geomertical parameters.
obstime = datetime.strptime('2022-09-05T16:30:00', '%Y-%m-%dT%H:%M:%S')
center = SkyCoord(135*u.degree, 0*u.degree, 5*u.R_sun, obstime=obstime,
                  observer='earth', frame=frames.HeliographicStonyhurst)
model_shock = ellipsoid(center, 10*u.R_sun, 10*u.R_sun, 10*u.R_sun, 0 * u.degree)

# %%
# Extract the ellipsoid mesh coordinates and convert them to HEE frame.
ellipsoid_coordinates = model_shock.coordinates
x, y = coord_to_heexy(ellipsoid_coordinates, parse_time(obstime))

# %%
# Get the the center coordinates of the ellipsoid in the HEE frame.
x_c, y_c = coord_to_heexy(model_shock.center, obstime)

# %%
# Create a figure and plot the planets, spacecraft, and the ellipsoid.
fig, axis = plt.subplots(figsize=(6, 6), dpi=200)
fig, axis = orbit_plot(obstime, fig=fig, ax=axis)
axis.plot(x, y, '-', color='tab:red')
arrow_properties = dict(facecolor='black', edgecolor='black', arrowstyle='->', shrinkA=0, linewidth=1)
axis.annotate('', xy=(10*x_c, 10*y_c), xytext=(0, 0), arrowprops=arrow_properties)
axis.text(10*x_c - 0.01, 10*y_c, 'Shock', color='black', clip_on=True)

fig.tight_layout()

plt.show()
# %%
