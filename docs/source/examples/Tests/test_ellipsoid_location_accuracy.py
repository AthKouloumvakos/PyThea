r"""
Ellipsoid location on AIA images
--------------------------------

This is an example of how the accuracy of the ellipsoid possition is tested in PyThea.

We use an EUV image when Venus transited in front of the Sun and was observed by SDO/AIA and we oveplot the ellipsoid.
"""

# %%
import astropy.units as u
import matplotlib.pyplot as plt
import sunpy.map
from astropy.coordinates import SkyCoord, solar_system_ephemeris
from sunpy.coordinates import get_body_heliographic_stonyhurst

from PyThea.data.sample_data import aia_sample_data
from PyThea.geometrical_models import ellipsoid
from PyThea.sunpy_dev.map.maputils import prepare_maps

# %%
# Use high-precision ephemeris information from jplephem to have a more accurate position of Venus than the ephemeris provided by astropy.
solar_system_ephemeris.set('de432s')

# %%
# Load the SDO/AIA image from PyThea's sample data and prepare the image correcting for pointing and observer location.
aia_fits_venus = aia_sample_data.fetch('aia_lev1_1600a_2012_06_06t04_07_29_12z_image_lev1.fits')
aiamap = sunpy.map.Map(aia_fits_venus, sequence=True)
aiamap = prepare_maps(aiamap)[0]

# %%
# Define Venus radious and get the sky coordinates at the time of observation.

venus_radius = 6051.8 * u.km  # From here: https://nssdc.gsfc.nasa.gov/planetary/factsheet/venusfact.html
venus = get_body_heliographic_stonyhurst('venus', aiamap.date_average, observer=aiamap.observer_coordinate)
center = SkyCoord(venus)


# %%
# Construct the ellipsoid model. The model is centred at Venus with a radius equal to the planet's radius.

model_shock = ellipsoid(center, venus_radius, venus_radius, venus_radius, 0 * u.degree)

# %%
# Create a submap with Venus at its center.
venus_hpc = venus.transform_to(aiamap.coordinate_frame)
fov = 200 * u.arcsec
bottom_left = SkyCoord(venus_hpc.Tx - fov/2, venus_hpc.Ty - fov/2, frame=aiamap.coordinate_frame)
smap = aiamap.submap(bottom_left, width=fov, height=fov)

# %%
# Make a figure and plot the map and the ellipsoid. The ellipsoid should appear on top of Venus.
fig = plt.figure(dpi=200)
ax = fig.add_subplot(projection=smap)
smap.plot(axes=ax)
smap.draw_limb(axes=ax)
ax.grid(False)
model_shock.plot(ax, mode='Full')
ax.set_xlim([0, smap.data.shape[0]])
ax.set_ylim([0, smap.data.shape[1]])

plt.show()
