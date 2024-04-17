import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

omega = 360. * u.degree / (25.38 * 24 * 60 * 60 * u.second)  # rot-angle in rad/sec, sidereal period


def spiral(pos, sw_speed, date):
    r = np.arange(1, pos.radius.to_value(u.R_sun), 0.1) * u.R_sun
    alpha = omega * ((pos.radius - r) / sw_speed)
    hg_coord = SkyCoord(pos.lon+alpha, pos.lat, r,
                        frame='heliographic_carrington',
                        observer='earth',
                        obstime=date)
    return hg_coord


def footpoint(pos, sw_speed, date, r=1*u.R_sun):
    alpha = omega * ((pos.radius - r) / sw_speed)
    # The location of the footpoint:
    hg_coord = SkyCoord(pos.lon+alpha, pos.lat, r,
                        frame='heliographic_carrington',
                        observer='earth', obstime=date)
    return hg_coord
