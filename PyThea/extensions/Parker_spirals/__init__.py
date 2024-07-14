import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord


def spiral(pos, sw_speed, date, omega=360. * u.degree / (25.38 * 24 * 60 * 60 * u.second)):
    """
    Calculates the Parker spiral coordinates based on the input position, solar wind speed, and date.

    Parameters
    ----------
    pos : astropy.coordinates.SkyCoord
        The initial position in heliographic Carrington coordinates.
    sw_speed : astropy.units.Quantity
        The solar wind speed (should be in velocity units, e.g., km/s).
    date : astropy.time.Time
        The observation time.
    omega : astropy.units.Quantity, optional
        The solar rotation rate in radians per second (default is sidereal period).

    Returns
    -------
    astropy.coordinates.SkyCoord
        The calculated spiral coordinates in heliographic Carrington coordinates.
    """
    r = np.arange(1, pos.radius.to_value(u.R_sun), 0.1)
    r = np.append(r, pos.radius.to_value(u.R_sun)) * u.R_sun
    alpha = omega * ((pos.radius - r) / sw_speed)
    hg_coord = SkyCoord(pos.lon + alpha, pos.lat, r,
                        frame='heliographic_carrington',
                        observer='earth',
                        obstime=date)
    return hg_coord


def footpoint(pos, sw_speed, date, r=1*u.R_sun, omega=360. * u.degree / (25.38 * 24 * 60 * 60 * u.second)):
    """
    Calculates the footpoint location based on the input position, solar wind speed, and date.

    Parameters
    ----------
    pos : astropy.coordinates.SkyCoord
        The initial position in heliographic Carrington coordinates.
    sw_speed : astropy.units.Quantity
        The solar wind speed (should be in velocity units, e.g., km/s).
    date : astropy.time.Time
        The observation time.
    r : astropy.units.Quantity, optional
        The radial distance for the footpoint (default is 1 solar radius).
    omega : astropy.units.Quantity, optional
        The solar rotation rate in radians per second (default is 2.7e-6 rad/s).

    Returns
    -------
    astropy.coordinates.SkyCoord
        The calculated footpoint coordinates in heliographic Carrington coordinates.
    """
    alpha = omega * ((pos.radius - r) / sw_speed)
    hg_coord = SkyCoord(pos.lon + alpha, pos.lat, r,
                        frame='heliographic_carrington',
                        observer='earth', obstime=date)
    return hg_coord
