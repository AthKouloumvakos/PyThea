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

import operator
import re

import astropy.units as u
import numpy as np
import pandas as pd
import seaborn as sns
from astropy.coordinates import (CartesianRepresentation, Distance, SkyCoord,
                                 SphericalRepresentation, concatenate)
from numpy.linalg import norm
from scipy.spatial.transform import Rotation
from sunpy.coordinates import frames


def sphere(n=20):
    """
    Returns the (x,y,z) coordinates of a unit sphere centred at the origin.
    The coordinates are (n+1)-by-(n+1) matrices.
    """

    theta = np.linspace(-np.pi, np.pi, n+1)
    phi = np.linspace(-1*np.pi/2, 1 * np.pi/2, n+1)
    cosphi = np.cos(phi)
    sintheta = np.sin(theta)

    x =  np.outer(cosphi, np.cos(theta))
    y =  np.outer(cosphi, sintheta)
    z =  np.outer(np.sin(phi), np.ones_like(theta))

    return [x, y, z]


class spheroid:
    '''
    A class for the spheroid model.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The coordinates of the spheroid center.
    rcenter : `~astropy.units.Quantity`
        The radial distance of the spheroid center from the solar center.
    radaxis : `~astropy.units.Quantity`
        The spheroid's first (radial) semi-axis length.
    orthoaxis1 : `~astropy.units.Quantity`
        The spheroid's second semi-axis length (orthogonal axis wrt the radial axis).
    height : `~astropy.units.Quantity`
        The radial distance of the spheroid apex from the solar center (e.g., rcenter+radaxis).
    kappa : `~astropy.units.Quantity`
        Spheroid self-similar constant, defined as the ratio of the apex height to the length
        of the orthoaxis1. (This is proportional to the aspect ratio of the spheroid)
    epsilon : `~astropy.units.Quantity`
        The eccentricity of the spheroid.

    Notes
    -----
    Details about spheroids can be found here: https://en.wikipedia.org/wiki/Spheroid
    The second and third semi-axis are equal in the spheroid. If all the axis are equal we have a sphere.
    Tilt defaults to zero, since any rotation along the radial axis is trivial for the spheroid.
    '''
    @u.quantity_input
    def __init__(self, center, radaxis: (u.R_sun), orthoaxis1: (u.R_sun), n=40):
        self.center = center.transform_to(frames.HeliographicStonyhurst)  # This is the SkyCoord of the ellipsoid center
        self.rcenter = center.radius
        self.radaxis = radaxis
        self.orthoaxis1 = orthoaxis1

        self.height = radaxis + self.rcenter
        self.kappa = orthoaxis1 / (self.height-1.*u.R_sun)

        if radaxis < orthoaxis1:
            self.epsilon = -1. * np.sqrt(1. - (radaxis/orthoaxis1)**2)
        elif radaxis > orthoaxis1:
            self.epsilon = 1. * np.sqrt(1. - (orthoaxis1/radaxis)**2)
        elif radaxis == orthoaxis1:
            self.epsilon = 0.

        self.tilt = 0 * u.degree
        self.n = n

    @property
    def coordinates(self):
        """
        Returns the coordinates of the spheroid cloud of points as `~astropy.coordinates.SkyCoord`.
        """
        [x__, y__, z__] = sphere(self.n)

        x_ = self.radaxis * x__ + self.rcenter
        y_ = self.orthoaxis1 * y__
        z_ = self.orthoaxis1 * z__

        x, y, z = self.rotate(x_, y_, z_)

        return SkyCoord(CartesianRepresentation(x, y, z),
                        frame=frames.HeliographicStonyhurst,
                        observer=self.center.observer,
                        obstime=self.center.obstime)

    def intersecting_curve(self):
        """
        Returns the coordinates of intersection of the spheroid with a unit sphere
         centred at the origin as `~astropy.coordinates.SkyCoord`.
         More information here: https://mathworld.wolfram.com/Sphere-SphereIntersection.html
        """
        d = self.rcenter
        r = self.radaxis
        R = 1 * u.R_sun

        step = 40
        theta = np.linspace(0, 2*np.pi, step)
        x_ = np.ones_like(theta) * (d**2 - r**2 + R**2) / (2 * d)
        # y_**2 - z_**2 = R**2 - x_**2
        alpha = np.sqrt(4 * d**2 * R**2 - (d**2 - r**2 + R**2)**2) / (2 * d)
        y_ = alpha * np.sin(theta)
        z_ = alpha * np.cos(theta)
        x, y, z = self.rotate(x_, y_, z_)

        return SkyCoord(CartesianRepresentation(x, y, z),
                        frame=frames.HeliographicStonyhurst,
                        observer=self.center.observer,
                        obstime=self.center.obstime)

    @staticmethod
    @u.quantity_input
    def hek_from_rab(rc: (u.R_sun), radaxis: (u.R_sun), orthoaxis1: (u.R_sun)):
        '''
        Returns the semi-axis length from the height, epsilon, and kappa parameters
        '''
        height = rc + radaxis
        kappa = orthoaxis1 / (height - 1 * u.R_sun)

        if radaxis < orthoaxis1:
            epsilon = -1. * np.sqrt(1 - (radaxis / orthoaxis1)**2)
        elif radaxis > orthoaxis1:
            epsilon = np.sqrt(1. - (orthoaxis1 / radaxis)**2)
        elif radaxis == orthoaxis1:
            epsilon = 0.
        else:
            height, kappa, epsilon = np.nan * u.R_sun, np.nan, np.nan

        return height, epsilon, kappa

    @staticmethod
    @u.quantity_input
    def rab_from_hek(height: (u.R_sun), epsilon, kappa):
        '''
        Returns the height, epsilon, and kappa parameters from the semi-axis length
        '''
        s = kappa * (height - 1. * u.R_sun)
        if epsilon < 0:
            d = s * np.sqrt(1. - epsilon**2)
        elif epsilon > 0:
            d = s / np.sqrt(1. - epsilon**2)
        elif epsilon == 0:
            d = s
        else:
            rc, radaxis, orthoaxis1 = np.nan * u.R_sun, np.nan * u.R_sun, np.nan * u.R_sun
            return rc, radaxis, orthoaxis1
        radaxis = d
        orthoaxis1 = s
        rc = height - radaxis

        return rc, radaxis, orthoaxis1

    @property
    def apex(self):
        """
        Returns the coordinates of the spheroid apex as `~astropy.coordinates.SkyCoord`.
        """
        rapex = self.center.radius + self.radaxis
        Spher_rep = SphericalRepresentation(self.center.lon, self.center.lat,
                                            Distance(np.abs(rapex)))
        return SkyCoord(Spher_rep.to_cartesian(),
                        frame=frames.HeliographicStonyhurst,
                        obstime=self.center.obstime,
                        observer=self.center.observer)

    @property
    def base(self):
        """
        Returns the coordinates of the spheroid base as `~astropy.coordinates.SkyCoord`.
        """
        rbase = (-self.center.radius + self.radaxis)
        Spher_rep = SphericalRepresentation(self.center.lon, self.center.lat,
                                            Distance(np.abs(rbase)))
        return SkyCoord(-np.sign(rbase) * Spher_rep.to_cartesian(),
                        frame=frames.HeliographicStonyhurst,
                        obstime=self.center.obstime,
                        observer=self.center.observer)

    def rotate(self, x_, y_, z_):
        """
        Rotates the (x,y,z) coordinates from the geometrical model's coordinate system
         to the Heliospheric coordinates.

        Returns the coordinates as `~astropy.coordinates.SkyCoord`.
        """
        Longc = self.center.lon.to_value(u.rad)
        Latc  = self.center.lat.to_value(u.rad)
        tilt = self.tilt.to_value(u.rad)

        v = np.transpose([x_.flatten(), y_.flatten(), z_.flatten()])
        v = Rotation.from_euler('xyz', [tilt, -Latc, Longc]).apply(v)

        x = np.reshape(v[:, 0], x_.shape) * x_.unit
        y = np.reshape(v[:, 1], y_.shape) * y_.unit
        z = np.reshape(v[:, 2], z_.shape) * z_.unit
        return x, y, z

    def plot(self, axis, mode='Skeleton', only_surface=False):
        '''
        Plots the spheroid model at the provided axes
        '''
        if mode == 'Full':
            axis.plot_coord(self.coordinates, color='red', linestyle='-', linewidth=0.2)
        elif (mode == 'Skeleton') and (only_surface is False):
            lw = 1.
            n = self.n
            palete = sns.color_palette('colorblind')

            # Plot the part of the spheroid mesh but filter lines bellow the solar surface
            axis.plot_coord(self.coordinates[:, int(n/2)], color=palete[0], linestyle='-', linewidth=lw)
            my_plot_coord(self.coordinates[:, 0], axis, color=palete[6], linestyle='-', linewidth=lw)
            axis.plot_coord(self.coordinates[0:int(n/2)+1, int(n/4)], color=palete[1], linestyle='-', linewidth=lw)
            axis.plot_coord(self.coordinates[int(n/2):n+1, int(n/4)], color=palete[1], linestyle='-', linewidth=lw)
            axis.plot_coord(self.coordinates[0:int(n/2)+1, int(3*n/4)], color=palete[1], linestyle='-', linewidth=lw)
            axis.plot_coord(self.coordinates[int(n/2):n+1, int(3*n/4)], color=palete[1], linestyle='-', linewidth=lw)
            my_plot_coord(self.coordinates[int(n/2), 0:int(n/4)+1], axis, color=palete[2], linestyle='-', linewidth=lw)
            my_plot_coord(self.coordinates[int(n/2), int(3*n/4):n+1], axis, color=palete[2], linestyle='-', linewidth=lw)
            axis.plot_coord(self.coordinates[int(n/2), int(n/4):int(3*n/4)+1], color=palete[3], linestyle='-', linewidth=lw)

        if (mode == 'Skeleton') or (only_surface is True):
            # Plot the part of the intersecting curve at the solar surface
            d_val = self.rcenter - self.radaxis
            if d_val > -1*u.R_sun and d_val < 1*u.R_sun:
                coords = self.intersecting_curve()
                axes_frame = axis._transform_pixel2world.frame_out
                coord_in_axes = coords.transform_to(axes_frame)
                rsun = getattr(axes_frame, 'rsun', u.R_sun)
                reference_distance = np.sqrt(axes_frame.observer.radius**2 - rsun**2)
                is_visible = coord_in_axes.spherical.distance <= reference_distance
                if np.any(is_visible):
                    axis.plot_coord(coords[is_visible], linestyle='', marker='+', markersize=3, color=(0.9, 0.24, 0.38, 1))
                    # coord_=coords[is_visible]
                    # s=''.join('X' if p else 'O' for p in np.diff(np.diff(coord_.lat))<5*u.deg)
                    # clist = list([m.span()[0], abs(operator.sub(*m.span()))] for m in re.finditer('X+', s))
                    # for c in clist:
                    #   axis.plot_coord(coord_[c[0]:(c[0]+c[1])], linestyle='-', linewidth=0.8, color=(0, 0, 0, 1))
                if np.any(~is_visible):
                    axis.plot_coord(coords[~is_visible], linestyle='', marker='+', markersize=3, color=(0, 0, 0, 1))
                    # coord_ = coords[~is_visible]
                    # s = ''.join('X' if p else 'O' for p in np.diff(np.diff(coord_.lat))<5*u.deg)
                    # clist = list([m.span()[0], abs(operator.sub(*m.span()))] for m in re.finditer('X+', s))
                    # for c in clist:
                    #    axis.plot_coord(coord_[c[0]:(c[0]+c[1])], linestyle='--', linewidth=0.8, color=(0, 0, 0, 1))

        axis.plot_coord(self.apex, marker='o', color=(0, 0, 0, 1))
        axis.plot_coord(self.base, marker='x', color=(0, 0, 0, 1))
        axis.plot_coord(self.center, marker='+', color=(0, 0, 0, 1))
        axis.plot_coord(concatenate((self.apex, self.base)), linestyle='-', linewidth=0.5, color=(0, 0, 0, 1))

    def to_dataframe(self):
        """
        Returns the spheroid parameters as `~pandas.DataFrame`.
        """
        center_ = self.center.transform_to(frames.HeliographicCarrington)
        data_dict = {
            'hgln': self.center.lon.to_value(u.degree),
            'hglt': self.center.lat.to_value(u.degree),
            'crln': center_.lon.to_value(u.degree),
            'crlt': center_.lat.to_value(u.degree),
            'rcenter': self.rcenter.to_value(u.R_sun),
            'radaxis': self.radaxis.to_value(u.R_sun),
            'orthoaxis1': self.orthoaxis1.to_value(u.R_sun),
            'height': self.height.to_value(u.R_sun),
            'kappa': self.kappa,
            'epsilon': self.epsilon,
        }

        return pd.DataFrame(data_dict, index=[self.center.obstime.to_datetime()])

    def __str__(self):
        '''
        Returns the Spheroid object and parameters in a printable sting array.
        '''
        center_ = self.center.transform_to(frames.HeliographicCarrington)
        output = '<Spheroid object \n'
        output += 'HGLN = %3.2f degrees \n'%self.center.lon.to_value(u.degree)
        output += 'HGLT = %3.2f degrees \n'%self.center.lat.to_value(u.degree)
        output += 'CRLN = %3.2f degrees \n'%center_.lon.to_value(u.degree)
        output += 'CRLT = %3.2f degrees \n'%center_.lat.to_value(u.degree)
        output += 'rcenter = %3.2f Rsun \n'%self.rcenter.to_value(u.R_sun)
        output += 'radaxis = %3.2f Rsun \n'%self.radaxis.to_value(u.R_sun)
        output += 'height = %3.2f Rsun \n'%self.radaxis.to_value(u.R_sun)
        output += 'kappa = %3.2f \n'%self.kappa
        output += 'epsilon = %3.2f '%self.epsilon
        output += '>'
        return output


class ellipsoid(spheroid):
    '''
    A class for the ellipsoid model.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The coordinates of the ellipsoid center.
    rcenter : `~astropy.units.Quantity`
        The radial distance of the ellipsoid center from the solar center.
    radaxis : `~astropy.units.Quantity`
        The ellipsoid's first (radial) semi-axis length.
    orthoaxis1 : `~astropy.units.Quantity`
        The ellipsoid's second semi-axis length (orthogonal axis wrt the radial axis).
    orthoaxis2 : `~astropy.units.Quantity`
        The ellipsoid's third semi-axis axis length (orthogonal axis wrt the radial axis).
    height : `~astropy.units.Quantity`
        The radial distance of the ellipsoid apex from the solar center (e.g., rcenter+radaxis).
    kappa : `~astropy.units.Quantity`
        Ellipsoid self-similar constant, defined as the ratio of the apex height to the length
        of the orthoaxis1. (This is proportional to the first aspect ratio of the ellipsoid)
    alpha : `~astropy.units.Quantity`
        The ellipsoid's second aspect ratio, namely the ratio between the second and the third
         semi-axis.
    epsilon : `~astropy.units.Quantity`
        The eccentricity of the ellipsoid.
    tilt : `~astropy.units.Quantity`
        The tilt angle of the ellipsoid. The axis of rotation is the first (radial) semi-axis.

    Notes
    -----
    Details about ellipsoids can be found here: https://en.wikipedia.org/wiki/Ellipsoid
    '''
    def __init__(self, center, radaxis: (u.R_sun), orthoaxis1: (u.R_sun), orthoaxis2: (u.R_sun), tilt: (u.degree), n=40):
        super().__init__(center, radaxis, orthoaxis1, n)

        self.orthoaxis2 = orthoaxis2
        self.tilt = tilt
        self.alpha = orthoaxis1 / orthoaxis2

    def intersecting_curve(self):
        """
        Returns the coordinates of intersection of the ellipsoid with a unit sphere
         centred at the origin as `~astropy.coordinates.SkyCoord`.
        We use vedo package for this calculation.
        """
        from vedo import Ellipsoid, Sphere
        sph = Sphere(pos=(0, 0, 0), r=1)
        ell = Ellipsoid(pos=(self.rcenter.to_value(1*u.R_sun), 0, 0),
                        axis1=(2*self.radaxis.to_value(1*u.R_sun), 0, 0),
                        axis2=(0, 2*self.orthoaxis1.to_value(1*u.R_sun), 0),
                        axis3=(0, 0, 2*self.orthoaxis2.to_value(1*u.R_sun)))
        sic = sph.intersectWith(ell)
        poi = sic.points()
        x, y, z = self.rotate(poi[:, 0]*u.R_sun, poi[:, 1]*u.R_sun, poi[:, 2]*u.R_sun)

        return SkyCoord(CartesianRepresentation(x, y, z),
                        frame=frames.HeliographicStonyhurst,
                        observer=self.center.observer,
                        obstime=self.center.obstime)

    @property
    def coordinates(self):
        """
        Returns the coordinates of the ellipsoid cloud of points as `~astropy.coordinates.SkyCoord`.
        """
        [x__, y__, z__] = sphere(self.n)

        x_ = self.radaxis * x__ + self.rcenter
        y_ = self.orthoaxis1 * y__
        z_ = self.orthoaxis2 * z__

        x, y, z = self.rotate(x_, y_, z_)

        return SkyCoord(CartesianRepresentation(x, y, z),
                        frame=frames.HeliographicStonyhurst,
                        observer=self.center.observer,
                        obstime=self.center.obstime)

    def to_dataframe(self):
        """
        Returns the ellipsoid parameters as `~pandas.DataFrame`.
        """
        center_ = self.center.transform_to(frames.HeliographicCarrington)
        data_dict = {
            'hgln': self.center.lon.to_value(u.degree),
            'hglt': self.center.lat.to_value(u.degree),
            'crln': center_.lon.to_value(u.degree),
            'crlt': center_.lat.to_value(u.degree),
            'rcenter': self.rcenter.to_value(u.R_sun),
            'radaxis': self.radaxis.to_value(u.R_sun),
            'orthoaxis1': self.orthoaxis1.to_value(u.R_sun),
            'orthoaxis2': self.orthoaxis2.to_value(u.R_sun),
            'tilt': self.tilt.to_value(u.degree),
            'height': self.height.to_value(u.R_sun),
            'kappa': self.kappa,
            'epsilon': self.epsilon,
            'alpha': self.alpha,
        }

        return pd.DataFrame(data_dict, index=[self.center.obstime.to_datetime()])

    @staticmethod
    @u.quantity_input
    def heka_from_rabc(rc: (u.R_sun), radaxis: (u.R_sun), orthoaxis1: (u.R_sun), orthoaxis2: (u.R_sun)):
        '''
        Returns the height, epsilon, kappa, and alpha parameters from the semi-axis length
        '''
        height = rc + radaxis
        kappa = orthoaxis1 / (height - 1 * u.R_sun)
        alpha = orthoaxis1 / orthoaxis2

        if radaxis < orthoaxis1:
            epsilon = -1. * np.sqrt(1 - (radaxis/orthoaxis1)**2)
        elif radaxis > orthoaxis1:
            epsilon = np.sqrt(1. - (orthoaxis1/radaxis)**2)
        elif radaxis == orthoaxis1:
            epsilon = 0.
        else:
            height, kappa, epsilon, alpha = np.nan * u.R_sun, np.nan, np.nan, np.nan

        return height, epsilon, kappa, alpha

    @staticmethod
    @u.quantity_input
    def rabc_from_heka(height: (u.R_sun), epsilon, kappa, alpha):
        '''
        Returns the semi-axis length from the height, epsilon, kappa, and alpha parameters
        '''
        s = kappa * (height - 1. * u.R_sun)
        if epsilon < 0:
            d = s * np.sqrt(1. - epsilon**2)
        elif epsilon > 0:
            d = s / np.sqrt(1. - epsilon**2)
        elif epsilon == 0:
            d = s
        else:
            rc, radaxis, orthoaxis1, orthoaxis2 = np.nan * u.R_sun, np.nan * u.R_sun, np.nan * u.R_sun, np.nan * u.R_sun
            return rc, radaxis, orthoaxis1, orthoaxis2
        radaxis = d
        orthoaxis1 = s
        orthoaxis2 = s / alpha
        rc = height - radaxis

        return rc, radaxis, orthoaxis1, orthoaxis2

    def __str__(self):
        '''
        Returns the Ellipsoid object and parameters in a printable sting array.
        '''
        center_ = self.center.transform_to(frames.HeliographicCarrington)
        output = '<Ellipsoid object \n'
        output += 'HGLN = %3.2f degrees \n'%self.center.lon.to_value(u.degree)
        output += 'HGLT = %3.2f degrees \n'%self.center.lat.to_value(u.degree)
        output += 'CRLN = %3.2f degrees \n'%center_.lon.to_value(u.degree)
        output += 'CRLT = %3.2f degrees \n'%center_.lat.to_value(u.degree)
        output += 'rcenter = %3.2f Rsun \n'%self.rcenter.to_value(u.R_sun)
        output += 'radaxis = %3.2f Rsun \n'%self.radaxis.to_value(u.R_sun)
        output += 'orthoaxis1 = %3.2f Rsun \n'%self.orthoaxis1.to_value(u.R_sun)
        output += 'orthoaxis2 = %3.2f Rsun \n'%self.orthoaxis2.to_value(u.R_sun)
        output += 'height = %3.2f Rsun \n'%self.radaxis.to_value(u.R_sun)
        output += 'kappa = %3.2f \n'%self.kappa
        output += 'alpha = %3.2f \n'%self.alpha
        output += 'epsilon = %3.2f '%self.epsilon
        output += '>'
        return output


class gcs():
    """
    A class of the GCS CME model.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The coordinates of the GCS center (defined at the apex circle center).
    rcenter : `~astropy.units.Quantity`
        The radial distance of the GCS center from the solar center.
    height : `~astropy.units.Quantity`
        The radial distance of the GCS apex from the solar center (e.g., rcenter+rapex).
    kappa : `~astropy.units.Quantity`
        Aspect ratio of the GCS.
    alpha : `~astropy.units.Quantity`
        Half angular width of the GCS model.
    tilt : `~astropy.units.Quantity`
        The tilt angle of the GCS.

    Notes
    -----
    This is based on IDL script in here: https://hesperia.gsfc.nasa.gov/ssw/stereo/secchi/idl/scraytrace/shellskeleton.pro
    and with some additions from here: https://gitlab.physik.uni-kiel.de/ET/gcs_python/-/blob/master/gcs/geometry.py
    A full description of the GCS model is given in Thernisien et al. (2006) and Thernisien (2011).
    """
    def __init__(self, center, height: (u.R_sun), alpha: (u.degree), kappa, tilt: (u.degree),
                 nbvertsl=10, nbvertcirc=20, nbvertcircshell=90):
        self.center = center.transform_to(frames.HeliographicStonyhurst)  # This is the SkyCoord of the center of circle at the apex
        self.rcenter = center.radius  # rcenter = height-rapex
        self.height = height  # GCS height in apex
        self.alpha = alpha  # GCS CME width (in half angle)
        self.kappa = kappa  # Aspect ratio
        self.tilt = tilt  # Tilt angle

        self.h = self.h_()  # Height of the cone (CME legs)
        self.rapex = self.rapex_()  # Cross-section radius at the apex

        self.nbvertsl = nbvertsl
        self.nbvertcirc = nbvertcirc
        self.nbvertcircshell = nbvertcircshell

    @staticmethod
    @u.quantity_input
    def rcenter_(height: (u.R_sun), alpha: (u.degree), kappa):
        """
        Returns the distance of the GCS center from the solar center (OC1 in T2011, see eq. 20).
        """
        h = height * (1. - kappa) * np.cos(alpha) / (1. + np.sin(alpha))
        rapex = kappa * (h / np.cos(alpha) +  h * np.tan(alpha)) / (1 - kappa ** 2)
        return height-rapex

    def h_(self):
        """
        Returns the full height (h) of the cone (OD distance in T2011, see eq. 2).
        We use eq. 3 in T2006 to derive the h from height, kappa, and alpha parameters
        """
        h = self.height * (1. - self.kappa) * np.cos(self.alpha) / (1. + np.sin(self.alpha))
        return h

    def rapex_(self):
        """
        Returns the cross-section radius at the apex (see eq. 29 in T2011).
        """
        rapex = self.kappa * (self.h_() / np.cos(self.alpha) +  self.h_() * np.tan(self.alpha)) / (1 - self.kappa ** 2)
        return rapex

    @property
    def coordinates(self):
        """
        Returns the coordinates of the ellipsoid cloud of points as `~astropy.coordinates.SkyCoord`.
        """
        x_, y_, z_ = self.cloud()
        x, y, z = self.rotate(x_ * u.R_sun, y_ * u.R_sun, z_ * u.R_sun)

        return SkyCoord(CartesianRepresentation(x, y, z),
                        frame=frames.HeliographicStonyhurst,
                        observer=self.center.observer,
                        obstime=self.center.obstime)

    def rotate(self, x_, y_, z_):
        """
        Rotates the (x,y,z) coordinates from the geometrical model's coordinate system
        to the Heliospheric coordinates.

        Returns the coordinates as `~astropy.coordinates.SkyCoord`.
        """
        Longc = self.center.lon.to_value(u.rad)
        Latc = self.center.lat.to_value(u.rad)
        tilt = self.tilt.to_value(u.rad)

        v = np.transpose([x_.flatten(), y_.flatten(), z_.flatten()])
        v = Rotation.from_euler('zyz', [tilt, -Latc + np.pi/2, Longc]).apply(v)

        x = np.reshape(v[:, 0], x_.shape) * x_.unit
        y = np.reshape(v[:, 1], y_.shape) * y_.unit
        z = np.reshape(v[:, 2], z_.shape) * z_.unit

        return x, y, z

    def shell_skeleton(self):
        """
        Computes the skeleton axis of the GCS CME model.
        Converted and modified from the IDL script shellskeleton.pro
        https://hesperia.gsfc.nasa.gov/ssw/stereo/secchi/idl/scraytrace/shellskeleton.pro
        """
        self.height.to_value(u.R_sun)
        alpha = self.alpha.to_value(u.rad)
        kappa = self.kappa
        h = self.h.to_value(u.R_sun)
        gamma = np.arcsin(self.kappa)

        # calculate the points of the straight line part
        pRstart = np.array([0, np.sin(alpha), np.cos(alpha)])
        pLstart = np.array([0, -np.sin(alpha), np.cos(alpha)])
        vertsl = np.linspace(0, h, self.nbvertsl)
        pslR = np.outer(vertsl, pRstart)
        pslL = np.outer(vertsl, pLstart)
        rsl = np.tan(gamma) * norm(pslR, axis=1)
        casl = np.full(self.nbvertsl, -alpha)

        # calculate the points of the circular part
        beta = np.linspace(-alpha, np.pi/2, self.nbvertcirc)
        hf = h
        h = hf / np.cos(alpha)
        rho = hf * np.tan(alpha)

        X0 = (rho + h * kappa ** 2 * np.sin(beta)) / (1 - kappa ** 2)
        rc = np.sqrt((h ** 2 * kappa ** 2 - rho ** 2) / (1 - kappa ** 2) + X0 ** 2)
        cac = beta

        pcR = np.array([np.zeros(beta.shape), X0 * np.cos(beta), h + X0 * np.sin(beta)]).T
        pcL = np.array([np.zeros(beta.shape), -X0 * np.cos(beta), h + X0 * np.sin(beta)]).T

        # This part is from here
        # https://gitlab.physik.uni-kiel.de/ET/gcs_python/-/blob/master/gcs/geometry.py
        r = np.concatenate((rsl, rc[1:], np.flipud(rc)[1:], np.flipud(rsl)[1:]))
        ca = np.concatenate((casl, cac[1:], np.pi-np.flipud(cac)[1:], np.pi-np.flipud(casl)[1:]))
        p = np.concatenate((pslR, pcR[1:], np.flipud(pcL)[1:], np.flipud(pslL)[1:]))

        return p, r, ca

    def cloud(self):
        """
        Returns the cloud of points for the GCS model.

        Notes
        -----
        Converted and modified from the IDL script cmecloud.pro.
        https://hesperia.gsfc.nasa.gov/ssw/stereo/secchi/idl/scraytrace/cmecloud.pro
        https://gitlab.physik.uni-kiel.de/ET/gcs_python/-/blob/master/gcs/geometry.py
        """

        # Compute the shell's skeleton axis
        p, r, ca = self.shell_skeleton()

        # This part is from here
        # https://gitlab.physik.uni-kiel.de/ET/gcs_python/-/blob/master/gcs/geometry.py
        ptheta = np.linspace(0, 2*np.pi, self.nbvertcircshell)
        pphi = np.arange(0, 2*(self.nbvertcirc + self.nbvertsl) - 3)
        uv, vv = np.meshgrid(ptheta, pphi)
        uv, vv = uv.flatten(), vv.flatten()
        mesh = r[vv, np.newaxis] * np.array([np.cos(uv), np.sin(uv) * np.cos(ca[vv]), np.sin(uv) * np.sin(ca[vv])]).T + p[vv]

        return mesh[:, 0], mesh[:, 1], mesh[:, 2]

    def plot(self, axis, mode=False, only_surface=False):
        '''
        Plots the GCS model at the provided axes
        '''
        axis.plot_coord(self.coordinates, color='green', linestyle='-', linewidth=0.4)

    def to_dataframe(self):
        """
        Returns the GCS parameters as `~pandas.DataFrame`.
        """
        center_ = self.center.transform_to(frames.HeliographicCarrington)
        data_dict = {
            'hgln': self.center.lon.to_value(u.degree),
            'hglt': self.center.lat.to_value(u.degree),
            'crln': center_.lon.to_value(u.degree),
            'crlt': center_.lat.to_value(u.degree),
            'rcenter': self.rcenter.to_value(u.R_sun),
            'height': self.height.to_value(u.R_sun),
            'rapex': self.rapex.to_value(u.R_sun),
            'alpha': self.alpha.to_value(u.degree),
            'kappa': self.kappa,
            'tilt': self.tilt.to_value(u.degree),
        }

        return pd.DataFrame(data_dict, index=[self.center.obstime.to_datetime()])

    def __str__(self):
        '''
        Returns the GCS object and parameters in a printable sting array.
        '''
        center_ = self.center.transform_to(frames.HeliographicCarrington)
        output = '<GCS object \n'
        output += 'HGLN = %3.2f degrees \n'%self.center.lon.to_value(u.degree)
        output += 'HGLT = %3.2f degrees \n'%self.center.lat.to_value(u.degree)
        output += 'CRLN = %3.2f degrees \n'%center_.lon.to_value(u.degree)
        output += 'CRLT = %3.2f degrees \n'%center_.lat.to_value(u.degree)
        output += 'rcenter = %3.2f Rsun \n'%self.rcenter.to_value(u.R_sun)
        output += 'height = %3.2f Rsun \n'%self.radaxis.to_value(u.R_sun)
        output += 'alpha = %3.2f \n'%self.alpha
        output += 'kappa = %3.2f \n'%self.kappa
        output += 'tilt = %3.2f '%self.tilt
        output += '>'
        return output


def my_plot_coord(coord, axis, **kargs):

    ass = coord.spherical.distance > 1*u.R_sun  # above solar surface
    s = ''.join('X' if p else 'O' for p in ass)
    clist = list([m.span()[0], abs(operator.sub(*m.span()))] for m in re.finditer('X+', s))
    for c in clist:
        axis.plot_coord(coord[c[0]:(c[0]+c[1])], **kargs)
