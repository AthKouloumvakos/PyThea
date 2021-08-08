import numpy as np
import pandas as pd
import datetime
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R  
import astropy.units as u
from astropy.coordinates import (
    SkyCoord,
    Distance,
    SphericalRepresentation,
    CartesianRepresentation,
    concatenate
)
from sunpy.coordinates import frames
from scipy.spatial.transform import Rotation
from numpy.linalg import norm
import seaborn as sns
import re
import operator

def sphere(n= 20, plot = 0):
    """
    @summary: generates three (n+1)-by-((N+1) matrices representing the coordinates of a unit sphere
    @param n: number of points in x and y
    @return: x,y,z coordinates of the unit sphere
        
    """

    theta = np.linspace( -np.pi, np.pi, n+1)
    phi = np.linspace(-1*np.pi/2, 1 * np.pi/2, n+1)
    cosphi = np.cos(phi)
    #cosphi[0]=0
    #cosphi[n]=0;
    sintheta = np.sin(theta)
    #sintheta[0] = 0
    #sintheta[n]=0

    x =  np.outer(cosphi, np.cos(theta))
    y =  np.outer(cosphi, sintheta)
    z =  np.outer(np.sin(phi), np.ones_like(theta))
    
    return [x, y, z]

class spheroid:
    @u.quantity_input
    def __init__(self, center, radaxis: (u.R_sun), orthoaxis1: (u.R_sun), n=40):
        self.center = center.transform_to(frames.HeliographicStonyhurst) # This is the SkyCoord of the ellipsoid center
        self.rcenter = center.radius
        self.radaxis = radaxis
        self.orthoaxis1 = orthoaxis1

        self.height = radaxis + self.rcenter
        self.kappa = orthoaxis1 / (self.height-1.*u.R_sun)
        
        if radaxis<orthoaxis1:
            self.epsilon = -1. * np.sqrt(1. - (radaxis/orthoaxis1)**2)
        elif radaxis>orthoaxis1:
            self.epsilon = 1. * np.sqrt(1. - (orthoaxis1/radaxis)**2)
        elif radaxis==orthoaxis1:
            self.epsilon = 0.

        self.tilt = 0 * u.degree
        self.n = n
    
    @property
    def coordinates(self):
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
        # e.g., see https://mathworld.wolfram.com/Sphere-SphereIntersection.html
        d = self.rcenter
        r = self.radaxis
        R = 1 * u.R_sun

        step = 40
        theta = np.linspace(0, 2*np.pi, step+1)[:-1]
        x_ = np.ones_like(theta) * (d**2 - r**2 + R**2) / (2 * d)
        # y_**2 - z_**2 = R**2 - x_**2
        alpha = np.sqrt(4 * d**2 * R**2 - (d**2 - r**2 + R**2)**2 ) / (2 * d)
        y_ = alpha * np.sin(theta)
        z_ = alpha * np.cos(theta)
        x, y, z = self.rotate(x_, y_, z_)

        return SkyCoord(CartesianRepresentation(x, y, z),
                        frame=frames.HeliographicStonyhurst,
                        observer=self.center.observer,
                        obstime=self.center.obstime)

    @staticmethod
    @u.quantity_input
    def heka_to_rabc(rc: (u.R_sun), radaxis: (u.R_sun), orthoaxis1: (u.R_sun)):
        height = rc + radaxis
        kappa = orthoaxis1 / (height - 1* u.R_sun)

        if radaxis<orthoaxis1:
            epsilon = -1. * np.sqrt(1 - (radaxis/orthoaxis1)**2)
        elif radaxis>orthoaxis1:
            epsilon = np.sqrt(1. - (orthoaxis1/radaxis)**2)
        elif radaxis==orthoaxis1:
            epsilon = 0.
        else:
            height, kappa, epsilon = np.nan * u.R_sun, np.nan, np.nan

        return height, kappa, epsilon

    @staticmethod
    @u.quantity_input
    def rabc_to_heka(height: (u.R_sun), kappa, epsilon):
        s = kappa * (height - 1. * u.R_sun)
        if epsilon<0:
            d = s * np.sqrt(1. - epsilon**2)
        elif epsilon>0:
            d = s / np.sqrt(1. - epsilon**2)
        elif epsilon==0:
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
        return SkyCoord(SphericalRepresentation(self.center.lon, self.center.lat, 
                                                Distance(self.center.radius + self.radaxis, allow_negative=True)),
                             frame=frames.HeliographicStonyhurst,
                             obstime=self.center.obstime,
                             observer=self.center.observer)
    @property
    def base(self):
        c = np.sign(-self.center.radius + self.radaxis) # TODO: A hack because negative distance values are not allowed, wait for astropy-4.3
        return SkyCoord(SphericalRepresentation(-c*self.center.lon, -c*self.center.lat, 
                                                Distance( c*(-self.center.radius + self.radaxis), allow_negative=True)),
                             frame=frames.HeliographicStonyhurst,
                             obstime=self.center.obstime,
                             observer=self.center.observer)
    def rotate(self, x_, y_, z_):
        Longc = self.center.lon.to_value(u.rad)
        Latc  = self.center.lat.to_value(u.rad)
        tilt = self.tilt.to_value(u.rad)

        v = np.transpose([x_.flatten(), y_.flatten(), z_.flatten()]) 
        v = Rotation.from_euler('xyz', [tilt, -Latc, Longc]).apply(v)
        
        x = np.reshape(v[:,0], x_.shape) * x_.unit
        y = np.reshape(v[:,1], y_.shape) * y_.unit
        z = np.reshape(v[:,2], z_.shape) * z_.unit
        return x, y, z

    def plot(self, axis, redused=False):
        if redused is False:
            axis.plot_coord(self.coordinates, color='red', linestyle='-', linewidth=0.2)
        elif redused is True:
            lw = 1.
            n = self.n
            palete = sns.color_palette("colorblind")
            
            # Plot the part of the intersecting curve at the solar surface
            coords = self.intersecting_curve()
            axes_frame = axis._transform_pixel2world.frame_out
            coord_in_axes = coords.transform_to(axes_frame)
            rsun = getattr(axes_frame, 'rsun', u.R_sun)
            reference_distance = np.sqrt(axes_frame.observer.radius**2 - rsun**2)
            is_visible = coord_in_axes.spherical.distance <= reference_distance
            if np.any(is_visible):
                axis.plot_coord(coords[is_visible], linestyle='-', linewidth=0.8, color=(0,0,0,1))
                #coord_ = coords[is_visible]
                #s = ''.join('X' if p else 'O' for p in np.diff(np.diff(coord_.lat))<5*u.deg)
                #clist = list([m.span()[0], abs(operator.sub(*m.span()))] for m in re.finditer('X+', s))
                #for c in clist:
                #   axis.plot_coord(coord_[c[0]:(c[0]+c[1])], linestyle='-', linewidth=0.8, color=(0,0,0,1))
            if np.any(~is_visible):
                axis.plot_coord(coords[~is_visible], linestyle='--', linewidth=0.8, color=(0,0,0,1))
                #coord_ = coords[~is_visible]
                #s = ''.join('X' if p else 'O' for p in np.diff(np.diff(coord_.lat))<5*u.deg)
                #clist = list([m.span()[0], abs(operator.sub(*m.span()))] for m in re.finditer('X+', s))
                #for c in clist:
                #    axis.plot_coord(coord_[c[0]:(c[0]+c[1])], linestyle='--', linewidth=0.8, color=(0,0,0,1))

            # Plot the part of the spheroid mesh but filter lines bellow the solar surface
            axis.plot_coord(self.coordinates[:,int(n/2)], color=palete[0],linestyle='-', linewidth=lw)
            my_plot_coord(self.coordinates[:,0], axis, color=palete[6], linestyle='-', linewidth=lw)
            axis.plot_coord(self.coordinates[0:int(n/2)+1,int(n/4)], color=palete[1], linestyle='-', linewidth=lw)
            axis.plot_coord(self.coordinates[int(n/2):n+1,int(n/4)], color=palete[1], linestyle='-', linewidth=lw)
            axis.plot_coord(self.coordinates[0:int(n/2)+1,int(3*n/4)], color=palete[1], linestyle='-', linewidth=lw)
            axis.plot_coord(self.coordinates[int(n/2):n+1,int(3*n/4)], color=palete[1], linestyle='-', linewidth=lw)
            my_plot_coord(self.coordinates[int(n/2),0:int(n/4)+1], axis, color=palete[2],linestyle='-', linewidth=lw)
            my_plot_coord(self.coordinates[int(n/2),int(3*n/4):n+1], axis, color=palete[2],linestyle='-', linewidth=lw)
            axis.plot_coord(self.coordinates[int(n/2),int(n/4):int(3*n/4)+1], color=palete[3],linestyle='-', linewidth=lw)
            axis.plot_coord(self.apex, marker = 'o',color=(0,0,0,1))
            axis.plot_coord(self.base, marker = 'x',color=(0,0,0,1))
            axis.plot_coord(concatenate((self.apex, self.base)),linestyle='-', linewidth=0.5, color=(0,0,0,1))

    def to_dataframe(self):
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

class gcs():
    """
    A class of the GCS CME model.
    This is based on IDL script in here: https://hesperia.gsfc.nasa.gov/ssw/stereo/secchi/idl/scraytrace/shellskeleton.pro
    and with some additions from here: https://gitlab.physik.uni-kiel.de/ET/gcs_python/-/blob/master/gcs/geometry.py
    """
    def __init__(self, center, height: (u.R_sun), alpha: (u.degree), kappa, tilt: (u.degree),
                 nbvertsl=10, nbvertcirc=20, nbvertcircshell=90):
        self.center = center.transform_to(frames.HeliographicStonyhurst) # This is the SkyCoord of the center of cyrcle at the apex
        self.rcenter = center.radius # rcenter = height-rappex
        self.height = height # GCS height in appex
        self.alpha = alpha # GCS CME width (in half angle)
        self.kappa = kappa # ratio
        self.tilt = tilt

        self.distjunc = height * (1-kappa) * np.cos(alpha) / (1.+np.sin(alpha))
        self.rappex = self.kappa * (self.distjunc / np.cos(self.alpha) + \
                      self.distjunc * np.tan(self.alpha)) / (1 - self.kappa ** 2)

        self.nbvertsl = nbvertsl
        self.nbvertcirc = nbvertcirc
        self.nbvertcircshell = nbvertcircshell

    @staticmethod
    @u.quantity_input
    def rcenter_(height: (u.R_sun), alpha: (u.degree), kappa):
        h = height * (1. - kappa) * np.cos(alpha) / (1. + np.sin(alpha))
        rappex = kappa * (h / np.cos(alpha) +  h * np.tan(alpha)) / (1 - kappa ** 2)
        return height-rappex

    @property
    def coordinates(self):
        x_, y_, z_ = self.cloud()
        x, y, z = self.rotate(x_* u.R_sun, y_* u.R_sun, z_* u.R_sun)

        return SkyCoord(CartesianRepresentation(x, y, z),
                        frame=frames.HeliographicStonyhurst,
                        observer=self.center.observer,
                        obstime=self.center.obstime)

    def rotate(self, x_, y_, z_):
        Longc = self.center.lon.to_value(u.rad)
        Latc = self.center.lat.to_value(u.rad)
        tilt = self.tilt.to_value(u.rad)

        v = np.transpose([x_.flatten(), y_.flatten(), z_.flatten()]) 
        v = Rotation.from_euler('zyz', [tilt, -Latc+np.pi/2, Longc]).apply(v)

        x = np.reshape(v[:,0], x_.shape) * x_.unit
        y = np.reshape(v[:,1], y_.shape) * y_.unit
        z = np.reshape(v[:,2], z_.shape) * z_.unit

        return x, y, z

    def shell_skeleton(self):
        """
        Computes the skeleton axis of the GCS CME model.
        Converted and modified from the IDL script shellskeleton.pro
        https://hesperia.gsfc.nasa.gov/ssw/stereo/secchi/idl/scraytrace/shellskeleton.pro
        """
        height = self.height.to_value(u.R_sun)
        alpha = self.alpha.to_value(u.rad)
        kappa = self.kappa
        distjunc = self.distjunc.to_value(u.R_sun)
        gamma = np.arcsin(self.kappa)

        # calculate the points of the straight line part
        pRstart = np.array([0, np.sin(alpha), np.cos(alpha)])
        pLstart = np.array([0, -np.sin(alpha), np.cos(alpha)])
        vertsl = np.linspace(0, distjunc, self.nbvertsl)
        pslR = np.outer(vertsl, pRstart)
        pslL = np.outer(vertsl, pLstart)
        rsl = np.tan(gamma) * norm(pslR, axis=1)
        casl = np.full(self.nbvertsl, -alpha)

        # calculate the points of the circular part
        beta = np.linspace(-alpha, np.pi/2, self.nbvertcirc)
        hf = distjunc
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
        Build a cloud of the GCS model.
        Converted and modified from the IDL script cmecloud.pro.
        https://hesperia.gsfc.nasa.gov/ssw/stereo/secchi/idl/scraytrace/cmecloud.pro
        https://gitlab.physik.uni-kiel.de/ET/gcs_python/-/blob/master/gcs/geometry.py
        """
        height = self.height.to_value(u.R_sun)
        alpha = self.alpha.to_value(u.rad)
        kappa = self.kappa
        distjunc = self.distjunc.to_value(u.R_sun)
        
        # Compute the shell's skeleton axis
        p, r, ca = self.shell_skeleton()
        
        # This part is from here
        # https://gitlab.physik.uni-kiel.de/ET/gcs_python/-/blob/master/gcs/geometry.py
        ptheta = np.linspace(0, 2*np.pi, self.nbvertcircshell)
        pphi = np.arange(0, 2*(self.nbvertcirc + self.nbvertsl) - 3)
        uv, vv = np.meshgrid(ptheta, pphi)
        uv, vv = uv.flatten(), vv.flatten()
        mesh = r[vv, np.newaxis] * np.array([np.cos(uv), np.sin(uv) * np.cos(ca[vv]), np.sin(uv) * np.sin(ca[vv])]).T + p[vv]

        return mesh[:,0], mesh[:,1], mesh[:,2]

    def plot(self, axis, redused=False):
        axis.plot_coord(self.coordinates, color='red', linestyle='-', linewidth=0.4)
        
    def to_dataframe(self):
        center_ = self.center.transform_to(frames.HeliographicCarrington)
        data_dict = {
                     'hgln': self.center.lon.to_value(u.degree),
                     'hglt': self.center.lat.to_value(u.degree),
                     'crln': center_.lon.to_value(u.degree),
                     'crlt': center_.lat.to_value(u.degree),
                     'rcenter': self.rcenter.to_value(u.R_sun),
                     'height': self.height.to_value(u.R_sun),
                     'rappex': self.rappex.to_value(u.R_sun),
                     'alpha': self.alpha.to_value(u.degree),
                     'kappa': self.kappa,
                     'tilt': self.tilt.to_value(u.degree),
                    }
        
        return pd.DataFrame(data_dict, index=[self.center.obstime.to_datetime()])

def my_plot_coord(coord, axis, **kargs):
    ass = coord.spherical.distance > 1*u.R_sun # above solar surface
    s = ''.join('X' if p else 'O' for p in ass)
    clist = list([m.span()[0], abs(operator.sub(*m.span()))] for m in re.finditer('X+', s))
    for c in clist:
        axis.plot_coord(coord[c[0]:(c[0]+c[1])], **kargs)
