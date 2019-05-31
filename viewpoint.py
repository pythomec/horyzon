"""
Viewpoint class
"""

import numbers

import numpy as np
import pylab as plt

from . import elevation_angle as elevation
from . import visibility as vis

direction2degrees = {dir: deg for dir, deg in zip(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'],
                                                  np.arange(0, 360, 45))}


class Viewpoint:
    """Class Viewpoint represents view from a particular point and handles plotting of panorama and horizon.

    """

    def __init__(self, altitude, lon, lat, above_ground=5, dtheta=360 / 3600 * 2, dr=0.001, maxr=1):
        """Constructor of Viewpoint class

        :param altitude: DataArray. Altitude vs longitude and latitude.
        :param lon: longitude of the viewpoint
        :param lat: latitude of the viewpoint
        :param above_ground: height of the observer above ground [m] (default: 5)
        :param dtheta: step in the polar angle [°], default: 360/3600*2
        :param dr: step in radial coordinate
        :param maxr: maximum distance where to look for peaks
        """

        self.altitude = altitude
        elevation_angle = elevation.compute_elevation(altitude, lat=lat, lon=lon, above_ground=above_ground)
        self.elevation_angle_polar = vis.data2polar(elevation_angle, dtheta=dtheta, dr=dr, maxr=maxr)
        self.mask, self.ridges = vis.get_mask_and_ridges(self.elevation_angle_polar)
        self.horizon = vis.compute_horizon(self.mask, self.elevation_angle_polar)

    @property
    def lon(self):
        """Longitude of the viewpoint

        """
        return self.elevation_angle_polar.attrs['lon']

    @property
    def lat(self):
        """Latitude of the viewpoint

        """
        return self.elevation_angle_polar.attrs['lat']

    @staticmethod
    def plot_panorama_scatter(elevation_angles_polar, mask, rotate=0, y_in_degrees=False, **kwargs):
        """Plot panorama based on elevations in polar coordinates and a pixel mask

        :param elevation_angles_polar: DataArray. Viewing angles in polar coordinates
        :param mask: DataArray. Pixel mask, can be either all visible pixels or just ridges
        :param rotate: rotate azimuth [°]
        :param y_in_degrees: y axis shows elevation angle (True) or projection to a vertical plane? Default: False
        :param kwargs: kwargs passed to plotting function (scatter)
        :return: None
        """

        # get horizontal and vertical angles of all visible (non-masked) points
        thetamesh, rmesh = np.meshgrid(mask.theta, mask.r)
        azimuth = thetamesh.ravel()[mask.values.ravel() == 1]
        angles = elevation_angles_polar.values.ravel()[mask.values.ravel() == 1]

        # rotate the view horizontally and convert y to [m] if requested
        azimuth = (azimuth + rotate) % (360)
        y = angles if y_in_degrees else np.arctan(np.deg2rad(angles))

        # plotting
        plt.scatter(azimuth, y, s=1, **kwargs)

    def plot_panorama(self, direction='S', y_in_degrees=False, figsize=(10, 2.5), newfig=True,
                      xlabel='', ylabel=''):
        """Plot panoramic view of the surroundings

        :param direction: center of the panorama aims at this direction ('N','NE', ... or number [°]), default: 'S'
        :param y_in_degrees: y axis shows elevation angle (True) or projection to a vertical plane? Default: False
        :param figsize: (dx,dy), default (10, 2.5)
        :param newfig: open new figure? default: True
        :return: None
        """

        # find proper horizontal rotation
        direction = direction2degrees.get(direction, direction)
        if not isinstance(direction, numbers.Real):
            raise ValueError('unknown direction %s' % (direction,))

        rotate = (180 - direction) % 360

        if newfig:
            plt.figure(figsize=figsize)

        # plot all visible points
        self.plot_panorama_scatter(self.elevation_angle_polar, self.mask, rotate=rotate, y_in_degrees=y_in_degrees,
                                   c='gray', alpha=0.2)
        # highlight edges
        self.plot_panorama_scatter(self.elevation_angle_polar, self.ridges, rotate=rotate, y_in_degrees=y_in_degrees,
                                   c='k')

        # anotate the plot
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if y_in_degrees:
            plt.ylabel('[°]')
        else:
            plt.yticks([])
        plt.title('lon =%.2f°, lat =%.2f°' % (self.lon, self.lat))
        plt.xlim([0, 360])
        labels = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
        plt.xticks(ticks=(np.array([direction2degrees[l] for l in labels]) + rotate) % 360,
                   labels=labels)
        plt.tight_layout()

    def plot_horizon_lonlat(self, *args, **kwargs):
        """Plot horizon in lon, lat coordinates

        :param args: args passed to matplotlib.plot
        :param kwargs: kwargs passed to matplotlib.plot
        :return: None
        """

        plt.plot(self.horizon.lon, self.horizon.lat, *args, **kwargs)
