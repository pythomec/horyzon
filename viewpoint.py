"""
Viewpoint class
"""

import numbers
from functools import lru_cache
import datetime
from collections.abc import Iterable

import numpy as np
import pylab as plt
import xarray as xr

from . import elevation_angle as elevation
from . import visibility as vis
from . import sun
from . import timezones as tzs

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

    def sun_position(self, when):
        """Compute sun position at particular date and time

        :param when: date or datetime object or array of these, must contain time zone, date is expanded to the whole
                     sun trajectory from midnight to midnight
        :return: dataset with sun eleveation vs azimuth and is_visible flag
        """

        if not isinstance(when, Iterable):
            when = [when]

        elevations = [self._sun_trajectory(w) if isinstance(w, datetime.date) else
                      sun.get_sun_position(lon=self.lon, lat=self.lat,
                                           when=tzs.fill_default_tz(w, lon=self.lon, lat=self.lat)) for w in when]

        elevations = xr.concat(elevations, dim='datetime')
        return elevations

        horizon = self.horizon
        shift_by_sun_radius = sun.sun_radius_deg if is_edge_visible else 0
        is_visible = horizon.elevation.interp(theta=sun_elevation.azimuth) < sun_elevation + shift_by_sun_radius

        sun_trajectory = xr.Dataset(data_vars={'elevation': ('datetime', sun_elevation),
                                               'is_visible': ('datetime', is_visible)},
                                    coords=sun_elevation.coords)

        return sun_trajectory


    def sun_trajectory(self, when, is_edge_visible = True):
        """Return sun trajectory and visibility mask for given date

        :param when: date or datetime object, must contain time zone
        :return: dataset with sun elevation vs azimuth/time and is_visible mask
        """

        when = datetime.datetime(when.year, when.month, when.day, tzinfo=when.tzinfo)
        sun_elevation = self._sun_trajectory(when)

        horizon = self.horizon
        shift_by_sun_radius = sun.sun_radius_deg if is_edge_visible else 0
        is_visible = horizon.elevation.interp(theta=sun_elevation.azimuth) < sun_elevation + shift_by_sun_radius

        sun_trajectory = xr.Dataset(data_vars={'elevation': ('datetime', sun_elevation),
                                               'is_visible': ('datetime', is_visible)},
                                    coords=sun_elevation.coords)

        return sun_trajectory

    @lru_cache(maxsize=16) # TODO: set reasonable cache size
    def _sun_trajectory(self, when):
        """Return sun trajectory and visibility mask for given date

        The position is corrected on refraction. Is_visible mask is true even if the center of the Sun is below
        horizon but part of the sun disk is still visible.

        :param when: datetime object, must contain time zone
        :return: edataset with sun elevation vs azimuth/time and is_visible mask
        """

        return sun.get_sun_trajectory(self.lat, self.lon, when)

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

    def plot_sun_trajectory(self, when, y_in_degrees=False):
        trajectory = self.sun_trajectory(when, is_edge_visible=False)

        sun_trajectory_x = trajectory.azimuth.where(trajectory.is_visible)

        sun_trajectory_y = trajectory.elevation.where(trajectory.is_visible)
        if not y_in_degrees:
            sun_trajectory_y = np.arctan(np.deg2rad(sun_trajectory_y))

        plt.plot(sun_trajectory_x, sun_trajectory_y, 'y')

    def plot_panorama(self, direction='S', y_in_degrees=False, figsize=(10, 2.5), newfig=True,
                      xlabel='', ylabel='', when=None):
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

        # plot the Sun
        if when is not None:
            self.plot_sun_trajectory(when, y_in_degrees=y_in_degrees)

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
