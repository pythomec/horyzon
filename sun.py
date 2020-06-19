"""
Functions related to sun position and surface irradiation
"""

import numpy as np
import xarray as xr

from pysolar import solar # TODO: make it an optional import
import datetime

from . import timezones as tzs

sun_radius_deg = 0.26667

def get_sun_position(lat, lon, when):
    """Get sun elevation and azimuth for given position(s) and datetime

    Includes correction on refraction.

    :param lat: latitude, float or array
    :param lon: longitude, float or array
    :param when: datetime, must contain time zone
    :return: elevation, azimuth, floats or arrays
    """

    # TODO: find meaning of the elevation parameter below (altitude of the observer?) and use correct value instead of 0
    azimuth, elevation = solar.get_position(latitude_deg=lat, longitude_deg=lon, when=when, elevation=0)

    return xr.Dataset(data_vars={'elevation': np.array(elevation).reshape(-1,1,1),
                                 'azimuth': np.array(azimuth)[:, None, None]}, coords={'datetime': when, 'lon':lon,
                                                                                      'lat': lat})


def get_sun_trajectory(lat, lon, when, step_minutes=1):
    """Compute sun trajectory from midnight to midnight as viewed from a particular point

    :param lat: latitude of the observer, float or array
    :param lon: longitude of the observer, float or array
    :param when: date or datetime object, must contain timezone info
    :return: DataArray elevations, coords: time, azimuth
    """

    when = tzs.fill_default_tz(when)

    dates = [when + datetime.timedelta(minutes=i*step_minutes) for i in range(60//step_minutes*24)]
    trajectory = [get_sun_position(lat, lon, d) for d in dates]
    trajectory = xr.concat(trajectory, dim='datetime')

    # TODO: add lat, lon as attributes or datarray dimensions

    return trajectory

