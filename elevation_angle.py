import numpy as np
import xarray as xr


def geo_to_cart(lat, lon, alt):
    """Convert geographic coordinates to cartesian

    :param lat: array of latitude coordinates
    :param lon: array of longtitude coordinates
    :param alt: array of altitudes
    :return: array of cartesian coordinates
    """
    Rz = 6.378e6
    x = (Rz + alt) * np.sin(lon) * np.cos(lat)
    y = (Rz + alt) * np.cos(lon) * np.cos(lat)
    z = (Rz + alt) * np.sin(lat)

    return np.dstack((x, y, z))


def unit_vector(vector):
    """Return the unit vector of the vector.

    :param vector: vector to be normalised
    :return: normalised vector
    """
    return vector / np.dstack([np.linalg.norm(vector, axis=2)] * 3)


def angle_between(v1, v2):
    """Calculate angle between two vectors

    :param v1: first vector
    :param v2: second vector
    :return: angle between the vectors
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.nan_to_num(np.clip(np.dot(v1_u, v2_u[0, 0, :]), -1.0, 1.0)))


def compute_elevation(d, lat, lon, above_ground=5):
    """Compute elevation angle between observer and data. Take into account curvature of the Earth.
    
    :param d: DataArray with altitudes
    :param lat: latitude of the viewpoint
    :param lon: longitude of the viewpoint
    :param above_ground: height of the observer above terrain [m] (default: 5 m)
    :return: DataArray with elevation angles from given observation point

    """
    latp = np.radians(d.lat)
    lonp = np.radians(d.lon)
    la, lo = np.meshgrid(latp, lonp)

    lato = np.radians(lat)
    lono = np.radians(lon)
    alto = d.data[np.argmin(np.abs(latp - lato)), np.argmin(np.abs(lonp - lono))] + above_ground

    # convert coordinates of the observer and of all points on the map to cartesian geocentric coordinates
    Z = geo_to_cart(la.T, lo.T, d)
    O = geo_to_cart(lato, lono, alto)

    # compute elevation angle
    OZ = Z - O
    A = angle_between(OZ, O)
    A = -np.degrees(A) + 90

    an = xr.DataArray(A.T, coords={'lon': d.lon, 'lat': d.lat}, dims=('lon', 'lat'), name='elevation',
                      attrs={**d.attrs, 'lon': lon, 'lat': lat, 'above_ground': above_ground})

    return an
