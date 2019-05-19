import numpy as np

import matplotlib.pyplot as plt
import scipy
import xarray
from pathlib import *


def geo_to_cart(lat, lon, alt):
    """
    Conversion of geographic coordinates to cartesian
    :param lat: array of latitude coordinates
    :param lon: array of longtitude coordinates
    :param alt: array of altitudes
    :return: array of cartesian coordinates
    """
    Rz=6.378e6
    x = (Rz + alt) * np.sin(lon) * np.cos(lat)
    y = (Rz + alt) * np.cos(lon) * np.cos(lat)
    z = (Rz + alt) * np.sin(lat)
    #print('shape z', np.shape(z))
    return np.dstack((x, y, z))


def unit_vector(vector):
    """ Returns the unit vector of the vector.  
    :param vector: vector to be normalised
    :return: normalised vector
    """
    return vector / np.dstack([np.linalg.norm(vector, axis=2)] * 3)


def angle_between(v1, v2):
    """
    Calculates angle between two vectors
    :param v1: first vector
    :param v2: second vector
    :return: angle between the vectors
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u[0, 0, :]), -1.0, 1.0))


def vis_ang_xr(d, observer, m_above = 5):
    """
    
    :param d: DataArray with altitudes
    :param observer: tuple (latitude,longtitude) of observer
    :return: DataArray with elavation angles from give obsarvation point
    """
    latp = np.radians(d.lat)
    lonp = np.radians(d.lon)
    la, lo = np.meshgrid(latp, lonp)

    lato = np.radians(observer[0])
    lono = np.radians(observer[1])
    alto = d.data[np.argmin(np.abs(latp - lato)), np.argmin(np.abs(lonp - lono))] + m_above

    altitudes = d.data  # .reshape(len(latp)*len(lonp)
    #print(np.shape(altitudes), np.shape(la), np.shape(lo))

    Z = geo_to_cart(la.T, lo.T, altitudes)
    O = geo_to_cart(lato, lono, alto)
    OZ = Z - O
    #print(np.shape(O))
    #print(np.shape(OZ))

    #print(np.shape(O))
    OZ = Z - O
    A = angle_between(OZ, O)
    an = d.copy()
    an.data = - np.degrees(A) + 180
    an.name = 'observation angle'

    an.attrs['observer'] = observer

    return an


##test call

# d = load_grt('CR_SR.grd')
# nejvyssi = d.where(d == d.max(), drop=True).squeeze()
# print(gerlach.lat, gerlach.lon)
#
# nejnizsi = d.where(d == d.min(), drop=True).squeeze()
# print(nizko.lat, nizko.lon)
#
# A, an = vis_ang_xr(d, (50, 15))
#
# print(an)
# fig, ax = plt.subplots(nrows=2, sharex=True, sharey=True)
# s = ax[0].imshow(-np.degrees(A) + 180, origin='bottom')
# plt.colorbar(s, ax=ax[0])
# r = ax[1].imshow(d.data, origin='bottom')
# plt.colorbar(r, ax=ax[1])
#print(A)