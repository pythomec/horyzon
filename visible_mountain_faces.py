import numpy as np
import matplotlib.pyplot as plt
import scipy
import xarray as xr


Rz=6.378e6


def direction_to_px(d, observer):
    """
    Calculates direction from the observer to all pixels in the grid
    :param d: DataArray with altitudes as a function of lat and lon
    :param observer: tuple of observer coordinates (lat,lon)
    :return: DataArray with direction vectors as a function of coordinates
    """
    latp = np.radians(d.lat)
    lonp = np.radians(d.lon)
    la, lo = np.meshgrid(latp, lonp)
    lato = np.radians(observer[0])
    lono = np.radians(observer[1])
    d = np.dstack((la.T - lato, lo.T - lono))
    dire = d / np.dstack([np.linalg.norm(d, axis=2)] * 2)
    directions = xr.DataArray(dire, coords={'lat': latp, 'lon': lonp, 'comp': [1, 2]}, dims=('lat', 'lon', 'comp'))
    print(np.shape(d), np.shape(dire), np.shape(d))
    return directions


def solar_shades(d, sunphi, sunelev):
    """
    WIP: Sinmplified calculation of illuminated mountain faces - projection of surface graidient to the direction of sunrayss

    :param d: DataArray with altitudes as a function of lat and lon
    :param sunphi: azimuth of the sun 0 is the local meridian
    :param sunelev: elevation angle of the sun
    :return: 
    """""
    u, v = np.gradient(d.data)
    grad = np.sqrt(u ** 2 + v ** 2)
    gradx = u / np.nanmax(grad)
    grady = v / np.nanmax(grad)
    sunx = gradx * np.sin(sunelev) * np.cos(sunphi)  # is the sun elevation correctly implemented?
    # - scalar product would be better
    suny = grady * np.sin(sunelev) * np.sin(sunphi)
    sun = d.copy()
    sun.data = sunx + suny
    return sun


def gradient_vision_direction(d, observer):
    """
    Calculatiats the projection of the surface gradients to the direction of the vision of the observer
    flat Earth aprox ...

    :param d: DataArray with altitudes as a function of lat and lon
    :param observer: tuple of observer coordinates (lat,lon)
    :return: DataArray with gradiant value in the direction
    """

    u, v = np.gradient(d.data)
    grad = np.sqrt(u ** 2 + v ** 2)
    gradx = u / np.nanmax(grad)
    grady = v / np.nanmax(grad)
    dirangle_xa = direction_to_px(d, observer)
    print(dirangle_xa)
    dirangle = dirangle_xa.data
    visgrad = d.copy()
    df = gradx * dirangle[:, :, 0] + grady * dirangle[:, :, 1]
    visgrad.data = df
    # plt.imshow(df)
    return visgrad


# test calls
# d = load_grt('CR_SR.grd')
# visgrad = gradient_vision_direction(d, (49.1, 20))
# limits = [d.lon[0], d.lon[-1], d.lat[0], d.lat[-1]]
# phis = np.linspace(-np.pi / 2.0, np.pi / 2.0, 2)
#
# fig2, ax2 = plt.subplots()
# ax2.imshow(visgrad.data, origin='bottom', vmin=0, vmax=np.nanmax(visgrad.data), extent=limits, aspect='equal')
# ax2.scatter([20], [49.1], s=10, color='r')
#
# i = 0
# for phi in phis:
#     i = i + 1
#     limits = [19, 21, 48.5, 49.8]
#     sun = solar_shades(d.loc[48.5:49.8, 19:21], phi, np.cos(phi))
#     figs, ax = plt.subplots(figsize=(15, 10))
#     ax.imshow(sun.data, origin='bottom', vmin=0, vmax=np.nanmax(sun.data), extent=limits, cmap=plt.get_cmap('jet'),
#               aspect='equal')  # ,alpha=0.6)
#     plt.tight_layout()
#     plt.savefig('Movie/JET_sun_phi_{}_Tatry.png'.format(i))
