import numpy as np
from scipy.interpolate import RegularGridInterpolator, RectBivariateSpline
import xarray as xr

def data2polar(data, lon=None, lat=None, dtheta= 360/3600, dr=0.001, maxr=1):
    """Interpolate elevation angle from geographic coordinates to polar coordinates around the observer.

    :param data: DataArray with elevation angles
    :param lon: longitude of the viewpoint, if None use data.attrs['lon']
    :param lat: latitude of the viewpoint, if None use data.attrs['lat']
    :param dtheta: step in polar coordinate [°], default: 360/3600
    :param dr: step in radial coordinate
    :param maxr: maximum distance where to look for peaks
    :return: DataArray with elevation angles in polar coordinates

    """
    # TODO: dr should be real distance along Earth surface, independent of the original geographic coordinates
    # TODO: default dtheta and dr should be adjusted to grid size (see e.g. generate_polar_grid)
    # TODO: allow data in different coordinate system than lon/lat

    lon = lon or data.attrs['lon']
    lat = lat or data.attrs['lat']

    # create polar grid
    theta = np.arange(0, 360, dtheta)
    r = np.arange(dr, maxr + dr/2, dr)
    mesh_theta, mesh_r = np.meshgrid(theta, r)

    grid_x = lon + mesh_r * np.cos(np.deg2rad(mesh_theta))
    grid_y = lat + mesh_r * np.sin(np.deg2rad(mesh_theta))

    # interpolate to the grid
    d1, d2 = data.dims
    interp_angle = RectBivariateSpline(data[d1], data[d2], data)
    ang_in_polar = interp_angle(grid_x, grid_y, grid=False)
    ang_in_polar = xr.DataArray(ang_in_polar, coords={'r': r, 'theta': theta, d1:(('r','theta'),grid_x),
                                                      d2:(('r','theta'),grid_y)},
                                dims=['r', 'theta'])

    # set attributes
    ang_in_polar.attrs['lon'] = lon
    ang_in_polar.attrs['lat'] = lat

    return ang_in_polar


def get_mask_and_ridges(ang_in_polar):
    """Find visible points and ridges in polar coordinates

    :param ang_in_polar: DataArray with viewing angles in polar coordinates
    :param return_horizon: bool. Compute and return horizon? Default: False
    :return: (mask, edges), DataArrays, masks marking visible points and ridges

    """

    # compute visibility mask
    cummax = np.maximum.accumulate(ang_in_polar, axis=0)
    mask = ang_in_polar >= cummax

    # compute ridges
    ridges = (ang_in_polar.diff('r').values <= 0) & mask.values[:-1, :]
    r, theta = ang_in_polar.r, ang_in_polar.theta
    ridges = xr.DataArray(np.vstack((ridges, [0] * len(theta))), coords={'r': r, 'theta': theta},
                         dims=['r', 'theta'])

    return mask, ridges

def compute_horizon(mask, data):
    """Compute horizon and evaluate its coordinates and data along it

    :param mask: DataArray. Visibility mask
    :param data: DataArray. Data to be evaluated along horizon. Horizontal dimension should be 'theta'.
    :return: DataSet. Values (data, r, lon, lat) along the horizon versus polar angle.
    """

    # for each angle find index of the last non-zero point
    ind_last_visible = mask.shape[0] - np.argmax(mask.values[::-1,:], axis=0) - 1

    # evaluate data and its coordinates along the horizon
    horizon = {data.name:('theta',data.values[ind_last_visible, np.arange(data.shape[1])])}
    for c, val in data.coords.items():
        if c == 'theta':
            continue

        if val.ndim == 1:
            horizon[c] = ('theta', val.values[ind_last_visible])
        elif val.ndim == 2:
            horizon[c] = ('theta', val.values[ind_last_visible, np.arange(val.shape[1])])

    horizon = xr.Dataset(data_vars = horizon, coords={'theta': data.theta})

    return horizon

###############################################

# TODO: generate_polar_grid is obsolete -> remove after moving grid size estimates to ang2polar
def generate_polar_grid(x0, y0, x_grid, y_grid, dx=None, dy=None, dr=None, size='standard'):
    """
    Generate a polar grid around (x0,y0) point, covering the original cartesian grid.

    Angle coordinates are in degrees. Angle resolution is chosen so as to not skip any pixels
    of the original grid even at the largest radius.

    :param x0, y0: origin of the polar coordinates
    :param x_grid, y_grid: M-sized resp. N-sized array defining the rectangular grid
    :param dx, dy: rectangular grid resolution (auto-detected if omitted)
    :param dr: by default chosen to match the smaller of (dx, dy)
    :param size: {'outside', 'standard', 'inside'} defines how much the original grid is covered.
    :return: r_polar (R-sized), alpha_polar (A-sized), x_polar (RxA-sized), y_polar (RxA-sized)
    """

    if dx is None:
        dx = np.mean(np.diff(x_grid))
    if dy is None:
        dy = np.mean(np.diff(y_grid))

    # By default, dr is the same as the smaller grid size
    if dr is None:
        dr = min(dx, dy)

    # Horizontal and vertical distances to the original grid corners
    x_corners_dist = abs(x_grid[[0, -1]] - x0)
    y_corners_dist = abs(y_grid[[0, -1]] - y0)

    # The angular grid size is governed by the fartherst points from the observer
    d_alpha = min(np.arctan(dy / max(x_corners_dist)), np.arctan(dx / max(y_corners_dist))) / np.pi * 180.

    if size == 'outside':
        #  a) maximum variant - include even the farthest conrner
        max_r = np.sqrt(max(x_corners_dist) ** 2 + max(y_corners_dist) ** 2)
        # max_r = np.max(np.sqrt((xx-x0)**2 + (yy-y0)**2)) # using xx,yy from meshgrid
    elif size == 'standard':
        #  b) 'reasonable' variant - cut out some corners but include edges
        max_r = max(max(x_corners_dist), max(y_corners_dist))
    elif size == 'inside':
        #  c) minimum variant - fit within the grid borders
        max_r = min(min(x_corners_dist), min(y_corners_dist))
        # Here, the angular grid size is governed by the grid size at the closest border
        d_alpha = min(np.arctan(dy / min(x_corners_dist)), np.arctan(dx / min(y_corners_dist))) / np.pi * 180.

    # PRODUCE THE GRID
    r_polar = np.arange(0, max_r, dr)
    # Divide 360 degrees evenly, avoiding overlap of the last element
    alpha_polar = np.linspace(0, 360, int(360 / d_alpha))[:-1]

    aa, rr = np.meshgrid(alpha_polar, r_polar)

    # Calculate the X,Y locations of the grid
    x_polar = x0 + rr * np.sin(aa / 180. * np.pi)
    y_polar = y0 + rr * np.cos(aa / 180. * np.pi)

    return r_polar, alpha_polar, x_polar, y_polar


def interpolate_from_grid_to_any_points(x_grid, y_grid, data, x_out, y_out):
    """
    Interpolate `data` from an rectangular M x N grid to a set of points.

    :param x_grid, y_grid: M-sized resp. N-sized array defining the rectangular grid
    :param data: M x N array of data to be interpolated from
    :param x_out, y_out: target locations for the interpolation
    :return: data_polar ... the interpolated data of the same shape as x_out and y_out
    """

    data_interp = RegularGridInterpolator(points=(x_grid, y_grid), values=data, method='nearest', bounds_error=False)
    data_polar = data_interp((x_out, y_out))

    return data_polar   