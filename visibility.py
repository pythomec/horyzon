import numpy as np
from scipy.interpolate import RegularGridInterpolator


def find_visible_area(Z_polar):
    """
    The boolean mask of the visible area and the radial indices of this area boundary.

    Syntax:
    visible_mask, max_i_r = find_visible_area(Z_polar)
    """
    len_r, len_alpha = np.shape(Z_polar)
    current_max_Z = np.zeros(len_alpha)
    max_i_r = np.zeros(len_alpha, dtype=int)
    # visible_part = []
    visible_mask = np.empty_like(Z_polar, dtype=bool)
    for i_r, Z_circle in enumerate(Z_polar):
        visible_mask[i_r, :] = Z_circle >= current_max_Z
        if any(visible_mask[i_r, :]):
            #         visible_part.append((i_r, visible_mask[i_r, :]))
            current_max_Z[visible_mask[i_r, :]] = Z_circle[visible_mask[i_r, :]]
            max_i_r[visible_mask[i_r, :]] = i_r
    return visible_mask, max_i_r


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
    r_polar = np.arange(1, max_r, dr)
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