import numpy as np
import xarray as xr
from scipy.io import netcdf
from scipy.interpolate import interp2d
from scipy.ndimage import maximum_filter

def load_grt(fname, name = 'altitude', dlat = None, dlon = None):
    '''Load topographic data from grt file.

    :param fname: Path to grt file
    :param name: name of the returned DataArray
    :param dlat: requested resolution in latitude (None = do not resample)
    :param dlon: requested resolution in longitude (None = do not resample)
    :return: xarray dataset with lon, lat as coordinates
    '''

    vars = netcdf.netcdf_file(fname, mmap=False).variables

    z = vars['z'][:].copy().reshape(vars['dimension'][::-1])[::-1,:]
    lon = np.linspace(*vars['x_range'][:].copy(), vars['dimension'][0])
    lat = np.linspace(*vars['y_range'][:].copy(), vars['dimension'][1])

    # TODO: check if dlon and dlat is consistent with vars['spacing']

    d = xr.DataArray(z, coords={'lon': lon, 'lat': lat}, dims=['lat', 'lon'],
                     name=name)

    d = resample(d, dlat, dlon)

    return d

def resample(data, dx = None, dy = None):
    ''' Resample data using maximum filter

    :param data: xr.DataArray
    :param dx: resolution of the first dim
    :param dy: resolution of the second dim
    :return: resampled DataArray
    '''

    if dx is None and dy is None:
        return data

    dimx, dimy = data.dims

    data_dx = data.coords[dimx][1] - data.coords[dimx][0]
    data_dy = data.coords[dimy][1] - data.coords[dimy][0]

    if dx is None or dx < 2*data_dx:
        # no resampling in x
        filt_size_x = 1
        resample_dx = data_dx
    else:
        filt_size_x = int(np.ceil(dx / data_dx))
        filt_size_x += 1 - (filt_size_x % 2) # make it odd
        resample_dx = dx / filt_size_x

    if dy is None or dy < 2*data_dy:
        # no resampling in y
        filt_size_y = 1
        resample_dy = data_dy
    else:
        filt_size_y = int(np.ceil(dy / data_dy))
        filt_size_y += 1 - (filt_size_y % 2)  # make it odd
        resample_dy = dy / filt_size_y

    # resample
    if resample_dx != data_dx or resample_dy != data_dy:
        new_x = np.arange(data.coords[dimx][0], data.coords[dimx][-1]+resample_dx/2, resample_dx)
        new_y = np.arange(data.coords[dimy][0], data.coords[dimy][-1]+resample_dy/2, resample_dy)

        # TODO: replace with DataArray.interp after upgrading to newer xarrays
        z = interp2d(data.coords[dimy], data.coords[dimx], data)
        data = xr.DataArray(z(new_y, new_x), coords={dimx:new_x, dimy: new_y},
                            dims=[dimx, dimy], attrs = data.attrs)

    # use 2d maximum filter and downsample
    if filt_size_x !=1 or filt_size_y != 1:
        new_data = maximum_filter(data, size=(filt_size_x,filt_size_y),
                                  origin=(int((filt_size_x-1)/2), int((filt_size_y-1)/2)))[::filt_size_x, ::filt_size_y]
        data = xr.DataArray(new_data, coords={dimx:data.coords[dimx][::filt_size_x],
                                              dimy:data.coords[dimy][::filt_size_y]},
                            dims=[dimx, dimy], attrs=data.attrs)

    return data