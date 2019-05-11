import numpy as np
import xarray as xr
from scipy.io import netcdf

def load_grt(fname, name = 'altitude'):
    '''Load topographic data from grt file.

    :param fname: Path to grt file
    :return: xarray dataset with lon, lat as coordinates
    '''

    vars = netcdf.netcdf_file(fname, mmap=False).variables

    z = vars['z'][:].copy().reshape(vars['dimension'][::-1])[::-1,:]
    lon = np.linspace(*vars['x_range'][:].copy(), vars['dimension'][0])
    lat = np.linspace(*vars['y_range'][:].copy(), vars['dimension'][1])

    # TODO: dlon and dlat is consistent with vars['spacing']

    d = xr.DataArray(z, coords={'lon': lon, 'lat': lat}, dims=['lat', 'lon'],
                     name=name)

    return d