# Horyzon

Horyzon is a topographic library which can compute field of view for any place on Earth and make panoramatic plots such 
as this view on Tatra Mountains in Slovakia: 

![Example output](docs/Slovakia_Low_and_High_Tatras.png) 

Horyzon currently uses Global Multi-Resolution Topographic data (GMRT) from https://www.gmrt.org/ .

Horyzon is in an early stage of development, its API is unstable and the features are not yet optimized for speed. 

## Getting Started


### Prerequisites

Horyzon depends on common libraries such as pandas, xarray, requests. Optionally, 
[pyresample](https://github.com/pytroll/pyresample) can be used to transform data to a different geographic projection. 

### Installing

Clone Horyzon from gitlab

```
git clone https://github.com/pythomec/horyzon.git
```

and use GMRT MapTool https://www.gmrt.org/GMRTMapTool/ to download data in your region of interest. 
Save the grid data in GMT v3 Compatible NetCDF format. Automatic data download is envisaged for future versions.

### Usage

1. Download topographic data from https://www.gmrt.org/GMRTMapTool/
2. Load the data as xarray.DataArray

 ```
 from horyzon import load_grt, Viewpoint  

 alt = load_grt(path_to_the_data)
 ```
 
3. Create Viewpoint object (it may take some time)

 ```
 lon, lat = 20, 49
 view = Viewpoint(altitude=alt, lon=lon, lat=lat)
 ```
 
4. Plot panorama looking towards North

```
view.plot_panorama(direction='N') 
```

5. And voilà ...

![High Tatras in the distance](docs/Slovakia_High_Tatras_in_distance.png)

6. You can also plot field of view, e.g. on the altitude map

```
import pylab as plt

plt.figure()
view.altitude.plot()
view.plot_horizon_lonlat('r') 
```