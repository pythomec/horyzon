# Horyzon

Horyzon computes field of view for any place on Earth and annotates objects and peaks visible on the horizon. 
Horyzon currently uses Global Multi-Resolution Topographic data (GMRT) from https://www.gmrt.org/ . 

## Getting Started


### Prerequisites

Horyzon depends on common libraries such as pandas, xarray, requests. Optionally, 
[pyresample](https://github.com/pytroll/pyresample) can be used to transform data to a different geographic projection. 

### Installing

Clone Horyzon from gitlab

```
git clone https://gitlab.com/jakub.seidl/horyzon 
```

and use GMRT MapTool https://www.gmrt.org/GMRTMapTool/ to download data in your region of interest. 
Save the grid data in GMT v3 Compatible NetCDF format. Automatic data download is envisaged for future versions.  
