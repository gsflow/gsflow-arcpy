GSFlow ArcPy
============

Series of Python/ArcPy (ArcGIS) scripts for developing inputs for a GSFLOW model.

## Script Execution Order
- fishnet_generator.py
- hru_parameters.py
- dem_parameters.py
- veg_parameters.py
- soil_raster_prep.py
- soil_parameters.py
- impervious_parameters.py
- prism_4km_normals.py / prism_800m_normals.py
- ppt_ratio_parameters.py
- *Iterate to define the stream network*
  - dem_2_streams.py
  - crt_fill_parameters.py
- stream_parameters.py
- prms_template_fill.py

## Ancillary Data

Almost all of the following data can be downloaded for a study area using the [USGS Geospatial Data Gateway](http://datagateway.nrcs.usda.gov/).  They could probably also be downloaded using the [National Map Viewer](http://viewer.nationalmap.gov/viewer/), but this has not been tested.  Specific download instructions are provided for each dataset below.

#### Elevation

Elevation data is set using the 10m (1/3 arc-second) or 30m (1 arc-second) National Elevation Dataset (NED) rasters.  These can be easily downloaded in 1x1 degree tiles for the CONUS from the [USGS FTP](rockyftp.cr.usgs.gov) in the folder vdelivery/Datasets/Staged/Elevation.

#### LANDFIRE

Vegetation type and cover percent is set using [LANDFIRE](http://www.landfire.gov/).  The data can be downloaded in [tiles](http://www.landfire.gov/viewer/) or for the entire [CONUS](http://www.landfire.gov/lf_mosaics.php) (but only version 1.3.0?).  The default remap files were developed and tested for version 1.2.0 LANDFIRE 2010 version 1.2.0.

#### Soils

Available water capacity (AWC), percent sand, percent clay, and saturated hydraulic conductivity (Ksat) can be set using SSURGO or STATSGO.  The easiest way to acquire these data are through the [USGS Geospatial Data Gateway](http://datagateway.nrcs.usda.gov/).  Shapefiles of the soil properties can be extracted using the [NRCS Soil Data Viewer](http://www.nrcs.usda.gov/wps/portal/nrcs/detailfull/soils/home/?cid=nrcs142p2_053620).  The shapefiles then need to be converted to raster.

#### PRISM

PRISM precipitation, minimum temperature, and maximum temperature 30 year normals for the CONUS can be downloaded from the [PRISM site](http://www.prism.oregonstate.edu/normals/).

#### CRT

User must have [Cascade Routing Tool](http://water.usgs.gov/ogw/CRT/) (CRT) version 1.3.1

#### Remap files

Example ASCII remap files are provided, although it may be necessary to modify these to include new or missing LANDFIRE vegetation types.  In versions of ArcGIS before 10.2, you could have comments after the values (indicated by a /*) but this was removed in 10.2.  Now, comments must be on a separate line and begin with the "#" symbol.  The convert_remap_10p2.py script will convert the ArcGIS 10.1 style ASCII remap files to the ArcGIS 10.2 style.

## Requirements

+ Python 2.7
+ ArcPy