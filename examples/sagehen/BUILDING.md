# Building the Sagehen Example

The following are the command line calls that were used to build the Sagehen example datasets.  Many of these processes could be done in ArcGIS and much of the data can be acquired more easily using the [USDA Geospatial Data Gateway](https://datagateway.nrcs.usda.gov/).

Many of these steps involve using the [GDAL](http://www.gdal.org/gdal_utilities.html) and [OGR](http://www.gdal.org/ogr_utilities.html) utilities which are available with the GDAL Python module.

Some of the larger CONUS scale datasets (LANDFIRE, NLCD impervious, etc.) can be downloaded once and then referenced as needed.

## Sagehen example folder

All of the steps in this document will assume that files are downloaded directly to the Sagehen example folder and that any command line calls will be executed from within the Sagehen folder.

## Sagehen Gage

The gage shapefile can be generated from the gage.geojson text file.  The coordinates for the gage location were defined based on the [USGS NWIS site description](https://waterdata.usgs.gov/nwis/inventory/?site_no=10343500).  The contents of the gage GeoJSON file are shown below.

```
{"type": "FeatureCollection", "name": "gage", "features": [ { "type": "Feature", "geometry": { "type": "Point", "coordinates": [ -120.23799763, 39.43156897, 0.0 ] } } ] }
```
To build the shapefile from the GeoJSON using ogr2ogr (the destination is first and the source is second when using ogr2ogr).
```
> ogr2ogr shapefiles\gage.shp shapefiles\gage.geojson
```

## NHD

The [NHDPlus high resolution datasets](https://nhd.usgs.gov/data.html) can be downloaded by HUC4.  The HUC4 file for the Sagehen area.  https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU8/HighResolution/Shape/NHD_H_16050102_Shape.zip

## Digital Elevation Model (DEM)

The example DEM can be derived from the 10m National Elevation Dataset (NED) by first downloading the 1x1 degree tile that covers the study area.  The 30m NED files can be acquired by changing the "13" to "1" in the URL.
+ [n40w121.zip](ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/Elevation/13/IMG/n40w121.zip).  After downloading the file, extract the contents into the "dem" folder.

For an area that intersect multiple NED 1x1 degree tiles, the tiles could be merged using the following command
```
gdalmerge something
```

The following command will clip (but not project) the DEM to the Sagehen fishnet:
```
gdalwarp -r "bilinear" -te -120.325 39.40 -120.23 39.468 -overwrite -ot Float32 -srcnodata None -dstnodata -3.4028234663852886e+38 -of HFA -co COMPRESSED=YES dem\imgn40w121_13.img dem\ned10m_nad83.img
```

## LANDFIRE

The [LANDFIRE](https://www.landfire.gov) existing vegetation cover and existing vegetation type rasters can be downloaded in 1x1 degree tiles using the [LANDFIRE Viewer](https://www.landfire.gov/viewer/) or the [LANDFIRE Data Access Tool](https://www.landfire.gov/datatool.php).  The full CONUS mosaics (3.3 GB) can also be downloaded directly from the following links.
+ [Existing Vegetation Cover (EVC) LF2014](https://www.landfire.gov/bulk/downloadfile.php?FNAME=US_140_mosaic-US_140EVC_12052016.zip&TYPE=landfire)
+ [Existing Vegetation Type (EVT) LF2014](https://www.landfire.gov/bulk/downloadfile.php?FNAME=US_140_mosaic-US_140EVT_04252017.zip&TYPE=landfire)

The example EVC and EVT images were clipped using the following command (assuming the downloaded images are in a "veg" sub-folder):
```
gdal_translate -r "nearest" -projwin -120.344 39.461 -120.210 39.406 -projwin_srs "EPSG:4326" -of HFA veg\US_140EVC_12052016\Grid\us_140evc veg\us_140evc.img
gdal_translate -r "nearest" -projwin -120.344 39.461 -120.210 39.406 -projwin_srs "EPSG:4326" -of HFA veg\US_140EVT_04252017\Grid\us_140evt veg\us_140evt.img
```

## NLCD Impervious

The [NLCD 2011](https://www.mrlc.gov/nlcd11_data.php) full CONUS mosaic (800 MB) can be downloaded from the following link:
+ [NLCD 2011 Percent Developed Imperviousness](http://www.landfire.gov/bulk/downloadfile.php?TYPE=nlcd2011&FNAME=nlcd_2011_impervious_2011_edition_2014_10_10.zip)

```
gdal_translate -r "nearest" -projwin -120.344 39.461 -120.210 39.406 -projwin_srs "EPSG:4326" -of HFA impervious\nlcd_2011_impervious_2011_edition_2014_10_10\nlcd_2011_impervious_2011_edition_2014_10_10.img impervious\nlcd2011_imp.img
```

## PRISM

http://prism.oregonstate.edu/normals/

Maximum temperature:
```
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmax_30yr_normal_800mM2_01_bil.bil prism\tmax\PRISM_tmax_30yr_normal_800mM2_01_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmax_30yr_normal_800mM2_02_bil.bil prism\tmax\PRISM_tmax_30yr_normal_800mM2_02_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmax_30yr_normal_800mM2_03_bil.bil prism\tmax\PRISM_tmax_30yr_normal_800mM2_03_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmax_30yr_normal_800mM2_04_bil.bil prism\tmax\PRISM_tmax_30yr_normal_800mM2_04_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmax_30yr_normal_800mM2_05_bil.bil prism\tmax\PRISM_tmax_30yr_normal_800mM2_05_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmax_30yr_normal_800mM2_06_bil.bil prism\tmax\PRISM_tmax_30yr_normal_800mM2_06_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmax_30yr_normal_800mM2_07_bil.bil prism\tmax\PRISM_tmax_30yr_normal_800mM2_07_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmax_30yr_normal_800mM2_08_bil.bil prism\tmax\PRISM_tmax_30yr_normal_800mM2_08_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmax_30yr_normal_800mM2_09_bil.bil prism\tmax\PRISM_tmax_30yr_normal_800mM2_09_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmax_30yr_normal_800mM2_10_bil.bil prism\tmax\PRISM_tmax_30yr_normal_800mM2_10_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmax_30yr_normal_800mM2_11_bil.bil prism\tmax\PRISM_tmax_30yr_normal_800mM2_11_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmax_30yr_normal_800mM2_12_bil.bil prism\tmax\PRISM_tmax_30yr_normal_800mM2_12_bil.bil
```

Minimum temperature:
```
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmin_30yr_normal_800mM2_01_bil.bil prism\tmin\PRISM_tmin_30yr_normal_800mM2_01_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmin_30yr_normal_800mM2_02_bil.bil prism\tmin\PRISM_tmin_30yr_normal_800mM2_02_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmin_30yr_normal_800mM2_03_bil.bil prism\tmin\PRISM_tmin_30yr_normal_800mM2_03_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmin_30yr_normal_800mM2_04_bil.bil prism\tmin\PRISM_tmin_30yr_normal_800mM2_04_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmin_30yr_normal_800mM2_05_bil.bil prism\tmin\PRISM_tmin_30yr_normal_800mM2_05_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmin_30yr_normal_800mM2_06_bil.bil prism\tmin\PRISM_tmin_30yr_normal_800mM2_06_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmin_30yr_normal_800mM2_07_bil.bil prism\tmin\PRISM_tmin_30yr_normal_800mM2_07_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmin_30yr_normal_800mM2_08_bil.bil prism\tmin\PRISM_tmin_30yr_normal_800mM2_08_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmin_30yr_normal_800mM2_09_bil.bil prism\tmin\PRISM_tmin_30yr_normal_800mM2_09_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmin_30yr_normal_800mM2_10_bil.bil prism\tmin\PRISM_tmin_30yr_normal_800mM2_10_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmin_30yr_normal_800mM2_11_bil.bil prism\tmin\PRISM_tmin_30yr_normal_800mM2_11_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_tmin_30yr_normal_800mM2_12_bil.bil prism\tmin\PRISM_tmin_30yr_normal_800mM2_12_bil.bil
```

Precipitation:
```
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_ppt_30yr_normal_800mM2_01_bil.bil prism\ppt\PRISM_ppt_30yr_normal_800mM2_01_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_ppt_30yr_normal_800mM2_02_bil.bil prism\ppt\PRISM_ppt_30yr_normal_800mM2_02_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_ppt_30yr_normal_800mM2_03_bil.bil prism\ppt\PRISM_ppt_30yr_normal_800mM2_03_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_ppt_30yr_normal_800mM2_04_bil.bil prism\ppt\PRISM_ppt_30yr_normal_800mM2_04_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_ppt_30yr_normal_800mM2_05_bil.bil prism\ppt\PRISM_ppt_30yr_normal_800mM2_05_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_ppt_30yr_normal_800mM2_06_bil.bil prism\ppt\PRISM_ppt_30yr_normal_800mM2_06_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_ppt_30yr_normal_800mM2_07_bil.bil prism\ppt\PRISM_ppt_30yr_normal_800mM2_07_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_ppt_30yr_normal_800mM2_08_bil.bil prism\ppt\PRISM_ppt_30yr_normal_800mM2_08_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_ppt_30yr_normal_800mM2_09_bil.bil prism\ppt\PRISM_ppt_30yr_normal_800mM2_09_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_ppt_30yr_normal_800mM2_10_bil.bil prism\ppt\PRISM_ppt_30yr_normal_800mM2_10_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_ppt_30yr_normal_800mM2_11_bil.bil prism\ppt\PRISM_ppt_30yr_normal_800mM2_11_bil.bil
gdal_translate -r "nearest" -projwin -120.328 39.47 -120.222 39.397 -projwin_srs "EPSG:4326" -of HFA prism\PRISM_ppt_30yr_normal_800mM2_12_bil.bil prism\ppt\PRISM_ppt_30yr_normal_800mM2_12_bil.bil
```
