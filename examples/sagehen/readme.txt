Gsflow-Arcpy example model readme

Included with the Sagehen example problem:

DEM - National elevation dataset 10m in NAD83 projection
Boundary shapefile - Sagehen watershed polygon shapefile - watershed.shp
Outlet point shapefile - Sagehen outlet point shapefile - model_points.shp
Flowlines - NHD stream lines clipped to Sagehen model boundary 
Vegetation data - LANDFIRE existing vegetation cover and type - us_140evc.img and us_140evt.img
Soil data - Rasters of available water capacity, clay percent, sand percent, and saturated hydraulic conductivity - awc.img, clay.img, sand.img, ksat.img
Climate data - PRISM 30 year monthly normals for precipitation and max/min temps
Impervious dataset - National land coverage database, impervious cover 2011 - nlcd2011_imp.img
Configuration file - Configuration file specific to Sagehen example model with all folders and fileneames set
Batch file to run all Gsflow-Arcpy scripts - runscripts.bat

Upon running the batch file, the Sagehen example model will run through Gsflow-Arcpy and output parameter files to the inputs folder (located in the prms folder) 
A control file is also pesent in the prms folder, and a data file is provided in the inputs folder
Control file - sagehen.control 
Data file - sagehen_datafile_1sta.data
GSFLOW executable should be downloaded from the USGS webpage - https://water.usgs.gov/ogw/gsflow/
