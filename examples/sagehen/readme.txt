Gsflow-Arcpy example model readme 

Included with the Sagehen example problem:

DEM - National elevation dataset 10m in NAD83 projection
Boundary shapefile - Sagehen watershed polygon shapefile - watershed.shp
Outlet point shapefile - Sagehen outlet point shapefile - model_points.shp
Flowlines - NHD stream lines clipped to Sagehen model boundary 
CRT executable (CRT1.4_beta.exe)
Vegetation data - LANDFIRE existing vegetation cover and type - us_140evc.img and us_140evt.img
Soil data - Rasters of available water capacity, clay percent, sand percent, and saturated hydraulic conductivity - awc.img, clay.img, sand.img, ksat.img
Climate data - PRISM 30 year monthly normals for precipitation and max/min temps
Impervious dataset - National land coverage database, impervious cover 2011 - nlcd2011_imp.img
Configuration file - Configuration file specific to Sagehen example model with all folders and fileneames 
Batch file to run all Gsflow-Arcpy scripts - runscripts.bat

Upon running the batch file, the Sagehen example model will run through Gsflow-Arcpy and output parameter files in the hru_params folder 
A control file is provided in the prms folder, and a data file is provided in the inputs folder, where the parameter files will need to be transferred to run GSFLOW
Control file - sagehen.control 
Data file - sagehen_datafile_1sta.data
Output data will be stored in the outputs folder, and users can compare their results to the original outputs found in the output_compare folder
GSFLOW executable (gsflow.exe) 

