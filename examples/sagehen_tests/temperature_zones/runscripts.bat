python ..\..\..\scripts\fishnet_generator.py -i temp_zones_parameters.ini

python ..\..\..\scripts\hru_parameters.py -i temp_zones_parameters.ini
python ..\..\..\scripts\dem_parameters.py -i temp_zones_parameters.ini
python ..\..\..\scripts\dem_2_streams.py -i temp_zones_parameters.ini
python ..\..\..\scripts\crt_fill_parameters.py -i temp_zones_parameters.ini
python ..\..\..\scripts\dem_2_streams.py -i temp_zones_parameters.ini
python ..\..\..\scripts\crt_fill_parameters.py -i temp_zones_parameters.ini
python ..\..\..\scripts\stream_parameters.py -i temp_zones_parameters.ini
python ..\..\..\scripts\veg_parameters.py -i temp_zones_parameters.ini
python ..\..\..\scripts\soil_raster_prep.py -i temp_zones_parameters.ini
python ..\..\..\scripts\soil_parameters.py -i temp_zones_parameters.ini
python ..\..\..\scripts\prism_800m_normals.py -i temp_zones_parameters.ini
python ..\..\..\scripts\temp_adjust_parameters.py -i temp_zones_parameters.ini
python ..\..\..\scripts\impervious_parameters.py -i temp_zones_parameters.ini
python ..\..\..\scripts\prms_template_fill.py -i temp_zones_parameters.ini
