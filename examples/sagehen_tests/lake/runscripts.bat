python ..\..\..\scripts\fishnet_generator.py -i test_parameters_lakes.ini --overwrite

python ..\..\..\scripts\hru_parameters.py -i test_parameters_lakes.ini
python ..\..\..\scripts\dem_parameters.py -i test_parameters_lakes.ini
python ..\..\..\scripts\dem_2_streams.py -i test_parameters_lakes.ini
python ..\..\..\scripts\crt_fill_parameters.py -i test_parameters_lakes.ini
python ..\..\..\scripts\dem_2_streams.py -i test_parameters_lakes.ini
python ..\..\..\scripts\crt_fill_parameters.py -i test_parameters_lakes.ini
python ..\..\..\scripts\stream_parameters.py -i test_parameters_lakes.ini
python ..\..\..\scripts\veg_parameters.py -i test_parameters_lakes.ini
python ..\..\..\scripts\soil_raster_prep.py -i test_parameters_lakes.ini
python ..\..\..\scripts\soil_parameters.py -i test_parameters_lakes.ini
python ..\..\..\scripts\prism_800m_normals.py -i test_parameters_lakes.ini
python ..\..\..\scripts\ppt_ratio_parameters.py -i test_parameters_lakes.ini
python ..\..\..\scripts\impervious_parameters.py -i test_parameters_lakes.ini
python ..\..\..\scripts\prms_template_fill.py -i test_parameters_lakes.ini