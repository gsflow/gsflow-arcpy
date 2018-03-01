#--------------------------------
# Name:         hru_parameters.py
# Purpose:      GSFLOW HRU parameters
# Notes:        ArcGIS 10.2 Version
# Python:       2.7
#--------------------------------

import argparse
import ConfigParser
import datetime as dt
import logging
import os
import sys

import arcpy
from arcpy import env

import support_functions as support


def hru_parameters(config_path, overwrite_flag=False, debug_flag=False):
    """Calculate GSFLOW HRU Parameters

    Args:
        config_file (str): Project config file path
        ovewrite_flag (bool): if True, overwrite existing files
        debug_flag (bool): if True, enable debug level logging

    Returns:
        None
    """

    # Initialize hru parameters class
    hru = support.HRUParameters(config_path)

    # Open input parameter config file
    inputs_cfg = ConfigParser.ConfigParser()
    try:
        inputs_cfg.readfp(open(config_path))
    except Exception as e:
        logging.error(
            '\nERROR: Config file could not be read, '
            'is not an input file, or does not exist\n'
            '  config_file = {}\n'
            '  Exception: {}\n'.format(config_path, e))
        sys.exit()

    # Log DEBUG to file
    log_file_name = 'hru_parameters_log.txt'
    log_console = logging.FileHandler(
        filename=os.path.join(hru.log_ws, log_file_name), mode='w')
    log_console.setLevel(logging.DEBUG)
    log_console.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger('').addHandler(log_console)
    logging.info('\nGSFLOW HRU Parameters')

    # Read parameters from config file
    study_area_orig_path = inputs_cfg.get('INPUTS', 'study_area_path')
    try:
        set_lake_flag = inputs_cfg.getboolean('INPUTS', 'set_lake_flag')
    except ConfigParser.NoOptionError:
        set_lake_flag = False
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'set_lake_flag', set_lake_flag))

    if set_lake_flag:
        lake_orig_path = inputs_cfg.get('INPUTS', 'lake_path')
        lake_zone_field = inputs_cfg.get('INPUTS', 'lake_zone_field')
        lake_area_pct = inputs_cfg.getfloat('INPUTS', 'lake_area_pct')

    # Model points
    model_inputs_path = inputs_cfg.get('INPUTS', 'model_points_path')
    try:
        model_points_zone_field = inputs_cfg.get(
            'INPUTS', 'model_points_zone_field')
    except:
        model_points_zone_field = 'FID'
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'model_points_zone_field', model_points_zone_field))
    try:
        model_points_type_field = inputs_cfg.get(
            'INPUTS', 'model_points_type_field')
    except:
        model_points_type_field = 'TYPE'
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'model_points_type_field', model_points_type_field))

    # Control flags
    try:
        calc_flow_acc_dem_flag = inputs_cfg.getboolean(
            'INPUTS', 'calc_flow_acc_dem_flag')
    except ConfigParser.NoOptionError:
        calc_flow_acc_dem_flag = False
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'calc_flow_acc_dem_flag', calc_flow_acc_dem_flag))

    try:
        calc_topo_index_flag = inputs_cfg.getboolean(
            'INPUTS', 'calc_topo_index_flag')
    except ConfigParser.NoOptionError:
        calc_topo_index_flag = False
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'calc_topo_index_flag', calc_topo_index_flag))

    # try:
    #     set_ppt_zones_flag = inputs_cfg.getboolean(
    #         'INPUTS', 'set_ppt_zones_flag')
    # except ConfigParser.NoOptionError:
    #     set_ppt_zones_flag = False
    #     logging.info(
    #         '  Missing INI parameter, setting {} = {}'.format(
    #             'set_ppt_zones_flag', set_ppt_zones_flag))

    # Check input paths
    if not arcpy.Exists(hru.polygon_path):
        logging.error(
            '\nERROR: Fishnet ({}) does not exist'.format(
                hru.polygon_path))
        sys.exit()

    if set_lake_flag:
        if not arcpy.Exists(lake_orig_path):
            logging.error(
                '\nERROR: Lake layer ({}) does not exist'.format(
                    lake_orig_path))
            sys.exit()
        # lake_path must be a polygon shapefile
        if arcpy.Describe(lake_orig_path).datasetType != 'FeatureClass':
            logging.error(
                '\nERROR: lake_path must be a polygon shapefile')
            sys.exit()
        # Check lake_zone_field
        if lake_zone_field.upper() in ['', 'FID', 'NONE']:
            lake_zone_field = arcpy.Describe(lake_orig_path).OIDFieldName
            logging.warning(
                '\n  NOTE: Using {} to set {}\n'.format(
                    lake_zone_field, hru.lake_id_field))
        elif not arcpy.ListFields(lake_orig_path, lake_zone_field):
            logging.error(
                '\nERROR: lake_zone_field field {} does not exist\n'.format(
                    lake_zone_field))
            sys.exit()
        # Need to check that lake_zone_field is an int type
        elif not [f.type for f in arcpy.Describe(lake_orig_path).fields
                  if (f.name == lake_zone_field and
                      f.type in ['SmallInteger', 'Integer'])]:
            logging.error(
                '\nERROR: lake_zone_field field {} must be an '
                'integer type\n'.format(lake_zone_field))
            sys.exit()

    # Check model points
    if not os.path.isfile(model_inputs_path):
        logging.error(
            '\nERROR: Model points shapefiles does not exist'
            '\nERROR:   {}'.format(model_inputs_path))
        sys.exit()
    # model_points_path must be a point shapefile
    elif arcpy.Describe(model_inputs_path).datasetType != 'FeatureClass':
        logging.error(
            '\nERROR: model_points_path must be a point shapefile')
        sys.exit()

    # For now, study area has to be a polygon
    if arcpy.Describe(study_area_orig_path).datasetType != 'FeatureClass':
        logging.error(
            '\nERROR: For now, study area must be a polygon shapefile')
        sys.exit()


    # Build output folder if necessary
    hru_temp_ws = os.path.join(hru.param_ws, 'hru_temp')
    if not os.path.isdir(hru_temp_ws):
        os.mkdir(hru_temp_ws)
    # Output paths
    study_area_path = os.path.join(hru_temp_ws, 'study_area.shp')
    lake_path = os.path.join(hru_temp_ws, 'lakes.shp')
    lake_clip_path = os.path.join(hru_temp_ws, 'lake_clip.shp')
    model_points_path = os.path.join(hru_temp_ws, 'model_points.shp')


    # Set ArcGIS environment variables
    arcpy.CheckOutExtension('Spatial')
    env.overwriteOutput = True
    env.pyramid = 'PYRAMIDS -1'
    # env.pyramid = 'PYRAMIDS 0'
    env.workspace = hru.param_ws
    env.scratchWorkspace = hru.scratch_ws

    # Create HRU points at polygon centroids
    if not arcpy.Exists(hru.point_path):
        logging.info('\n  Building HRU point shapefile')
        # FeatureToPoint will copy all fields in hru.polygon_path
        # arcpy.FeatureToPoint_management(
        #    hru.polygon_path, hru.point_path)
        # Build point_path directly
        arcpy.CreateFeatureclass_management(
            os.path.dirname(hru.point_path),
            os.path.basename(hru.point_path), 'POINT')
        arcpy.DefineProjection_management(hru.point_path, hru.sr)
        arcpy.AddField_management(
            hru.point_path, hru.fid_field, 'LONG')
        hru_centroid_list = [
            row for row in arcpy.da.SearchCursor(
                hru.polygon_path, ['OID@', 'SHAPE@XY'])]
        with arcpy.da.InsertCursor(
                hru.point_path,
                ['OID@', 'SHAPE@XY', hru.fid_field]) as update_c:
            for hru_centroid in hru_centroid_list:
                update_c.insertRow(
                    [hru_centroid[0], hru_centroid[1], hru_centroid[0]])
        del hru_centroid_list
    # Check existing HRU points
    else:
        # Remove any extra fields
        field_remove_list = [
            f.name for f in arcpy.ListFields(hru.point_path)
            if f.name not in ['FID', 'Shape', hru.fid_field]]
        # Skip if there is only one field in the shapefile
        if field_remove_list and len(field_remove_list) > 1:
            logging.info('\n  Removing HRU point fields')
            for field in field_remove_list:
                logging.debug('    {}'.format(field))
                try:
                    arcpy.DeleteField_management(hru.point_path, field)
                except Exception as e:
                    logging.debug('    Unhandled exception: {}'.format(e))
                    continue
        # Save original FID
        if len(arcpy.ListFields(hru.point_path, hru.fid_field)) == 0:
            arcpy.AddField_management(
                hru.point_path, hru.fid_field, 'LONG')
        arcpy.CalculateField_management(
            hru.point_path, hru.fid_field, '!FID!', 'PYTHON')
        if len(arcpy.ListFields(hru.point_path, 'Id')) > 0:
            arcpy.DeleteField_management(hru.point_path, 'Id')
        del field_remove_list

    # Add all output fields
    logging.info('\nAdding fields if necessary')
    logging.info(
        '  Note: You may see duplicate field names when writing to a network '
        'drive')

    # HRU/DEM Fields
    support.add_field_func(hru.polygon_path, hru.fid_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.id_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.type_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.dem_mean_field, 'DOUBLE')
    #support.add_field_func(hru.polygon_path, hru.dem_median_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.dem_min_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.dem_max_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.dem_adj_field, 'DOUBLE')
    if calc_flow_acc_dem_flag:
        support.add_field_func(hru.polygon_path, hru.dem_flowacc_field, 'DOUBLE')
        support.add_field_func(hru.polygon_path, hru.dem_sum_field, 'DOUBLE')
        support.add_field_func(hru.polygon_path, hru.dem_count_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.dem_sink_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.crt_elev_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.crt_fill_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.dem_aspect_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.dem_slope_deg_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.dem_slope_rad_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.dem_slope_pct_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.area_field, 'DOUBLE')
    if calc_topo_index_flag:
        support.add_field_func(hru.polygon_path, hru.topo_index_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.row_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.col_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.x_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.y_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.lat_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.lon_field, 'DOUBLE')

    # Lake fields
    support.add_field_func(hru.polygon_path, hru.lake_id_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.lake_area_field, 'DOUBLE')

    # Stream fields
    support.add_field_func(hru.polygon_path, hru.iseg_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.irunbound_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.flow_dir_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.krch_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.irch_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.jrch_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.iseg_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.reach_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.rchlen_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.maxreach_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.outseg_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.iupseg_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.subbasin_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.segbasin_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.outflow_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.strm_top_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.strm_slope_field, 'DOUBLE')

    # Sink field
    support.add_field_func(hru.polygon_path, hru.hru_sink_field, 'LONG')

    # PPT Zone fields
    # if set_ppt_zones_flag:
    support.add_field_func(hru.polygon_path, hru.ppt_zone_id_field, 'SHORT')
    support.add_field_func(hru.polygon_path, hru.hru_psta_field, 'SHORT')

    # DEM based
    support.add_field_func(hru.polygon_path, hru.jh_tmax_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.jh_tmin_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.jh_coef_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.snarea_thresh_field, 'DOUBLE')

    # Aspect based
    support.add_field_func(hru.polygon_path, hru.tmax_adj_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.tmin_adj_field, 'DOUBLE')

    # Vegetation fields
    support.add_field_func(hru.polygon_path, hru.cov_type_field, 'SHORT')
    support.add_field_func(hru.polygon_path, hru.covden_sum_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.covden_win_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.rad_trncf_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.snow_intcp_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.srain_intcp_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.wrain_intcp_field, 'DOUBLE')

    # Soil fields
    support.add_field_func(hru.polygon_path, hru.awc_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.clay_pct_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.sand_pct_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.ksat_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.soil_type_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.soil_root_max_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.moist_init_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.moist_max_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.rechr_init_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.rechr_max_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.ssr2gw_rate_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.ssr2gw_k_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.slowcoef_lin_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.slowcoef_sq_field, 'DOUBLE')

    # Impervious fields
    support.add_field_func(hru.polygon_path, hru.imperv_pct_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.carea_max_field, 'DOUBLE')

    # PRISM mean monthly fields
    month_list = ['{:02d}'.format(m) for m in range(1, 13)]
    # month_list.extend(['14'])
    for prism_data_name in ['PPT', 'TMAX', 'TMIN']:
        for month in month_list:
            support.add_field_func(
                hru.polygon_path,
                '{}_{}'.format(prism_data_name, month), 'DOUBLE')
    # PRISM mean monthly PPT ratio fields
    for month in month_list:
        if month == '14':
            continue
        support.add_field_func(
            hru.polygon_path, 'PPT_RT_{}'.format(month), 'DOUBLE')

    # Id field is added by default to new fishnets
    if arcpy.ListFields(hru.polygon_path, 'Id'):
        arcpy.DeleteField_management(hru.polygon_path, 'Id')

    logging.info('\nCalculating parameters')
    # Keep original FID for subsetting in zonal stats
    logging.info('  Saving original HRU FID to {}'.format(
        hru.fid_field))
    arcpy.CalculateField_management(
        hru.polygon_path, hru.fid_field, '!FID!', 'PYTHON')

    # Cell X/Y
    logging.info('  Calculating cell X/Y')
    cell_xy_func(hru.polygon_path, hru.x_field, hru.y_field)

    # Create unique ID, start at top left corner, work down rows
    # Row/Col numbered from top left corner (1's based numbering)
    logging.info('  Calculating cell ID/row/col')
    cell_id_col_row_func(
        hru.polygon_path, hru.id_field, hru.col_field, hru.row_field,
        hru.extent, hru.cs)

    # Cell Lat/Lon
    logging.info('  Calculating cell lat/lon')
    cell_lat_lon_func(
        hru.polygon_path, hru.lat_field, hru.lon_field, hru.sr.GCS)

    # Cell Area
    logging.info('  Calculating cell area (acres)')
    arcpy.CalculateField_management(
        hru.polygon_path, hru.area_field, '!SHAPE.AREA@acres!', 'PYTHON')

    # Reset HRU_TYPE
    logging.info('\nResetting {} to 0'.format(hru.type_field))
    arcpy.CalculateField_management(
        hru.polygon_path, hru.type_field, 0, 'PYTHON')
    # Reset LAKE_ID
    if set_lake_flag:
        logging.info('Resetting {} to 0'.format(hru.lake_id_field))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.lake_id_field, 0, 'PYTHON')
    if set_lake_flag:
        logging.info('Resetting {} to 0'.format(hru.lake_area_field))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.lake_area_field, 0, 'PYTHON')

    # Calculate HRU Type
    logging.info('\nCalculating cell HRU Type')
    study_area_desc = arcpy.Describe(study_area_orig_path)
    study_area_sr = study_area_desc.spatialReference
    logging.debug('  Study area: {}'.format(study_area_orig_path))
    logging.debug('  Study area spat. ref.:  {}'.format(
        study_area_sr.name))
    logging.debug('  Study area GCS:         {}'.format(
        study_area_sr.GCS.name))
    # If study area spat_ref doesn't match hru_param spat_ref
    # Project study area to hru_param spat ref
    # Otherwise, read study_area directly
    if hru.sr.name != study_area_sr.name:
        logging.info('  Projecting study area...')
        # Set preferred transforms
        transform_str = support.transform_func(hru.sr, study_area_sr)
        logging.debug('    Transform: {}'.format(transform_str))
        # Project study area shapefile
        arcpy.Project_management(
            study_area_orig_path, study_area_path, hru.sr,
            transform_str, study_area_sr)
        del transform_str
    else:
        arcpy.Copy_management(study_area_orig_path, study_area_path)
    support.zone_by_centroid_func(
        study_area_path, hru.type_field, 1,
        hru.polygon_path, hru.point_path, hru)

    # Calculate HRU Type for lakes (HRU_TYPE = 2)
    if set_lake_flag:
        logging.info('\nCalculating cell HRU Type & ID for lakes')
        lake_layer = 'lake_layer'
        lake_desc = arcpy.Describe(lake_orig_path)
        lake_sr = lake_desc.spatialReference
        logging.debug('  Lakes: {}'.format(lake_orig_path))
        logging.debug('  Lakes spat. ref.:  {}'.format(lake_sr.name))
        logging.debug('  Lakes GCS:         {}'.format(lake_sr.GCS.name))

        # If lakes spat_ref doesn't match hru_param spat_ref
        # Project lakes to hru_param spat ref
        # Otherwise, read lakes directly
        if hru.sr.name != lake_sr.name:
            logging.info('  Projecting lakes...')
            # Set preferred transforms
            transform_str = support.transform_func(hru.sr, lake_sr)
            logging.debug('    Transform: {}'.format(transform_str))
            # Project lakes shapefile
            arcpy.Project_management(
                lake_orig_path, lake_path, hru.sr, transform_str, lake_sr)
            arcpy.MakeFeatureLayer_management(lake_path, lake_layer)
            del lake_path, transform_str
        else:
            arcpy.MakeFeatureLayer_management(
                lake_orig_path, lake_layer)

        # Clip lakes by study area after projecting lakes
        logging.info('  Clipping lakes...')
        arcpy.Clip_analysis(lake_layer, study_area_path, lake_clip_path)
        # Remove all unnecesary fields
        for field in arcpy.ListFields(lake_clip_path):
            if field.name not in [lake_zone_field, 'Shape']:
                try:
                    arcpy.DeleteField_management(lake_clip_path, field.name)
                except Exception as e:
                    logging.debug('    Unhandled exception: {}'.format(e))
                    continue

        # Set lake HRU_TYPE
        logging.info('  Setting lake {}'.format(hru.type_field))
        support.zone_by_area_func(
            lake_clip_path, hru.type_field, 2,
            hru.polygon_path, hru, hru.area_field,
            hru.lake_area_field, lake_area_pct)
        # Set lake ID
        logging.info('  Setting {}'.format(hru.lake_id_field))
        support.zone_by_area_func(
            lake_clip_path, hru.lake_id_field, lake_zone_field,
            hru.polygon_path, hru, hru.area_field,
            hru.lake_area_field, lake_area_pct)
        # Cleanup
        del lake_layer, lake_desc, lake_sr


    # Read in model points shapefile
    logging.info('\nChecking model points shapefile')
    model_points_desc = arcpy.Describe(model_inputs_path)
    model_points_sr = model_points_desc.spatialReference
    logging.debug('  Points: {}'.format(model_inputs_path))
    logging.debug('  Points spat. ref.:  {}'.format(model_points_sr.name))
    logging.debug('  Points GCS:         {}'.format(model_points_sr.GCS.name))

    # If model points spat_ref doesn't match hru_param spat_ref
    # Project model points to hru_param spat ref
    # Otherwise, read model points directly
    if hru.sr.name != model_points_sr.name:
        logging.info(
            '  Model points projection does not match fishnet.\n'
            '  Projecting model points.\n')
        # Set preferred transforms
        transform_str = support.transform_func(hru.sr, model_points_sr)
        logging.debug('    Transform: {}'.format(transform_str))
        arcpy.Project_management(
            model_inputs_path, model_points_path,
            hru.sr, transform_str, model_points_sr)
    else:
        arcpy.Copy_management(model_inputs_path, model_points_path)
    model_points_lyr = 'model_points_lyr'
    arcpy.MakeFeatureLayer_management(model_points_path, model_points_lyr)

    # Check model point types
    logging.info('  Checking model point types')
    model_point_types = [str(r[0]).upper() for r in arcpy.da.SearchCursor(
        model_points_path, [model_points_type_field])]
    if not set(model_point_types).issubset(set(['OUTLET', 'SUBBASIN', 'SWALE'])):
        logging.error('\nERROR: Unsupported model point type(s) found, exiting')
        logging.error('\n  Model point types: {}\n'.format(model_point_types))
        sys.exit()
    # elif 'OUTLET' not in model_point_types and 'SWALE' not in model_point_types:
    elif not set(model_point_types).issubset(set(['OUTLET', 'SWALE'])):
        logging.error(
            '\nERROR: At least one model point must be an OUTLET or SWALE, '
            'exiting\n')
        sys.exit()
    else:
        logging.debug('  {}'.format(', '.join(model_point_types)))

    if 'SWALE' in model_point_types:
        arcpy.SelectLayerByAttribute_management(
            model_points_lyr, 'NEW_SELECTION', '"TYPE" = \'SWALE\'')

        logging.info('  Setting swale (sink) cells to {}=3'.format(
            hru.type_field))
        hru_polygon_lyr = 'hru_polygon_lyr'
        arcpy.MakeFeatureLayer_management(hru.polygon_path, hru_polygon_lyr)
        arcpy.SelectLayerByAttribute_management(hru_polygon_lyr, 'CLEAR_SELECTION')
        arcpy.SelectLayerByLocation_management(
            hru_polygon_lyr, 'INTERSECT', model_points_lyr)
        arcpy.CalculateField_management(
            hru_polygon_lyr, hru.type_field, 3, 'PYTHON')
        arcpy.SelectLayerByAttribute_management(hru_polygon_lyr, 'CLEAR_SELECTION')
        arcpy.SelectLayerByAttribute_management(model_points_lyr, 'CLEAR_SELECTION')
        arcpy.Delete_management(hru_polygon_lyr)
        del hru_polygon_lyr
    arcpy.Delete_management(model_points_lyr)
    del model_points_lyr


    # Setting HRU_PSTA to default value of 1
    if all([row[0] == 0 for row in arcpy.da.SearchCursor(
            hru.polygon_path, [hru.hru_psta_field])]):
        logging.info('Setting {} to 1'.format(
            hru.hru_psta_field))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.hru_psta_field, '1', 'PYTHON')

    # Cleanup
    del study_area_desc, study_area_sr


def cell_xy_func(hru_param_path, x_field, y_field):
    """"""
    fields = ('SHAPE@XY', x_field, y_field)
    with arcpy.da.UpdateCursor(hru_param_path, fields) as u_cursor:
        for row in u_cursor:
            row[1], row[2] = row[0]
            u_cursor.updateRow(row)
            del row


def cell_lat_lon_func(hru_param_path, lat_field, lon_field, gcs_sr):
    """"""
    fields = ('SHAPE@XY', lon_field, lat_field)
    with arcpy.da.UpdateCursor(hru_param_path, fields, '', gcs_sr) as u_cursor:
        for row in u_cursor:
            row[1], row[2] = row[0]
            u_cursor.updateRow(row)
            del row


def cell_id_col_row_func(hru_param_path, id_field, col_field, row_field,
                         extent, cs):
    """"""
    fields = ('SHAPE@XY', col_field, row_field, id_field)
    with arcpy.da.UpdateCursor(hru_param_path, fields) as u_cursor:
        num_cols = (extent.XMax - extent.XMin) / cs
        for row in u_cursor:
            #  Row/Col are 1's based indices
            row[1] = ((row[0][0] - extent.XMin) // cs) + 1
            row[2] = ((extent.YMax - row[0][1]) // cs) + 1
            #  Create unique ID, start at top left corner, work down rows
            row[3] = row[1] + (row[2] - 1) * num_cols
            u_cursor.updateRow(row)
            del row


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='HRU Parameters',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i', '--ini', required=True,
        help='Project input file', metavar='PATH')
    parser.add_argument(
        '-o', '--overwrite', default=False, action='store_true',
        help='Force overwrite of existing files')
    parser.add_argument(
        '-d', '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action='store_const', dest='loglevel')
    args = parser.parse_args()

    # Convert input file to an absolute path
    if os.path.isfile(os.path.abspath(args.ini)):
        args.ini = os.path.abspath(args.ini)
    return args


if __name__ == '__main__':
    args = arg_parse()

    logging.basicConfig(level=args.loglevel, format='%(message)s')
    logging.info('\n{}'.format('#' * 80))
    log_f = '{:<20s} {}'
    logging.info(log_f.format(
        'Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', os.getcwd()))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))

    # Calculate GSFLOW HRU Parameters
    hru_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
