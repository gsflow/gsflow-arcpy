#--------------------------------
# Name:         dem_2_streams.py
# Purpose:      GSFLOW Flow Parameters
# Notes:        ArcGIS 10.2+ Version
# Python:       2.7
#--------------------------------

import argparse
from collections import defaultdict
import ConfigParser
import datetime as dt
import logging
import math
import os
import sys
import time

import arcpy
from arcpy import env
import numpy as np

import support_functions as support


def flow_parameters(config_path, overwrite_flag=False, debug_flag=False):
    """Calculate GSFLOW Flow Parameters

    Args:
        config_file (str): Project config file path
        ovewrite_flag (bool): if True, overwrite existing files
        debug_flag (bool): if True, enable debug level logging

    Returns:
        None
    """

    # Initialize hru_parameters class
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
    log_file_name = 'dem_2_stream_log.txt'
    log_console = logging.FileHandler(
        filename=os.path.join(hru.log_ws, log_file_name), mode='w')
    log_console.setLevel(logging.DEBUG)
    log_console.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger('').addHandler(log_console)
    logging.info('\nGSFLOW DEM To Streams')

    # Check whether lake parameters should be calculated
    try:
        set_lake_flag = inputs_cfg.getboolean('INPUTS', 'set_lake_flag')
    except ConfigParser.NoOptionError:
        set_lake_flag = False
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'set_lake_flag', set_lake_flag))

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

    # Flow parameters
    flow_acc_threshold = inputs_cfg.getint('INPUTS', 'flow_acc_threshold')
    flow_length_threshold = inputs_cfg.getint('INPUTS', 'flow_length_threshold')
    try:
        calc_flow_dir_points_flag = inputs_cfg.getboolean(
            'INPUTS', 'calc_flow_dir_points_flag')
    except ConfigParser.NoOptionError:
        calc_flow_dir_points_flag = False
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'calc_flow_dir_points_flag', calc_flow_dir_points_flag))
    try:
        lake_seg_offset = inputs_cfg.getint('INPUTS', 'lake_seg_offset')
    except ConfigParser.NoOptionError:
        lake_seg_offset = 0
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'lake_seg_offset', lake_seg_offset))
    if lake_seg_offset < 0:
        logging.error(
            '\nERROR: lake_seg_offset must be an integer greater than 0')
        sys.exit()

    # Check input paths
    dem_temp_ws = os.path.join(hru.param_ws, 'dem_rasters')
    dem_path = os.path.join(dem_temp_ws, 'dem.img')
    if not arcpy.Exists(dem_path):
        logging.error(
            '\nERROR: Projected/clipped DEM ({}) does not exist'
            '\nERROR: Try rerunning dem_parameters.py'.format(dem_path))
        sys.exit()
    if not arcpy.Exists(hru.polygon_path):
        logging.error(
            '\nERROR: Fishnet ({}) does not exist'.format(hru.polygon_path))
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

    # Build output folder if necessary
    flow_temp_ws = os.path.join(hru.param_ws, 'flow_rasters')
    if not os.path.isdir(flow_temp_ws):
        os.mkdir(flow_temp_ws)

    # Output paths
    hru_type_path = os.path.join(flow_temp_ws, 'hru_type.img')
    dem_adj_path = os.path.join(flow_temp_ws, 'dem_adj.img')
    lake_id_path = os.path.join(flow_temp_ws, 'lake_id.img')
    dem_sink_path = os.path.join(flow_temp_ws, 'dem_sink.img')
    dem_fill_path = os.path.join(flow_temp_ws, 'dem_fill.img')
    flow_dir_path = os.path.join(flow_temp_ws, 'flow_dir.img')
    flow_dir_points = os.path.join(flow_temp_ws, 'flow_dir_points.shp')
    flow_acc_full_path = os.path.join(flow_temp_ws, 'flow_acc_full.img')
    flow_acc_sub_path = os.path.join(flow_temp_ws, 'flow_acc_sub.img')
    flow_mask_path = os.path.join(flow_temp_ws, 'flow_mask.img')
    stream_link_path = os.path.join(flow_temp_ws, 'stream_link.img')
    stream_link_a_path = os.path.join(flow_temp_ws, 'stream_link_a.img')
    stream_link_b_path = os.path.join(flow_temp_ws, 'stream_link_b.img')
    stream_order_path = os.path.join(flow_temp_ws, 'stream_order.img')
    stream_length_path = os.path.join(flow_temp_ws, 'stream_length.img')
    watersheds_path = os.path.join(flow_temp_ws, 'watersheds.img')
    outlet_path = os.path.join(flow_temp_ws, 'outlet.img')
    swale_path = os.path.join(flow_temp_ws, 'swale.img')
    subbasin_path = os.path.join(flow_temp_ws, 'subbasin.img')
    basin_path = os.path.join(flow_temp_ws, 'basin.img')
    streams_path = os.path.join(flow_temp_ws, 'streams.shp')
    model_points_path = os.path.join(flow_temp_ws, 'model_points.shp')

    # Set ArcGIS environment variables
    arcpy.CheckOutExtension('Spatial')
    env.overwriteOutput = True
    # env.pyramid = 'PYRAMIDS -1'
    env.pyramid = 'PYRAMIDS 0'
    env.workspace = flow_temp_ws
    env.scratchWorkspace = hru.scratch_ws

    # Set environment parameters
    env.extent = hru.extent
    env.cellsize = hru.cs
    env.outputCoordinateSystem = hru.sr

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

    # Check model_points_zone_field
    if model_points_zone_field.upper() in ['', 'FID', 'NONE']:
        model_points_fid_field = arcpy.Describe(model_points_path).OIDFieldName
        logging.warning(
            '  NOTE: Using {}+1 to set {}'.format(
                model_points_fid_field, hru.subbasin_field))
        model_points_zone_field = 'ZONE_VALUE'
        if not arcpy.ListFields(model_points_path, model_points_zone_field):
            arcpy.AddField_management(
                model_points_path, model_points_zone_field, 'LONG')
        arcpy.CalculateField_management(
            model_points_path, model_points_zone_field,
            '!{}! + 1'.format(model_points_fid_field), 'PYTHON')
    elif not arcpy.ListFields(model_points_path, model_points_zone_field):
        logging.error(
            '\nERROR: model_points_zone_field {} does not exist\n'.format(
                model_points_zone_field))
        sys.exit()
    # Need to check that model_points_zone_field is an int type
    elif not [f.type for f in arcpy.Describe(model_points_path).fields
              if (f.name == model_points_zone_field and
                  f.type in ['SmallInteger', 'Integer'])]:
        logging.error(
            '\nERROR: model_points_zone_field {} must be an integer type\n'.format(
                model_points_zone_field))
        sys.exit()

    # Need to check that model_points_zone_field is all positive values
    if min([row[0] for row in arcpy.da.SearchCursor(
            model_points_path, [model_points_zone_field])]) <= 0:
        logging.error(
            '\nERROR: model_points_zone_field values must be positive\n'.format(
                model_points_zone_field))
        sys.exit()

    # Check that subbasin values increment from 1 to nsub
    logging.info('  Checking subbasin numbering')
    subbasin_id_list = sorted(list(set(
        [row[0] for row in arcpy.da.SearchCursor(
            model_points_path, [model_points_zone_field])])))
    if subbasin_id_list != range(1, len(subbasin_id_list) + 1):
        logging.error(
            '\nERROR: SUB_BASINs must be sequential starting from 1'
            '\nERROR:   {}'.format(subbasin_id_list))
        sys.exit()
    subbasin_input_count = len(subbasin_id_list)
    logging.debug('    {} subbasins'.format(subbasin_input_count))

    # Check model point types
    logging.info('  Checking model point types')
    model_point_types = [str(r[0]).upper() for r in arcpy.da.SearchCursor(
        model_points_path, [model_points_type_field])]
    if not set(model_point_types).issubset(set(['OUTLET', 'SUBBASIN', 'SWALE'])):
        logging.error('\nERROR: Unsupported model point type(s) found, exiting')
        logging.error('\n  Model point types: {}\n'.format(model_point_types))
        sys.exit()
    elif not set(model_point_types).intersection(set(['OUTLET', 'SWALE'])):
        logging.error(
            '\nERROR: At least one model point must be an OUTLET or SWALE, '
            'exiting\n')
        sys.exit()
    else:
        logging.debug('  {}'.format(', '.join(model_point_types)))

    # Check DEM field
    logging.info('\nAdding DEM fields if necessary')
    support.add_field_func(hru.polygon_path, hru.iseg_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.irunbound_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.flow_dir_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.dem_sink_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.outflow_field, 'DOUBLE')


    if set_lake_flag:
        # Check lake cell elevations
        logging.info('\nChecking lake cell {}'.format(hru.dem_adj_field))
        lake_elev_dict = defaultdict(list)
        fields = [
            hru.type_field, hru.lake_id_field,
            hru.dem_adj_field, hru.id_field]
        for row in arcpy.da.SearchCursor(hru.polygon_path, fields):
            if int(row[0]) != 2:
                continue
            lake_elev_dict[int(row[1])].append(float(row[2]))
        del fields
        logging.info('  {:>7} {:>12} {:>12} {:>12} {:>12}'.format(
            'Lake ID', 'Minimum', 'Mean', 'Maximum', 'Std. Dev.'))
        for lake_id, lake_elev_list in lake_elev_dict.items():
            lake_elev_array = np.array(lake_elev_list)
            logging.info('  {:7} {:12f} {:12f} {:12f} {:12f}'.format(
                lake_id, np.min(lake_elev_array), np.mean(lake_elev_array),
                np.max(lake_elev_array), np.std(lake_elev_array)))
            if np.std(lake_elev_array) > 1:
                logging.warning(
                    '  Please check the lake cell elevations\n'
                    '  They may need to be manually adjusted'.format(lake_id))
                raw_input('  Press ENTER to continue')
            del lake_elev_array

        # Build Lake raster
        logging.debug('  LAKE_ID')
        arcpy.PolygonToRaster_conversion(
            hru.polygon_path, hru.lake_id_field, lake_id_path,
            'CELL_CENTER', '', hru.cs)
        lake_id_obj = arcpy.sa.Raster(lake_id_path)


    logging.info('\nExporting HRU polygon parameters to raster')
    logging.debug('  HRU_TYPE')
    arcpy.PolygonToRaster_conversion(
        hru.polygon_path, hru.type_field, hru_type_path,
        'CELL_CENTER', '', hru.cs)
    hru_type_obj = arcpy.sa.Raster(hru_type_path)

    # Convert DEM_ADJ to raster
    logging.debug('  DEM_ADJ')
    arcpy.PolygonToRaster_conversion(
        hru.polygon_path, hru.dem_adj_field, dem_adj_path,
        'CELL_CENTER', '', hru.cs)
    dem_adj_obj = arcpy.sa.Raster(dem_adj_path)
    # dem_adj_obj = arcpy.sa.Float(arcpy.sa.Raster(dem_adj_path))

    hru_polygon_lyr = 'hru_polygon_lyr'
    arcpy.MakeFeatureLayer_management(hru.polygon_path, hru_polygon_lyr)
    arcpy.SelectLayerByAttribute_management(hru_polygon_lyr, 'CLEAR_SELECTION')
    arcpy.CalculateField_management(
        hru_polygon_lyr, hru.outflow_field, 0, 'PYTHON')

    if 'OUTLET' in model_point_types:
        arcpy.SelectLayerByAttribute_management(
            model_points_lyr, 'NEW_SELECTION', '"TYPE" = \'OUTLET\'')

        arcpy.SelectLayerByLocation_management(
            hru_polygon_lyr, 'INTERSECT', model_points_lyr)
        arcpy.CalculateField_management(
            hru_polygon_lyr, hru.outflow_field, 1, 'PYTHON')

        # The point of all of this code is to determine the flow direction
        #   at the outlet points since it won't be computed.
        # It might be easier to compute fill and flow dir. on the full raster
        logging.info('  Computing OUTLET point flow direction')

        # Get HRU values at outlet points
        outlet_points = [
            (int(r[0]), int(r[1]))
            for r in arcpy.da.SearchCursor(
                hru_polygon_lyr, [hru.col_field, hru.row_field])]

        # Get elevations and type of neighboring cells
        # Multiplying the cellsize by 1.5 is needed to get all possible
        #   neighbors but it can return extra cells that will need to be skipped
        # It might be easier to use the Select tool directly
        arcpy.SelectLayerByLocation_management(
            hru_polygon_lyr, 'WITHIN_A_DISTANCE', model_points_lyr, 1.5 * hru.cs)
        elev_dict = dict()
        hru_type_dict = dict()
        fields = [hru.col_field, hru.row_field, hru.dem_adj_field, hru.type_field]
        for row in arcpy.da.SearchCursor(hru_polygon_lyr, fields):
            elev_dict[(int(row[0]), int(row[1]))] = float(row[2])
            hru_type_dict[(int(row[0]), int(row[1]))] = int(row[3])

        # For each outlet cell, cycle through flow directions and find ?.
        # Outlet cells should exit to an inactive cell or out of the grid.
        outlet_flowdir = {}
        for outlet_pt in outlet_points:
            logging.debug('    Outlet Point: {}'.format(outlet_pt))
            outlet_slopes = []
            # Search non-diagonals first.
            for fd in [1, 4, 16, 64, 2, 8, 32, 128]:
                if support.next_row_col(fd, outlet_pt) not in elev_dict.keys():
                    # Don't compute other slopes if next cell is outside the grid
                    outlet_slopes.append([-9999, fd])
                    break
                elif hru_type_dict[support.next_row_col(fd, outlet_pt)] != 0:
                    # Only compute slope to inactive cells
                    continue
                else:
                    # Compute slope to next cell
                    slope = (elev_dict[support.next_row_col(fd, outlet_pt)] -
                             elev_dict[outlet_pt])
                    if fd in [2, 8, 32, 128]:
                        # For diagonals, adjust slope
                        # I think Arc approximates root(2) to 1.5
                        slope /= 1.5
                    outlet_slopes.append([slope, fd])
                logging.debug('    {:>3d} {}'.format(fd, slope))

            if not outlet_slopes:
                logging.error(
                    '\nERROR: The OUTLET model point is not at the '
                    'edge of the study area or model grid.\n'
                    '  Col: {0} Rol: {1}'.format(*outlet_pt))
                sys.exit()

            # Assign the flow direction with the steepest (positive) slope
            outlet_slope, outlet_fd = min(outlet_slopes)
            outlet_flowdir[outlet_pt] = outlet_fd
            if outlet_slope > 0:
                logging.warning(
                    '\n  WARNING: The OUTLET model point flow direction may '
                    'be invalid')
            logging.debug('    Flow Direction: {}'.format(outlet_fd))

        logging.info('  Building OUTLET point raster')
        outlet_array = np.zeros((hru.rows, hru.cols)).astype(np.uint8)
        for outlet_pt in outlet_points:
            outlet_array[outlet_pt[1] - 1, outlet_pt[0] - 1] = outlet_flowdir[outlet_pt]
        support.array_to_raster(
            outlet_array, outlet_path,
            arcpy.Point(hru.extent.XMin, hru.extent.YMin, 0),
            hru.cs, outlet_array)
        outlet_obj = arcpy.sa.Raster(outlet_path)

    if 'SWALE' in model_point_types:
        logging.info('  Building SWALE point raster')
        arcpy.SelectLayerByAttribute_management(
            model_points_lyr, 'NEW_SELECTION', '"TYPE" = \'SWALE\'')

        # DEADBEEF - Should SWALE points be written to OUTFLOWHRU.TXT?
        arcpy.SelectLayerByLocation_management(
            hru_polygon_lyr, 'INTERSECT', model_points_lyr)
        arcpy.CalculateField_management(
            hru_polygon_lyr, hru.outflow_field, 1, 'PYTHON')

        arcpy.PointToRaster_conversion(
            model_points_lyr, model_points_type_field, swale_path,
            "", "", hru.cs)
        swale_obj = arcpy.sa.Raster(swale_path)
        arcpy.SelectLayerByAttribute_management(
            model_points_lyr, 'CLEAR_SELECTION')

    arcpy.Delete_management(hru_polygon_lyr)



    logging.info('\nCalculating flow direction')
    # This will force all active cells to flow to an outlet
    logging.debug('  Setting DEM_ADJ values to 20000 for inactivate cells')
    dem_mod_obj = arcpy.sa.Con(hru_type_obj > 0, dem_adj_obj, 20000.0)
    if 'OUTLET' in model_point_types:
        logging.debug('  Setting DEM_ADJ values to NoData for OUTLET cells')
        dem_mod_obj = arcpy.sa.Con(arcpy.sa.IsNull(outlet_obj), dem_mod_obj)
    if 'SWALE' in model_point_types:
        logging.debug('  Setting DEM_ADJ values to NoData for SWALE cells')
        dem_mod_obj = arcpy.sa.Con(arcpy.sa.IsNull(swale_obj), dem_mod_obj)

    logging.info('  Filling DEM_ADJ (8-way)')
    dem_fill_obj = arcpy.sa.Fill(dem_mod_obj)
    del dem_mod_obj

    if 'OUTLET' in model_point_types:
        logging.debug('  Resetting OUTLET cell values')
        dem_fill_obj = arcpy.sa.Con(
            arcpy.sa.IsNull(outlet_obj), dem_fill_obj, dem_adj_obj)

    logging.info('  Calculating sinks (8-way)')
    # Threshold of 0.001 is needed to avoid noise from 32/64 bit conversion
    dem_sink_obj = arcpy.sa.Con(hru_type_obj > 0, dem_fill_obj - dem_adj_obj)
    dem_sink_obj = arcpy.sa.Con(dem_sink_obj > 0.001, dem_sink_obj)

    logging.info('  Calculating flow direction')
    flow_dir_obj = arcpy.sa.FlowDirection(dem_fill_obj, False)

    logging.debug('  Setting flow direction to NoData for inactive cells')
    flow_dir_obj = arcpy.sa.SetNull(hru_type_obj == 0, flow_dir_obj)

    if 'OUTLET' in model_point_types:
        logging.debug('  Resetting OUTLET cell flow direction')
        flow_dir_obj = arcpy.sa.Con(
            ~arcpy.sa.IsNull(outlet_obj), outlet_obj, flow_dir_obj)
        del outlet_obj
    if 'SWALE' in model_point_types:
        logging.debug('  Resetting SWALE cell flow direction')
        flow_dir_obj = arcpy.sa.Con(
            ~arcpy.sa.IsNull(swale_obj), 1, flow_dir_obj)
        del swale_obj

    logging.debug('  Resetting DEM_ADJ values for inactive cell')
    dem_fill_obj = arcpy.sa.Con(hru_type_obj == 0, dem_adj_obj, dem_fill_obj)

    flow_dir_obj.save(flow_dir_path)
    dem_fill_obj.save(dem_fill_path)
    dem_sink_obj.save(dem_sink_path)
    del dem_sink_obj


    # Save flow direction as points
    if calc_flow_dir_points_flag:
        logging.info('\nFlow direction points')
        # ArcGIS fails for raster_to_x conversions on a network path
        # You have to go through an in_memory file first
        flow_dir_temp = os.path.join('in_memory', 'flow_dir')
        arcpy.RasterToPoint_conversion(flow_dir_path, flow_dir_temp)
        try:
            arcpy.CopyFeatures_management(flow_dir_temp, flow_dir_points)
        except:
            time.sleep(1)
            logging.warning('Copy feature failed')
        arcpy.Delete_management(flow_dir_temp)
        del flow_dir_temp

        # Reclassify flow directions to angles, assuming 1 is 0
        remap_cb = (
            'def Reclass(value):\n'
            '    if value == 1: return 0\n'
            '    elif value == 2: return 45\n'
            '    elif value == 4: return 90\n'
            '    elif value == 8: return 135\n'
            '    elif value == 16: return 180\n'
            '    elif value == 32: return 225\n'
            '    elif value == 64: return 270\n'
            '    elif value == 128: return 315\n')
        arcpy.CalculateField_management(
            flow_dir_points, 'grid_code',
            'Reclass(!{}!)'.format('grid_code'), 'PYTHON', remap_cb)

    # Write flow direction to hru_polygon
    logging.debug('  Extracting flow direction at points')
    vt_list = [[flow_dir_path, hru.flow_dir_field]]
    mem_point_path = os.path.join('in_memory', 'hru_point')
    arcpy.CopyFeatures_management(hru.point_path, mem_point_path)
    arcpy.sa.ExtractMultiValuesToPoints(mem_point_path, vt_list, 'NONE')
    logging.debug('  Reading flow direction values at point')
    data_dict = defaultdict(dict)
    fields = [hru.flow_dir_field, hru.fid_field]
    with arcpy.da.SearchCursor(mem_point_path, fields) as s_cursor:
        for row in s_cursor:
            # Set nodata cells to 0
            if row[0] is not None and row[1] is not None:
                data_dict[int(row[1])][hru.flow_dir_field] = int(row[0])
            del row
    logging.debug('  Writing flow direction values to polygon')
    fields = [hru.flow_dir_field, hru.fid_field]
    with arcpy.da.UpdateCursor(hru.polygon_path, fields) as u_cursor:
        for row in u_cursor:
            row_dict = data_dict.get(int(row[-1]), None)
            for i, field in enumerate(fields[:-1]):
                if row_dict:
                    row[i] = row_dict[field]
                else:
                    row[i] = 0
            u_cursor.updateRow(row)
            del row_dict, row


    # DEADBEEF - This whole section seems to only be needed if the outflows
    #   are not specified by the user.
    # # Subbasins
    # # Select the HRU cells that intersect the subbasin point cells
    # logging.debug('  Reading input subbasin points')
    # hru_polygon_lyr = 'hru_polygon_lyr'
    # arcpy.MakeFeatureLayer_management(hru.polygon_path, hru_polygon_lyr)
    # arcpy.SelectLayerByLocation_management(
    #     hru_polygon_lyr, 'intersect', model_points_path)
    # input_xy_dict = dict()
    # fields = [hru.col_field, hru.row_field, hru.x_field, hru.y_field]
    # for row in arcpy.da.SearchCursor(hru_polygon_lyr, fields):
    #     input_xy_dict[(int(row[0]), int(row[1]))] = (int(row[2]), int(row[3]))
    # arcpy.Delete_management(hru_polygon_lyr)
    # del hru_polygon_lyr
    # # for k,v in input_xy_dict.items():
    # #    logging.debug('    {} {}'.format(k,v))

    # logging.info('\nBuilding all subbasin points')
    # # First calculate downstream cell for all cells
    # logging.debug('  Calculating downstream cells')
    # out_cell_dict = dict()
    # hru_type_dict = dict()
    # cell_xy_dict = dict()
    # fields = [
    #     hru.type_field, hru.flow_dir_field, hru.id_field,
    #     hru.col_field, hru.row_field, hru.x_field, hru.y_field]
    # for row in arcpy.da.SearchCursor(hru.polygon_path, fields):
    #     cell = (int(row[3]), int(row[4]))
    #     out_cell_dict[cell] = support.next_row_col(int(row[1]), cell)
    #     hru_type_dict[cell] = int(row[0])
    #     cell_xy_dict[cell] = (int(row[5]), int(row[6]))

    # # Identify all active/lake cells that exit the model
    # #   or flow to an inactive cell
    # logging.debug('  Identifying active cells that exit the model')
    # out_cell_xy_list = []
    # for cell, cell_xy in sorted(cell_xy_dict.items()):
    #     #  DEADBEEF - This is finding exit cells that aren't already gauges
    #     # if cell in input_xy_dict.keys():
    #     #    continue
    #     # elif cell not in hru_type_dict.keys():
    #     if cell not in hru_type_dict.keys():
    #         continue
    #     elif hru_type_dict[cell] not in [1, 2]:
    #         continue
    #     elif cell not in out_cell_dict.keys():
    #         continue
    #     elif out_cell_dict[cell] not in hru_type_dict.keys():
    #         out_cell_xy_list.append(cell_xy)
    #     elif (out_cell_dict[cell] in hru_type_dict.keys() and
    #           hru_type_dict[out_cell_dict[cell]] not in [1, 2]):
    #         out_cell_xy_list.append(cell_xy)

    # # Outflow cells exit the model to inactive cells or out of the domain
    # # These cells will be used to set the OUTFLOW_HRU.DAT for CRT
    # #   in crt_fill_parameters and stream_parameters
    # logging.info('  Flag outflow cells')
    # fields = [hru.type_field, hru.x_field, hru.y_field, hru.outflow_field]
    # with arcpy.da.UpdateCursor(hru.polygon_path, fields) as u_cursor:
    #     for row in u_cursor:
    #         cell_xy = (row[1], row[2])
    #         # Inactive cells can't be outflow cells
    #         if int(row[0]) == 0:
    #             continue
    #         elif out_cell_xy_list and cell_xy in out_cell_xy_list:
    #             row[3] = 1
    #         else:
    #             row[3] = 0
    #         u_cursor.updateRow(row)
    # del out_cell_dict, hru_type_dict, cell_xy_dict


    # DEADBEEF - This was added for sinks or ocean so that there would be
    #   subbasin points along the edge?
    # fields = ['SHAPE@XY', model_points_zone_field]
    # with arcpy.da.InsertCursor(model_points_path, fields) as insert_c:
    #     for out_cell_xy in sorted(out_cell_xy_list):
    #         insert_c.insertRow([out_cell_xy, subbasin_input_count + 1])
    # del fields
    # del out_cell_xy_list


    # Flow Accumulation
    logging.info('\nCalculating initial flow accumulation')
    flow_acc_full_obj = arcpy.sa.FlowAccumulation(flow_dir_obj)
    logging.info('  Only keeping flow_acc >= {}'.format(flow_acc_threshold))
    flow_acc_full_obj = arcpy.sa.Con(
        flow_acc_full_obj >= flow_acc_threshold, flow_acc_full_obj)
    flow_acc_full_obj.save(flow_acc_full_path)

    # Flow accumulation and stream link with lakes
    logging.info('\nCalculating flow accumulation & stream link (w/ lakes)')
    flow_acc_mask_obj = arcpy.sa.Con(
        (hru_type_obj >= 1) & (hru_type_obj <= 3) & (flow_acc_full_obj > 0), 1)
    stream_link_obj = arcpy.sa.StreamLink(flow_acc_mask_obj, flow_dir_obj)
    stream_link_obj.save(stream_link_a_path)
    del flow_acc_mask_obj, stream_link_obj

    # Flow accumulation and stream link without lakes
    logging.info('Calculating flow accumulation & stream link (w/o lakes)')
    flow_acc_mask_obj = arcpy.sa.Con(
        (hru_type_obj == 1) | (hru_type_obj == 3) & (flow_acc_full_obj > 0), 1)
    # flow_acc_obj.save(flow_acc_sub_path)
    stream_link_obj = arcpy.sa.StreamLink(flow_acc_mask_obj, flow_dir_obj)
    stream_link_obj.save(stream_link_b_path)
    del flow_acc_mask_obj, stream_link_obj

    # Initial Stream Link
    # logging.info('\nCalculating initial stream link')
    # stream_link_obj = StreamLink(flow_acc_obj, flow_dir_obj)
    # stream_link_obj.save(stream_link_path)
    # Calculate stream link with and without lakes
    # Initial Stream Order (w/ lakes)
    logging.info('Calculating stream order (w/ lakes)')
    logging.debug(
        '  Using SHREVE ordering so after 1st order are removed, '
        '2nd order will only be dangles')
    stream_order_obj = arcpy.sa.StreamOrder(
        stream_link_a_path, flow_dir_obj, 'SHREVE')
    stream_order_obj.save(stream_order_path)

    # Stream Length (cell count w/o lakes)
    logging.info('Calculating stream length (cell count w/o lakes)')
    stream_length_obj = arcpy.sa.Lookup(stream_link_b_path, 'Count')
    stream_length_obj.save(stream_length_path)

    # Filter 1st order segments
    logging.info(
         '\nFilter all 1st order streams with length < {}'
         '\nKeep all higher order streams'.format(flow_length_threshold))
    # Stream length is nodata for lakes, so put lakes back in
    # This is needed to remove short 1st order streams off of lakes
    flow_acc_mask_obj = (
        (hru_type_obj == 3) | (hru_type_obj == 2) | (stream_order_obj >= 2) |
        ((stream_order_obj == 1) &
         (stream_length_obj >= flow_length_threshold)))
    flow_acc_mask_obj.save(flow_mask_path)
    flow_acc_sub_obj = arcpy.sa.Con(flow_acc_mask_obj, flow_acc_full_obj)
    flow_acc_sub_obj.save(flow_acc_sub_path)
    del flow_acc_mask_obj, stream_order_obj, stream_length_obj

    # Final Stream Link
    logging.info('\nCalculating final stream link')
    flow_acc_mask_obj = arcpy.sa.Con((flow_acc_sub_obj >= 1), 1)
    stream_link_obj = arcpy.sa.StreamLink(flow_acc_mask_obj, flow_dir_obj)
    # Get count of streams for automatically setting lake_seg_offset
    if not lake_seg_offset:
        lake_seg_count = int(stream_link_obj.maximum)
        # NOTE: This call fails in 10.6.1
        # lake_seg_count = int(
        #     arcpy.GetCount_management(stream_link_obj).getOutput(0))
        n = 10 ** math.floor(math.log10(lake_seg_count))
        lake_seg_offset = int(math.ceil((lake_seg_count + 1) / n)) * int(n)
        logging.info(
             '  lake_segment_offset was not set in the input file\n'
             '  Using automatic lake segment offset: {}'.format(
                 lake_seg_offset))
    elif set_lake_flag:
        logging.info(
             '  Using manual lake segment offset: {}'.format(lake_seg_offset))

    # Include lake cells into 'stream_link' before calculating watersheds
    # Watershed function doesn't work for negative values
    # Convert lakes to large positive numbers for Watershed
    # ISEG needs to be negative values though
    if set_lake_flag:
        logging.info(
             '  Including lakes as {0} + {1}\n'
             '  This will allow for a watershed/subbasin for the lakes\n'
             '  {2} will be save as negative of {0} though'.format(
                 hru.lake_id_field, lake_seg_offset, hru.iseg_field))
        stream_link_obj = arcpy.sa.Con(
            (hru_type_obj == 2),
            (lake_id_obj + lake_seg_offset), stream_link_obj)
    stream_link_obj.save(stream_link_path)
    del flow_acc_mask_obj

    # Watersheds
    logging.info('Calculating watersheds')
    watersheds_obj = arcpy.sa.Watershed(flow_dir_obj, stream_link_obj)
    watersheds_obj.save(watersheds_path)
    del stream_link_obj, watersheds_obj

    # Subbasins
    logging.info('Calculating subbasins')
    subbasin_obj = arcpy.sa.Watershed(
        flow_dir_obj, model_points_path, model_points_zone_field)
    subbasin_obj.save(subbasin_path)
    del subbasin_obj

    # Basins
    logging.info('Calculating basins')
    basin_obj = arcpy.sa.Basin(flow_dir_obj)
    basin_obj.save(basin_path)
    del basin_obj

    # Clear subbasin value if HRU_TYPE is 0
    logging.info('Clearing subbasin ID for inactive cells')
    subbasin_obj = arcpy.sa.SetNull(
        hru_type_obj, arcpy.sa.Raster(subbasin_path), 'VALUE=0')
    subbasin_obj.save(subbasin_path)
    del subbasin_obj
    del hru_type_obj


    # Stream polylines
    logging.info('Calculating stream polylines')
    # ArcGIS fails for raster_to_x conversions on a network path
    # You have to go through an in_memory file first
    streams_temp = os.path.join('in_memory', 'streams')
    arcpy.sa.StreamToFeature(
        stream_link_path, flow_dir_obj, streams_temp, 'NO_SIMPLIFY')
    arcpy.CopyFeatures_management(streams_temp, streams_path)
    arcpy.Delete_management(streams_temp)
    del streams_temp


    # Write values to hru_polygon
    logging.info('\nExtracting stream parameters')
    vt_list = [
        [watersheds_path, hru.irunbound_field],
        [stream_link_path, hru.iseg_field],
        # [flow_dir_path, hru.flow_dir_field],
        [subbasin_path, hru.subbasin_field],
        [hru_type_path, hru.type_field]]
    mem_point_path = os.path.join('in_memory', 'hru_point')
    arcpy.CopyFeatures_management(hru.point_path, mem_point_path)
    arcpy.sa.ExtractMultiValuesToPoints(mem_point_path, vt_list, 'NONE')
    del vt_list

    # Read values from points
    logging.info('  Reading cell values')
    data_dict = defaultdict(dict)
    fields = [
        hru.irunbound_field, hru.iseg_field,
        hru.subbasin_field, hru.type_field, hru.fid_field]
    # fields = [
    #    hru.irunbound_field, hru.iseg_field, hru.flow_dir_field,
    #    hru.subbasin_field, hru.type_field, hru.fid_field]
    with arcpy.da.SearchCursor(mem_point_path, fields) as s_cursor:
        for row in s_cursor:
            for i, field in enumerate(fields[:-1]):
                # Set nodata or inactive cells to 0
                if row[i] is None or (int(row[-2]) == 0):
                    data_dict[int(row[-1])][field] = 0
                else:
                    data_dict[int(row[-1])][field] = int(row[i])
            del row
    del fields

    # ISEG for lake cells must be -1 * LAKE_ID, not LAKE_ID + OFFSET
    for k in data_dict.keys():
        irunbound = data_dict[k][hru.irunbound_field]
        iseg = data_dict[k][hru.iseg_field]
        if irunbound > lake_seg_offset:
            data_dict[k][hru.irunbound_field] = lake_seg_offset - irunbound
        if iseg > lake_seg_offset:
            data_dict[k][hru.iseg_field] = lake_seg_offset - iseg

    # data_dict = dict([(k,v) for k,v in data_dict.items()])
    # Write values to polygon
    logging.info('  Writing values to polygons')
    fields = [
        hru.irunbound_field, hru.iseg_field,
        hru.subbasin_field, hru.type_field, hru.fid_field]
    # fields = [
    #    hru.irunbound_field, hru.iseg_field, hru.flow_dir_field,
    #    hru.subbasin_field, hru.type_field, hru.fid_field]
    with arcpy.da.UpdateCursor(hru.polygon_path, fields) as u_cursor:
        for row in u_cursor:
            row_dict = data_dict.get(int(row[-1]), None)
            for i, field in enumerate(fields[:-1]):
                if row_dict:
                    row[i] = row_dict[field]
                else:
                    row[i] = 0
            u_cursor.updateRow(row)
            del row_dict, row
    del fields


    # Write sink values to hru_polygon
    vt_list = []
    if arcpy.Exists(dem_sink_path):
        vt_list.append([dem_sink_path, hru.dem_sink_field])
    if vt_list:
        logging.info('\nExtracting sink values')
        for vt_item in vt_list:
            logging.debug('  {}: {}'.format(
                vt_item[1], os.path.basename(vt_item[0])))
        mem_point_path = os.path.join('in_memory', 'hru_point')
        arcpy.CopyFeatures_management(hru.point_path, mem_point_path)
        arcpy.sa.ExtractMultiValuesToPoints(mem_point_path, vt_list, 'NONE')

        # Read sink values from points
        logging.info('  Reading sink values')
        data_dict = defaultdict(dict)
        fields = [field for path, field in vt_list] + [hru.fid_field]
        with arcpy.da.SearchCursor(mem_point_path, fields) as s_cursor:
            for row in s_cursor:
                for i, field in enumerate(fields[:-1]):
                    # Set nodata or inactive cells to 0
                    if row[i] is None:
                        data_dict[int(row[-1])][field] = 0
                    else:
                        data_dict[int(row[-1])][field] = float(row[i])
                del row

        # Write sink values to polygon
        logging.info('  Writing sink values to polygons')
        fields = [field for path, field in vt_list] + [hru.fid_field]
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as u_cursor:
            for row in u_cursor:
                row_dict = data_dict.get(int(row[-1]), None)
                for i, field in enumerate(fields[:-1]):
                    if row_dict:
                        row[i] = row_dict[field]
                    else:
                        row[i] = 0
                u_cursor.updateRow(row)
                del row_dict, row

    # Cleanup
    arcpy.Delete_management(mem_point_path)
    del mem_point_path, vt_list, data_dict, field

    # Re-Calculate HRU_ELEV
    # logging.info('Calculating HRU_ELEV from DEM_ADJ')
    # logging.info('  Converting from meters to feet')
    # arcpy.CalculateField_management(
    #    hru.polygon_path, hru_elev_field,
    #    # Convert meters to feet
    #    '!{}! * 3.28084'.format(dem_adj_field), 'PYTHON')

    # Cleanup
    del dem_fill_obj
    if set_lake_flag:
        del lake_id_obj
    del flow_dir_obj
    del flow_acc_full_obj
    del flow_acc_sub_obj


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='DEM To Streams',
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

    # Calculate GSFLOW Flow Parameters
    flow_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
