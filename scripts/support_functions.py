#--------------------------------
# Name:         support_functions.py
# Purpose:      GSFLOW parameter support functions
# Notes:        ArcGIS 10.2+ Version
# Python:       2.7
#--------------------------------

from collections import defaultdict
import ConfigParser
import itertools
import logging
import math
from operator import itemgetter
import os
import re
import sys
from time import sleep

import numpy as np

import arcpy
from arcpy import env


class HRUParameters():
    """"""

    def __init__(self, config_path):
        # Open input parameter config file
        inputs_cfg = ConfigParser.ConfigParser()
        try:
            inputs_cfg.readfp(open(config_path))
        except IOError:
            logging.error(
                '\nERROR: Config file does not exist\n'
                '  {}\n'.format(config_path))
            sys.exit()
        except ConfigParser.MissingSectionHeaderError:
            logging.error(
                '\nERROR: Config file is missing a section header\n'
                '    Please make sure the following line is at the '
                'beginning of the file\n[INPUTS]\n')
            sys.exit()
        except Exception as e:
            logging.error(
                '\nERROR: Config file could not be read\n'
                '  {}\n'.format(config_path))

        # Open field list config file
        # Use script directory (from sys.argv[0]) in case script is a
        #   relative path (i.e. called from a project folder)
        field_list_path = os.path.join(
            os.path.dirname(sys.argv[0]), 'field_list.ini')
        # #field_list_path =  inputs_cfg.get('INPUTS', 'field_list_path')

        logging.debug('\nReading Field List File')
        fields_cfg = ConfigParser.ConfigParser()
        try:
            fields_cfg.readfp(open(field_list_path))
        except IOError:
            logging.error(
                '\nERROR: Field list file does not exist\n'
                '  {}\n'.format(field_list_path))
            sys.exit()
        except ConfigParser.MissingSectionHeaderError:
            logging.error(
                '\nERROR: Field list file is missing a section header\n'
                '    Please make sure the following line is at the '
                'beginning of the file\n[FIELDS]\n')
            sys.exit()
        except Exception as e:
            logging.error(
                '\nERROR: Field list file could not be read\n'
                '  {}\n'.format(field_list_path))

        # Read parameters from config file
        logging.debug('\nReading Input File')
        logging.debug('  {}'.format(os.path.basename(config_path)))
        self.polygon_path = inputs_cfg.get('INPUTS', 'hru_fishnet_path')
        self.point_path = inputs_cfg.get('INPUTS', 'hru_centroid_path')
        self.cs = inputs_cfg.getfloat('INPUTS', 'hru_cellsize')
        self.fid_field = inputs_cfg.get('INPUTS', 'orig_fid_field')
        self.type_field = fields_cfg.get('FIELDS', 'type_field')

        # Check inputs
        if self.cs <= 0:
            logging.error('\nERROR: Fishnet cellsize must be greater than 0')
            sys.exit()

        #
        self.param_ws = inputs_cfg.get('INPUTS', 'parameter_folder')
        if not os.path.isdir(self.param_ws):
            os.mkdir(self.param_ws)

        # Log workspace
        self.log_ws = os.path.join(self.param_ws, 'logs')
        if not os.path.isdir(self.log_ws):
            os.mkdir(self.log_ws)

        # Scratch workspace
        try:
            scratch_name = inputs_cfg.get('INPUTS', 'scratch_name')
        except ConfigParser.NoOptionError:
            scratch_name = 'in_memory'
            logging.info(
                '  Missing INI parameter, setting {} = {}'.format(
                    'scratch_name', scratch_name))
        if scratch_name == 'in_memory':
            self.scratch_ws = scratch_name
        else:
            scratch_ws = os.path.join(self.param_ws, scratch_name)
            if not os.path.isdir(scratch_ws):
                os.mkdir(scratch_ws)
            self.scratch_ws = scratch_ws

        # Set spatial reference of hru shapefile
        if arcpy.Exists(self.polygon_path):
            hru_desc = arcpy.Describe(self.polygon_path)
            self.sr = hru_desc.spatialReference
            logging.debug('  Fishnet cellsize:   {}'.format(self.cs))
            logging.debug('  Fishnet spat. ref.: {}'.format(self.sr.name))
            logging.debug('  Fishnet GCS:        {}'.format(self.sr.GCS.name))
            logging.debug('  Fishnet extent:     {}'.format(
                extent_string(hru_desc.extent)))
            self.extent = round_extent(hru_desc.extent, 6)
            logging.debug('  Fishnet extent:     {}'.format(
                extent_string(self.extent)))

            # DEADBEEF - I'm not sure why I would adjust the extent
            # If the extent doesn't match the reference point, the script
            #   should probably terminate
            # self.extent = adjust_extent_to_snap(
            #    hru_desc.extent, self.ref_pnt, self.cs, method='ROUND')
            # self.extent = adjust_extent_to_snap(
            #    hru_param_desc.extent, snap_pnt, self.cs, method='ROUND')

            # Check that the fishnet extent is close to an even multiple of the cellsize
            self.cols = abs(self.extent.XMax - self.extent.XMin) / self.cs
            self.rows = abs(self.extent.YMax - self.extent.YMin) / self.cs
            if ((self.cols - round(self.cols)) > 0.001 or
                    (self.cols - round(self.cols)) > 0.001):
                logging.error(
                    '\nWARNING: {} does not appear to have the {} cellsize '
                    'specified in the INI file\n  This may be a rounding '
                    'issue.'.format(os.path.basename(self.polygon_path, self.cs)))
                logging.debug('  Cols: {}\m  Rows: {}'.format(self.cols, self.rows))
                raw_input('Press ENTER to continue')
            # Round to the nearest integer
            self.cols = int(round(self.cols))
            self.rows = int(round(self.rows))
            logging.debug('  Fishnet shape: {} {}'.format(self.cols, self.rows))

            self.ref_x = self.extent.XMin
            self.ref_y = self.extent.YMin
            logging.debug('  Fishnet snap: {} {}'.format(self.ref_x, self.ref_y))

        else:
            logging.debug('  Fishnet cellsize:   {}'.format(self.cs))


        # Some fields are dependent on the control flags
        try:
            set_lake_flag = inputs_cfg.getboolean('INPUTS', 'set_lake_flag')
        except ConfigParser.NoOptionError:
            set_lake_flag = False
        try:
            calc_flow_acc_dem_flag = inputs_cfg.getboolean(
                'INPUTS', 'calc_flow_acc_dem_flag')
        except ConfigParser.NoOptionError:
            calc_flow_acc_dem_flag = False
        try:
            calc_topo_index_flag = inputs_cfg.getboolean(
                'INPUTS', 'calc_topo_index_flag')
        except ConfigParser.NoOptionError:
            calc_topo_index_flag = False
        try:
            clip_root_depth_flag = inputs_cfg.getboolean(
                'INPUTS', 'clip_root_depth_flag')
        except ConfigParser.NoOptionError:
            clip_root_depth_flag = False
        # try:
        #     set_ppt_zones_flag = inputs_cfg.getboolean(
        #         'INPUTS', 'set_ppt_zones_flag')
        # except ConfigParser.NoOptionError:
        #     set_lake_flag = False

        # Read in all field names
        self.id_field = fields_cfg.get('FIELDS', 'id_field')
        self.type_field = fields_cfg.get('FIELDS', 'type_field')
        self.dem_mean_field = fields_cfg.get('FIELDS', 'dem_mean_field')
        self.dem_max_field = fields_cfg.get('FIELDS', 'dem_max_field')
        self.dem_min_field = fields_cfg.get('FIELDS', 'dem_min_field')
        self.dem_adj_field = fields_cfg.get('FIELDS', 'dem_adj_field')
        if calc_flow_acc_dem_flag:
            self.dem_sum_field = fields_cfg.get('FIELDS', 'dem_sum_field')
            self.dem_count_field = fields_cfg.get('FIELDS', 'dem_count_field')
            self.dem_flowacc_field = fields_cfg.get('FIELDS', 'dem_flowacc_field')
        else:
            self.dem_sum_field = 'DEM_SUM'
            self.dem_count_field = 'DEM_COUNT'
            self.dem_flowacc_field = 'DEM_FLOW_AC'
        self.dem_sink_field = fields_cfg.get('FIELDS', 'dem_sink_field')
        self.crt_elev_field = fields_cfg.get('FIELDS', 'crt_elev_field')
        self.crt_fill_field = fields_cfg.get('FIELDS', 'crt_fill_field')
        self.dem_aspect_field = fields_cfg.get('FIELDS', 'dem_aspect_field')
        self.dem_slope_deg_field = fields_cfg.get('FIELDS', 'dem_slope_deg_field')
        self.dem_slope_rad_field = fields_cfg.get('FIELDS', 'dem_slope_rad_field')
        self.dem_slope_pct_field = fields_cfg.get('FIELDS', 'dem_slope_pct_field')
        self.area_field = fields_cfg.get('FIELDS', 'area_field')
        self.topo_index_field = fields_cfg.get('FIELDS', 'topo_index_field')
        self.row_field = fields_cfg.get('FIELDS', 'row_field')
        self.col_field = fields_cfg.get('FIELDS', 'col_field')
        self.x_field = fields_cfg.get('FIELDS', 'x_field')
        self.y_field = fields_cfg.get('FIELDS', 'y_field')
        self.lat_field = fields_cfg.get('FIELDS', 'lat_field')
        self.lon_field = fields_cfg.get('FIELDS', 'lon_field')
        self.hru_sink_field = fields_cfg.get('FIELDS', 'hru_sink_field')

        if set_lake_flag:
            self.lake_id_field = fields_cfg.get('FIELDS', 'lake_id_field')
            self.lake_area_field = fields_cfg.get('FIELDS', 'lake_area_field')
        else:
            self.lake_id_field = 'LAKE_ID'
            self.lake_area_field = 'LAKE_AREA'

        # DEM based
        self.jh_tmax_field = fields_cfg.get('FIELDS', 'jh_tmax_field')
        self.jh_tmin_field = fields_cfg.get('FIELDS', 'jh_tmin_field')
        self.jh_coef_field = fields_cfg.get('FIELDS', 'jh_coef_field')
        self.snarea_thresh_field = fields_cfg.get('FIELDS', 'snarea_thresh_field')
        self.tmax_adj_field = fields_cfg.get('FIELDS', 'tmax_adj_field')
        self.tmin_adj_field = fields_cfg.get('FIELDS', 'tmin_adj_field')

        # Vegetation
        self.cov_type_field = fields_cfg.get('FIELDS', 'cov_type_field')
        self.covden_sum_field = fields_cfg.get('FIELDS', 'covden_sum_field')
        self.covden_win_field = fields_cfg.get('FIELDS', 'covden_win_field')
        self.snow_intcp_field = fields_cfg.get('FIELDS', 'snow_intcp_field')
        self.wrain_intcp_field = fields_cfg.get('FIELDS', 'wrain_intcp_field')
        self.srain_intcp_field = fields_cfg.get('FIELDS', 'srain_intcp_field')
        self.rad_trncf_field = fields_cfg.get('FIELDS', 'rad_trncf_field')

        # Soil
        self.awc_field = fields_cfg.get('FIELDS', 'awc_field')
        self.clay_pct_field = fields_cfg.get('FIELDS', 'clay_pct_field')
        self.sand_pct_field = fields_cfg.get('FIELDS', 'sand_pct_field')
        self.ksat_field = fields_cfg.get('FIELDS', 'ksat_field')
        self.soil_type_field = fields_cfg.get('FIELDS', 'soil_type_field')
        self.soil_root_max_field = fields_cfg.get('FIELDS', 'soil_root_max_field')
        self.moist_init_field = fields_cfg.get('FIELDS', 'moist_init_field')
        self.moist_max_field = fields_cfg.get('FIELDS', 'moist_max_field')
        self.rechr_init_field = fields_cfg.get('FIELDS', 'rechr_init_field')
        self.rechr_max_field = fields_cfg.get('FIELDS', 'rechr_max_field')
        self.ssr2gw_rate_field = fields_cfg.get('FIELDS', 'ssr2gw_rate_field')
        self.ssr2gw_k_field = fields_cfg.get('FIELDS', 'ssr2gw_k_field')
        self.slowcoef_lin_field = fields_cfg.get('FIELDS', 'slowcoef_lin_field')
        self.slowcoef_sq_field = fields_cfg.get('FIELDS', 'slowcoef_sq_field')

        # Impervious Parameter Fields
        self.imperv_pct_field = fields_cfg.get('FIELDS', 'imperv_pct_field')
        self.carea_max_field = fields_cfg.get('FIELDS', 'carea_max_field')

        # Streams
        self.irunbound_field = fields_cfg.get('FIELDS', 'irunbound_field')
        self.iseg_field = fields_cfg.get('FIELDS', 'iseg_field')
        self.flow_dir_field = fields_cfg.get('FIELDS', 'flow_dir_field')
        self.krch_field = fields_cfg.get('FIELDS', 'krch_field')
        self.irch_field = fields_cfg.get('FIELDS', 'irch_field')
        self.jrch_field = fields_cfg.get('FIELDS', 'jrch_field')
        self.reach_field = fields_cfg.get('FIELDS', 'reach_field')
        self.rchlen_field = fields_cfg.get('FIELDS', 'rchlen_field')
        self.maxreach_field = fields_cfg.get('FIELDS', 'maxreach_field')
        self.outseg_field = fields_cfg.get('FIELDS', 'outseg_field')
        self.iupseg_field = fields_cfg.get('FIELDS', 'iupseg_field')
        self.strm_top_field = fields_cfg.get('FIELDS', 'strm_top_field')
        self.strm_slope_field = fields_cfg.get('FIELDS', 'strm_slope_field')
        self.subbasin_field = fields_cfg.get('FIELDS', 'subbasin_field')
        self.segbasin_field = fields_cfg.get('FIELDS', 'segbasin_field')
        self.outflow_field = fields_cfg.get('FIELDS', 'outflow_field')

        # if set_ppt_zones_flag:
        self.ppt_zone_id_field = fields_cfg.get('FIELDS', 'ppt_zone_id_field')
        self.hru_psta_field = fields_cfg.get('FIELDS', 'hru_psta_field')


def next_row_col(flow_dir, cell):
    """"""
    i_next, j_next = cell
    # Upper left cell is 0,0
    if flow_dir in [1, 2, 128]:
        i_next += 1
    elif flow_dir in [8, 16, 32]:
        i_next -= 1
    if flow_dir in [2, 4, 8]:
        j_next += 1
    elif flow_dir in [32, 64, 128]:
        j_next -= 1
    return i_next, j_next


def field_stat_func(input_path, value_field, stat='MAXIMUM'):
    """"""
    value_list = []
    with arcpy.da.SearchCursor(input_path, value_field) as s_cursor:
        for row in s_cursor:
            value_list.append(row[0])
    if stat.upper() in ['MAXIMUM', 'MAX']:
        return max(value_list)
    elif stat.upper() in ['MINIMUM', 'MIN']:
        return min(value_list)
    else:
        return float(sum(value_list)) / len(value_list)


def add_field_func(hru_param_path, field_name, field_type='DOUBLE'):
    """"""
    if len(field_name) >= 11:
        logging.error(
            'Field: {}\n  Shapefile field names cannot be longer than 10 '
            'characters, exiting'.format(field_name))
        sys.exit()

    # Try adding the field a few times
    # This was needed when add fields to a shapefile on a network drive
    i = 0
    while not arcpy.ListFields(hru_param_path, field_name):
        logging.info('  Field: {}'.format(field_name))
        try:
            arcpy.AddField_management(
                hru_param_path, field_name, field_type)
        except Exception as e:
            logging.debug('    Exception: {}'.format(e))

        if i > 20:
            logging.error('  The field could not be added')
            sys.exit()
        else:
            sleep(1.0)
            i += 1


def transform_func(spat_ref_a, spat_ref_b):
    """"""
    # Set preferred transforms
    if ((spat_ref_a.GCS.name == 'GCS_WGS_1984' and
         spat_ref_b.GCS.name == 'GCS_North_American_1983') or
        (spat_ref_a.GCS.name == 'GCS_North_American_1983' and
         spat_ref_b.GCS.name == 'GCS_WGS_1984')):
        return 'NAD_1983_To_WGS_1984_5'
    else:
        return None
        # return '#'


def valid_raster_func(raster_path, raster_name, hru_param, cs=10):
    """This will check for matching spat. ref., snap point, and cellsize

    Assume raster is not valid unless it passes all checks
    """
    if not arcpy.Exists(raster_path):
        return False
    logging.debug('\nReading existing {} raster'.format(raster_name))
    raster_obj = arcpy.sa.Raster(raster_path)
    raster_sr = raster_obj.spatialReference
    raster_extent = raster_obj.extent
    raster_cs = raster_obj.meanCellWidth
    logging.debug('  {} spat. ref.: {}'.format(
        raster_name, raster_sr.name))
    logging.debug('  {} GCS:        {}'.format(
        raster_name, raster_sr.GCS.name))
    logging.debug('  {} extent:     {}'.format(
        raster_name, extent_string(raster_extent)))
    logging.debug('  {} cellsize:   {}'.format(raster_name, raster_cs))
    ref_pnt = arcpy.Point(hru_param.ref_x, hru_param.ref_y)
    if raster_sr.name != hru_param.sr.name:
        logging.error(
            '\nERROR: The {} spatial reference does not match '
            'the hru_param spatial reference'.format(raster_name))
        return False
    elif not snapped(raster_extent, ref_pnt, hru_param.cs):
        logging.error(
            '\nWARNING: The {} is not snapped to the hru_param'.format(
                raster_name))
        return False
    elif raster_cs != cs:
        logging.error(
            '\nERROR: The {} needs to have a {}m cellsize'.format(
                raster_name, cs))
        return False
    elif not raster_extent.contains(hru_param.extent):
        logging.error(
            '\nERROR: The {} extent is too small'.format(raster_name))
        return False
    else:
        return True


def zonal_stats_func(zs_dict, polygon_path, point_path, hru_param,
                     nodata_value=-999, default_value=0):
    """"""
    for zs_field, (raster_path, zs_stat) in sorted(zs_dict.items()):
        logging.info('  {}: {}'.format(zs_field, zs_stat))
        logging.info('    {}'.format(raster_path))
        # Check inputs
        zs_stat_list = ['MEAN', 'MINIMUM', 'MAXIMUM', 'MAJORITY', 'SUM']
        zs_field_list = arcpy.ListFields(polygon_path, zs_field)
        if zs_stat not in zs_stat_list:
            sys.exit()
        elif len(zs_field_list) == 0:
            logging.error('\nERROR: Zonal stats field {} doesn\'t exist'.format(
                zs_field))
            sys.exit()

    # Check that the shapefiles have a spatial reference
    if arcpy.Describe(polygon_path).spatialReference.name == 'Unknown':
        logging.error(
            '\nERROR: HRU centroids  is not projected '
            '(i.e. does not have a prj file)')
        sys.exit()
    if arcpy.Describe(point_path).spatialReference.name == 'Unknown':
        logging.error(
            '\nERROR: HRU centroids does not appear to be projected '
            '(or does not have a prj file)'
            '\nERROR: Try deleting the centroids (i.e. "_label.shp") and '
            'rerunning hru_parameters.py\n')
        sys.exit()

    # Check that ORIG_FID is in point_path (HRU centroids)
    if len(arcpy.ListFields(point_path, hru_param.fid_field)) == 0:
        logging.error(
             '\nERROR: HRU centroids does not have the field: {}'
             '\nERROR: Try deleting the centroids (i.e. "_label.shp") and '
             'rerunning hru_parameters.py\n'.format(hru_param.fid_field))
        sys.exit()

    # Check for duplicate ORIG_FID values
    hru_param_count = int(arcpy.GetCount_management(point_path).getOutput(0))
    if field_duplicate_check(point_path, hru_param.fid_field, hru_param_count):
        logging.error('\nERROR: There are duplicate {} values\n'.format(
            hru_param.fid_field))
        sys.exit()
    # DEADBEEF - remove once field_duplicate_check() is full developed
    # fid_list = [r[0] for r in arcpy.da.SearchCursor(point_path, [hru_param.fid_field])]
    # if len(fid_list) != len(set(fid_list)):
    #    logging.error(
    #        ('\nERROR: There are duplicate {} values\n').format(hru_param.fid_field))
    #    sys.exit()

    # Create memory objects
    point_subset_path = os.path.join('in_memory', 'point_subset')
    hru_raster_path = os.path.join('in_memory', 'hru_raster')
    # point_subset_path = os.path.join(env.scratchWorkspace, 'point_subset.shp')
    # hru_raster_path = os.path.join(env.scratchWorkspace, 'hru_raster.img')
    # Set environment parameters for polygon to raster conversion
    env.extent = hru_param.extent
    env.outputCoordinateSystem = polygon_path
    # env.cellSize = hru_param.cs

    # Only ~65536 objects can be processed by zonal stats
    block_size = 65000
    for i, x in enumerate(xrange(0, hru_param_count, block_size)):
        logging.info('  FIDS: {}-{}'.format(x, x + block_size))
        # Select a subset of the cell centroids
        logging.debug('    Selecting FID subset')
        subset_str = '"{0}" >= {1} AND "{0}" < {2}'.format(
            hru_param.fid_field, x, x + block_size)
        arcpy.Select_analysis(
            point_path, point_subset_path, subset_str)
        # Convert points subset to raster
        logging.debug('    Converting shapefile to raster')
        arcpy.FeatureToRaster_conversion(
            point_subset_path, hru_param.fid_field,
            hru_raster_path, hru_param.cs)

        # Zonal stats
        logging.debug('    Calculating zonal stats')
        data_dict = defaultdict(dict)
        for zs_field, (raster_path, zs_stat) in sorted(zs_dict.items()):
            zs_name = '{}_{}'.format(zs_field.upper(), i)
            logging.info('    {}: {}'.format(zs_stat.upper(), zs_name))
            # For some reason with 10.2, ZS doesn't work with cs at HRU cs
            env.cellSize = arcpy.sa.Raster(raster_path).meanCellWidth
            # Calculate zonal statistics
            zs_table = os.path.join('in_memory', zs_name)
            # zs_table = os.path.join(env.scratchWorkspace, zs_name+'.dbf')
            zs_obj = arcpy.sa.ZonalStatisticsAsTable(
                hru_raster_path, 'Value', raster_path,
                zs_table, 'DATA', zs_stat.upper())

            # Read values from points
            logging.debug('    Reading values from zs table')
            # Fields 1 & 4 are the 'Value' (ORIG_FID) and the stat (SUM, MEAN, etc)
            fields = [
                f.name for f_i, f in enumerate(arcpy.ListFields(zs_table))
                if f_i in [1, 4]]
            logging.debug('    Fields: {}'.format(', '.join(fields)))
            for row in arcpy.da.SearchCursor(zs_table, fields):
                # Set NoData value for cells that are entirely NoData
                if row[1] is None:
                    data_dict[int(row[0])][zs_field] = nodata_value
                else:
                    data_dict[int(row[0])][zs_field] = float(row[1])
            try:
                arcpy.Delete_management(zs_obj)
            except Exception as e:
                pass
            try:
                arcpy.Delete_management(zs_table)
            except Exception as e:
                pass
            del zs_table, zs_obj, fields

        # Write values to polygon
        logging.info('    Writing values to polygons')
        zs_fields = sorted(zs_dict.keys())
        fields = zs_fields + [hru_param.fid_field]
        with arcpy.da.UpdateCursor(polygon_path, fields, subset_str) as u_cursor:
            for row in u_cursor:
                # Create an empty dictionary if FID does not exist
                # Missing FIDs did not have zonal stats calculated
                row_dict = data_dict.get(int(row[-1]), None)
                for i, zs_field in enumerate(zs_fields):
                    # If stats were calculated for only some parameters,
                    #   then set missing parameter value to nodata value (-999)
                    if row_dict:
                        try:
                            row[i] = row_dict[zs_field]
                        except KeyError:
                            row[i] = nodata_value
                    # Otherwise, if no stats were calculated,
                    #   reset value to 0 (shapefile default)
                    else:
                        row[i] = default_value
                u_cursor.updateRow(row)

        # Cleanup
        del data_dict
        if arcpy.Exists(point_subset_path):
            arcpy.Delete_management(point_subset_path)
        if arcpy.Exists(hru_raster_path):
            arcpy.Delete_management(hru_raster_path)

    arcpy.ClearEnvironment('extent')
    arcpy.ClearEnvironment('outputCoordinateSystem')
    arcpy.ClearEnvironment('cellSize')


def field_duplicate_check(table_path, field_name, n=None):
    """Check if there are duplicate values in a shapefile field

    For now assume table_path is actually a shapefile that can be read
        with arcpy.da.SearchCursor()

    Args:
        table_path (str): File path of the table to search
        field_name (str): Field/column name to search

    Returns:
        bool: True if there are duplicate values in the field, False otherwise
    """

    # Eventually check that field is in table
    field_obj = arcpy.ListFields(table_path, field_name)[0]

    if n is None:
        n = int(arcpy.GetCount_management(table_path).getOutput(0))
    n32_max = 3000000
    logging.debug('\n  Testing for duplicate values')
    logging.debug('    field:    {}'.format(field_name))
    logging.debug('    features: {}'.format(n))
    logging.debug('    n32_max:  {}'.format(n32_max))
    logging.debug('    2**32:    {}'.format(2 ** 32))
    logging.debug('    maxsize:  {}'.format(sys.maxsize))

    if sys.maxsize > 2**32 or n < n32_max:
        logging.debug('    Reading values')
        # If 64-bit or row count is low, read all values into memory
        fid_list = [
            r[0] for r in arcpy.da.SearchCursor(table_path, [field_name])]
        if len(fid_list) != len(set(fid_list)):
            return True
        else:
            logging.debug('    No duplicates')
            return False
    elif field_obj.type in ['Integer', 'SmallInteger']:
        # This approach will only work with integers
        block_size = 500000
        fid_ranges = []
        for i, x in enumerate(xrange(0, n, block_size)):
            logging.debug('    FIDS: {}-{}'.format(x, x + block_size))
            subset_str = '"{0}" >= {1} AND "{0}" < {2}'.format(
                arcpy.Describe(table_path).OIDFieldName, x, x + block_size)

            # Don't sort here since it gets sorted in group_ranges()
            fid_list = [r[0] for r in arcpy.da.SearchCursor(
                table_path, [field_name], subset_str)]

            # Return True if there are duplicates in a subset
            if len(fid_list) != len(set(fid_list)):
                return True

            # Group consecutive values into ranges
            fid_new_ranges = list(group_ranges(fid_list))

            # Check if any of the new ranges overlap each other
            # Skip for now since subset duplicates were checked above
            # if ranges_overlap(fid_new_ranges):
            #    return True

            # Check if any of the new ranges overlap the existing ranges
            if ranges_overlap(fid_ranges + fid_new_ranges):
                return True

            # Merge subset ranges into main range list
            fid_ranges = list(merge_ranges(fid_ranges + fid_new_ranges))
            del fid_new_ranges

        logging.debug('    FID ranges: {}'.format(fid_ranges))
        logging.debug('    No duplicates')
        return False
    else:
        # For now, assume there are not duplicates if the values can't
        #   all be read in
        logging.debug(
            '    Assuming no duplicates since field type is not integer, \n'
            '      table contains more than {} rows, and 32-bit Python is \n'
            '      limited to 2 GB of memory.')
        return False
    # DEADBEEF - I'm not sure this will actually work on a really large table
    # else:
    #    duplicate_flag = False
    #    fid_prev = None
    #    cursor = arcpy.SearchCursor(
    #        table_path, fields=field_name, sort_fields=field_name+' A')
    #    for row in cursor:
    #        fid = row.getValue(field_name)
    #        if fid == fid_prev and fid_prev is not None:
    #            duplicate_flag = True
    #            break
    #        else:
    #            fid_prev = fid
    #    del cursor, row
    #    logging.debug('    No duplicates')
    #    return duplicate_flag


def group_ranges(input_list):
    """Group

    Copied from:
    http://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list

    Args:
        input_list (list): list of numbers to group into ranges

    Yields
         tuple: pairs (min, max)
    """
    for k, g in itertools.groupby(enumerate(sorted(input_list)), lambda (i, x): i-x):
        group = map(itemgetter(1), g)
        yield group[0], group[-1]


def merge_ranges(ranges):
    """Merge overlapping and adjacent integer ranges

    Yield the merged ranges in order
    The argument must be an iterable of pairs (start, stop).

    Copied from:
    http://codereview.stackexchange.com/questions/21307/consolidate-list-of-ranges-that-overlap

    >>> list(merge_ranges([(5,7), (3,5), (-1,3)]))
    [(-1, 7)]
    >>> list(merge_ranges([(5,6), (3,4), (1,2)]))
    [(1, 2), (3, 4), (5, 6)]
    >>> list(merge_ranges([]))
    []

    Args:
        ranges (list): Iterable of pairs (min, max)

    Yields:
        tuple: pairs (min, max)
    """
    ranges = iter(sorted(ranges))
    current_start, current_stop = next(ranges)
    for start, stop in ranges:
        if start > (current_stop + 1):
            # Gap between segments: output current segment and start a new one.
            yield current_start, current_stop
            current_start, current_stop = start, stop
        else:
            # Segments adjacent or overlapping: merge.
            current_stop = max(current_stop, stop)
    yield current_start, current_stop


def ranges_overlap(ranges):
    """Test if ranges overlap

    Args:
        ranges (list): Iterable of pairs (min, max)
    Returns:
         bool: True if ranges overlap each other, False otherwise
    """
    for r1, r2 in itertools.combinations(ranges, 2):
        if r1[1] > r2[0] and r1[0] < r2[1]:
            return True
    return False


def extent_string(extent_obj):
    """"""
    return ' '.join(str(extent_obj).split()[:4])
    # return ' '.join(['{0:.4f}'.format(s) for s in str(extent_obj).split()[:4]])


def round_extent(extent_obj, n=10):
    """"""
    return arcpy.Extent(
        round(extent_obj.XMin, n), round(extent_obj.YMin, n),
        round(extent_obj.XMax, n), round(extent_obj.YMax, n))


# This adjusts one extent to a snap point
# This is similar to the GDAL implementation
def adjust_extent_to_snap(extent_obj, snap_pnt, cs, method='EXPAND', digits=0):
    """"""
    if method.upper() == 'ROUND':
        extent_xmin = math.floor(
            (extent_obj.XMin - snap_pnt.X) / cs + 0.5) * cs + snap_pnt.X
        extent_ymin = math.floor(
            (extent_obj.YMin - snap_pnt.Y) / cs + 0.5) * cs + snap_pnt.Y
        extent_xmax = math.floor(
            (extent_obj.XMax - snap_pnt.X) / cs + 0.5) * cs + snap_pnt.X
        extent_ymax = math.floor(
            (extent_obj.YMax - snap_pnt.Y) / cs + 0.5) * cs + snap_pnt.Y
    elif method.upper() == 'EXPAND':
        extent_xmin = math.floor(
            (extent_obj.XMin - snap_pnt.X) / cs) * cs + snap_pnt.X
        extent_ymin = math.floor(
            (extent_obj.YMin - snap_pnt.Y) / cs) * cs + snap_pnt.Y
        extent_xmax = math.ceil(
            (extent_obj.XMax - snap_pnt.X) / cs) * cs + snap_pnt.X
        extent_ymax = math.ceil(
            (extent_obj.YMax - snap_pnt.Y) / cs) * cs + snap_pnt.Y
    elif method.upper() == 'SHRINK':
        extent_xmin = math.ceil(
            (extent_obj.XMin - snap_pnt.X) / cs) * cs + snap_pnt.X
        extent_ymin = math.ceil(
            (extent_obj.YMin - snap_pnt.Y) / cs) * cs + snap_pnt.Y
        extent_xmax = math.floor(
            (extent_obj.XMax - snap_pnt.X) / cs) * cs + snap_pnt.X
        extent_ymax = math.floor(
            (extent_obj.YMax - snap_pnt.Y) / cs) * cs + snap_pnt.Y

    extent_list = [extent_xmin, extent_ymin, extent_xmax, extent_ymax]
    extent_list = [round(float(x), digits) for x in extent_list]
    if digits == 0:
        extent_list = [int(x) for x in extent_list]
    return arcpy.Extent(*extent_list)


def buffer_extent_func(extent_obj, extent_buffer):
    """"""
    # extent_obj = arcpy.Describe(extent_feature).extent
    extent_xmin = extent_obj.XMin - extent_buffer
    extent_ymin = extent_obj.YMin - extent_buffer
    extent_xmax = extent_obj.XMax + extent_buffer
    extent_ymax = extent_obj.YMax + extent_buffer
    return arcpy.Extent(
        extent_xmin, extent_ymin, extent_xmax, extent_ymax)


# Check if rasters are aligned to snap_raster
# Check if rasters have same cellsize as snap_raster
def snapped(extent_obj, snap_pnt, cs, tol=0.001):
    """Compute difference between value and value rounded to the nearest even
    mulitple of the cell size"""
    def diff(value, cs):
        return abs(value - math.floor(value / cs + 0.5) * cs)
    if ((diff(snap_pnt.X - extent_obj.XMin, cs) <= tol) and
            (diff(snap_pnt.X - extent_obj.XMax, cs) <= tol) and
            (diff(snap_pnt.Y - extent_obj.YMin, cs) <= tol) and
            (diff(snap_pnt.Y - extent_obj.YMax, cs) <= tol)):
        return True
    else:
        return False
    # if (((snap_pnt.X - extent_obj.XMin) % cs == 0) and
    #         ((snap_pnt.X - extent_obj.XMax) % cs == 0) and
    #         ((snap_pnt.Y - extent_obj.YMin) % cs == 0) and
    #         ((snap_pnt.Y - extent_obj.YMax) % cs == 0)):
    #     return True


def get_ini_file(workspace, ini_re, function_str='function'):
    """"""
    # Get ini file name
    ini_file_list = build_file_list(workspace, ini_re)
    # Filter field list ini file
    ini_file_list = [
        item for item in ini_file_list if '_field_list.ini' not in item]
    if len(ini_file_list) == 1:
        config_filepath = ini_file_list[0]
    elif len(ini_file_list) > 1:
        ini_file_len_max = max([len(item) for item in ini_file_list])
        # An ini file was not passed as an arguement to the script
        # Look for ini files in the working directory
        # If only one, use it
        # If more, let user pick from list)
        # If none, error out
        print('\nThere is more than one INI file present in the folder')
        print('  {0:2s}  {1}'.format('# ', 'INI File'))
        print('  {0:2s}  {1}'.format('==', '=' * ini_file_len_max))
        for i, ini_file in enumerate(ini_file_list):
            print('  {0:2d}  {1}'.format(i, ini_file))
        config_filepath = None
        while not config_filepath:
            usr_input = raw_input('\nPlease select an INI file to use: ')
            try:
                ini_file_index = int(usr_input)
                config_filepath = ini_file_list[ini_file_index]
            except (ValueError, IndexError):
                pass
        print('  Using {}\n'.format(config_filepath))
        del ini_file_len_max, usr_input
    else:
        print(
            '\nERROR: No suitable ini files were found\n'
            'ERROR: Please set input file when calling {}\n'
            'ERROR: For example: test.py test.ini\n'.format(function_str))
        sys.exit()
    config_filename = os.path.basename(config_filepath)
    print('{0:<20s} {1}'.format('INI File Name:', config_filename))
    return config_filepath


def get_param(param_str, param_default, config, section='INPUTS'):
    """"""
    param_type = type(param_default)
    try:
        if param_type is float:
            param_value = config.getfloat('INPUTS', param_str)
        elif param_type is int:
            param_value = config.getint('INPUTS', param_str)
        elif param_type is bool:
            param_value = config.getboolean('INPUTS', param_str)
        elif param_type is list or param_type is tuple:
            param_value = [
                # i for i in re.split('\W+', config.get('INPUTS', param_str)) if i]
                i.strip() for i in config.get('INPUTS', param_str).split(',')
                if i.strip()]
        elif param_type is str or param_default is None:
            param_value = config.get('INPUTS', param_str)
            if param_value.upper() == 'NONE':
                param_value = None
        else:
            logging.error('ERROR: Unknown Input Type: {}'.format(param_type))
            sys.exit()
    except Exception as e:
        param_value = param_default
        if param_type is str and param_value.upper() == 'NONE':
            param_value = None
        logging.warning('  NOTE: {} = {}'.format(param_str, param_value))
    return param_value


def build_file_list(ws, test_re, test_other_re=None):
    """"""
    if test_other_re is None:
        test_other_re = re.compile('a^')
    if os.path.isdir(ws):
        return sorted([
            os.path.join(ws, item) for item in os.listdir(ws)
            if (os.path.isfile(os.path.join(ws, item)) and
                test_re.match(item) or test_other_re.match(item))])
    else:
        return []


def get_prism_data_name():
    """"""
    #  Get PRISM data name
    data_name_dict = dict()
    data_name_dict[1] = 'PPT'
    data_name_dict[2] = 'TMAX'
    data_name_dict[3] = 'TMIN'
    data_name_dict[4] = 'ALL'
    print('\nPlease select which PRISM product(s) to calculate')
    print('  {0:2s}  {1}'.format('# ', 'PRISM Data'))
    print('  {0:2s}  {1}'.format('==', '=========='))
    for i, data_name in sorted(data_name_dict.items()):
        print('  {0:2d}  {1}'.format(i, data_name))
    data_name = None
    while not data_name:
        usr_input = raw_input(
            '\nPlease select a PRISM data product to calculate: ')
        try:
            data_name_index = int(usr_input)
            data_name = data_name_dict[data_name_index]
        except (ValueError, IndexError):
            pass
    print('  Using {}\n'.format(data_name))
    return data_name


def project_hru_extent_func(hru_extent, hru_cs, hru_sr,
                            target_extent, target_cs, target_sr):
    """"""
    logging.debug('\n  Projecting extent')
    logging.debug('  HRU Extent:   {}'.format(extent_string(hru_extent)))
    logging.debug('  HRU cellsize: {}'.format(hru_cs))
    logging.debug('  HRU spatref:  {}'.format(hru_sr.name))
    logging.debug('  Target snap:     {}'.format(target_extent.lowerLeft))
    logging.debug('  Target cellsize: {}'.format(target_cs))
    logging.debug('  Target spatref:  {}'.format(target_sr.name))

    # DEADBEEF - Arc10.2 ProjectRaster does not honor extent
    # Project the HRU extent to the raster spatial reference
    hru_corners = [
        [hru_extent.XMin, hru_extent.YMax],
        [hru_extent.XMax, hru_extent.YMax],
        [hru_extent.XMax, hru_extent.YMin],
        [hru_extent.XMin, hru_extent.YMin],
        [hru_extent.XMin, hru_extent.YMax]]

    # Add points between corners
    hru_points = []
    for point_a, point_b in zip(hru_corners[:-1], hru_corners[1:]):
        steps = float(max(
            abs(point_b[0] - point_a[0]),
            abs(point_b[1] - point_a[1]))) / hru_cs
        for x, y in zip(np.linspace(point_a[0], point_b[0], steps + 1),
                        np.linspace(point_a[1], point_b[1], steps + 1)):
            hru_points.append(arcpy.Point(x, y))

    # Project all points to output spatial reference and get projected extent
    transform = transform_func(hru_sr, target_sr)
    if transform:
        projected_extent = arcpy.Polygon(
            arcpy.Array(hru_points), hru_sr).projectAs(
                target_sr, transform).extent
    else:
        projected_extent = arcpy.Polygon(
            arcpy.Array(hru_points), hru_sr).projectAs(target_sr).extent
    logging.debug('  Projected Extent: {}'.format(
        extent_string(projected_extent)))

    # Adjust extent to match snap
    projected_extent = adjust_extent_to_snap(
        projected_extent, target_extent.lowerLeft, target_cs,
        method='EXPAND', digits=6)
    logging.debug('  Snapped Extent:   {}'.format(
        extent_string(projected_extent)))

    # Buffer extent 4 input cells
    projected_extent = buffer_extent_func(projected_extent, 4 * target_cs)
    # This will cause problems when target cellsize is in decimal degrees
    # projected_extent = buffer_extent_func(
    #     projected_extent, 4 * max(target_cs, hru_cs))
    logging.debug('  Buffered Extent:  {}'.format(
        extent_string(projected_extent)))
    return projected_extent


def project_raster_func(input_raster, output_raster, output_sr,
                        proj_method, output_cs, transform_str,
                        reg_point, input_sr, hru_param, in_memory=True):
    """"""
    # Input raster can be a raster object or a raster path
    # print isinstance(input_raster, Raster), isinstance(input_raster, str)
    # cellsize is the "actual" input cellsize
    #   and is needed to get the snapping
    # This could be passed as an input to the function
    try:
        input_extent = arcpy.sa.Raster(input_raster).extent
        input_cs = arcpy.sa.Raster(input_raster).meanCellWidth
    except Exception as e:
        input_extent = input_raster.extent
        input_cs = input_raster.meanCellWidth

    # DEADBEEF - ArcGIS 10.2+ ProjectRaster function does not honor extent
    # Clip the input raster with the projected HRU extent first
    # Project extent from "output" to "input" to get clipping extent
    proj_extent = project_hru_extent_func(
        hru_param.extent, hru_param.cs, output_sr,
        input_extent, input_cs, input_sr)

    if in_memory:
        clip_path = os.path.join('in_memory', 'clip_raster')
    else:
        clip_path = output_raster.replace('.img', '_clip.img')

    env.extent = proj_extent
    arcpy.Clip_management(
        input_raster, ' '.join(str(proj_extent).split()[:4]), clip_path)
    arcpy.ClearEnvironment('extent')

    # Then project the clipped raster
    arcpy.ProjectRaster_management(
        clip_path, output_raster, output_sr, proj_method.upper(), output_cs,
        transform_str, reg_point, input_sr)

    # Cleanup
    arcpy.Delete_management(clip_path)


def cell_area_func(hru_param_path, area_field):
    """"""
    arcpy.CalculateField_management(
        hru_param_path, area_field, '!SHAPE.AREA@acres!', 'PYTHON')


def zone_by_area_func(zone_path, zone_field, zone_value, hru_param_path,
                      hru_param, hru_area_field='HRU_AREA',
                      zone_area_field=None, area_pct=50):
    """Flag cells that are inside a feature based on an area weighting

    Set values that are in zone, but don't reset values that are out of zone

    Args:
        zone_path (str):
        zone_field (str):
        zone_value (int):
        hru_param_path (str):
        hru_param: class:`HRUParameters`
        hru_area_field (str):
        zone_area_field (str):
        area_pct ():

    Returns:
        None
    """
    zone_value_field = 'ZONE_VALUE'
    int_area_field = 'INT_AREA'
    int_pct_field = 'INT_PCT'

    # Need to set zone value into a field before intersect
    # If zone_value is FID, add 1 so that only non-lake cells are 0
    arcpy.AddField_management(zone_path, zone_value_field, 'LONG')
    arcpy.AddField_management(zone_path, int_area_field, 'DOUBLE')
    arcpy.AddField_management(zone_path, int_pct_field, 'DOUBLE')
    if zone_value == arcpy.Describe(zone_path).OIDFieldName:
        arcpy.CalculateField_management(
            zone_path, zone_value_field,
            '!{}! + 1'.format(zone_value), 'PYTHON')

    # If zone value is an INT, save it into a field first
    elif type(zone_value) is int:
        # zone_value = int(zone_value)
        arcpy.CalculateField_management(
            zone_path, zone_value_field, zone_value, 'PYTHON')
    # Use zone_value field directly
    else:
        # zone_value_field = zone_value
        arcpy.CalculateField_management(
            zone_path, zone_value_field,
            '!{}!'.format(zone_value), 'PYTHON')

    # Calculate area of HRU cell if necessary
    # if not arcpy.ListFields(zone_path, area_field):
    #    arcpy.AddField_management(zone_path, area_field, 'DOUBLE')
    #    cell_area_func(zone_path, area_field)

    # Intersect the zone layer with the fishnet
    # zone_int_path = os.path.join('in_memory', 'hru_lakes')
    zone_int_path = zone_path.replace('.shp', '_intersect.shp')
    arcpy.Intersect_analysis(
        (hru_param_path, zone_path), zone_int_path)

    # Calculate using cell_area_func to force units to match
    cell_area_func(zone_int_path, int_area_field)

    n = int(arcpy.GetCount_management(zone_int_path).getOutput(0))
    block_size = 200000
    for i, x in enumerate(range(0, n, block_size)):
        logging.debug('    FIDS: {}-{}'.format(x, x + block_size))
        subset_str = '"{0}" >= {1} AND "{0}" < {2}'.format(
            hru_param.fid_field, x, x + block_size)
            # arcpy.Describe(zone_int_path).OIDFieldName, x, x + block_size)

        # Read in FID of selected cells
        hru_cell_dict = dict()
        fields = [
            hru_param.fid_field, hru_area_field,
            int_area_field, zone_value_field]
        with arcpy.da.SearchCursor(zone_int_path, fields, subset_str) as s_cursor:
            for row in s_cursor:
                if (100 * float(row[2]) / float(row[1])) >= area_pct:
                    hru_cell_dict[int(row[0])] = [float(row[3]), float(row[2])]

        # Set value of selected HRU cells
        fields = [hru_param.fid_field, zone_field]
        if zone_area_field:
            fields.append(zone_area_field)
        with arcpy.da.UpdateCursor(hru_param_path, fields, subset_str) as u_cursor:
            for row in u_cursor:
                # Remove items to speed up subsequent searches
                try:
                    if len(fields) == 3:
                        row[1], row[2] = hru_cell_dict.pop(int(row[0]))
                    elif len(fields) == 2:
                        row[1] = hru_cell_dict.pop(int(row[0]))
                    u_cursor.updateRow(row)
                except KeyError:
                    pass
        del hru_cell_dict


def zone_by_centroid_func(zone_path, zone_field, zone_value,
                          hru_param_path, hru_point_path, hru_param):
    """Flag cells that are inside a feature based on the centroid location

    Set values that are in zone, but don't reset values that are out of zone

    Args:
        zone_path (str):
        zone_field (str):
        zone_value (int):
        hru_param_path (str):
        hru_point_path (str):
        hru_param: class:`HRUParameters`

    Returns:
        None
    """
    logging.debug('\nzone_by_centroid_func')
    logging.debug('  {}'.format(zone_path))
    # Need to set zone value into a field before intersect
    # If zone_value is FID, add 1 so that only zone cells are 0
    zone_value_field = 'ZONE_VALUE'
    arcpy.AddField_management(zone_path, zone_value_field, 'LONG')
    if zone_value == arcpy.Describe(zone_path).OIDFieldName:
        arcpy.CalculateField_management(
            zone_path, zone_value_field,
            '!{}! + 1'.format(zone_value), 'PYTHON')
    # Save zone value into a field first
    elif type(zone_value) is int:
        arcpy.CalculateField_management(
            zone_path, zone_value_field, zone_value, 'PYTHON')
        # DEADBEEF
        # Use zone_value field directly
        # zone_value = int(zone_value)
    else:
        arcpy.CalculateField_management(
            zone_path, zone_value_field,
            '!{}!'.format(zone_value), 'PYTHON')
        # DEADBEEF
        # Use zone_value field directly
        # zone_value_field = zone_value

    # Intersect the zone layer with the fishnet
    # zone_int_path = os.path.join('in_memory', 'hru_ppt_zones')
    zone_int_path = zone_path.replace('.shp', '_intersect.shp')
    arcpy.Intersect_analysis(
        (hru_point_path, zone_path), zone_int_path)

    # DEADBEEF - Why do I make a layer and select all features?
    zone_int_layer = 'zone_int_layer'
    arcpy.MakeFeatureLayer_management(zone_int_path, zone_int_layer)
    arcpy.SelectLayerByAttribute_management(zone_int_layer, 'CLEAR_SELECTION')
    arcpy.SelectLayerByAttribute_management(zone_int_layer, 'SWITCH_SELECTION')

    n = int(arcpy.GetCount_management(zone_int_layer).getOutput(0))
    block_size = 200000
    for i, x in enumerate(range(0, n, block_size)):
        logging.debug('    FIDS: {}-{}'.format(x, x + block_size))
        subset_str = '"{0}" >= {1} AND "{0}" < {2}'.format(
            hru_param.fid_field, x, x + block_size)
            # arcpy.Describe(zone_int_layer).OIDFieldName, x, x + block_size)

        # Read in FID of selected cells
        hru_point_dict = dict()
        fields = (hru_param.fid_field, zone_value_field)
        with arcpy.da.SearchCursor(zone_int_layer, fields, subset_str) as s_cursor:
            for row in s_cursor:
                hru_point_dict[int(row[0])] = row[1]

        # Set value of selected HRU cells
        fields = (hru_param.fid_field, zone_field)
        with arcpy.da.UpdateCursor(hru_param_path, fields, subset_str) as u_cursor:
            for row in u_cursor:
                # Remove items to speed up subsequent searches
                try:
                    row[1] = hru_point_dict.pop(int(row[0]))
                    u_cursor.updateRow(row)
                except KeyError:
                    pass
        del hru_point_dict

    # Cleanup
    arcpy.Delete_management(zone_int_layer)


def jensen_haise_func(hru_param_path, jh_coef_field, dem_field,
                      jh_tmin_field, jh_tmax_field, unit_scalar):
    """"""
    jh_cb = (
        'def ea(temp_c):\n'
        '    return 6.1078 * math.exp((17.269 * temp_c) / (temp_c + 237.3))\n'
        'def jensen_haise(elev, t_low, t_high, unit_scalar):\n'
        '    return 27.5 - 0.25 * (ea(t_high) - ea(t_low)) - (elev * unit_scalar / 1000)\n')
    arcpy.CalculateField_management(
        hru_param_path, jh_coef_field,
        'jensen_haise(!{}!, !{}!, !{}!, {})'.format(
            dem_field, jh_tmin_field, jh_tmax_field, unit_scalar),
        'PYTHON', jh_cb)


def remap_check(remap_path):
    """"""
    # Check that the file exists
    if not os.path.isfile(remap_path):
        logging.error(
            '\nERROR: ASCII remap file ({}) does not exist\n'.format(
                os.path.basename(remap_path)))
        sys.exit()

    # Read in the remap
    with open(remap_path, 'r') as remap_f:
        remap_lines = remap_f.readlines()
        line_count = len(remap_lines)

    # Problems reading/applying ASCII remap files can be caused by:
    #   Blank lines
    #   Empty lines at the end (from a final newline character)
    #   ArcGIS 10.2 - Old style in line comments (/*)
    #   ArcGIS 10.2 - Comments can't be in first line (?)
    #   ArcGIS 10.2 - Comments can't be longer than 80 characters (?)
    # If either of these are present in the file, resave the filtered lines

    # Check for old style comments (/*) in ASCII remap files
    # This could be changed to save the comments at the end of the file
    if arcpy.GetInstallInfo()['Version'].startswith('10.2'):
        if any([l for l in remap_lines if '/*' in l]):
            logging.error((
                '\nERROR: ASCII remap file ({}) has pre-ArcGIS 10.2 '
                'comments (\*)\n  Try running the "convert_remap_arc10p2.py"'
                'script\n').format(os.path.basename(remap_path)))
            sys.exit()

    # First check for final newline character
    save_flag = False
    if remap_lines and remap_lines[-1] and remap_lines[-1].endswith('\n'):
        logging.debug('  Final newline character')
        save_flag = True

    # Then remove empty lines and strip white space and newline characters
    # If lines were removed, resave the filtered remap file
    remap_lines = [l.strip() for l in remap_lines]
    remap_lines = [l for l in remap_lines if l]
    if len(remap_lines) != line_count:
        logging.debug('  Whitespace or empty lines')
        save_flag = True

    # Trim comments longer than 80 characters
    if arcpy.GetInstallInfo()['Version'].startswith('10.2'):
        if any([len(l) > 80 for l in remap_lines if "#" in l]):
            remap_lines = [l[:79] if "#" in l else l for l in remap_lines]
            logging.debug('  Lines longer than 80 characters')
            save_flag = True

    # If lines were removed, resave the filtered remap file
    if save_flag:
        logging.warning(
            '  The ASCII remap file ({}) will be overwritten'.format(
                os.path.basename(remap_path)))
        with open(remap_path, 'w') as remap_f:
            for i, line in enumerate(remap_lines):
                # Don't write newline character on last line
                # This causes an error in ArcGIS 10.2.2
                if (i + 1) < len(remap_lines):
                    remap_f.write(line + '\n')
                else:
                    remap_f.write(line)
    return True


def remap_code_block(remap_path):
    """"""
    with open(remap_path) as remap_f:
        lines = remap_f.readlines()
    remap_cb = ''
    for line in lines:
        # Skip comment lines
        if '#' in line:
            continue
        # Remove remap description
        line = line.strip().split('/*')[0]
        # Split line on spaces and semi-colon
        l_split = [item.strip() for item in re.split('[ :]+', line)]
        # Remap as a range if a min, max and value are all present
        if len(l_split) == 3:
            range_remap_flag = True
        # Otherwise remap directly
        elif len(l_split) == 2:
            range_remap_flag = False
        # Skip lines that don't match format
        else:
            continue
        # Write remap code block
        if not range_remap_flag:
            if not remap_cb:
                remap_cb = '    if value == {}: return {}\n'.format(*l_split)
            else:
                remap_cb += ('    elif value == {}: '
                             'return {}\n'.format(*l_split))
        else:
            if not remap_cb:
                remap_cb = (
                    '    if (value >= {} and value <= {}): '
                    'return {}\n').format(*l_split)
            else:
                remap_cb += (
                    '    elif (value > {} and value <= {}): '
                    'return {}\n').format(*l_split)
    remap_cb = 'def Reclass(value):\n' + remap_cb
    return remap_cb


# def reclass_ascii_float_func(raster_path, remap_path):
#    # Read remap file into memory
#    with open(remap_path) as remap_f: lines = remap_f.readlines()
#    remap_f.close()
#    first_line = True
#    raster_obj = arcpy.sa.Raster(raster_path)
#    for l in lines:
#        # Skip comment lines
#        if '#' in l: continue
#        # Remove remap description
#        l = l.split('/*')[0]
#        # Split line on spaces and semi-colon
#        l_split = map(float, [item for item in re.split('[ :]+', l) if item])
#        # Remap as a range if a min, max and value are all present
#        if len(l_split) == 3: range_remap_flag = True
#        # Otherwise remap directly
#        elif len(l_split) == 2: range_remap_flag = False
#        # Skip lines that don't match format
#        else: continue
#        # Write remap code block
#        if not range_remap_flag:
#            raster_obj = Con(raster_obj == l_split[0], l_split[1], raster_obj)
#        elif first_line:
#            raster_obj = Con(
#                ((raster_obj >= l_split[0]) & (raster_obj <= l_split[1])),
#                l_split[2], raster_obj)
#        else:
#            raster_obj = Con(
#                ((raster_obj > l_split[0]) & (raster_obj <= l_split[1])),
#                l_split[2], raster_obj)
#        first_line = False
#    return raster_obj


def is_number(s):
    """"""
    try:
        float(s)
        return True
    except ValueError:
        return False


def raster_path_to_array(input_path, mask_extent=None, return_nodata=False):
    """"""
    return raster_obj_to_array(
        arcpy.sa.Raster(input_path), mask_extent, return_nodata)


def raster_obj_to_array(input_obj, mask_extent=None, return_nodata=False):
    """"""
    input_nodata = input_obj.noDataValue
    input_cs = input_obj.meanCellHeight
    input_rows, input_cols = input_obj.height, input_obj.width
    input_extent = input_obj.extent
    output_array = arcpy.RasterToNumPyArray(input_obj)
    # DEADBEEF - get_extent_intersection() and extent_shape() are not present
    #   Commenting out for now
    # if mask_extent:
    #     int_extent = get_extent_intersection([input_extent, mask_extent])
    #     int_pnt = arcpy.Point()
    #     int_pnt.X = int_extent.XMin
    #     int_pnt.Y = int_extent.YMin
    #     int_rows, int_cols = extent_shape(int_extent, input_cs)
    #     output_array = arcpy.RasterToNumPyArray(
    #         input_obj, int_pnt, int_cols, int_rows)
    # else:
    #     output_array = arcpy.RasterToNumPyArray(input_obj)
    # Integer type raster can't have NaN values, will only set floats to NaN
    if (output_array.dtype == np.float32 or
        output_array.dtype == np.float64):
        output_array[output_array == input_nodata] = np.NaN
        output_nodata = np.NaN
    else:
        output_nodata = int(input_nodata)
    if return_nodata:
        return output_array, output_nodata
    else:
        return output_array


def array_to_raster(input_array, output_path, pnt, cs, mask_array=None):
    """"""
    output_array = np.copy(input_array)

    # Float arrays have to have nodata set to some value (-9999)
    if (output_array.dtype == np.float32 or
            output_array.dtype == np.float64):
        output_nodata = -9999
        output_array[np.isnan(output_array)] = output_nodata
    # Boolean arrays need to be converted unsigned ints
    elif output_array.dtype == np.bool:
        output_array = output_array.astype(np.uint8)
        output_nodata = 255
    elif output_array.dtype == np.uint8:
        output_nodata = 255

    # If a mask array is give, assume all 0 values are nodata
    if np.any(mask_array):
        output_array[mask_array == 0] = output_nodata

    output_obj = arcpy.NumPyArrayToRaster(
        output_array, pnt, cs, cs, output_nodata)
    output_obj.save(output_path)
    del output_obj
    arcpy.DefineProjection_management(
        output_path, env.outputCoordinateSystem)
    arcpy.CalculateStatistics_management(output_path)


# def flood_fill(test_array, four_way_flag=True, edge_flt=None):
#     """Flood fill algorithm"""
#     input_array = np.copy(test_array)
#     input_rows, input_cols = input_array.shape
#     h_max = np.nanmax(input_array) + 10.0
#     logging.debug("  Hmax: %s" % (h_max / 2.0))
#
#     # Since ArcGIS doesn't ship with SciPy (only numpy), don't use ndimage module
#     if four_way_flag:
#         el = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]]).astype(np.bool)
#     else:
#         el = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]]).astype(np.bool)
#     # if four_way_flag:
#     #    el = ndimage.generate_binary_structure(2,1).astype(np.int)
#     # else:
#     #    el = ndimage.generate_binary_structure(2,2).astype(np.int)
#
#     # Build data/inside/edge masks
#     data_mask = ~np.isnan(input_array)
#
#     # Since ArcGIS doesn't ship with SciPy (only numpy), don't use ndimage module
#     inside_mask = np_binary_erosion(data_mask, structure=el)
#     # inside_mask = ndimage.binary_erosion(data_mask, structure=el)
#
#     edge_mask = (data_mask & ~inside_mask)
#     # Initialize output array as max value test_array except edges
#     output_array = np.copy(input_array)
#     output_array[inside_mask] = h_max
#     # Set edge pixels less than edge_flt to edge_flt
#     if edge_flt:
#         output_array[edge_mask & (output_array <= edge_flt)] = edge_flt
#
#     # Build priority queue and place edge pixels into queue
#     put = heapq.heappush
#     get = heapq.heappop
#     fill_heap = [
#         (output_array[t_row, t_col], int(t_row), int(t_col), 1)
#         for t_row, t_col in np.transpose(np.where(edge_mask))]
#     heapq.heapify(fill_heap)
#     # logging.info("    Queue Size: %s" % len(fill_heap))
#
#     # Cleanup
#     del data_mask, edge_mask, el
#     # logging.info("    Prep Time: %s" % (clock()-start_total))
#
#     while True:
#         try:
#             h_crt, t_row, t_col, edge_flag = get(fill_heap)
#         except IndexError:
#             break
#         for n_row, n_col in [
#             ((t_row-1), t_col), ((t_row+1), t_col),
#             (t_row, (t_col-1)), (t_row, (t_col+1))]:
#             # Skip cell if outside array edges
#             if edge_flag:
#                 try:
#                     if not inside_mask[n_row, n_col]:
#                         continue
#                 except IndexError:
#                     continue
#             if output_array[n_row, n_col] == h_max:
#                 output_array[n_row, n_col] = max(
#                     h_crt, input_array[n_row, n_col])
#                 put(fill_heap,
#                     (output_array[n_row, n_col], n_row, n_col, 0))
#     return output_array
#
#
# def np_binary_erosion(input_array,
#                       structure=np.ones((3, 3)).astype(np.bool)):
#     """NumPy binary erosion function
#
#     No error checking on input array (type)
#     No error checking on structure element (# of dimensions, shape, type, etc.)
#
#     Args:
#         input_array: Binary NumPy array to be eroded. Non-zero (True) elements
#             form the subset to be eroded
#         structure: Structuring element used for the erosion. Non-zero elements
#             are considered True. If no structuring element is provided, an
#             element is generated with a square connectivity equal to one.
#
#     Returns:
#         binary_erosion: Erosion of the input by the stucturing element
#     """
#     rows, cols = input_array.shape
#
#     # Pad output array (binary_erosion) with extra cells around the edge
#     # so that structuring element will fit without wrapping.
#     # A 3x3 structure, will need 1 additional cell around the edge
#     # A 5x5 structure, will need 2 additional cells around the edge
#     output_shape = tuple(
#         ss + dd - 1 for ss, dd in zip(input_array.shape, structure.shape))
#     input_pad_array = np.zeros(output_shape).astype(np.bool)
#     input_pad_array[1: rows+1, 1: cols+1] = input_array
#     binary_erosion = np.zeros(output_shape).astype(np.bool)
#
#     # Cast structure element to boolean
#     struc_mask = structure.astype(np.bool)
#
#     # Iterate over each cell
#     for row in xrange(rows):
#         for col in xrange(cols):
#             # The value of the output pixel is the minimum value of all the
#             #   pixels in the input pixel's neighborhood.
#             binary_erosion[row+1, col+1] = np.min(
#                 input_pad_array[row: row+3, col: col+3][struc_mask])
#     return binary_erosion[1: rows+1, 1: cols+1]
