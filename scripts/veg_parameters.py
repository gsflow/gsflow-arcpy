#--------------------------------
# Name:         veg_parameters
# Purpose:      GSFLOW vegetation parameters
# Notes:        ArcGIS 10.2+ Version
# Created       2017-08-25
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


def veg_parameters(config_path, overwrite_flag=False, debug_flag=False):
    """Calculate GSFLOW Vegetation Parameters

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
    log_file_name = 'veg_parameters_log.txt'
    log_console = logging.FileHandler(
        filename=os.path.join(hru.log_ws, log_file_name), mode='w')
    log_console.setLevel(logging.DEBUG)
    log_console.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger('').addHandler(log_console)
    logging.info('\nGSFLOW Vegetation Parameters')

    # Landfire Vegetation Type
    veg_type_orig_path = inputs_cfg.get('INPUTS', 'veg_type_orig_path')
    veg_type_cs = inputs_cfg.getint('INPUTS', 'veg_type_cellsize')
    try:
        veg_type_field = inputs_cfg.get('INPUTS', 'veg_type_field')
    except ConfigParser.NoOptionError:
        veg_type_field = None
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'veg_type_field', veg_type_field))

    # Landfire Vegetation Cover
    veg_cover_orig_path = inputs_cfg.get('INPUTS', 'veg_cover_orig_path')
    veg_cover_cs = inputs_cfg.getint('INPUTS', 'veg_cover_cellsize')

    # Remap
    remap_ws = inputs_cfg.get('INPUTS', 'remap_folder')
    cov_type_remap_name = inputs_cfg.get('INPUTS', 'cov_type_remap')
    covden_sum_remap_name = inputs_cfg.get('INPUTS', 'covden_sum_remap')
    covden_win_remap_name = inputs_cfg.get('INPUTS', 'covden_win_remap')
    snow_intcp_remap_name = inputs_cfg.get('INPUTS', 'snow_intcp_remap')
    srain_intcp_remap_name = inputs_cfg.get('INPUTS', 'srain_intcp_remap')
    wrain_intcp_remap_name = inputs_cfg.get('INPUTS', 'wrain_intcp_remap')
    root_depth_remap_name = inputs_cfg.get('INPUTS', 'root_depth_remap')

    # Get remap conversion factors
    try:
        snow_intcp_remap_factor = inputs_cfg.getfloat(
            'INPUTS', 'snow_intcp_remap_factor')
    except ConfigParser.NoOptionError:
        snow_intcp_remap_factor = 0.01
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'snow_intcp_remap_factor', snow_intcp_remap_factor))
    try:
        wrain_intcp_remap_factor = inputs_cfg.getfloat(
            'INPUTS', 'wrain_intcp_remap_factor')
    except ConfigParser.NoOptionError:
        wrain_intcp_remap_factor = 0.01
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'wrain_intcp_remap_factor', wrain_intcp_remap_factor))
    try:
        srain_intcp_remap_factor = inputs_cfg.getfloat(
            'INPUTS', 'srain_intcp_remap_factor')
    except ConfigParser.NoOptionError:
        srain_intcp_remap_factor = 0.01
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'srain_intcp_remap_factor', srain_intcp_remap_factor))

    # Check input paths
    if not arcpy.Exists(hru.polygon_path):
        logging.error(
            '\nERROR: Fishnet ({}) does not exist'.format(
                hru.polygon_path))
        sys.exit()
    # Check that either the original vegetation raster exist
    if not arcpy.Exists(veg_cover_orig_path):
        logging.error(
            '\nERROR: Vegetation cover raster does not exist')
        sys.exit()
    if not arcpy.Exists(veg_type_orig_path):
        logging.error(
            '\nERROR: Vegetation type raster does not exist')
        sys.exit()
    # Vegetation cover can be set from another field in the raster
    # This is mostly for US_120EVT
    if not veg_type_field:
        logging.info('\n  Using VALUE field to set vegetation type')
        veg_type_field = 'VALUE'
    elif len(arcpy.ListFields(veg_type_orig_path, veg_type_field)) == 0:
        logging.info(
            '  veg_type_field {} does not exist\n  Using VALUE '
            'field to set vegetation type'.format(veg_type_field))
        veg_type_field = 'VALUE'
    elif arcpy.ListFields(veg_type_orig_path, veg_type_field)[0].type not in ['Integer', 'SmallInteger']:
        logging.info(
            '  veg_type_field {} is not an integer type\n  Using VALUE '
            'field to set vegetation type'.format(veg_type_field))
        veg_type_field = 'VALUE'

    # Check that remap folder is valid
    if not os.path.isdir(remap_ws):
        logging.error('\nERROR: Remap folder does not exist')
        sys.exit()
    # Check that remap files exist
    # Check remap files comment style
    cov_type_remap_path = os.path.join(remap_ws, cov_type_remap_name)
    covden_sum_remap_path = os.path.join(remap_ws, covden_sum_remap_name)
    covden_win_remap_path = os.path.join(remap_ws, covden_win_remap_name)
    snow_intcp_remap_path = os.path.join(remap_ws, snow_intcp_remap_name)
    srain_intcp_remap_path = os.path.join(remap_ws, srain_intcp_remap_name)
    wrain_intcp_remap_path = os.path.join(remap_ws, wrain_intcp_remap_name)
    root_depth_remap_path = os.path.join(remap_ws, root_depth_remap_name)
    remap_path_list = [
        cov_type_remap_path, covden_sum_remap_path, covden_win_remap_path,
        snow_intcp_remap_path, srain_intcp_remap_path,
        wrain_intcp_remap_path, root_depth_remap_path]
    for remap_path in remap_path_list:
        support.remap_check(remap_path)

    # Check other inputs
    if veg_type_cs <= 0:
        logging.error('\nERROR: Veg. type cellsize must be greater than 0')
        sys.exit()
    if veg_cover_cs <= 0:
        logging.error('\nERROR: Veg. cover cellsize must be greater than 0')
        sys.exit()

    # Build output folders if necesssary
    veg_temp_ws = os.path.join(hru.param_ws, 'veg_rasters')
    if not os.path.isdir(veg_temp_ws):
        os.mkdir(veg_temp_ws)
    # Output paths
    veg_cover_path = os.path.join(veg_temp_ws, 'veg_cover.img')
    veg_type_path = os.path.join(veg_temp_ws, 'veg_type.img')
    cov_type_path = os.path.join(veg_temp_ws, 'cov_type.img')
    covden_sum_path = os.path.join(veg_temp_ws, 'covden_sum.img')
    covden_win_path = os.path.join(veg_temp_ws, 'covden_win.img')
    snow_intcp_path = os.path.join(veg_temp_ws, 'snow_intcp.img')
    wrain_intcp_path = os.path.join(veg_temp_ws, 'wrain_intcp.img')
    srain_intcp_path = os.path.join(veg_temp_ws, 'srain_intcp.img')
    root_depth_path = os.path.join(veg_temp_ws, 'root_depth.img')
    rad_trncf_path = os.path.join(veg_temp_ws, 'rad_trncf.img')

    # Set ArcGIS environment variables
    arcpy.CheckOutExtension('Spatial')
    env.overwriteOutput = True
    env.pyramid = 'PYRAMIDS -1'
    # env.pyramid = 'PYRAMIDS 0'
    env.workspace = veg_temp_ws
    env.scratchWorkspace = hru.scratch_ws

    # Check fields
    logging.info('\nAdding vegetation fields if necessary')
    support.add_field_func(hru.polygon_path, hru.cov_type_field, 'SHORT')
    support.add_field_func(hru.polygon_path, hru.covden_sum_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.covden_win_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.rad_trncf_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.snow_intcp_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.srain_intcp_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.wrain_intcp_field, 'DOUBLE')
    # support.add_field_func(hru.polygon_path, hru.root_depth_field, 'DOUBLE')

    # Check that remaps have all necessary values
    logging.info(
        '\nChecking remap tables against all raster cells'
        '  (i.e. even those outside the study area)')
    check_remap_keys(cov_type_remap_path, veg_type_orig_path)
    check_remap_keys(covden_sum_remap_path, veg_cover_orig_path)
    check_remap_keys(root_depth_remap_path, veg_type_orig_path)

    # Assume all vegetation rasters will need to be rebuilt
    # Check veg cover and veg type rasters
    # This will check for matching spat. ref., snap point, and cellsize

    # Project/clip veg cover to match HRU
    logging.info('\nProjecting/clipping vegetation cover raster')
    veg_cover_orig_sr = arcpy.sa.Raster(veg_cover_orig_path).spatialReference
    # Remove existing clipped/projected veg cover raster
    if arcpy.Exists(veg_cover_path):
        arcpy.Delete_management(veg_cover_path)
    # Set preferred transforms
    transform_str = support.transform_func(hru.sr, veg_cover_orig_sr)
    logging.debug('  Transform: {}'.format(transform_str))
    logging.debug('  Projection method: NEAREST')

    # Project veg cover
    # DEADBEEF - Arc10.2 ProjectRaster does not extent
    support.project_raster_func(
        veg_cover_orig_path, veg_cover_path, hru.sr,
        'NEAREST', veg_cover_cs, transform_str,
        '{} {}'.format(hru.ref_x, hru.ref_y), veg_cover_orig_sr, hru)
    # env.extent = hru.extent
    # arcpy.ProjectRaster_management(
    #    veg_cover_orig_path, veg_cover_path, hru.sr,
    #    'NEAREST', veg_cover_cs, transform_str,
    #    '{} {}'.format(hru.ref_x, hru.ref_y),
    #    veg_cover_orig_sr)
    # arcpy.ClearEnvironment('extent')
    del transform_str, veg_cover_orig_sr

    # Project/clip veg type to match HRU
    logging.info('Projecting/clipping vegetation type raster')
    veg_type_orig_sr = arcpy.sa.Raster(veg_type_orig_path).spatialReference
    # Remove existing clipped/projected veg type raster
    if arcpy.Exists(veg_type_path):
        arcpy.Delete_management(veg_type_path)
    # Set preferred transforms
    transform_str = support.transform_func(hru.sr, veg_type_orig_sr)
    logging.debug('  Transform: {}'.format(transform_str))
    logging.debug('  Projection method: NEAREST')
    # Use a different field to calculate vegetation type
    if veg_type_field != 'VALUE':
        logging.info(
            '  Calculating vegetation type from {} field'.format(
                veg_type_field))
        veg_type_obj = arcpy.sa.Lookup(veg_type_orig_path, veg_type_field)
    else:
        veg_type_obj = arcpy.sa.Raster(veg_type_orig_path)

    # Project veg type
    # DEADBEEF - Arc10.2 ProjectRaster does not honor extent
    support.project_raster_func(
        veg_type_obj, veg_type_path, hru.sr,
        'NEAREST', veg_type_cs, transform_str,
        '{} {}'.format(hru.ref_x, hru.ref_y), veg_type_orig_sr, hru)
    # env.extent = hru.extent
    # arcpy.ProjectRaster_management(
    #    veg_type_obj, veg_type_path, hru.sr,
    #    'NEAREST', veg_type_cs, transform_str,
    #    '{} {}'.format(hru.ref_x, hru.ref_y),
    #    veg_type_orig_sr)
    # arcpy.ClearEnvironment('extent')
    del transform_str, veg_type_orig_sr, veg_type_obj

    # Reclassifying vegetation cover type
    logging.info('\nCalculating COV_TYPE')
    logging.debug('  Reclassifying: {}'.format(cov_type_remap_path))
    cov_type_obj = arcpy.sa.ReclassByASCIIFile(
        veg_type_path, cov_type_remap_path)
    cov_type_obj.save(cov_type_path)
    del cov_type_obj

    # Summer cover density
    logging.info('Calculating COVDEN_SUM')
    logging.debug('  Reclassifying: {}'.format(covden_sum_remap_path))
    covden_sum_obj = arcpy.sa.ReclassByASCIIFile(
        veg_cover_path, covden_sum_remap_path)
    covden_sum_obj *= 0.01
    covden_sum_obj.save(covden_sum_path)
    del covden_sum_obj

    # Winter cover density
    logging.info('Calculating COVDEN_WIN')
    logging.debug('  Reclassifying: {}'.format(covden_win_remap_path))
    covden_win_obj = arcpy.sa.ReclassByASCIIFile(
        cov_type_path, covden_win_remap_path)
    covden_win_obj *= 0.01
    covden_win_obj *= arcpy.sa.Raster(covden_sum_path)
    covden_win_obj.save(covden_win_path)
    del covden_win_obj

    # Snow interception storage capacity
    logging.info('Calculating SNOW_INTCP')
    logging.debug('  Reclassifying: {}'.format(snow_intcp_remap_path))
    snow_intcp_obj = arcpy.sa.ReclassByASCIIFile(
        cov_type_path, snow_intcp_remap_path)
    snow_intcp_obj *= snow_intcp_remap_factor
    snow_intcp_obj.save(snow_intcp_path)
    del snow_intcp_obj

    # Winter rain interception storage capacity
    logging.info('Calculating WRAIN_INTCP')
    logging.debug('  Reclassifying: {}'.format(wrain_intcp_remap_path))
    wrain_intcp_obj = arcpy.sa.ReclassByASCIIFile(
        cov_type_path, wrain_intcp_remap_path)
    wrain_intcp_obj *= wrain_intcp_remap_factor
    wrain_intcp_obj.save(wrain_intcp_path)
    del wrain_intcp_obj

    # Summer rain interception storage capacity
    logging.info('Calculating SRAIN_INTCP')
    logging.debug('  Reclassifying: {}'.format(srain_intcp_remap_path))
    srain_intcp_obj = arcpy.sa.ReclassByASCIIFile(
        cov_type_path, srain_intcp_remap_path)
    srain_intcp_obj *= srain_intcp_remap_factor
    srain_intcp_obj.save(srain_intcp_path)
    del srain_intcp_obj

    # Root depth
    logging.info('Calculating ROOT_DEPTH')
    logging.debug('  Reclassifying: {}'.format(root_depth_remap_path))
    root_depth_obj = arcpy.sa.ReclassByASCIIFile(
        veg_type_path, root_depth_remap_path)
    root_depth_obj.save(root_depth_path)
    del root_depth_obj

    # Short-wave radiation transmission coefficent
    logging.info('Calculating {}'.format(hru.rad_trncf_field))
    rad_trncf_obj = 0.9917 * arcpy.sa.Exp(
        -2.7557 * arcpy.sa.Raster(covden_win_path))
    rad_trncf_obj.save(rad_trncf_path)
    del rad_trncf_obj

    # List of rasters, fields, and stats for zonal statistics
    zs_veg_dict = dict()
    zs_veg_dict[hru.cov_type_field] = [cov_type_path, 'MAJORITY']
    zs_veg_dict[hru.covden_sum_field] = [covden_sum_path, 'MEAN']
    zs_veg_dict[hru.covden_win_field] = [covden_win_path, 'MEAN']
    zs_veg_dict[hru.snow_intcp_field] = [snow_intcp_path, 'MEAN']
    zs_veg_dict[hru.srain_intcp_field] = [srain_intcp_path, 'MEAN']
    zs_veg_dict[hru.wrain_intcp_field] = [wrain_intcp_path, 'MEAN']
    # zs_veg_dict[hru.root_depth_field] = [root_depth_path, 'MEAN']
    zs_veg_dict[hru.rad_trncf_field] = [rad_trncf_path, 'MEAN']

    # Calculate zonal statistics
    logging.info('\nCalculating vegetation zonal statistics')
    support.zonal_stats_func(
        zs_veg_dict, hru.polygon_path, hru.point_path, hru)

    # Short-wave radiation transmission coefficient
    # logging.info('\nCalculating {}'.format(hru.rad_trncf_field))
    # arcpy.CalculateField_management(
    #    hru.polygon_path, hru.rad_trncf_field,
    #    '0.9917 * math.exp(-2.7557 * !{}!)'.format(hru.covden_win_field),
    #    'PYTHON')

    # Clear COV_TYPE values for lake cells (HRU_TYPE == 2)
    if True:
        logging.info('\nClearing lake nodata vegetation parameters')
        # logging.info(
        #     '\nClearing vegetation parameters for lake and inactive cells')
        hru_polygon_layer = "hru_polygon_layer"
        arcpy.MakeFeatureLayer_management(
            hru.polygon_path, hru_polygon_layer)
        arcpy.SelectLayerByAttribute_management(
            hru_polygon_layer, "NEW_SELECTION",
            '"{0}" = 2 OR ("{0}" = 0 AND "{1}" = 0)'.format(
                hru.type_field, hru.dem_adj_field))
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.cov_type_field, 0, 'PYTHON')
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.covden_sum_field, 0, 'PYTHON')
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.covden_win_field, 0, 'PYTHON')
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.snow_intcp_field, 0, 'PYTHON')
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.srain_intcp_field, 0, 'PYTHON')
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.wrain_intcp_field, 0, 'PYTHON')
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.rad_trncf_field, 0, 'PYTHON')
        arcpy.Delete_management(hru_polygon_layer)
        del hru_polygon_layer


def get_remap_keys(remap_path):
    """"""
    with open(remap_path) as remap_f:
        lines = remap_f.readlines()
    remap_f.close()
    return [int(l.split(':')[0].strip()) for l in lines if l and '#' not in l]


def get_raster_values(raster_path):
    """"""
    return [
        int(row[0])
        for row in arcpy.da.SearchCursor(raster_path, ['Value'])]


def check_remap_keys(remap_path, raster_path):
    """"""
    logging.info('  {} - {}'.format(
        os.path.basename(remap_path), os.path.basename(raster_path)))
    remap_keys = get_remap_keys(remap_path)
    raster_values = get_raster_values(raster_path)
    missing_keys = sorted(list(set(raster_values) - set(remap_keys)))
    for key in missing_keys:
        logging.warning(
            '    Raster value {} is not in the remap table'.format(key))


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Vegetation Parameters',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i', '--ini', required=True,
        help='Project input file', metavar='PATH')
    parser.add_argument(
        '-o', '--overwrite', default=False, action="store_true",
        help='Force overwrite of existing files')
    parser.add_argument(
        '-d', '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action="store_const", dest="loglevel")
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

    # Calculate GSFLOW Vegetation Parameters
    veg_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
