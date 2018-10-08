#--------------------------------
# Name:         soil_raster_prep.py
# Purpose:      GSFLOW soil raster prep
# Notes:        ArcGIS 10.2+ Version
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


def soil_raster_prep(config_path, overwrite_flag=False, debug_flag=False):
    """Prepare GSFLOW soil rasters

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
    log_file_name = 'soil_prep_log.txt'
    log_console = logging.FileHandler(
        filename=os.path.join(hru.log_ws, log_file_name), mode='w')
    log_console.setLevel(logging.DEBUG)
    log_console.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger('').addHandler(log_console)
    logging.info('\nPrepare GSFLOW Soil Rasters')

    soil_orig_ws = inputs_cfg.get('INPUTS', 'soil_orig_folder')
    awc_name = inputs_cfg.get('INPUTS', 'awc_name')
    clay_pct_name = inputs_cfg.get('INPUTS', 'clay_pct_name')
    sand_pct_name = inputs_cfg.get('INPUTS', 'sand_pct_name')
    soil_proj_method = 'NEAREST'
    soil_cs = inputs_cfg.getint('INPUTS', 'soil_cellsize')
    fill_soil_nodata_flag = inputs_cfg.getboolean(
        'INPUTS', 'fill_soil_nodata_flag')

    # Use Ksat to calculate ssr2gw_rate and slowcoef_lin
    ksat_name = inputs_cfg.get('INPUTS', 'ksat_name')

    # Read and apply soil depth raster
    # Otherwise soil depth will only be derived from rooting depth
    try:
        soil_depth_flag = inputs_cfg.getboolean('INPUTS', 'soil_depth_flag')
    except ConfigParser.NoOptionError:
        soil_depth_flag = False
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'soil_depth_flag', soil_depth_flag))
    if soil_depth_flag:
        soil_depth_name = inputs_cfg.get('INPUTS', 'soil_depth_name')

    # Use geology based multipliers to adjust ssr2gw_rate
    # Otherwise default value set in config file will be used
    try:
        ssr2gw_mult_flag = inputs_cfg.getboolean('INPUTS', 'ssr2gw_mult_flag')
    except ConfigParser.NoOptionError:
        ssr2gw_mult_flag = False
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'ssr2gw_flag', ssr2gw_flag))
    if ssr2gw_mult_flag:
        ssr2gw_mult_name = inputs_cfg.get('INPUTS', 'ssr2gw_mult_name')

    # Check input paths
    if not arcpy.Exists(hru.polygon_path):
        logging.error(
            '\nERROR: Fishnet ({}) does not exist'.format(
                hru.polygon_path))
        sys.exit()

    # All of the soil rasters must exist
    awc_orig_path = os.path.join(soil_orig_ws, awc_name)
    clay_pct_orig_path = os.path.join(soil_orig_ws, clay_pct_name)
    sand_pct_orig_path = os.path.join(soil_orig_ws, sand_pct_name)
    ksat_orig_path = os.path.join(soil_orig_ws, ksat_name)
    if soil_depth_flag:
        soil_depth_orig_path = os.path.join(soil_orig_ws, soil_depth_name)
    if ssr2gw_mult_flag:
        ssr2gw_mult_orig_path = os.path.join(soil_orig_ws, ssr2gw_mult_name)

    # Check that either the original or projected/clipped raster exists
    if not arcpy.Exists(awc_orig_path):
        logging.error('\nERROR: AWC raster does not exist')
        sys.exit()
    if not arcpy.Exists(clay_pct_orig_path):
        logging.error('\nERROR: Clay raster does not exist')
        sys.exit()
    if not arcpy.Exists(sand_pct_orig_path):
        logging.error('\nERROR: Sand raster does not exist')
        sys.exit()
    if not arcpy.Exists(ksat_orig_path):
        logging.error('\nERROR: Ksat raster does not exist')
        sys.exit()
    if soil_depth_flag and not arcpy.Exists(soil_depth_orig_path):
        logging.error('\nERROR: Soil depth raster does not exist')
        sys.exit()
    if ssr2gw_mult_flag and not arcpy.Exists(ssr2gw_mult_orig_path):
        logging.error('\nERROR: Geology based raster for ssr2gw multiplier '
                      'does not exist')
        sys.exit()

    # Check other inputs
    if soil_cs <= 0:
        logging.error('\nERROR: soil cellsize must be greater than 0')
        sys.exit()
    soil_proj_method_list = ['BILINEAR', 'CUBIC', 'NEAREST']
    if soil_proj_method.upper() not in soil_proj_method_list:
        logging.error('\nERROR: Soil projection method must be: {}'.format(
            ', '.join(soil_proj_method_list)))
        sys.exit()

    # Build output folder if necessary
    soil_temp_ws = os.path.join(hru.param_ws, 'soil_rasters')
    if not os.path.isdir(soil_temp_ws):
        os.mkdir(soil_temp_ws)
    # Output paths
    awc_path = os.path.join(soil_temp_ws, 'awc.img')
    clay_pct_path = os.path.join(soil_temp_ws, 'clay_pct.img')
    sand_pct_path = os.path.join(soil_temp_ws, 'sand_pct.img')
    ksat_path = os.path.join(soil_temp_ws, 'ksat.img')
    soil_depth_path = os.path.join(soil_temp_ws, 'soil_depth.img')
    ssr2gw_mult_path = os.path.join(soil_temp_ws, 'ssr2gw_mult.img')

    # Set ArcGIS environment variables
    arcpy.CheckOutExtension('Spatial')
    env.overwriteOutput = True
    env.pyramid = 'PYRAMIDS -1'
    # env.pyramid = 'PYRAMIDS 0'
    env.workspace = soil_temp_ws
    env.scratchWorkspace = hru.scratch_ws

    # Available Water Capacity (AWC)
    logging.info('\nProjecting/clipping AWC raster')
    soil_orig_sr = arcpy.sa.Raster(awc_orig_path).spatialReference
    logging.debug('  AWC GCS:  {}'.format(
        soil_orig_sr.GCS.name))
    # Remove existing projected raster
    if arcpy.Exists(awc_path):
        arcpy.Delete_management(awc_path)
    # Set preferred transforms
    transform_str = support.transform_func(hru.sr, soil_orig_sr)
    logging.debug('  Transform: {}'.format(transform_str))
    logging.debug('  Projection method: NEAREST')
    # Project soil raster
    # DEADBEEF - Arc10.2 ProjectRaster does not honor extent
    support.project_raster_func(
        awc_orig_path, awc_path, hru.sr,
        soil_proj_method, soil_cs, transform_str,
        '{} {}'.format(hru.ref_x, hru.ref_y), soil_orig_sr, hru)
    # env.extent = hru.extent
    # arcpy.ProjectRaster_management(
    #    awc_orig_path, awc_path, hru.sr,
    #    soil_proj_method, soil_cs, transform_str,
    #    '{} {}'.format(hru.ref_x, hru.ref_y),
    #    soil_orig_sr)
    # arcpy.ClearEnvironment('extent')

    # Percent clay
    logging.info('Projecting/clipping clay raster')
    soil_orig_sr = arcpy.sa.Raster(clay_pct_orig_path).spatialReference
    logging.debug('  Clay GCS: {}'.format(
        soil_orig_sr.GCS.name))
    # Remove existing projected raster
    if arcpy.Exists(clay_pct_path):
        arcpy.Delete_management(clay_pct_path)
    # Set preferred transforms
    transform_str = support.transform_func(hru.sr, soil_orig_sr)
    logging.debug('  Transform: {}'.format(transform_str))
    logging.debug('  Projection method: NEAREST')
    # Project soil raster
    # DEADBEEF - Arc10.2 ProjectRaster does not extent
    support.project_raster_func(
        clay_pct_orig_path, clay_pct_path, hru.sr,
        soil_proj_method, soil_cs, transform_str,
        '{} {}'.format(hru.ref_x, hru.ref_y), soil_orig_sr, hru)
    # env.extent = hru.extent
    # arcpy.ProjectRaster_management(
    #    clay_pct_orig_path, clay_pct_path, hru.sr,
    #    soil_proj_method, soil_cs, transform_str,
    #    '{} {}'.format(hru.ref_x, hru.ref_y),
    #    soil_orig_sr)
    # arcpy.ClearEnvironment('extent')

    # Percent sand
    logging.info('Projecting/clipping sand raster')
    soil_orig_sr = arcpy.sa.Raster(sand_pct_orig_path).spatialReference
    logging.debug('  Sand GCS: {}'.format(
        soil_orig_sr.GCS.name))
    # Remove existing projected raster
    if arcpy.Exists(sand_pct_path):
        arcpy.Delete_management(sand_pct_path)
    # Set preferred transforms
    transform_str = support.transform_func(hru.sr, soil_orig_sr)
    logging.debug('  Transform: {}'.format(transform_str))
    logging.debug('  Projection method: NEAREST')
    # Project soil raster
    # DEADBEEF - Arc10.2 ProjectRaster does not honor extent
    support.project_raster_func(
        sand_pct_orig_path, sand_pct_path, hru.sr,
        soil_proj_method, soil_cs, transform_str,
        '{} {}'.format(hru.ref_x, hru.ref_y), soil_orig_sr, hru)
    # env.extent = hru.extent
    # arcpy.ProjectRaster_management(
    #    sand_pct_orig_path, sand_pct_path, hru.sr,
    #    soil_proj_method, soil_cs, transform_str,
    #    '{} {}'.format(hru.ref_x, hru.ref_y),
    #    soil_orig_sr)
    # arcpy.ClearEnvironment('extent')

    # Hydraulic conductivity
    logging.info('Projecting/clipping ksat raster')
    ksat_orig_sr = arcpy.sa.Raster(ksat_orig_path).spatialReference
    logging.debug('  Ksat GCS: {}'.format(
        soil_orig_sr.GCS.name))
    # Remove existing projected raster
    if arcpy.Exists(ksat_path):
        arcpy.Delete_management(ksat_path)
    # Set preferred transforms
    transform_str = support.transform_func(hru.sr, ksat_orig_sr)
    logging.debug('  Transform: {}'.format(transform_str))
    logging.debug('  Projection method: NEAREST')
    # Project ksat raster
    # DEADBEEF - Arc10.2 ProjectRaster does not honor extent
    support.project_raster_func(
        ksat_orig_path, ksat_path, hru.sr,
        soil_proj_method, soil_cs, transform_str,
        '{} {}'.format(hru.ref_x, hru.ref_y), soil_orig_sr, hru)
    # env.extent = hru.extent
    # arcpy.ProjectRaster_management(
    #    ksat_orig_path, ksat_path, hru.sr,
    #    soil_proj_method, soil_cs, transform_str,
    #    '{} {}'.format(hru.ref_x, hru.ref_y),
    #    soil_orig_sr)
    # arcpy.ClearEnvironment('extent')

    # Soil depth is only needed if clipping root depth
    if soil_depth_flag:
        logging.info('\nProjecting/clipping depth raster')
        soil_orig_sr = arcpy.sa.Raster(soil_depth_orig_path).spatialReference
        logging.debug('  Depth GCS: {}'.format(
            soil_orig_sr.GCS.name))
        # Remove existing projected raster
        if arcpy.Exists(soil_depth_path):
            arcpy.Delete_management(soil_depth_path)
        # Set preferred transforms
        transform_str = support.transform_func(hru.sr, soil_orig_sr)
        logging.debug('  Transform: {}'.format(transform_str))
        logging.debug('  Projection method: NEAREST')
        # Project soil raster
        # DEADBEEF - Arc10.2 ProjectRaster does not honor extent
        support.project_raster_func(
            soil_depth_orig_path, soil_depth_path, hru.sr,
            soil_proj_method, soil_cs, transform_str,
            '{} {}'.format(hru.ref_x, hru.ref_y), soil_orig_sr, hru)
        # env.extent = hru.extent
        # arcpy.ProjectRaster_management(
        #    soil_depth_orig_path, soil_depth_path, hru.sr,
        #    soil_proj_method, soil_cs, transform_str,
        #    '{} {}'.format(hru.ref_x, hru.ref_y),
        #    soil_orig_sr)
        # arcpy.ClearEnvironment('extent')

    # Geology based multiplier for gravity drainage (ssr2gw multiplier)
    if ssr2gw_mult_flag:
        logging.info('\nProjecting/clipping ssr2gw multiplier raster')
        soil_orig_sr = arcpy.sa.Raster(ssr2gw_mult_orig_path).spatialReference
        logging.debug('  Depth GCS: {}'.format(
            soil_orig_sr.GCS.name))
        # Remove existing projected raster
        if arcpy.Exists(ssr2gw_mult_path):
            arcpy.Delete_management(ssr2gw_mult_path)
        # Set preferred transforms
        transform_str = support.transform_func(hru.sr, soil_orig_sr)
        logging.debug('  Transform: {}'.format(transform_str))
        logging.debug('  Projection method: NEAREST')
        # Project soil raster
        # DEADBEEF - Arc10.2 ProjectRaster does not honor extent
        support.project_raster_func(
            ssr2gw_mult_orig_path, ssr2gw_mult_path, hru.sr,
            soil_proj_method, soil_cs, transform_str,
            '{} {}'.format(hru.ref_x, hru.ref_y), soil_orig_sr, hru)    

    # Fill soil nodata values using nibble
    if fill_soil_nodata_flag:
        logging.info('\nFilling soil nodata values using Nibble')
        soil_raster_list = [
            awc_path, clay_pct_path, sand_pct_path, ksat_path]
        if soil_depth_flag:
            soil_raster_list.append(soil_depth_path)
        for soil_raster_path in soil_raster_list:
            logging.info('  {}'.format(soil_raster_path))
            # DEADBEEF - Check if there is any nodata to be filled first?
            mask_obj = arcpy.sa.Int(1000 * arcpy.sa.SetNull(
                arcpy.sa.Raster(soil_raster_path) < 0,
                arcpy.sa.Raster(soil_raster_path)))
            input_obj = arcpy.sa.Con(arcpy.sa.IsNull(mask_obj), 0, mask_obj)
            nibble_obj = 0.001 * arcpy.sa.Nibble(input_obj, mask_obj, 'ALL_VALUES')
            nibble_obj.save(soil_raster_path)
            arcpy.BuildPyramids_management(soil_raster_path)


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Soil Raster Prep',
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

    # Prepare GSFLOW Soil Rasters
    soil_raster_prep(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
