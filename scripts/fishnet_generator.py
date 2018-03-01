#--------------------------------
# Name:         fishnet_generator.py
# Purpose:      GSFLOW fishnet generator
# Notes:        ArcGIS 10.2 Version
# Python:       2.7
#--------------------------------

import argparse
import ConfigParser
import datetime as dt
from decimal import Decimal
import logging
import os
import sys

import arcpy
from arcpy import env

import support_functions as support


def fishnet_func(config_path, overwrite_flag=False, debug_flag=False):
    """GSFLOW Fishnet Generator

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
    log_file_name = 'fishnet_generator_log.txt'
    log_console = logging.FileHandler(
        filename=os.path.join(hru.log_ws, log_file_name), mode='w')
    log_console.setLevel(logging.DEBUG)
    log_console.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger('').addHandler(log_console)
    logging.info('\nGSFLOW Fishnet Generator')

    # Warn the user if the fishnet already exists
    # It might be better to not allow the user to do this at all and force them
    # to manually remove the file.
    if arcpy.Exists(hru.polygon_path) and not overwrite_flag:
        logging.warning('\nWARNING: The existing fishnet/grid will be '
                        'over written\n  {}'.format(hru.polygon_path))
        raw_input('Press ENTER to continue')

    # Check input paths
    study_area_path = inputs_cfg.get('INPUTS', 'study_area_path')
    if not arcpy.Exists(study_area_path):
        logging.error(
            '\nERROR: Study area ({}) does not exist'.format(
                study_area_path))
        sys.exit()

    # For now, study area has to be a polygon
    if arcpy.Describe(study_area_path).datasetType != 'FeatureClass':
        logging.error(
            '\nERROR: For now, study area must be a polygon shapefile')
        sys.exit()

    # Read Fishnet specific parameters from INI
    # If ref_x and ref_y are not specified, get from the study area extent
    try:
        hru.ref_x = inputs_cfg.getfloat('INPUTS', 'hru_ref_x')
    except:
        hru.ref_x = arcpy.Describe(study_area_path).extent.XMin
        logging.info(
        '  {0} parameter not set in INI, setting {0} = {1}'.format(
            'ref_x', hru.ref_x))
    try:
        hru.ref_y = inputs_cfg.getfloat('INPUTS', 'hru_ref_y')
    except:
        hru.ref_y = arcpy.Describe(study_area_path).extent.YMin
        logging.info(
            '  {0} parameter not set in INI, setting {0} = {1}'.format(
                'ref_y', hru.ref_y))
    try:
        buffer_cells = inputs_cfg.getint('INPUTS', 'hru_buffer_cells')
    except:
        buffer_cells = 2
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'buffer_cells', buffer_cells))
    try:
        snap_method = inputs_cfg.get('INPUTS', 'hru_param_snap_method')
    except:
        snap_method = 'EXPAND'
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'snap_method', snap_method))
    snap_method_list = ['EXPAND', 'ROUND', 'SHRINK']
    if snap_method not in snap_method_list:
        logging.error('\nERROR: {} must be: {}'.format(
            'snap_method', ', '.join(snap_method_list)))
        sys.exit()

    # Log input hru parameters
    logging.info('\nFishnet Parameters')
    logging.info('  Cellsize:      {}'.format(hru.cs))
    logging.info('  Snap point:    {} {}'.format(hru.ref_x, hru.ref_y))
    logging.debug('  Buffer cells:  {}'.format(buffer_cells))

    # Read reference point as string for determining number of digits
    try:
        digits = abs(min(
            Decimal(inputs_cfg.get('INPUTS', 'hru_ref_x')).as_tuple().exponent,
            Decimal(inputs_cfg.get('INPUTS', 'hru_ref_y')).as_tuple().exponent))
    except ConfigParser.NoOptionError:
        digits = 10
    logging.debug('  Extent digits: {}'.format(digits))

    # Check inputs
    if buffer_cells < 0:
        logging.error('\nERROR: Buffer cells must be greater than 0')
        sys.exit()

    # Build output folder if necessary
    fishnet_temp_ws = os.path.join(hru.param_ws, 'fishnet_temp')
    if not os.path.isdir(fishnet_temp_ws):
        os.mkdir(fishnet_temp_ws)
    # Output paths
    study_area_proj_path = os.path.join(
        fishnet_temp_ws, 'projected_study_area.shp')

    # Set ArcGIS environment variables
    arcpy.CheckOutExtension('Spatial')
    env.overwriteOutput = True
    env.pyramid = 'PYRAMIDS -1'
    # env.pyramid = 'PYRAMIDS 0'
    env.workspace = hru.param_ws
    env.scratchWorkspace = hru.scratch_ws

    # Get spatial reference of study_area
    hru.sr = arcpy.Describe(study_area_path).spatialReference

    # If study area spat_ref doesn't match hru_param spat_ref
    # Project study are to hru_param and get projected extent
    # Otherwise, read study_area extent directly
    study_area_extent = arcpy.Describe(study_area_path).extent
    logging.debug('\n  Study area: {}'.format(study_area_path))
    logging.debug('  Study area spat. ref.: {}'.format(hru.sr.name))
    logging.debug('  Study area GCS:        {}'.format(hru.sr.GCS.name))
    logging.info('  Study Area extent: {}'.format(
        support.extent_string(study_area_extent)))

    # Check if the study area shapefile is projeted
    if (hru.sr.name in ['GCS_North_American_1983', 'GCS_WGS_1984'] or
            hru.sr.GCS.name == hru.sr.name):
        logging.warning(
            '\nWARNING: The study area shapefile does not appear to be projected.'
            '\n  This will likely cause problems or not work at all.'
            '\n  Projection: {}'.format(hru.sr.name))
        raw_input('Press ENTER to continue\n')

    # Buffer extent
    buffer_extent = support.buffer_extent_func(
        study_area_extent, buffer_cells * hru.cs)
    logging.info('  Buffered Extent:   {}'.format(
        support.extent_string(buffer_extent)))

    # Adjust study area extent to reference points
    # Set the number of digits of rounding based on the number digits
    #   int the reference points
    hru.ref_pnt = arcpy.Point(hru.ref_x, hru.ref_y)
    hru.extent = support.adjust_extent_to_snap(
        buffer_extent, hru.ref_pnt, hru.cs,
        method=snap_method, digits=digits)
    logging.info('  Snapped Extent:    {}'.format(
        support.extent_string(hru.extent)))

    # Build hru_param
    logging.info('\nBuilding HRU parameter fishnet')
    build_fishnet_func(
        hru.polygon_path, hru.point_path, hru.extent, hru.cs, hru.sr)

    # Write initial parameters to hru_param (X/Y, ROW/COL, Unique ID)
    # set_hru_id_func(hru.polygon_path, hru.extent, hru.cs)


def build_fishnet_func(hru_polygon_path, hru_point_path, extent, cs, sr):
    """"""
    # Remove existing
    if arcpy.Exists(hru_polygon_path):
        arcpy.Delete_management(hru_polygon_path)
    if arcpy.Exists(hru_point_path):
        arcpy.Delete_management(hru_point_path)
    # Calculate LL/UR corner points
    origin_pnt = (extent.XMin, extent.YMin)
    yaxis_pnt = (extent.XMin, extent.YMin + cs)
    corner_pnt = (extent.XMax, extent.YMax)
    origin_str = ' '.join(map(str, origin_pnt))
    yaxis_str = ' '.join(map(str, yaxis_pnt))
    corner_str = ' '.join(map(str, corner_pnt))
    logging.debug('  Origin: {}'.format(origin_str))
    logging.debug('  Y-Axis: {}'.format(yaxis_str))
    logging.debug('  Corner: {}'.format(corner_str))
    # Build fishnet & labels
    arcpy.CreateFishnet_management(
        hru_polygon_path, origin_str, yaxis_str, cs, cs,
        '0', '0', corner_str, 'LABELS', '#', 'POLYGON')
    arcpy.DefineProjection_management(hru_polygon_path, sr)
    arcpy.DefineProjection_management(hru_point_path, sr)


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Fishnet Generator',
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

    # Run GSFLOW Fishnet Generator
    fishnet_func(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
