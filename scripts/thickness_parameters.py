#--------------------------------
# Name:         thickness_parameters.py
# Purpose:      GSFLOW thickness parameters
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


def thickness_parameters(config_path, overwrite_flag=False, debug_flag=False):
    """Calculate GSFLOW Thickness Parameters

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
    log_file_name = 'thickness_parameters_log.txt'
    log_console = logging.FileHandler(
        filename=os.path.join(hru.log_ws, log_file_name), mode='w')
    log_console.setLevel(logging.DEBUG)
    log_console.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger('').addHandler(log_console)
    logging.info('\nGSFLOW Thickness Parameters')

    # Input folders
    # thickness_ws = os.path.join(parameter_ws, 'thickness_rasters')
    # if not os.path.isdir(thickness_ws):
    #    os.mkdir(thickness_ws)

    # Check input paths
    if not arcpy.Exists(hru.polygon_path):
        logging.error(
            '\nERROR: Fishnet ({}) does not exist'.format(
                hru.polygon_path))
        sys.exit()

    # Build output folders if necesssary
    thickness_temp_ws = os.path.join(hru.param_ws, 'thickness_rasters')
    if not os.path.isdir(thickness_temp_ws):
        os.mkdir(thickness_temp_ws)

    # Set ArcGIS environment variables
    arcpy.CheckOutExtension('Spatial')
    env.overwriteOutput = True
    env.pyramid = 'PYRAMIDS -1'
    # env.pyramid = 'PYRAMIDS 0'
    env.workspace = thickness_temp_ws
    env.scratchWorkspace = hru.scratch_ws

    # Check field
    logging.info('\nAdding thickness fields if necessary')
    support.add_field_func(hru.polygon_path, hru.alluv_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.alluv_thick_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.lay1_thick_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.lay2_thick_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.lay3_thick_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.lay4_thick_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.lay1_bottom_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.lay2_bottom_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.lay3_bottom_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.lay4_bottom_field, 'DOUBLE')

    # Make a fishnet layer for calculating fields
    hru_polygon_layer = "hru_polygon_layer"
    arcpy.MakeFeatureLayer_management(
        hru.polygon_path, hru_polygon_layer)
    arcpy.SelectLayerByAttribute_management(
        hru_polygon_layer, "NEW_SELECTION",
        '"{}" >= 0 '.format(hru.type_field))

    # Calculate layer bottom values
    logging.info('Calculating {}'.format(hru.lay1_bottom_field))
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.lay1_bottom_field,
        '!{}! - !{}!'.format(hru.dem_adj_field, hru.lay1_thick_field),
        'PYTHON')

    logging.info('Calculating {}'.format(hru.lay2_bottom_field))
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.lay2_bottom_field,
        '!{}! - !{}!'.format(hru.lay1_bottom_field, hru.lay2_thick_field),
        'PYTHON')

    logging.info('Calculating {}'.format(hru.lay3_bottom_field))
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.lay3_bottom_field,
        '!{}! - !{}!'.format(hru.lay2_bottom_field, hru.lay3_thick_field),
        'PYTHON')

    logging.info('Calculating {}'.format(hru.lay4_bottom_field))
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.lay4_bottom_field,
        '!{}! - !{}!'.format(hru.lay3_bottom_field, hru.lay4_thick_field),
        'PYTHON')

    # Cleanup
    arcpy.SelectLayerByAttribute_management(
        hru_polygon_layer, "CLEAR_SELECTION")
    arcpy.Delete_management(hru_polygon_layer)
    del hru_polygon_layer


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Alluvial Thickness Parameters',
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

    # Calculate GSFLOW Thickness Parameters
    thickness_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
