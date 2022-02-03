#--------------------------------
# Name:         prism_800m_normals.py
# Purpose:      GSFLOW PRISM parameters from default 800m normals
# Notes:        ArcGIS 10.2 Version
# Python:       2.7
#--------------------------------

import argparse
import ConfigParser
import datetime as dt
import logging
import os
import re
import sys

import arcpy
from arcpy import env

import support_functions as support


def prism_800m_parameters(config_path, data_name='ALL'):
    """Calculate GSFLOW PRISM Parameters

    Parameters
    ----------
    config_path : str
        Project configuration file (.ini) path.
    data_name : {'PPT', 'TMAX', 'TMIN', 'ALL'}
        DAYMET data type (the default is 'PPT').

    Returns
    -------
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
    log_file_name = 'prism_800m_normals_log.txt'
    log_console = logging.FileHandler(
        filename=os.path.join(hru.log_ws, log_file_name), mode='w')
    log_console.setLevel(logging.DEBUG)
    log_console.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger('').addHandler(log_console)
    logging.info('\nGSFLOW PRISM Parameters')

    # PRISM
    prism_ws = inputs_cfg.get('INPUTS', 'prism_folder')
    prism_proj_method = inputs_cfg.get('INPUTS', 'prism_projection_method')
    prism_cs = inputs_cfg.getint('INPUTS', 'prism_cellsize')
    calc_jh_coef_flag = inputs_cfg.getboolean(
        'INPUTS', 'calc_prism_jh_coef_flag')

    if calc_jh_coef_flag:
        # DEADBEEF - This could/should be moved to support_functions.py since it is
        #   in this script and in both PRISM scripts.
        # DEM Units
        dem_units = inputs_cfg.get('INPUTS', 'dem_units').lower()
        dem_unit_types = {
            'meters': 'meter', 'm': 'meter', 'meter': 'meter',
            'feet': 'feet', 'ft': 'meter', 'foot': 'meter',}
        try:
            dem_units = dem_unit_types[dem_units]
        except:
            logging.error(
                '\nERROR: DEM unit "{}" is not supported\n'.format(dem_units))
            sys.exit()
        # Many expressions are hardcoded to units of feet
        # If dem_units are in meters, scale DEM_ADJ to get to feet
        if dem_units == 'meter':
            dem_unit_scalar = 0.3048
        else:
            dem_unit_scalar = 1.0

    # Check input paths
    if not arcpy.Exists(hru.polygon_path):
        logging.error(
            '\nERROR: Fishnet ({}) does not exist'.format(
                hru.polygon_path))
        sys.exit()
    # Check that PRISM folder is valid
    if not os.path.isdir(prism_ws):
        logging.error(
            '\nERROR: PRISM folder ({}) does not exist'.format(prism_ws))
        sys.exit()
    proj_method_list = ['BILINEAR', 'CUBIC', 'NEAREST']
    if prism_proj_method.upper() not in proj_method_list:
        logging.error('\nERROR: PRISM projection method must be: {}'.format(
            ', '.join(proj_method_list)))
        sys.exit()
    logging.debug('  Projection method:    {}'.format(
        prism_proj_method.upper()))

    # Check other inputs
    if prism_cs <= 0:
        logging.error('\nERROR: PRISM cellsize must be greater than 0\n')
        sys.exit()
    elif prism_cs > hru.cs:
        logging.warning(
            '\nWARNING: The "prism_cellsize" parameter should generally be '
            'set less than or equal \nto the fishnet cellsize.\n  '
            'Larger values may result in cells not having PRISM values.')
        raw_input('Press ENTER to continue')

    # Set ArcGIS environment variables
    arcpy.CheckOutExtension('Spatial')
    env.overwriteOutput = True
    env.pyramid = 'PYRAMIDS 0'
    env.workspace = hru.param_ws
    env.scratchWorkspace = hru.scratch_ws

    # PRISM data names
    if data_name == 'ALL':
        data_name_list = ['PPT', 'TMAX', 'TMIN']
    else:
        data_name_list = [data_name]

    # Set month list
    month_list = ['{:02d}'.format(m) for m in range(1, 13)]
    # month_list.extend(['annual'])

    # Check fields
    logging.info('\nAdding PRISM fields if necessary')
    for data_name in data_name_list:
        for month in month_list:
            support.add_field_func(
                hru.polygon_path, '{}_{}'.format(data_name, month), 'DOUBLE')

    # Process each PRISM data type
    logging.info('\nProjecting/clipping PRISM mean monthly rasters')
    for data_name in data_name_list:
        logging.info('\n{}'.format(data_name))
        prism_normal_re = re.compile(
            'PRISM_(?P<type>%s)_30yr_normal_800mM3_(?P<month>\d{2})_bil.bil$' % data_name,
            re.IGNORECASE)

        # Search all files & subfolders in prism folder
        #   for images that match data type
        input_raster_dict = dict()
        for root, dirs, files in os.walk(prism_ws):
            for file_name in files:
                prism_normal_match = prism_normal_re.match(file_name)
                if prism_normal_match:
                    month_str = prism_normal_match.group('month')
                    input_raster_dict[month_str] = os.path.join(root, file_name)
        if not input_raster_dict:
            logging.error(
                '\nERROR: No PRISM rasters were found matching the '
                'following pattern:\n  {}\n\nDouble check that the script '
                'and folder are for the same resolution '
                '(800m vs 4km)\n\n'.format(prism_normal_re.pattern))
            logging.error()
            sys.exit()

        # PRISM input data workspace
        # input_ws = os.path.join(prism_ws, data_name.lower())
        # if not os.path.isdir(input_ws):
        #    logging.error(
        #        '\nERROR: The PRISM {} folder does not exist'.format(
        #            data_name.lower()))
        #    sys.exit()

        # PRISM output data workspace
        output_ws = os.path.join(
            hru.param_ws, data_name.lower() + '_rasters')
        if not os.path.isdir(output_ws):
            os.mkdir(output_ws)

        # Remove all non year/month rasters in PRISM temp folder
        logging.info('  Removing existing PRISM files')
        for item in os.listdir(output_ws):
            # if prism_normal_re.match(item) and overwrite_flag:
            if prism_normal_re.match(item):
                os.remove(os.path.join(output_ws, item))

        # Extract, project/resample, clip
        # Process images by month
        zs_prism_dict = dict()
        # env.extent = hru.extent
        for month in month_list:
            logging.info('  Month: {}'.format(month))

            # Projected/clipped PRISM raster
            input_raster = input_raster_dict[month]
            # input_name = 'PRISM_{}_30yr_normal_800mM2_{1}_bil.bil'.format(
            #    data_name.lower(), input_month)
            # input_raster = os.path.join(input_ws, input_name)
            output_name = 'PRISM_{}_30yr_normal_800mM3_{}.img'.format(
                data_name.lower(), month)
            output_raster = os.path.join(output_ws, output_name)

            # Set preferred transforms
            input_sr = arcpy.sa.Raster(input_raster).spatialReference
            transform_str = support.transform_func(hru.sr, input_sr)
            if transform_str:
                logging.debug('  Transform: {}'.format(transform_str))

            # Project PRISM rasters to HRU coordinate system
            # DEADBEEF - Arc10.2 ProjectRaster does not extent
            support.project_raster_func(
                input_raster, output_raster, hru.sr,
                prism_proj_method.upper(), prism_cs, transform_str,
                '{} {}'.format(hru.ref_x, hru.ref_y), input_sr, hru)
            # arcpy.ProjectRaster_management(
            #    input_raster, output_raster, hru.sr,
            #    prism_proj_method.upper(), prism_cs, transform_str,
            #    '{} {}'.format(hru.ref_x, hru.ref_y),
            #    input_sr)

            # Save parameters for calculating zonal stats
            zs_field = '{}_{}'.format(data_name, month)
            zs_prism_dict[zs_field] = [output_raster, 'MEAN']

            # Cleanup
            del input_raster, output_raster, output_name
            del input_sr, transform_str, zs_field

        # Cleanup
        # arcpy.ClearEnvironment('extent')

        # Calculate zonal statistics
        logging.info('\nCalculating PRISM zonal statistics')
        support.zonal_stats_func(
            zs_prism_dict, hru.polygon_path, hru.point_path, hru)
        del zs_prism_dict

    # Jensen-Haise Potential ET air temperature coefficient
    # Update Jensen-Haise PET estimate using PRISM air temperature
    # DEADBEEF - First need to figure out month with highest Tmax
    #            Then get Tmin for same month
    if calc_jh_coef_flag:
        logging.info('\nRe-Calculating JH_COEF_HRU')
        logging.info('  Using PRISM temperature values')
        tmax_field_list = ['!TMAX_{:02d}!'.format(m) for m in range(1, 13)]
        tmin_field_list = ['!TMIN_{:02d}!'.format(m) for m in range(1, 13)]
        tmax_expr = 'max([{}])'.format(','.join(tmax_field_list))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.jh_tmax_field, tmax_expr, 'PYTHON')
        # Sort TMAX and get TMIN for same month
        tmin_expr = 'max(zip([{}],[{}]))[1]'.format(
            ','.join(tmax_field_list), ','.join(tmin_field_list))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.jh_tmin_field, tmin_expr, 'PYTHON')

        # Pass unit scalar to convert DEM_ADJ to feet if necessary
        support.jensen_haise_func(
            hru.polygon_path, hru.jh_coef_field, hru.dem_adj_field,
            hru.jh_tmin_field, hru.jh_tmax_field, dem_unit_scalar)


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='PRISM 800m Normals',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i', '--ini', required=True,
        help='Project input file', metavar='PATH')
    parser.add_argument(
        '-t', '--type', default='ALL', choices=['TMAX', 'TMIN', 'PPT', 'ALL'],
        help='PRISM Data Type')
    parser.add_argument(
        '-d', '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action="store_const", dest="loglevel")
    args = parser.parse_args()

    # Convert relative paths to absolute paths
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

    prism_800m_parameters(config_path=args.ini, data_name=args.type)
