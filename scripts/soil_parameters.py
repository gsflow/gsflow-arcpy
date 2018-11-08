#--------------------------------
# Name:         soil_parameters.py
# Purpose:      GSFLOW soil parameters
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


def soil_parameters(config_path):
    """Calculate GSFLOW Soil Parameters

    Parameters
    ----------
    config_path : str
        Project configuration file (.ini) path.

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
    log_file_name = 'soil_parameters_log.txt'
    log_console = logging.FileHandler(
        filename=os.path.join(hru.log_ws, log_file_name), mode='w')
    log_console.setLevel(logging.DEBUG)
    log_console.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger('').addHandler(log_console)
    logging.info('\nGSFLOW Soil Parameters')

    # Input parameters
    try:
        soil_pct_flag = inputs_cfg.getboolean('INPUTS', 'soil_pct_flag')
    except ConfigParser.NoOptionError:
        soil_pct_flag = True
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'soil_pct_flag', soil_pct_flag))
    try:
        moist_init_ratio = inputs_cfg.getfloat('INPUTS', 'moist_init_ratio')
    except ConfigParser.NoOptionError:
        moist_init_ratio = 0.1
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'moist_init_ratio', moist_init_ratio))
    try:
        rechr_init_ratio = inputs_cfg.getfloat('INPUTS', 'rechr_init_ratio')
    except ConfigParser.NoOptionError:
        rechr_init_ratio = 0.1
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'rechr_init_ratio', rechr_init_ratio))

    # Read and apply ssr2gw multiplier raster
    # Otherwise default value will be used
    try:
        ssr2gw_mult_flag = inputs_cfg.getboolean('INPUTS', 'ssr2gw_mult_flag')
    except ConfigParser.NoOptionError:
        ssr2gw_mult_flag = False
    try:
        ssr2gw_k_default = inputs_cfg.getfloat('INPUTS', 'ssr2gw_k_default')
    except ConfigParser.NoOptionError:
        ssr2gw_k_default = 0.001
        logging.info(
        '  Missing INI parameter, setting {} = {}'.format(
            'ssr2gw_k_default', ssr2gw_k_default))

    # Read and apply soil depth raster
    # Otherwise soil depth will only be derived from rooting depth
    try:
        soil_depth_flag = inputs_cfg.getboolean('INPUTS', 'soil_depth_flag')
    except ConfigParser.NoOptionError:
        soil_depth_flag = False
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'soil_depth_flag', soil_depth_flag))

    # Input folders
    soil_temp_ws = os.path.join(hru.param_ws, 'soil_rasters')
    if not os.path.isdir(soil_temp_ws):
        os.mkdir(soil_temp_ws)

    # Input paths
    awc_path = os.path.join(soil_temp_ws, 'awc.img')
    clay_pct_path = os.path.join(soil_temp_ws, 'clay_pct.img')
    sand_pct_path = os.path.join(soil_temp_ws, 'sand_pct.img')
    ksat_path = os.path.join(soil_temp_ws, 'ksat.img')
    soil_depth_path = os.path.join(soil_temp_ws, 'soil_depth.img')
    soil_root_max_path = os.path.join(soil_temp_ws, 'soil_root_max.img')
    ssr2gw_mult_path = os.path.join(soil_temp_ws, 'ssr2gw_mult.img')

    # Check input paths
    if not arcpy.Exists(hru.polygon_path):
        logging.error(
            '\nERROR: Fishnet ({}) does not exist'.format(
                hru.polygon_path))
        sys.exit()
    # All of the soil rasters must exist
    # Check that the projected/clipped/filled raster exists
    if not arcpy.Exists(awc_path):
        logging.error('\nERROR: AWC raster does not exist')
        sys.exit()
    if not arcpy.Exists(clay_pct_path):
        logging.error('\nERROR: Clay raster does not exist')
        sys.exit()
    if not arcpy.Exists(sand_pct_path):
        logging.error('\nERROR: Sand raster does not exist')
        sys.exit()
    if not arcpy.Exists(ksat_path):
        logging.error('\nERROR: Ksat raster does not exist')
        sys.exit()
    if soil_depth_flag and not arcpy.Exists(soil_depth_path):
        logging.error('\nERROR: Soil depth raster does not exist')
        sys.exit()
    if ssr2gw_mult_flag and not arcpy.Exists(ssr2gw_mult_path):
        logging.error('\nERROR: SSR2GW multiplier raster does not exist')
        sys.exit()
    # Check soil init ratios
    if moist_init_ratio < 0 or moist_init_ratio > 1:
        logging.error('\nERROR: Soil moist_init_ratio must be between 0 & 1')
        sys.exit()
    if rechr_init_ratio < 0 or rechr_init_ratio > 1:
        logging.error('\nERROR: Soil rechr_init_ratio must be between 0 & 1')
        sys.exit()

    # DEM Slope is needed for SSR2GW_RATE
    dem_temp_ws = os.path.join(hru.param_ws, 'dem_rasters')
    dem_slope_path = os.path.join(dem_temp_ws, 'dem_slope.img')
    if not os.path.isdir(dem_temp_ws):
        logging.error(
            '\nERROR: DEM temp folder does not exist\n'
            '\nERROR: Try re-running dem_2_stream.py')
        sys.exit()
    if not os.path.isfile(dem_slope_path):
        logging.error(
            '\nERROR: Slope raster does not exist\n'
            '\nERROR: Try re-running dem_2_stream.py')
        sys.exit()

    # Output paths
    # soil_type_path = os.path.join(soil_temp_ws, 'soil_type.img')
    moist_max_path = os.path.join(soil_temp_ws, 'soil_moist_max.img')
    rechr_max_path = os.path.join(soil_temp_ws, 'soil_rechr_max.img')

    # Root depth is calculated by veg script
    veg_temp_ws = os.path.join(hru.param_ws, 'veg_rasters')
    root_depth_path = os.path.join(veg_temp_ws, 'root_depth.img')
    if not arcpy.Exists(root_depth_path):
        logging.error(
            '\nERROR: Root depth raster does not exists'
            '\nERROR: Try re-running veg_parameters script\n')
        sys.exit()


    # Set ArcGIS environment variables
    arcpy.CheckOutExtension('Spatial')
    env.overwriteOutput = True
    env.pyramid = 'PYRAMIDS -1'
    # env.pyramid = 'PYRAMIDS 0'
    env.workspace = soil_temp_ws
    env.scratchWorkspace = hru.scratch_ws

    # Check field
    logging.info('\nAdding soil fields if necessary')
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
    support.add_field_func(hru.polygon_path, hru.slowcoef_lin_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.slowcoef_sq_field, 'DOUBLE')
    support.add_field_func(hru.polygon_path, hru.ssr2gw_k_field, 'DOUBLE')


    # Compute soil depth as max of root and soil depth
    if soil_depth_flag:
        logging.info('\nComputing max soil depth from root and soil depth')
        soil_depth_obj = arcpy.sa.Con(
           arcpy.sa.Raster(root_depth_path) > arcpy.sa.Raster(soil_depth_path),
           arcpy.sa.Raster(root_depth_path), arcpy.sa.Raster(soil_depth_path))
        soil_depth_obj.save(soil_root_max_path)
    else:
        soil_depth_obj = arcpy.sa.Raster(root_depth_path)

    # Calculate maximum soil moisture
    # logging.info('\nCalculating soil {}'.format(hru.moist_max_field))
    # moist_max_obj = arcpy.sa.Raster(awc_path) * soil_depth_obj
    # moist_max_obj.save(moist_max_path)
    # del moist_max_obj

    # # Calculate soil recharge zone maximum
    # logging.info('Calculating soil {}'.format(hru.rechr_max_field))
    # # Minimum of rooting depth and 18 (inches?)
    # rechr_max_obj = arcpy.sa.Float(
    #     arcpy.sa.Con(soil_depth_obj < 18, soil_depth_obj, 18))
    # rechr_max_obj *= arcpy.sa.Raster(awc_path)
    # rechr_max_obj.save(rechr_max_path)
    # del rechr_max_obj

    # # Read in slope raster and convert to radians
    # dem_slope_obj = math.pi * arcpy.sa.Raster(dem_slope_path) / 180
    # porosity_obj = 0.475
    #
    # # Gravity drainage to groundwater reservoir linear coefficient
    # logging.info('\nCalculating SSR2GW_RATE')
    # logging.info('  Assuming slope is in degrees')
    # logging.info('  Porosity is currently fixed at: {}'.format(
    #     porosity_obj))
    # ssr2gw_rate_obj = (
    #     arcpy.sa.Raster(ksat_path) * porosity_obj * (1 - dem_slope_obj))
    # ssr2gw_rate_obj.save(ssr2gw_rate_path)
    # del ssr2gw_rate_obj
    #
    # # Gravity drainage to groundwater reservoir linear coefficient
    # logging.info('\nCalculating SLOWCOEF_L')
    # logging.info('  Assuming slope is in degrees')
    # logging.info('  Porosity is currently fixed at: {}'.format(
    # slowcoef_lin_obj = (
    #     arcpy.sa.Raster(ksat_path) * math.sin(dem_slope_obj) /
    #     (porosity_obj * hru_length_obj))
    # slowcoef_lin_obj.save(slowcoef_lin_path)
    # del slowcoef_lin_obj, hru_length_obj
    # del dem_slope_obj, porosity_obj
    # # This block ^^ could be used to perform operations on a raster level if wanted


    # List of rasters, fields, and stats for zonal statistics
    zs_soil_dict = dict()
    zs_soil_dict[hru.awc_field] = [awc_path, 'MEAN']
    zs_soil_dict[hru.clay_pct_field] = [clay_pct_path, 'MEAN']
    zs_soil_dict[hru.sand_pct_field] = [sand_pct_path, 'MEAN']
    zs_soil_dict[hru.ksat_field] = [ksat_path, 'MEAN']
    if soil_depth_flag:
        zs_soil_dict[hru.soil_root_max_field] = [soil_root_max_path, 'MEAN']
    else:
        zs_soil_dict[hru.soil_root_max_field] = [root_depth_path, 'MEAN']
    if ssr2gw_mult_flag:
        zs_soil_dict[hru.ssr2gw_k_field] = [ssr2gw_mult_path, 'MEAN']
    # zs_soil_dict[hru.moist_max_field] = [moist_max_path, 'MEAN']
    # zs_soil_dict[hru.rechr_max_field] = [rechr_max_path, 'MEAN']

    # Calculate zonal statistics
    logging.info('\nCalculating zonal statistics')
    support.zonal_stats_func(
        zs_soil_dict, hru.polygon_path, hru.point_path, hru)


    # Make a fishnet layer for calculating fields
    hru_polygon_layer = "hru_polygon_layer"
    arcpy.MakeFeatureLayer_management(
        hru.polygon_path, hru_polygon_layer)

    # Calculate maximum soil moisture
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.moist_max_field,
        '!{}! * !{}!'.format(hru.soil_root_max_field, hru.awc_field),
        'PYTHON')

    # Calculate soil recharge zone maximum
    logging.info('Calculating soil {}'.format(hru.rechr_max_field))
    # Minimum of rooting depth and 18 (inches)
    rech_max_cb = (
        'def rech_max_func(soil_root_max, awc):\n'
        '    if soil_root_max > 18: return 18*awc\n'
        '    else: return soil_root_max*awc\n')
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.rechr_max_field,
        'rech_max_func(!{}!, !{}!)'.format(
            hru.soil_root_max_field, hru.awc_field),
        'PYTHON', rech_max_cb)

    # Calculate SOIL_TYPE
    logging.info('\nCalculating {}'.format(hru.soil_type_field))
    arcpy.SelectLayerByAttribute_management(
        hru_polygon_layer, "NEW_SELECTION",
        '"{}" = 1'.format(hru.type_field))
    if soil_pct_flag:
        soil_type_pct = (50, 40)
    else:
        soil_type_pct = (0.50, 0.40)
    soil_type_cb = (
        'def soil_type_func(clay, sand):\n'
        '    if sand > {}: return 1\n'
        '    elif clay > {}: return 3\n'
        '    else: return 2\n').format(*soil_type_pct)
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.soil_type_field,
        'soil_type_func(!{}!, !{}!)'.format(
            hru.clay_pct_field, hru.sand_pct_field),
        'PYTHON', soil_type_cb)
    arcpy.SelectLayerByAttribute_management(
        hru_polygon_layer, "SWITCH_SELECTION")
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.soil_type_field, '0', 'PYTHON')

    # Calculate SOIL_MOIST_INIT & SOIL_RECHR_INIT from max values
    arcpy.SelectLayerByAttribute_management(
        hru_polygon_layer, "NEW_SELECTION",
        '"{}" = 1 AND "{}" >= 0'.format(
            hru.type_field, hru.moist_max_field))
    logging.info('\nCalculating {0} as {2} * {1}'.format(
        hru.moist_init_field, hru.moist_max_field, moist_init_ratio))
    arcpy.CalculateField_management(
        hru.polygon_path, hru.moist_init_field,
        '!{}! * {}'.format(hru.moist_max_field, moist_init_ratio), 'PYTHON')
    arcpy.SelectLayerByAttribute_management(
        hru_polygon_layer, "SWITCH_SELECTION")
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.moist_init_field, '0', 'PYTHON')
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.moist_max_field, '0', 'PYTHON')

    # Calculate SOIL_MOIST_INIT & SOIL_RECHR_INIT from max values
    arcpy.SelectLayerByAttribute_management(
        hru_polygon_layer, "NEW_SELECTION",
        '"{}" = 1 AND "{}" >= 0'.format(
            hru.type_field, hru.rechr_max_field))
    logging.info('Calculating {0} as {2} * {1}'.format(
        hru.rechr_init_field, hru.rechr_max_field, moist_init_ratio))
    arcpy.CalculateField_management(
        hru.polygon_path, hru.rechr_init_field,
        '!{}! * {}'.format(hru.rechr_max_field, moist_init_ratio), 'PYTHON')
    arcpy.SelectLayerByAttribute_management(
        hru_polygon_layer, "SWITCH_SELECTION",
        '"{}" != 1'.format(hru.type_field))
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.rechr_init_field, '0', 'PYTHON')
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.rechr_max_field, '0', 'PYTHON')

    # Calculate SSR2G_KFAC from ssr2gw_mult raster
    arcpy.SelectLayerByAttribute_management(
        hru_polygon_layer, "NEW_SELECTION",
        '"{}" = 1 AND "{}" >= 0'.format(
            hru.type_field, hru.ssr2gw_k_field))
    logging.info('Using {1} to calculate {0}'.format(
        hru.ssr2gw_k_field, ssr2gw_mult_path))
    arcpy.CalculateField_management(
        hru.polygon_path, hru.ssr2gw_k_field,
        '!{}!'.format(hru.ssr2gw_k_field), 'PYTHON')

    # Fill SSR2G_K multiplier value if field not set
    logging.info ('ssr2gw_k_default = {}'.format(ssr2gw_k_default))
    if (all([row[0] == 0 for row in arcpy.da.SearchCursor(
            hru.polygon_path, [hru.ssr2gw_k_field])])):
        logging.info('Filling {} from default value in config file'.format(
            hru.ssr2gw_k_field))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.ssr2gw_k_field,
            ssr2gw_k_default, 'PYTHON')
    else:
        logging.info(
            '{} appears to already have been set and '
            'will not be overwritten'.format(hru.ssr2gw_k_field))

    # Calculating ssr2gw_rate
    # Gravity drainage to groundwater reservoir linear coefficient
    # Default value is 0.1 (range 0-1)
    # Convert Ksat from um/s to in/day
    # ssr2gw_rate = ks / sat_threshold
    # sat_threshold = moist_max * (sand% / 100)
    logging.info('\nCalculating {}'.format(hru.ssr2gw_rate_field))
    logging.info('  assuming {} is in units of um/s'.format(hru.ksat_field))
    # porosity_flt = 0.475
    ssr2gw_exp = 1
    logging.debug('  using eqn: ssr2gw_rate = ks/sat threshold')
    # logging.debug('  default values: porosity_flt = 0.475')
    logging.debug('  default values: ssr2gw_exp = 1')
    arcpy.SelectLayerByAttribute_management(
        hru_polygon_layer, "NEW_SELECTION",
        '"{}" = 1 AND "{}" > 0 AND "{}" > 0'.format(
            hru.type_field, hru.moist_max_field, hru.sand_pct_field))
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.ssr2gw_rate_field,
        '(!{}! * (3600 * 24 / (2.54 * 10000))) * !{}! / (!{}! * (!{}!/ 100))'.format(
            hru.ksat_field, hru.ssr2gw_k_field, hru.moist_max_field, hru.sand_pct_field),
        'PYTHON')
    arcpy.SelectLayerByAttribute_management(
        hru_polygon_layer, "SWITCH_SELECTION")
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.ssr2gw_rate_field, '0', 'PYTHON')

    # Calculating slowceof_lin
    # Default value is 0.015 (range 0-1)
    # Convert Ksat from um/s to m/day
    logging.info('Calculating {}'.format(hru.slowcoef_lin_field))
    logging.info('  {} must be in um/s'.format(hru.ksat_field))
    arcpy.SelectLayerByAttribute_management(
        hru_polygon_layer, "NEW_SELECTION",
        '"{}" = 1'.format(hru.type_field))
    slowcoef_lin_cb = (
        'def slowcoef_lin(ksat, slope, cs):\n'
        '    return 0.1 * ksat * 0.0864 * math.sin(slope) / cs\n')
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.slowcoef_lin_field,
        'slowcoef_lin(!{0}!, !{1}!, {2})'.format(
            hru.ksat_field, hru.dem_slope_rad_field, hru.cs),
        'PYTHON', slowcoef_lin_cb)
    arcpy.SelectLayerByAttribute_management(
        hru_polygon_layer, "SWITCH_SELECTION")
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.slowcoef_lin_field, '0', 'PYTHON')

    # Calculating slowceof_sq
    # Default value is 0.015 (range 0-1)
    # Convert Ksat from um/s to m/day
    logging.info('Calculating {}'.format(hru.slowcoef_sq_field))
    logging.info('  {} must be in um/s'.format(hru.ksat_field))
    arcpy.SelectLayerByAttribute_management(
        hru_polygon_layer, "NEW_SELECTION",
        '"{}" = 1 AND "{}" > 0 AND "{}" > 0'.format(
            hru.type_field, hru.moist_max_field, hru.sand_pct_field))
    slowcoef_sq_cb = (
        'def slowcoef_sq(ksat, slope, moist_max, sand, cs):\n'
        '    return 0.9 * (ksat * 0.0864 * math.sin(slope) / '
        '(moist_max * (sand / 100) * cs))\n')
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.slowcoef_sq_field,
        'slowcoef_sq(!{0}!, !{1}!, !{2}!, !{3}!, {4})'.format(
            hru.ksat_field, hru.dem_slope_rad_field,
            hru.moist_max_field, hru.sand_pct_field, hru.cs),
        'PYTHON', slowcoef_sq_cb)
    arcpy.SelectLayerByAttribute_management(
        hru_polygon_layer, "SWITCH_SELECTION")
    arcpy.CalculateField_management(
        hru_polygon_layer, hru.slowcoef_sq_field, '0', 'PYTHON')

    # Cleanup
    arcpy.Delete_management(hru_polygon_layer)
    del hru_polygon_layer

    #  Reset soils values for lake cells (HRU_TYPE == 2)
    #  Also reset for ocean cells (HRU_TYPE == 0 and DEM_ADJ == 0)
    #  This will remove many NoData cells ( == -999)
    # if True:
    #    logging.info('\nClearing lake soil parameters')
    #    # logging.info(
    #    #    '\nClearing soil parameters for lake and inactive ocean cells')
    #    hru_polygon_layer = "hru_polygon_layer"
    #    arcpy.MakeFeatureLayer_management(
    #        hru.polygon_path, hru_polygon_layer)
    #    arcpy.SelectLayerByAttribute_management(
    #        hru_polygon_layer, "NEW_SELECTION",
    #        '"{0}" = 2 OR ("{0}" = 0 AND "{1}" = 0)'.format(
    #            hru.type_field, hru.dem_adj_field))
    #        # '"{0}" = 2 AND "{1}" = -999'.format(
    #        #    hru.type_field, hru.awc_field))
    #    arcpy.CalculateField_management(
    #        hru_polygon_layer, hru.awc_field, 0, 'PYTHON')
    #    arcpy.CalculateField_management(
    #        hru_polygon_layer, hru.clay_pct_field, 0, 'PYTHON')
    #    arcpy.CalculateField_management(
    #        hru_polygon_layer, hru.sand_pct_field, 0, 'PYTHON')
    #    arcpy.CalculateField_management(
    #        hru_polygon_layer, hru.ksat_field, 0, 'PYTHON')
    #    # DEADBEEF
    #    # Soil type must be 1-3, can i set it to 0?
    #    logging.info('  Setting default lake soil type to 2')
    #    arcpy.CalculateField_management(
    #        hru_polygon_layer, hru.soil_type_field, 2, 'PYTHON')
    #    arcpy.CalculateField_management(
    #        hru_polygon_layer, hru.moist_init_field, 0, 'PYTHON')
    #    arcpy.CalculateField_management(
    #        hru_polygon_layer, hru.moist_max_field, 0, 'PYTHON')
    #    arcpy.CalculateField_management(
    #        hru_polygon_layer, hru.rechr_init_field, 0, 'PYTHON')
    #    arcpy.CalculateField_management(
    #        hru_polygon_layer, hru.rechr_max_field, 0, 'PYTHON')
    #    arcpy.CalculateField_management(
    #        hru_polygon_layer, hru.ssr2gw_rate_field, 0, 'PYTHON')
    #    arcpy.CalculateField_management(
    #        hru_polygon_layer, hru.slowcoef_lin_field, 0, 'PYTHON')
    #    arcpy.Delete_management(hru_polygon_layer)
    #    del hru_polygon_layer


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Soil Parameters',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i', '--ini', required=True,
        help='Project input file', metavar='PATH')
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

    soil_parameters(config_path=args.ini)
