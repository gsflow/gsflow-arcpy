#--------------------------------
# Name:         ppt_ratio_parameters.py
# Purpose:      GSFLOW precipitation ratio parameters
# Notes:        ArcGIS 10.2+ Version
# Python:       2.7
#--------------------------------

import argparse
from collections import defaultdict
import ConfigParser
import datetime as dt
import logging
import os
import sys

import arcpy
from arcpy import env

import support_functions as support


def ppt_ratio_parameters(config_path, overwrite_flag=False, debug_flag=False):
    """Calculate GSFLOW Precipitation Ratio Parameters

    Args:
        config_file (str): Project config file path
        ovewrite_flag (bool): if True, overwrite existing files
        debug_flag (bool): if True, enable debug level logging

    Returns:
        None
    """

    # Hardcoded HRU field formats for now
    ppt_field_format = 'PPT_{:02d}'
    ratio_field_format = 'PPT_RT_{:02d}'

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
    log_file_name = 'ppt_ratio_parameters_log.txt'
    log_console = logging.FileHandler(
        filename=os.path.join(hru.log_ws, log_file_name), mode='w')
    log_console.setLevel(logging.DEBUG)
    log_console.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger('').addHandler(log_console)
    logging.info('\nGSFLOW Precipitation Ratio Parameters')

    # Units
    ppt_obs_units = support.get_param(
        'ppt_obs_units', 'mm', inputs_cfg).lower()
    ppt_units_list = ['mm', 'cm', 'm', 'in', 'ft']
    # Compare against the lower case of the values in the list
    #   but don't modify the acceptable units list
    if ppt_obs_units not in ppt_units_list:
        logging.error(
            '\nERROR: Invalid observed precipitation units ({})\n  '
            'Valid units are: {}'.format(
                ppt_obs_units, ', '.join(ppt_units_list)))
        sys.exit()

    # Convert units while reading obs values
    if ppt_obs_units == 'mm':
        units_factor = 1
    elif ppt_obs_units == 'cm':
        units_factor = 10
    elif ppt_obs_units == 'm':
        units_factor = 1000
    elif ppt_obs_units == 'in':
        units_factor = 25.4
    elif ppt_obs_units == 'ft':
        units_factor = 304.8
    else:
        units_factor = 1

    # Check input paths
    if not arcpy.Exists(hru.polygon_path):
        logging.error(
            '\nERROR: Fishnet ({}) does not exist'.format(
                hru.polygon_path))
        sys.exit()

    # PPT Zones
    set_ppt_zones_flag = inputs_cfg.getboolean('INPUTS', 'set_ppt_zones_flag')
    if set_ppt_zones_flag:
        ppt_zone_orig_path = inputs_cfg.get('INPUTS', 'ppt_zone_path')
        try:
            ppt_zone_field = inputs_cfg.get('INPUTS', 'ppt_zone_field')
        except:
            logging.error(
                '\nERROR: ppt_zone_field must be set in INI to apply '
                'zone specific ppt ratios\n')
            sys.exit()
        try:
            ppt_hru_id_field = inputs_cfg.get('INPUTS', 'ppt_hru_id_field')
        except:
            ppt_hru_id_field = None
            logging.warning(
                '  ppt_hru_id_field was not set in the INI file\n'
                '  PPT ratios will not be adjusted to match station '
                'values'.format(ppt_zone_field, hru.ppt_zone_id_field))

        # Field name for PSTA hard coded, but could be changed to be read from
        # config file like ppt_zone
        hru_psta_field = 'HRU_PSTA'

        try:
            ppt_obs_field_format = inputs_cfg.get(
                'INPUTS', 'ppt_obs_field_format')
        except:
            ppt_obs_field_format = 'PPT_{:02d}'
            logging.info('  Defaulting ppt_obs_field_format = {}'.format(
                ppt_obs_field_format))

        if not arcpy.Exists(ppt_zone_orig_path):
            logging.error(
                '\nERROR: PPT Zone ({}) does not exist'.format(
                    ppt_zone_orig_path))
            sys.exit()
        # ppt_zone_path must be a polygon shapefile
        if arcpy.Describe(ppt_zone_orig_path).datasetType != 'FeatureClass':
            logging.error(
                '\nERROR: ppt_zone_path must be a polygon shapefile')
            sys.exit()

        # Check ppt_zone_field
        if ppt_zone_field.upper() in ['FID', 'OID']:
            ppt_zone_field = arcpy.Describe(ppt_zone_orig_path).OIDFieldName
            logging.warning(
                '\n  NOTE: Using {} to set {}\n'.format(
                    ppt_zone_field, hru.ppt_zone_id_field))
        elif not arcpy.ListFields(ppt_zone_orig_path, ppt_zone_field):
            logging.error(
                '\nERROR: ppt_zone_field field {} does not exist\n'.format(
                    ppt_zone_field))
            sys.exit()
        # Need to check that field is an int type
        # Only check active cells (HRU_TYPE >0)?!
        elif not [f.type for f in arcpy.Describe(ppt_zone_orig_path).fields
                  if (f.name == ppt_zone_field and
                      f.type in ['SmallInteger', 'Integer'])]:
            logging.error(
                '\nERROR: ppt_zone_field field {} must be an integer type\n'.format(
                    ppt_zone_field))
            sys.exit()
        # Need to check that field values are all positive
        # Only check active cells (HRU_TYPE >0)?!
        elif min([row[0] for row in arcpy.da.SearchCursor(
                ppt_zone_orig_path, [ppt_zone_field])]) <= 0:
            logging.error(
                '\nERROR: ppt_zone_field values must be positive\n'.format(
                    ppt_zone_field))
            sys.exit()

        # Check hru_psta_field
        if not arcpy.ListFields(ppt_zone_orig_path, hru_psta_field):
            logging.error(
                '\nERROR: hru_psta_field field {} does not exist\n'.format(
                    hru_psta_field))
            sys.exit()
        # Need to check that field is an int type
        # Should we only check active cells (HRU_TYPE > 0)?
        elif not [f.type for f in arcpy.Describe(ppt_zone_orig_path).fields
                  if (f.name == hru_psta_field and
                      f.type in ['SmallInteger', 'Integer'])]:
            logging.error(
                '\nERROR: hru_psta_field field {} must be an integer type\n'.format(
                    hru_psta_field))
            sys.exit()
        # Need to check that field values are all positive
        # Should we only check active cells (HRU_TYPE > 0)?
        elif min([row[0] for row in arcpy.da.SearchCursor(
                ppt_zone_orig_path, [hru_psta_field])]) <= 0:
            logging.error(
                '\nERROR: hru_psta_field values must be positive\n'.format(
                    hru_psta_field))
            sys.exit()

        # Check ppt_hru_id_field
        # ppt_hru_id values are checked later
        if ppt_hru_id_field is not None:
            if not arcpy.ListFields(ppt_zone_orig_path, ppt_hru_id_field):
                logging.error(
                    '\nERROR: ppt_hru_id_field field {} does not exist\n'.format(
                        ppt_hru_id_field))
                sys.exit()
            # Need to check that field is an int type
            elif not [f.type for f in arcpy.Describe(ppt_zone_orig_path).fields
                      if (f.name == ppt_hru_id_field and
                          f.type in ['SmallInteger', 'Integer'])]:
                logging.error(
                    '\nERROR: ppt_hru_id_field field {} must be an integer type\n'.format(
                        ppt_hru_id_field))
                sys.exit()
            # Need to check that field values are not negative (0 is okay)
            elif min([row[0] for row in arcpy.da.SearchCursor(
                    ppt_zone_orig_path, [ppt_hru_id_field])]) < 0:
                logging.error(
                    '\nERROR: ppt_hru_id_field values cannot be negative\n'.format(
                        ppt_hru_id_field))
                sys.exit()
    else:
        # If a zone shapefile is not used, PPT must be set manually
        ppt_obs_list = inputs_cfg.get('INPUTS', 'ppt_obs_list')

        # Check that values are floats
        try:
            ppt_obs_list = map(float, ppt_obs_list.split(','))
        except ValueError:
            logging.error(
                '\nERROR: ppt_obs_list (mean monthly precipitation) '
                'values could not be parsed as floats\n')
            sys.exit()

        # Check that there are 12 values
        if len(ppt_obs_list) != 12:
            logging.error(
                '\nERROR: There must be exactly 12 mean monthly '
                'observed precipitation values based to ppt_obs_list\n')
            sys.exit()
        logging.info(
            '  Observed Mean Monthly PPT ({}):\n    {}\n'
            '    (Script will assume these are listed in month order, '
            'i.e. Jan, Feb, ...)'.format(
                ppt_obs_units, ', '.join(map(str, ppt_obs_list))))

        # Check if all the values are 0
        if ppt_obs_list == ([0.0] * 12):
            logging.error(
                '\nERROR: The observed precipitation values are all 0.\n'
                '  To compute PPT ratios, please set the ppt_obs_list '
                'parameter in the INI with\n  observed mean monthly PPT '
                'values (i.e. from a weather station)')
            sys.exit()

        # Get the PPT HRU ID
        try:
            ppt_hru_id = inputs_cfg.getint('INPUTS', 'ppt_hru_id')
        except:
            ppt_hru_id = 0

        # Check that the ppt_hru_id is a valid cell hru_id
        # If ppt_hru_id is 0, PPT ratios will not be adjusted
        if ppt_hru_id > 0:
            # Check that HRU_ID is valid
            logging.info('    PPT HRU_ID: {}'.format(ppt_hru_id))
            arcpy.MakeTableView_management(
                hru.polygon_path, "layer",
                "{} = {}".format(hru.id_field, ppt_hru_id))
            if (ppt_hru_id != 0 and
                    int(arcpy.GetCount_management("layer").getOutput(0)) == 0):
                logging.error(
                    '\nERROR: ppt_hru_id {0} is not a valid cell hru_id'
                    '\nERROR: ppt_ratios will NOT be forced to 1'
                    ' at cell {0}\n'.format(ppt_hru_id))
                ppt_hru_id = 0
            arcpy.Delete_management("layer")
        else:
            logging.info(
                '  PPT ratios will not be adjusted to match station values\n'
                '    (ppt_hru_id = 0)')

        # Could add a second check that HRU_PSTA has values >0

    # Build output folders if necessary
    ppt_ratio_temp_ws = os.path.join(hru.param_ws, 'ppt_ratio_temp')
    if not os.path.isdir(ppt_ratio_temp_ws):
        os.mkdir(ppt_ratio_temp_ws)
    ppt_zone_path = os.path.join(ppt_ratio_temp_ws, 'ppt_zone.shp')
    # ppt_zone_clip_path = os.path.join(ppt_ratio_temp_ws, 'ppt_zone_clip.shp')


    # Set ArcGIS environment variables
    arcpy.CheckOutExtension('Spatial')
    env.overwriteOutput = True
    # env.pyramid = 'PYRAMIDS -1'
    env.pyramid = 'PYRAMIDS 0'
    env.workspace = hru.param_ws
    env.scratchWorkspace = hru.scratch_ws

    # Set month list based on flags
    month_list = range(1, 13)
    ppt_field_list = [ppt_field_format.format(m) for m in month_list]
    ratio_field_list = [ratio_field_format.format(m) for m in month_list]

    # Check fields
    logging.info('\nAdding PPT ratio fields if necessary')
    # PPT zone fields
    support.add_field_func(
        hru.polygon_path, hru.ppt_zone_id_field, 'LONG')
    # PPT ratio fields
    for ratio_field in ratio_field_list:
        support.add_field_func(hru.polygon_path, ratio_field, 'DOUBLE')

    # Calculate PPT zone ID
    if set_ppt_zones_flag:
        logging.info('\nCalculating cell HRU Precipitation Zone ID')
        ppt_zone_desc = arcpy.Describe(ppt_zone_orig_path)
        ppt_zone_sr = ppt_zone_desc.spatialReference
        logging.debug('  Zones:      {}'.format(ppt_zone_orig_path))
        logging.debug('  Projection: {}'.format(ppt_zone_sr.name))
        logging.debug('  GCS:        {}'.format(ppt_zone_sr.GCS.name))

        # Reset PPT_ZONE_ID
        # if set_ppt_zones_flag:
        logging.info('  Resetting {} to 0'.format(hru.ppt_zone_id_field))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.ppt_zone_id_field, 0, 'PYTHON')

        # If ppt_zone spat_ref doesn't match hru_param spat_ref
        # Project ppt_zone to hru_param spat ref
        # Otherwise, read ppt_zone directly
        if hru.sr.name != ppt_zone_sr.name:
            logging.info('  Projecting precipitation zones...')
            # Set preferred transforms
            transform_str = support.transform_func(
                hru.sr, ppt_zone_sr)
            logging.debug('    Transform: {}'.format(transform_str))
            # Project ppt_zone shapefile
            arcpy.Project_management(
                ppt_zone_orig_path, ppt_zone_path, hru.sr,
                transform_str, ppt_zone_sr)
            del transform_str
        else:
            arcpy.Copy_management(ppt_zone_orig_path, ppt_zone_path)

        # # Remove all unnecessary fields
        # for field in arcpy.ListFields(ppt_zone_path):
        #     skip_field_list = ppt_obs_field_list + [ppt_zone_field, 'Shape']
        #     if field.name not in skip_field_list:
        #         try:
        #             arcpy.DeleteField_management(ppt_zone_path, field.name)
        #         except:
        #             pass

        # Set ppt zone ID
        logging.info('  Setting {}'.format(hru.ppt_zone_id_field))
        support.zone_by_centroid_func(
            ppt_zone_path, hru.ppt_zone_id_field, ppt_zone_field,
            hru.polygon_path, hru.point_path, hru)
        # support.zone_by_area_func(
        #    ppt_zone_layer, hru.ppt_zone_id_field, ppt_zone_field,
        #    hru.polygon_path, hru, hru_area_field, None, 50)

        # Set HRU_PSTA
        logging.info('  Setting {}'.format(hru.hru_psta_field))
        support.zone_by_centroid_func(
            ppt_zone_path, hru.hru_psta_field, hru_psta_field,
            hru.polygon_path, hru.point_path, hru)

        del ppt_zone_desc, ppt_zone_sr
    else:
        # Set all cells to zone 1
        arcpy.CalculateField_management(
            hru.polygon_path, hru.ppt_zone_id_field, 1, 'PYTHON')

    # Calculate ratios
    logging.info('\nCalculating mean monthly PPT ratios')
    if set_ppt_zones_flag:
        # Read mean monthly values for each zone
        ppt_obs_dict = dict()
        ppt_obs_field_list = [
            ppt_obs_field_format.format(m) for m in month_list]
        fields = [ppt_zone_field] + ppt_obs_field_list
        logging.debug('  Obs. Fields: {}'.format(', '.join(fields)))

        with arcpy.da.SearchCursor(ppt_zone_path, fields) as s_cursor:
            for row in s_cursor:
                ppt_obs_dict[int(row[0])] = map(float, row[1:13])

        # Convert values to mm if necessary to match PRISM
        if units_factor != 1:
            ppt_obs_dict = {z: p * units_factor for z, p in ppt_obs_dict}

        ppt_zone_list = sorted(ppt_obs_dict.keys())
        logging.debug('  PPT Zones: {}'.format(ppt_zone_list))

        # Print the observed PPT values
        logging.debug('  Observed PPT')
        for zone, ppt_obs in ppt_obs_dict.items():
            logging.debug('    {}: {}'.format(
                zone, ', '.join(['{:.2f}'.format(x) for x in ppt_obs])))

        # Default all zones to a ratio of 1
        ppt_ratio_dict = {z: [1] * 12 for z in ppt_zone_list}

        # Get list of HRU_IDs for each zone
        fields = [hru.ppt_zone_id_field, hru.id_field]
        zone_hru_id_dict = defaultdict(list)
        with arcpy.da.SearchCursor(hru.polygon_path, fields) as s_cursor:
            for row in s_cursor:
                zone_hru_id_dict[int(row[0])].append(int(row[1]))

        # Check that PPT_HRU_IDs are in the correct zone
        # Default all PPT Zone HRU IDs to 0
        ppt_hru_id_dict = {z: 0 for z in ppt_zone_list}
        if ppt_hru_id_field is not None:
            fields = [ppt_zone_field, ppt_hru_id_field]
            logging.debug('  PPT Zone ID field: {}'.format(ppt_zone_field))
            logging.debug('  PPT HRU ID field: {}'.format(ppt_hru_id_field))
            with arcpy.da.SearchCursor(ppt_zone_path, fields) as s_cursor:
                for row in s_cursor:
                    ppt_zone = int(row[0])
                    hru_id = int(row[1])
                    if hru_id == 0 or hru_id in zone_hru_id_dict[ppt_zone]:
                        ppt_hru_id_dict[ppt_zone] = hru_id
                        logging.debug('    {}: {}'.format(ppt_zone, hru_id))
                    else:
                        logging.error(
                            '\nERROR: HRU_ID {} is not in PPT ZONE {}'.format(
                                hru_id, ppt_hru_id_dict[ppt_zone]))
                        sys.exit()

        # Get gridded PPT values for each PPT_HRU_ID
        fields = [hru.ppt_zone_id_field, hru.id_field] + ppt_field_list
        # ppt_ratio_dict = dict()
        with arcpy.da.SearchCursor(hru.polygon_path, fields) as s_cursor:
            for row in s_cursor:
                ppt_zone = int(row[0])
                hru_id = int(row[1])
                if hru_id == 0:
                    pass
                elif hru_id in ppt_hru_id_dict.values():
                    ppt_gridded_list = map(float, row[2:14])
                    ppt_obs_list = ppt_obs_dict[ppt_zone]
                    ppt_ratio_list = [
                        float(o) / p if p > 0 else 0
                        for o, p in zip(ppt_obs_list, ppt_gridded_list)]
                    ppt_ratio_dict[ppt_zone] = ppt_ratio_list
        del ppt_hru_id_dict, zone_hru_id_dict, fields

        logging.debug('  PPT Ratio Adjustment Factors:')
        for k, v in ppt_ratio_dict.items():
            logging.debug('    {}: {}'.format(
                k, ', '.join(['{:.3f}'.format(x) for x in v])))

        # DEADBEEF - ZONE_VALUE is calculated in zone_by_centroid_func
        # There is probably a cleaner way of linking these two
        fields = [hru.ppt_zone_id_field] + ppt_field_list + ratio_field_list
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as u_cursor:
            for row in u_cursor:
                ppt_zone = int(row[0])
                for i, month in enumerate(month_list):
                    ppt_i = fields.index(ppt_field_format.format(month))
                    ratio_i = fields.index(ratio_field_format.format(month))

                    if ppt_zone in ppt_zone_list:
                        ppt_obs = ppt_obs_dict[ppt_zone][i]
                    else:
                        ppt_obs = 0

                    if ppt_obs > 0:
                        row[ratio_i] = (
                            ppt_ratio_dict[ppt_zone][i] * row[ppt_i] / ppt_obs)
                    else:
                        row[ratio_i] = 0
                u_cursor.updateRow(row)
            del row
    else:
        # Get gridded precip at PPT_HRU_ID
        fields = [hru.id_field] + ppt_field_list
        logging.debug('  Fields: {}'.format(', '.join(ppt_field_list)))

        # Convert values to mm if necessary to match PRISM
        if units_factor != 1:
            ppt_obs_list = [p * units_factor for p in ppt_obs_list]
            logging.debug(
                '\nConverted Mean Monthly PPT ({}):\n  {}'.format(
                    ppt_obs_units, ', '.join(map(str, ppt_obs_list))))

        # Scale all ratios so gridded PPT will match observed PPT at target cell
        if ppt_hru_id != 0:
            ppt_gridded_list = map(float, arcpy.da.SearchCursor(
                hru.polygon_path, fields,
                '"{}" = {}'.format(hru.id_field, ppt_hru_id)).next()[1:])
            logging.info('  Gridded PPT: {}'.format(
                ', '.join(['{:.2f}'.format(p) for p in ppt_gridded_list])))

            # Ratio of MEASURED or OBSERVED PPT to GRIDDED PPT
            # This will be multiplied by GRIDDED/OBSERVED below
            ppt_ratio_list = [
                float(o) / p if p > 0 else 0
                for o, p in zip(ppt_obs_list, ppt_gridded_list)]
            logging.info('  Obs./Gridded: {}'.format(
                ', '.join(['{:.3f}'.format(p) for p in ppt_ratio_list])))
        else:
            ppt_ratio_list = [1 for p in ppt_obs_list]

        # Use single mean monthly PPT for all cells
        # Assume ppt_obs_list is in month order
        fields = ppt_field_list + ratio_field_list
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as u_cursor:
            for row in u_cursor:
                for i, month in enumerate(month_list):
                    ppt_i = fields.index(ppt_field_format.format(month))
                    ratio_i = fields.index(ratio_field_format.format(month))

                    if ppt_obs_list[i] > 0:
                        row[ratio_i] = (
                            ppt_ratio_list[i] * row[ppt_i] / ppt_obs_list[i])
                    else:
                        row[ratio_i] = 0
                u_cursor.updateRow(row)
            del row


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Precipitation Ratio Parameters',
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

    # Convert relative paths to absolute path
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

    ppt_ratio_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
