#--------------------------------
# Name:         temp_adjust_parameters.py
# Purpose:      GSFLOW temperature adjustment parameters
# Notes:        ArcGIS 10.2+ Version
# Python:       2.7
#--------------------------------

import argparse
from collections import defaultdict
import ConfigParser
import datetime as dt
import logging
import os
import pprint
import sys

import arcpy
from arcpy import env

import support_functions as support


def temp_adjust_parameters(config_path):
    """Calculate GSFLOW Temperature Adjustment Parameters

    Parameters
    ----------
    config_path : str
        Project configuration file (.ini) path.

    Returns
    -------
    None

    """
    # Hardcoded HRU field formats for now
    tmax_field_fmt = 'TMAX_{:02d}'
    tmin_field_fmt = 'TMIN_{:02d}'
    tmax_adj_field_fmt = 'TMX_ADJ_{:02d}'
    tmin_adj_field_fmt = 'TMN_ADJ_{:02d}'

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
    log_file_name = 'temp_adjust_parameters_log.txt'
    log_console = logging.FileHandler(
        filename=os.path.join(hru.log_ws, log_file_name), mode='w')
    log_console.setLevel(logging.DEBUG)
    log_console.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger('').addHandler(log_console)
    logging.info('\nGSFLOW Temperature Adjustment Parameters')

    # Units
    temp_obs_units = support.get_param(
        'temp_obs_units', 'C', inputs_cfg).upper()
    temp_units_list = ['C', 'F', 'K']
    # Compare against the upper case of the values in the list
    #   but don't modify the acceptable units list
    if temp_obs_units not in temp_units_list:
        logging.error(
            '\nERROR: Invalid observed temperature units ({})\n  '
            'Valid units are: {}'.format(
                temp_obs_units, ', '.join(temp_units_list)))
        sys.exit()

    # Check input paths
    if not arcpy.Exists(hru.polygon_path):
        logging.error(
            '\nERROR: Fishnet ({}) does not exist'.format(
                hru.polygon_path))
        sys.exit()

    # Temperature calculation method
    try:
        temp_calc_method = inputs_cfg.get(
            'INPUTS', 'temperature_calc_method').upper()
    except:
        temp_calc_method = '1STA'
        logging.info('  Defaulting temperature_calc_method = {}'.format(
            temp_calc_method))

    temp_calc_options = ['ZONES', '1STA', 'LAPSE']
    if temp_calc_method not in temp_calc_options:
        logging.error(
            '\nERROR: Invalid temperature calculation method ({})\n  '
            'Valid methods are: {}'.format(
                temp_calc_method, ', '.join(temp_calc_options)))
        sys.exit()
    if temp_calc_method == 'LAPSE':
        logging.warning(
            '  If temperature calculation set to LAPSE,\n'
            '  it is not necessary to run the temp_adjust_parameters.py\n'
            '  Exiting\n')
        return False

    if temp_calc_method == 'ZONES':
        temp_zone_orig_path = inputs_cfg.get('INPUTS', 'temp_zone_path')
        try:
            temp_zone_field = inputs_cfg.get('INPUTS', 'temp_zone_field')
        except:
            logging.error(
                '\nERROR: temp_zone_field must be set in INI to apply '
                'zone specific temperature adjustments\n')
            sys.exit()
        try:
            temp_hru_id_field = inputs_cfg.get('INPUTS', 'temp_hru_id_field')
        except:
            temp_hru_id_field = None
            logging.warning(
                '  temp_hru_id_field was not set in the INI file\n'
                '  Temperature adjustments will not be changed to match station '
                'values'.format(temp_zone_field, hru.temp_zone_id_field))

        # Field name for TSTA hard coded, but could be changed to be read from
        # config file like temp_zone
        hru_tsta_field = 'HRU_TSTA'

        try:
            tmax_obs_field_fmt = inputs_cfg.get(
                'INPUTS', 'tmax_obs_field_format')
        except:
            tmax_obs_field_fmt = 'TMAX_{:02d}'
            logging.info('  Defaulting tmax_obs_field_format = {}'.format(
                tmax_obs_field_fmt))

        try:
            tmin_obs_field_fmt = inputs_cfg.get(
                'INPUTS', 'temp_obs_field_format')
        except:
            tmin_obs_field_fmt = 'TMIN_{:02d}'
            logging.info('  Defaulting tmin_obs_field_format = {}'.format(
                tmin_obs_field_fmt))

        if not arcpy.Exists(temp_zone_orig_path):
            logging.error(
                '\nERROR: Temperature Zone ({}) does not exist'.format(
                    temp_zone_orig_path))
            sys.exit()
        # temp_zone_path must be a polygon shapefile
        if arcpy.Describe(temp_zone_orig_path).datasetType != 'FeatureClass':
            logging.error(
                '\nERROR: temp_zone_path must be a polygon shapefile')
            sys.exit()

        # Check temp_zone_field
        if temp_zone_field.upper() in ['FID', 'OID']:
            temp_zone_field = arcpy.Describe(temp_zone_orig_path).OIDFieldName
            logging.warning(
                '\n  NOTE: Using {} to set {}\n'.format(
                    temp_zone_field, hru.temp_zone_id_field))
        elif not arcpy.ListFields(temp_zone_orig_path, temp_zone_field):
            logging.error(
                '\nERROR: temp_zone_field field {} does not exist\n'.format(
                    temp_zone_field))
            sys.exit()
        # Need to check that field is an int type
        # Should we only check active cells (HRU_TYPE > 0)?
        elif not [f.type for f in arcpy.Describe(temp_zone_orig_path).fields
                  if (f.name == temp_zone_field and
                      f.type in ['SmallInteger', 'Integer'])]:
            logging.error(
                '\nERROR: temp_zone_field field {} must be an integer type\n'.format(
                    temp_zone_field))
            sys.exit()
        # Need to check that field values are all positive
        # Should we only check active cells (HRU_TYPE > 0)?
        elif min([row[0] for row in arcpy.da.SearchCursor(
                temp_zone_orig_path, [temp_zone_field])]) <= 0:
            logging.error(
                '\nERROR: temp_zone_field values must be positive\n'.format(
                    temp_zone_field))
            sys.exit()

        # Check hru_tsta_field
        if not arcpy.ListFields(temp_zone_orig_path, hru_tsta_field):
            logging.error(
                '\nERROR: hru_tsta_field field {} does not exist\n'.format(
                    hru_tsta_field))
            sys.exit()
        # Need to check that field is an int type
        # Only check active cells (HRU_TYPE >0)?!
        elif not [f.type for f in arcpy.Describe(temp_zone_orig_path).fields
                  if (f.name == hru_tsta_field and
                      f.type in ['SmallInteger', 'Integer'])]:
            logging.error(
                '\nERROR: hru_tsta_field field {} must be an integer type\n'.format(
                    hru_tsta_field))
            sys.exit()
        # Need to check that field values are all positive
        # Only check active cells (HRU_TYPE >0)?!
        elif min([row[0] for row in arcpy.da.SearchCursor(
                temp_zone_orig_path, [hru_tsta_field])]) <= 0:
            logging.error(
                '\nERROR: hru_tsta_field values must be positive\n'.format(
                    hru_tsta_field))
            sys.exit()

        # Check temp_hru_id_field
        # temp_hru_id values are checked later
        if temp_hru_id_field is not None:
            if not arcpy.ListFields(temp_zone_orig_path, temp_hru_id_field):
                logging.error(
                    '\nERROR: temp_hru_id_field field {} does not exist\n'.format(
                        temp_hru_id_field))
                sys.exit()
            # Need to check that field is an int type
            elif not [f.type for f in arcpy.Describe(temp_zone_orig_path).fields
                      if (f.name == temp_hru_id_field and
                          f.type in ['SmallInteger', 'Integer'])]:
                logging.error(
                    '\nERROR: temp_hru_id_field field {} must be an integer type\n'.format(
                        temp_hru_id_field))
                sys.exit()
            # Need to check that field values are not negative (0 is okay)
            elif min([row[0] for row in arcpy.da.SearchCursor(
                    temp_zone_orig_path, [temp_hru_id_field])]) < 0:
                logging.error(
                    '\nERROR: temp_hru_id_field values cannot be negative\n'.format(
                        temp_hru_id_field))
                sys.exit()
    elif temp_calc_method == '1STA':
        # If a zone shapefile is not used, temperature must be set manually
        tmax_obs_list = inputs_cfg.get('INPUTS', 'tmax_obs_list')
        tmin_obs_list = inputs_cfg.get('INPUTS', 'tmin_obs_list')

        # Check that values are floats
        try:
            tmax_obs_list = map(float, tmax_obs_list.split(','))
        except ValueError:
            logging.error(
                '\nERROR: tmax_obs_list (mean monthly tmax) '
                'values could not be parsed as floats\n')
            sys.exit()
        try:
            tmin_obs_list = map(float, tmin_obs_list.split(','))
        except ValueError:
            logging.error(
                '\nERROR: tmin_obs_list (mean monthly tmin) '
                'values could not be parsed as floats\n')
            sys.exit()

        # Check that there are 12 values
        if len(tmax_obs_list) != 12:
            logging.error(
                '\nERROR: There must be exactly 12 mean monthly '
                'observed tmax values based to tmax_obs_list\n')
            sys.exit()
        logging.info(
            '  Observed mean monthly tmax ({}):\n    {}\n'
            '    (Script will assume these are listed in month order, '
            'i.e. Jan, Feb, ...)'.format(
                temp_obs_units, ', '.join(map(str, tmax_obs_list))))

        if len(tmin_obs_list) != 12:
            logging.error(
                '\nERROR: There must be exactly 12 mean monthly '
                'observed tmin values based to tmin_obs_list\n')
            sys.exit()
        logging.info(
            '  Observed mean monthly tmin ({}):\n    {}\n'
            '    (Script will assume these are listed in month order, '
            'i.e. Jan, Feb, ...)'.format(
                temp_obs_units, ', '.join(map(str, tmin_obs_list))))

        # Check if all the values are 0
        if tmax_obs_list == ([0.0] * 12):
            logging.error(
                '\nERROR: The observed tmax values are all 0.\n'
                '  To compute tmax adjustments, please set the tmax_obs_list '
                'parameter in the INI with\n  observed mean monthly tmax '
                'values (i.e. from a weather station)')
            sys.exit()
        if tmin_obs_list == ([0.0] * 12):
            logging.error(
                '\nERROR: The observed tmin values are all 0.\n'
                '  To compute tmin adjustments, please set the tmin_obs_list '
                'parameter in the INI with\n  observed mean monthly tmin '
                'values (i.e. from a weather station)')
            sys.exit()

        # Get the temperature HRU ID
        try:
            temp_hru_id = inputs_cfg.getint('INPUTS', 'temp_hru_id')
        except:
            temp_hru_id = 0

        # Check that the temp_hru_id is a valid cell hru_id
        # If temp_hru_id is 0, temperature adjustments will not be adjusted
        if temp_hru_id > 0:
            # Check that HRU_ID is valid
            logging.info('    Temperature HRU_ID: {}'.format(temp_hru_id))
            arcpy.MakeTableView_management(
                hru.polygon_path, "layer",
                "{} = {}".format(hru.id_field, temp_hru_id))
            if (temp_hru_id != 0 and
                    int(arcpy.GetCount_management("layer").getOutput(0)) == 0):
                logging.error(
                    '\nERROR: temp_hru_id {0} is not a valid cell hru_id'
                    '\nERROR: temp adjustments will NOT be forced to 1'
                    ' at cell {0}\n'.format(temp_hru_id))
                temp_hru_id = 0
            arcpy.Delete_management("layer")
        else:
            logging.info(
                '  Temperatures adjustments will not be adjusted to match '
                'station values\n    (temp_hru_id = 0)')

        # Could add a second check that HRU_TSTA has values >0

    # Build output folders if necessary
    temp_adj_temp_ws = os.path.join(hru.param_ws, 'temp_adjust')
    if not os.path.isdir(temp_adj_temp_ws):
        os.mkdir(temp_adj_temp_ws)
    temp_zone_path = os.path.join(temp_adj_temp_ws, 'temp_zone.shp')
    # temp_zone_clip_path = os.path.join(temp_adj_temp_ws, 'temp_zone_clip.shp')


    # Set ArcGIS environment variables
    arcpy.CheckOutExtension('Spatial')
    env.overwriteOutput = True
    # env.pyramid = 'PYRAMIDS -1'
    env.pyramid = 'PYRAMIDS 0'
    env.workspace = hru.param_ws
    env.scratchWorkspace = hru.scratch_ws

    # Set month list based on flags
    month_list = range(1, 13)
    tmax_field_list = [tmax_field_fmt.format(m) for m in month_list]
    tmin_field_list = [tmin_field_fmt.format(m) for m in month_list]
    tmax_adj_field_list = [tmax_adj_field_fmt.format(m) for m in month_list]
    tmin_adj_field_list = [tmin_adj_field_fmt.format(m) for m in month_list]

    # Check fields
    logging.info('\nAdding temperature adjust fields if necessary')
    # Temperature zone fields
    support.add_field_func(hru.polygon_path, hru.temp_zone_id_field, 'LONG')
    support.add_field_func(hru.polygon_path, hru.hru_tsta_field, 'SHORT')
    # Temperature adjustment fields
    for tmax_adj_field in tmax_adj_field_list:
        support.add_field_func(hru.polygon_path, tmax_adj_field, 'DOUBLE')
    for tmin_adj_field in tmin_adj_field_list:
        support.add_field_func(hru.polygon_path, tmin_adj_field, 'DOUBLE')

    # Calculate temperature zone ID
    if temp_calc_method == 'ZONES':
        logging.info('\nCalculating cell HRU Temperature Zone ID')
        temp_zone_desc = arcpy.Describe(temp_zone_orig_path)
        temp_zone_sr = temp_zone_desc.spatialReference
        logging.debug('  Zones:      {}'.format(temp_zone_orig_path))
        logging.debug('  Projection: {}'.format(temp_zone_sr.name))
        logging.debug('  GCS:        {}'.format(temp_zone_sr.GCS.name))

        # Reset temp_ZONE_ID
        logging.info('  Resetting {} to 0'.format(hru.temp_zone_id_field))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.temp_zone_id_field, 0, 'PYTHON')

        # If temp_zone spat_ref doesn't match hru_param spat_ref
        # Project temp_zone to hru_param spat ref
        # Otherwise, read temp_zone directly
        if hru.sr.name != temp_zone_sr.name:
            logging.info('  Projecting temperature zones...')
            # Set preferred transforms
            transform_str = support.transform_func(
                hru.sr, temp_zone_sr)
            logging.debug('    Transform: {}'.format(transform_str))
            # Project temp_zone shapefile
            arcpy.Project_management(
                temp_zone_orig_path, temp_zone_path, hru.sr,
                transform_str, temp_zone_sr)
            del transform_str
        else:
            arcpy.Copy_management(temp_zone_orig_path, temp_zone_path)

        # # Remove all unnecessary fields
        # for field in arcpy.ListFields(temp_zone_path):
        #     skip_field_list = temp_obs_field_list + [temp_zone_field, 'Shape']
        #     if field.name not in skip_field_list:
        #         try:
        #             arcpy.DeleteField_management(temp_zone_path, field.name)
        #         except:
        #             pass

        # Set temperature zone ID
        logging.info('  Setting {}'.format(hru.temp_zone_id_field))
        support.zone_by_centroid_func(
            temp_zone_path, hru.temp_zone_id_field, temp_zone_field,
            hru.polygon_path, hru.point_path, hru)
        # support.zone_by_area_func(
        #    temp_zone_layer, hru.temp_zone_id_field, temp_zone_field,
        #    hru.polygon_path, hru, hru_area_field, None, 50)

        # Set HRU_TSTA
        logging.info('  Setting {}'.format(hru.hru_tsta_field))
        support.zone_by_centroid_func(
            temp_zone_path, hru.hru_tsta_field, hru_tsta_field,
            hru.polygon_path, hru.point_path, hru)

        del temp_zone_desc, temp_zone_sr
    elif temp_calc_method == '1STA':
        # Set all cells to zone 1
        arcpy.CalculateField_management(
            hru.polygon_path, hru.temp_zone_id_field, 1, 'PYTHON')

    # Calculate adjustments
    logging.info('\nCalculating mean monthly temperature adjustments')
    if temp_calc_method == 'ZONES':
        # Read mean monthly values for each zone
        tmax_obs_dict = dict()
        tmin_obs_dict = dict()
        tmax_obs_field_list = [tmax_obs_field_fmt.format(m) for m in month_list]
        tmin_obs_field_list = [tmin_obs_field_fmt.format(m) for m in month_list]
        tmax_fields = [temp_zone_field] + tmax_obs_field_list
        tmin_fields = [temp_zone_field] + tmin_obs_field_list
        logging.debug('  Tmax Obs. Fields: {}'.format(', '.join(tmax_fields)))
        logging.debug('  Tmin Obs. Fields: {}'.format(', '.join(tmax_fields)))

        with arcpy.da.SearchCursor(temp_zone_path, tmax_fields) as s_cursor:
            for row in s_cursor:
                tmax_obs_dict[int(row[0])] = map(float, row[1:13])
        with arcpy.da.SearchCursor(temp_zone_path, tmin_fields) as s_cursor:
            for row in s_cursor:
                tmin_obs_dict[int(row[0])] = map(float, row[1:13])

        # Convert values to Celsius if necessary to match PRISM
        if temp_obs_units == 'F':
            tmax_obs_dict = {
                z: [(t - 32) * (5.0 / 9) for t in t_list]
                for z, t_list in tmax_obs_dict.items()}
            tmin_obs_dict = {
                z: [(t - 32) * (5.0 / 9) for t in t_list]
                for z, t_list in tmin_obs_dict.items()}
        elif temp_obs_units == 'K':
            tmax_obs_dict = {
                z: [(t - 273.15) for t in t_list]
                for z, t_list in tmax_obs_dict.items()}
            tmin_obs_dict = {
                z: [(t - 273.15) for t in t_list]
                for z, t_list in tmin_obs_dict.items()}

        tmax_zone_list = sorted(tmax_obs_dict.keys())
        tmin_zone_list = sorted(tmin_obs_dict.keys())
        logging.debug('  Tmax Zones: {}'.format(tmax_zone_list))
        logging.debug('  Tmin Zones: {}'.format(tmin_zone_list))

        # Print the observed temperature values
        logging.debug('  Observed Tmax')
        for zone, tmax_obs in tmax_obs_dict.items():
            logging.debug('    {}: {}'.format(
                zone, ', '.join(['{:.2f}'.format(x) for x in tmax_obs])))
        logging.debug('  Observed Tmin')
        for zone, tmin_obs in tmin_obs_dict.items():
            logging.debug('    {}: {}'.format(
                zone, ', '.join(['{:.2f}'.format(x) for x in tmin_obs])))

        # Default all zones to an adjustment of 0
        tmax_adj_dict = {z: [0] * 12 for z in tmax_zone_list}
        tmin_adj_dict = {z: [0] * 12 for z in tmin_zone_list}

        # Get list of HRU_IDs for each zone
        fields = [hru.temp_zone_id_field, hru.id_field]
        zone_hru_id_dict = defaultdict(list)
        with arcpy.da.SearchCursor(hru.polygon_path, fields) as s_cursor:
            for row in s_cursor:
                zone_hru_id_dict[int(row[0])].append(int(row[1]))

        # Check that TEMP_HRU_IDs are in the correct zone
        # Default all temperature zone HRU IDs to 0
        temp_hru_id_dict = {z: 0 for z in tmax_zone_list}
        if temp_hru_id_field is not None:
            fields = [temp_zone_field, temp_hru_id_field]
            logging.debug('  Temp Zone ID field: {}'.format(temp_zone_field))
            logging.debug('  Temp HRU ID field: {}'.format(temp_hru_id_field))
            with arcpy.da.SearchCursor(temp_zone_path, fields) as s_cursor:
                for row in s_cursor:
                    temp_zone = int(row[0])
                    hru_id = int(row[1])
                    if hru_id == 0 or hru_id in zone_hru_id_dict[temp_zone]:
                        temp_hru_id_dict[temp_zone] = hru_id
                        logging.debug('    {}: {}'.format(temp_zone, hru_id))
                    else:
                        logging.error(
                            '\nERROR: HRU_ID {} is not in temperature ZONE {}'.format(
                                hru_id, temp_hru_id_dict[temp_zone]))
                        sys.exit()

        # Get gridded tmax values for each TEMP_HRU_ID
        fields = [hru.temp_zone_id_field, hru.id_field] + tmax_field_list
        with arcpy.da.SearchCursor(hru.polygon_path, fields) as s_cursor:
            for row in s_cursor:
                temp_zone = int(row[0])
                hru_id = int(row[1])
                if hru_id == 0:
                    pass
                elif hru_id in temp_hru_id_dict.values():
                    tmax_gridded_list = map(float, row[2:14])
                    tmax_obs_list = tmax_obs_dict[temp_zone]
                    tmax_adj_list = [
                        float(o) - t
                        for o, t in zip(tmax_obs_list, tmax_gridded_list)]
                    tmax_adj_dict[temp_zone] = tmax_adj_list

        # Get gridded tmin values for each TEMP_HRU_ID
        fields = [hru.temp_zone_id_field, hru.id_field] + tmin_field_list
        with arcpy.da.SearchCursor(hru.polygon_path, fields) as s_cursor:
            for row in s_cursor:
                temp_zone = int(row[0])
                hru_id = int(row[1])
                if hru_id == 0:
                    pass
                elif hru_id in temp_hru_id_dict.values():
                    tmin_gridded_list = map(float, row[2:14])
                    tmin_obs_list = tmin_obs_dict[temp_zone]
                    tmin_adj_list = [
                        float(o) - t
                        for o, t in zip(tmin_obs_list, tmin_gridded_list)]
                    tmin_adj_dict[temp_zone] = tmin_adj_list
        del temp_hru_id_dict, zone_hru_id_dict, fields

        logging.debug('  Tmax Adjustment Factors:')
        for k, v in tmax_adj_dict.items():
            logging.debug('    {}: {}'.format(
                k, ', '.join(['{:.3f}'.format(x) for x in v])))

        logging.debug('  Tmin Adjustment Factors:')
        for k, v in tmin_adj_dict.items():
            logging.debug('    {}: {}'.format(
                k, ', '.join(['{:.3f}'.format(x) for x in v])))

        logging.debug('\nWriting adjustment values to hru_params')
        fields = [hru.temp_zone_id_field]
        fields.extend(tmax_field_list + tmax_adj_field_list)
        fields.extend(tmin_field_list + tmin_adj_field_list)
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as u_cursor:
            for row in u_cursor:
                zone = int(row[0])
                for i, month in enumerate(month_list):
                    tmax_i = fields.index(tmax_field_fmt.format(month))
                    tmax_adj_i = fields.index(tmax_adj_field_fmt.format(month))
                    row[tmax_adj_i] = (
                        row[tmax_i] - tmax_obs_dict[zone][i] +
                        tmax_adj_dict[zone][i])

                    tmin_i = fields.index(tmin_field_fmt.format(month))
                    tmin_adj_i = fields.index(tmin_adj_field_fmt.format(month))
                    row[tmin_adj_i] = (
                        row[tmin_i] - tmin_obs_dict[zone][i] +
                        tmin_adj_dict[zone][i])
                u_cursor.updateRow(row)
            del row

    elif temp_calc_method == '1STA':
        # Get gridded temperature at temp_HRU_ID
        tmax_fields = [hru.id_field] + tmax_field_list
        tmin_fields = [hru.id_field] + tmin_field_list
        logging.debug('  Tmax Fields: {}'.format(', '.join(tmax_field_list)))
        logging.debug('  Tmin Fields: {}'.format(', '.join(tmin_field_list)))

        # Convert values to Celsius if necessary to match PRISM
        if temp_obs_units == 'F':
            tmax_obs_list = [(t - 32) * (5.0 / 9) for t in tmax_obs_list]
            tmin_obs_list = [(t - 32) * (5.0 / 9) for t in tmin_obs_list]
        elif temp_obs_units == 'K':
            tmax_obs_list = [t - 273.15 for t in tmax_obs_list]
            tmin_obs_list = [t - 273.15 for t in tmin_obs_list]
        if temp_obs_units != 'C':
            logging.info(
                '\nConverted Mean Monthly Tmax ({}):\n  {}'.format(
                    temp_obs_units, ', '.join(map(str, tmax_obs_list))))
            logging.info(
                'Converted Mean Monthly Tmin ({}):\n  {}'.format(
                    temp_obs_units, ', '.join(map(str, tmin_obs_list))))

        # Scale all adjustments so gridded temperature will match observed
        # temperature at target cell
        if temp_hru_id != 0:
            tmax_gridded_list = map(float, arcpy.da.SearchCursor(
                hru.polygon_path, tmax_fields,
                '"{}" = {}'.format(hru.id_field, temp_hru_id)).next()[1:])
            logging.debug('  Gridded Tmax: {}'.format(
                ', '.join(['{:.2f}'.format(p) for p in tmax_gridded_list])))

            tmin_gridded_list = map(float, arcpy.da.SearchCursor(
                hru.polygon_path, tmin_fields,
                '"{}" = {}'.format(hru.id_field, temp_hru_id)).next()[1:])
            logging.debug('  Gridded Tmin: {}'.format(
                ', '.join(['{:.2f}'.format(p) for p in tmin_gridded_list])))

            # Difference of MEASURED or OBSERVED TEMP to GRIDDED TEMP
            tmax_adj_list = [
                float(o) - t
                for o, t in zip(tmax_obs_list, tmax_gridded_list)]
            logging.info('  Obs./Gridded: {}'.format(
                ', '.join(['{:.3f}'.format(p) for p in tmax_adj_list])))

            tmin_adj_list = [
                float(o) - t
                for o, t in zip(tmin_obs_list, tmin_gridded_list)]
            logging.info('  Obs./Gridded: {}'.format(
                ', '.join(['{:.3f}'.format(p) for p in tmin_adj_list])))
        else:
            tmax_adj_list = [0 for p in tmax_obs_list]
            tmin_adj_list = [0 for p in tmin_obs_list]

        # Use single mean monthly tmax for all cells
        # Assume tmax_obs_list is in month order
        fields = tmax_field_list + tmax_adj_field_list
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as u_cursor:
            for row in u_cursor:
                for i, month in enumerate(month_list):
                    tmax_i = fields.index(tmax_field_fmt.format(month))
                    tmax_adj_i = fields.index(tmax_adj_field_fmt.format(month))
                    row[tmax_adj_i] = (
                            row[tmax_i] - tmax_obs_list[i] + tmax_adj_list[i])
                u_cursor.updateRow(row)
            del row

        # Use single mean monthly tmax for all cells
        # Assume tmax_obs_list is in month order
        fields = tmin_field_list + tmin_adj_field_list
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as u_cursor:
            for row in u_cursor:
                for i, month in enumerate(month_list):
                    tmin_i = fields.index(tmin_field_fmt.format(month))
                    tmin_adj_i = fields.index(
                        tmin_adj_field_fmt.format(month))
                    row[tmin_adj_i] = (
                            row[tmin_i] - tmin_obs_list[i] + tmin_adj_list[i])
                u_cursor.updateRow(row)
            del row


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Temperature Adjustment Parameters',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i', '--ini', required=True,
        help='Project input file', metavar='PATH')
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

    temp_adjust_parameters(config_path=args.ini)
