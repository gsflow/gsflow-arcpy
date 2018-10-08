#--------------------------------
# Name:         prms_template_fill.py
# Purpose:      Fill PRMS Parameter File Template
# Notes:        ArcGIS 10.2+ Version
# Python:       2.7
#--------------------------------

import argparse
from collections import defaultdict
import ConfigParser
import datetime as dt
import logging
import operator
import os
import sys

import arcpy

import support_functions as support


def prms_template_fill(config_path, overwrite_flag=False, debug_flag=False):
    """Fill PRMS Parameter Template File

    Args:
        config_file (str): Project config file path
        ovewrite_flag (bool): if True, overwrite existing files
        debug_flag (bool): if True, enable debug level logging

    Returns:
        None
    """

    param_formats = {1: '{:d}', 2: '{:f}', 3: '{:f}', 4: '{}'}

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
    log_file_name = 'prms_template_log.txt'
    log_console = logging.FileHandler(
        filename=os.path.join(hru.log_ws, log_file_name), mode='w')
    log_console.setLevel(logging.DEBUG)
    log_console.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger('').addHandler(log_console)
    logging.info('\nFilling PRMS Parameter File Template')

    # Read parameters from config file
    hru.polygon_path = inputs_cfg.get('INPUTS', 'hru_fishnet_path')
    hru.fid_field = inputs_cfg.get('INPUTS', 'orig_fid_field')
    parameter_ws = inputs_cfg.get('INPUTS', 'parameter_folder')
    try:
        prms_parameter_ws = inputs_cfg.get('INPUTS', 'prms_parameter_folder')
    except ConfigParser.NoOptionError:
        prms_parameter_ws = inputs_cfg.get('INPUTS', 'parameter_folder')
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'prms_parameter_ws', prms_parameter_ws))
    prms_dimen_csv_path = inputs_cfg.get('INPUTS', 'prms_dimen_csv_path')
    prms_param_csv_path = inputs_cfg.get('INPUTS', 'prms_param_csv_path')

    # Get input DEM units and desired output HRU_ELEV units
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
    elev_units = inputs_cfg.getint('INPUTS', 'elev_units')
    elev_unit_types = {0: 'feet', 1: 'meter'}
    try:
        elev_units = elev_unit_types[elev_units]
    except:
        logging.error(
            '\nERROR: elev_units "{}" is not supported\n'.format(elev_units))
        sys.exit()
    if dem_units == 'feet' and elev_units == 'meter':
        elev_unit_scalar = 0.3048
    elif dem_units == 'meter' and elev_units == 'feet':
        elev_unit_scalar = (1.0 / 0.3048)
    else:
        elev_unit_scalar = 1.0

    # Write parameter/dimensions to separate files based on "PARAM_FILE"
    #   value in prms_parameters.csv and prms_dimensions.csv
    try:
        single_param_file_flag = inputs_cfg.getboolean(
            'INPUTS', 'single_param_file_flag')
    except ConfigParser.NoOptionError:
        single_param_file_flag = False
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'single_param_file_flag', single_param_file_flag))
    if single_param_file_flag:
        try:
            single_param_file_name = inputs_cfg.get(
                'INPUTS', 'single_param_file_name')
        except ConfigParser.NoOptionError:
            single_param_file_name = 'prms_inputs.param'
            logging.info(
                '  Missing INI parameter, setting {} = {}'.format(
                'single_param_file_name', single_param_file_name))

    # Write nhru gridded parameters as single column or array
    try:
        param_column_flag = inputs_cfg.getboolean(
            'INPUTS', 'param_column_flag')
    except ConfigParser.NoOptionError:
        param_column_flag = False
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'param_column_flag', param_column_flag))

    # Scratch workspace
    try:
        scratch_name = inputs_cfg.get('INPUTS', 'scratch_name')
    except ConfigParser.NoOptionError:
        scratch_name = 'in_memory'
        logging.info(
            '  Missing INI parameter, setting {} = {}'.format(
                'scratch_name', scratch_name))

    # Cascades
    crt_ws = os.path.join(parameter_ws, 'cascade_work')
    crt_dimension_path = os.path.join(crt_ws, 'parameter_dimensions.txt')
    crt_parameter_path = os.path.join(crt_ws, 'cascade.param')
    crt_gw_parameter_path = os.path.join(crt_ws, 'groundwater_cascade.param')
    # fill_ws = os.path.join(parameter_ws, 'fill_work')
    # crt_gw_parameter_path = os.path.join(
    #     fill_ws, 'groundwater_cascade.param')

    # Strings to search PRMS parameter file for
    # Newline character is required after title
    file_header_str = 'PRMS parameter file generated with gsflow-arcpy-tools version X\n'
    # file_header_str = 'Default file generated by model\nVersion: 1.7'
    dimen_header_str = '** Dimensions **'
    param_header_str = '** Parameters **'
    break_str = '####'

    # Check input paths
    if not arcpy.Exists(hru.polygon_path):
        logging.error(
            '\nERROR: The fishnet does not exist\n  {}'.format(
                hru.polygon_path))
        sys.exit()
    # if not os.path.isfile(prms_template_path):
    #    logging.error('\nERROR: The template parameter file does not exist\n')
    #    sys.exit()
    if not os.path.isfile(prms_dimen_csv_path):
        logging.error(
            '\nERROR: The dimensions CSV file does not exist\n  {}'.format(
                prms_dimen_csv_path))
        sys.exit()
    if not os.path.isfile(prms_param_csv_path):
        logging.error(
            '\nERROR: The parameters CSV file does not exist\n  {}'.format(
                prms_param_csv_path))
        sys.exit()
    if not os.path.isdir(crt_ws):
        logging.error(
            '\nERROR: Cascades folder does not exist'
            '\nERROR:   {}'
            '\nERROR: Try re-running CRT using stream_parameters.py\n'.format(
                crt_ws))
        sys.exit()
    if not os.path.isfile(crt_dimension_path):
        logging.error(
            '\nERROR: Cascades dimension file does not exist'
            '\nERROR:   {}'
            '\nERROR: Try re-running CRT using stream_parameters.py\n'.format(
                crt_dimension_path))
        sys.exit()
    if not os.path.isfile(crt_parameter_path):
        logging.error(
            '\nERROR: Cascades parameter file does not exist'
            '\nERROR:   {}'
            '\nERROR: Try re-running CRT using stream_parameters.py\n'.format(
                crt_parameter_path))
        sys.exit()
    if not os.path.isfile(crt_gw_parameter_path):
        logging.error(
            '\nERROR: Groundwater cascades parameter file does not exist'
            '\nERROR:   {}'
            '\nERROR: Try re-running CRT using crt_fill_work.py\n'.format(
                crt_gw_parameter_path))
        sys.exit()
    # if not os.path.isfile(crt_gw_parameter_path):
    #    logging.error(
    #        '\nERROR: Groundwater cascades parameter file does not exist'
    #        '\nERROR:   {}'
    #        '\nERROR: Try re-running CRT using stream_parameters\n'.format(
    #             crt_gw_parameter_path))
    #    sys.exit()


    # Get number of cells in fishnet
    fishnet_count = int(arcpy.GetCount_management(
        hru.polygon_path).getOutput(0))
    logging.info('  Fishnet cells: {}'.format(fishnet_count))


    # Read in dimensions from CSV
    logging.info('\nReading dimensions CSV')
    dimen_names = dict()
    dimen_files = dict()
    dimen_sizes = dict()
    with open(prms_dimen_csv_path, 'r') as input_f:
        dimen_lines = input_f.readlines()
    input_f.close()
    # Dimensions can be set to a value, a field, or not set
    dimen_lines = [l.strip().split(',') for l in dimen_lines]
    header = dimen_lines[0]
    for line in dimen_lines[1:]:
        dimen_name = line[header.index('NAME')]
        dimen_names[dimen_name] = dimen_name
        logging.debug('  {}'.format(dimen_name))

        # What should the default parameter file name be if not set?
        if single_param_file_flag:
            dimen_file = os.path.join(
                prms_parameter_ws, single_param_file_name)
        elif 'PARAM_FILE' not in header:
            dimen_file = os.path.join(prms_parameter_ws, 'prms_inputs.param')
            logging.info(
                '  PARAM_FILE field not in dimensions CSV\n'
                '  Defaulting to {}'.format(dimen_file))
        elif line[header.index('PARAM_FILE')] == '':
            dimen_file = os.path.join(prms_parameter_ws, 'prms_inputs.param')
            logging.info(
                '  PARAM_FILE value not set for dimension: {}\n'
                '  Defaulting to {}'.format(dimen_name, dimen_file))
        else:
            dimen_file = os.path.join(
                prms_parameter_ws, line[header.index('PARAM_FILE')] + '.param')
        dimen_files[dimen_name] = dimen_file

        dimen_size = line[header.index('SIZE')]
        if dimen_size.lower() in ['calculated', 'config_file']:
            dimen_sizes[dimen_name] = dimen_size
        elif not dimen_size:
            dimen_sizes[dimen_name] = ''
        else:
            # Don't force to integer type unless necessary since values are
            # written back out as strings
            dimen_sizes[dimen_name] = dimen_size
            # dimen_sizes[dimen_name] = int(dimen_size)
        del dimen_size

    # Set CALCULATED dimension values
    # These parameters equal the fishnet cell count
    for dimen_name in ['ngw', 'ngwcell', 'nhru', 'nhrucell', 'nssr']:
        if dimen_sizes[dimen_name].lower() == 'calculated':
            dimen_sizes[dimen_name] = fishnet_count
            logging.info('  {} = {}'.format(
                dimen_name, dimen_sizes[dimen_name]))

    # Getting number of lakes
    if dimen_sizes['nlake'].lower() == 'calculated':
        logging.info('\nCalculating number of lakes')
        #logging.info('  Lake cells are {} >= 0'.format(hru.lake_id_field))
        value_fields = (hru.id_field, hru.lake_id_field)
        with arcpy.da.SearchCursor(hru.polygon_path, value_fields) as s_cursor:
            dimen_sizes['nlake'] = max(list(
                [int(row[1]) for row in s_cursor if int(row[1]) > 0]))
        logging.info('  nlakes = {}'.format(dimen_sizes['nlake']))

    # Getting number of lake cells
    if dimen_sizes['nlake_hrus'].lower() == 'calculated':
        logging.info('\nCalculating number of lake cells')
        logging.info('  Lake cells are {} >= 0'.format(hru.lake_id_field))
        value_fields = (hru.id_field, hru.lake_id_field)
        with arcpy.da.SearchCursor(hru.polygon_path, value_fields) as s_cursor:
            dimen_sizes['nlake_hrus'] = len(list(
                [int(row[1]) for row in s_cursor if int(row[1]) > 0]))
        logging.info('  nlake cells = {}'.format(dimen_sizes['nlake_hrus']))

    # Getting number of stream cells
    if dimen_sizes['nreach'].lower() == 'calculated':
        logging.info('Calculating number of stream cells')
        logging.info('  Stream cells are {} >= 0'.format(hru.krch_field))
        value_fields = (hru.id_field, hru.krch_field)
        with arcpy.da.SearchCursor(hru.polygon_path, value_fields) as s_cursor:
            dimen_sizes['nreach'] = len(list(
                [int(row[1]) for row in s_cursor if int(row[1]) > 0]))
        logging.info('  nreach = {}'.format(dimen_sizes['nreach']))

    # Getting number of stream segments
    if dimen_sizes['nsegment'].lower() == 'calculated':
        logging.info('Calculating number of unique stream segments')
        logging.info('  Stream segments are {} >= 0'.format(hru.iseg_field))
        value_fields = (hru.id_field, hru.iseg_field)
        with arcpy.da.SearchCursor(hru.polygon_path, value_fields) as s_cursor:
            dimen_sizes['nsegment'] = len(list(set(
                [int(row[1]) for row in s_cursor if int(row[1]) > 0])))
        logging.info('  nsegment = {}'.format(dimen_sizes['nsegment']))

    # Getting number of subbasins
    if dimen_sizes['nsub'].lower() == 'calculated':
        logging.info('Calculating number of unique subbasins')
        logging.info('  Subbasins are {} >= 0'.format(hru.subbasin_field))
        value_fields = (hru.id_field, hru.subbasin_field)
        with arcpy.da.SearchCursor(hru.polygon_path, value_fields) as s_cursor:
            dimen_sizes['nsub'] = len(list(set(
                [int(row[1]) for row in s_cursor if int(row[1]) > 0])))
        logging.info('  nsub = {}'.format(dimen_sizes['nsub']))

    # Read in CRT dimensions
    if (dimen_sizes['ncascade'].lower() == 'calculated' or
            dimen_sizes['ncascdgw'].lower() == 'calculated'):
        logging.info('\nReading CRT dimensions')
        with open(crt_dimension_path, 'r') as input_f:
            crt_dimen_lines = [line.strip() for line in input_f.readlines()]
        input_f.close()
        crt_dimen_break_i_list = [
            i for i, x in enumerate(crt_dimen_lines) if x == break_str]
        for i in crt_dimen_break_i_list:
            logging.info('  {} = {}'.format(
                crt_dimen_lines[i + 1], crt_dimen_lines[i + 2]))
            dimen_sizes[crt_dimen_lines[i + 1]] = int(crt_dimen_lines[i + 2])
        del crt_dimen_lines, crt_dimen_break_i_list

    # Set CONFIG file dimension values
    config_file_dimensions = [
        d_name for d_name, d_size in sorted(dimen_sizes.items())
        if type(d_size) is str and d_size.lower() == 'config_file']
    if config_file_dimensions:
        logging.info('Reading configuration file dimensions')
        for dimen_name in config_file_dimensions:
            logging.info('  {}'.format(dimen_name))
            try:
                dimen_sizes[dimen_name] = inputs_cfg.getint(
                    'INPUTS', dimen_name)
            except ConfigParser.NoOptionError:
                logging.error(
                    '  Dimension set to "config_file" in {} but not found in '
                    'config file, exiting'.format(
                        os.path.basename(prms_dimen_csv_path)))


    # Link HRU fishnet field names to parameter names in '.param'
    param_names = dict()
    param_files = dict()
    param_dimen_counts = dict()
    param_dimen_names = dict()
    param_value_counts = dict()
    param_types = dict()
    param_defaults = dict()
    param_values = defaultdict(dict)

    # Read in parameters from CSV
    logging.info('\nReading parameters CSV')
    with open(prms_param_csv_path, 'r') as input_f:
        param_lines = input_f.readlines()
    input_f.close()
    param_lines = [l.strip().split(',') for l in param_lines]
    header = param_lines[0]
    for line in param_lines[1:]:
        # Get parameters from CSV line
        param_name = line[header.index('NAME')]
        logging.debug('  {}'.format(param_name))
        # This assumes multiple dimensions are separated by semicolon
        dimen_names = line[header.index('DIMENSION_NAMES')].split(';')

        # What should the default parameter file name be if not set?
        if single_param_file_flag:
            param_file = os.path.join(
                prms_parameter_ws, single_param_file_name)
        elif 'PARAM_FILE' not in header:
            param_file = os.path.join(prms_parameter_ws, 'prms_inputs.param')
            logging.info(
                '  PARAM_FILE field not in parameters CSV\n'
                '  Defaulting to {}'.format(param_file))
        elif line[header.index('PARAM_FILE')] == '':
            param_file = os.path.join(prms_parameter_ws, 'prms_inputs.param')
            logging.info(
                '  PARAM_FILE value not set for parameter: {}\n'
                '  Defaulting to {}'.format(param_name, param_file))
        else:
            param_file = os.path.join(
                prms_parameter_ws, line[header.index('PARAM_FILE')] + '.param')

        # Check that parameter type is 1, 2, 3, or 4
        param_type = int(line[header.index('TYPE')])
        if param_type not in [1, 2, 3, 4]:
            logging.error(
                '\nERROR: Parameter type {} is invalid'
                '\nERROR: {}'.format(param_type, line))
            sys.exit()

        # This will initially read defaults in as a list
        param_default = line[header.index('DEFAULT_VALUE'):]

        # Removing empty strings avoids checking ints/floats
        param_default = [l for l in param_default if l]

        # For empty lists, set to none
        if not param_default:
            param_default = None
        # For single value lists, get first value
        # Check that param_default is a number or field name
        elif len(param_default) == 1:
            param_default = param_default[0]
            if isfloat(param_default) and param_type == 1:
                param_default = int(param_default)
            elif isfloat(param_default) and param_type in [2, 3]:
                param_default = float(param_default)
            elif param_default.lower() in ['calculated', 'config_file', 'crt_file']:
                pass
            elif arcpy.ListFields(hru.polygon_path, param_default):
                pass
            else:
                logging.error(
                    '\nERROR: Default value {} was not parsed'
                    '\nERROR: {}'.format(param_default, line))
                sys.exit()
        # For multi-value lists, convert values to int/float
        elif len(param_default) >= 2:
            if param_type == 1:
                param_default = map(int, param_default)
            elif param_type in [2, 3]:
                param_default = map(float, param_default)
            else:
                logging.error(
                    '\nERROR: Default value {} was not parsed'
                    '\nERROR: {}'.format(param_default, line))
                sys.exit()

        # Check that dimension names are valid
        for dimen_name in dimen_names:
            if dimen_name not in dimen_sizes.keys():
                logging.error(
                    '\nERROR: The dimension {} is not set in the '
                    'dimension CSV file'.format(dimen_name))
                sys.exit()

        # Calculate number of dimensions
        dimen_count = str(len(dimen_names))

        # Calculate number of values
        values_count = prod(
            [int(dimen_sizes[dn]) for dn in dimen_names
             if dimen_sizes[dn]])

        # Write parameter to dictionaries
        param_names[param_name] = param_name
        param_files[param_name] = param_file
        param_dimen_counts[param_name] = dimen_count
        param_dimen_names[param_name] = dimen_names
        param_value_counts[param_name] = values_count
        param_types[param_name] = param_type
        param_defaults[param_name] = param_default

    # Apply default values to full dimension
    logging.info('\nSetting static parameters from defaults')
    for param_name, param_default in param_defaults.items():
        param_value_count = param_value_counts[param_name]
        # Skip if not set
        if param_default is None:
            continue
        # Skip if still a string (field names)
        elif type(param_default) is str:
            continue
        # For float/int, apply default across dimension size
        elif type(param_default) is float or type(param_default) is int:
            for i in range(param_value_count):
                param_values[param_name][i] = param_default
        # For lists of floats, match up one-to-one for now
        elif len(param_default) == param_value_count:
            for i in range(param_value_count):
                param_values[param_name][i] = param_default[i]
        else:
            logging.error(
                '\nERROR: The default value(s) ({0}) could not be '
                'broadcast to the dimension length ({1})'.format(
                    param_default, param_value_count))
            sys.exit()

    # Set CONFIG file parameter values
    config_file_parameters = [
        p_name for p_name, p_value in sorted(param_defaults.items())
        if type(p_value) is str and p_value.lower() == 'config_file']
    if config_file_parameters:
        logging.info('Reading configuration file parameters')
        for param_name in config_file_parameters:
            logging.info('  {}'.format(param_name))
            try:
                values = inputs_cfg.get('INPUTS', param_name)
            except ConfigParser.NoOptionError:
                logging.error(
                    '  Parameter set to "config_file" in {} but not found in '
                    'config file, exiting'.format(
                        os.path.basename(prms_dimen_csv_path)))

            # Convert comma separate strings to lists
            param_values[param_name] = {
                i: v for i, v in enumerate(values.split(','))}

            # Convert the strings to the appropriate type
            if param_types[param_name] == 1:
                param_values[param_name] = {
                    k: int(v) for k, v in param_values[param_name].items()}
            elif param_types[param_name] in [2, 3]:
                param_values[param_name] = {
                    k: float(v) for k, v in param_values[param_name].items()}

            # Try and honor dimension value from CSV
            # Repeat values if actual value count doesn't match expected count
            #   (from dimensions)
            # For now, only apply to INI parameters with a single value
            #   and dimensions greater than 1
            param_value_count = param_value_counts[param_name]
            if ((len(param_values[param_name]) != param_value_count) and
                    (len(param_values[param_name]) == 1) and
                    (param_value_count > 1)):
                value = param_values[param_name].copy()
                param_values[param_name] = {}
                for i in range(param_value_count):
                    param_values[param_name][i] = value[0]

    # Read in HRU parameter data from fishnet polygon
    logging.info('\nReading in variable parameters from fishnet')
    param_fields = {
        k: v for k, v in param_defaults.items()
        if (type(v) is str and
            v.lower() not in ['calculated', 'config_file', 'crt_file'])
    }
    value_fields = param_fields.values()

    # Use HRU_ID to uniquely identify each cell
    if hru.id_field not in value_fields:
        value_fields.append(hru.id_field)
    hru_id_i = value_fields.index(hru.id_field)

    # Read in each cell parameter value
    with arcpy.da.SearchCursor(hru.polygon_path, value_fields) as s_cursor:
        for row in s_cursor:
            for field_i, (param, field) in enumerate(param_fields.items()):
                if param_types[param] == 1:
                    param_values[param][row[hru_id_i]] = int(row[field_i])
                elif param_types[param] in [2, 3]:
                    param_values[param][row[hru_id_i]] = float(row[field_i])
                elif param_types[param] == 4:
                    param_values[param][row[hru_id_i]] = row[field_i]
                # param_values[param][row[hru_id_i]] = row[field_i]

    # Calculate number of columns
    with arcpy.da.SearchCursor(
            hru.polygon_path, (hru.id_field, hru.col_field)) as s_cursor:
        ncol = len(list(set([int(row[1]) for row in s_cursor])))

    # # DEADBEEF - Per Rich this is not needed anymore
    # # The following will override the parameter CSV values
    # # Calculate basin_area from active cells (land and lake)
    # logging.info('\nCalculating basin area')
    # param_names['basin_area'] = 'basin_area'
    # param_dimen_counts['basin_area'] = 1
    # param_dimen_names['basin_area'] = ['one']
    # param_value_counts['basin_area'] = dimen_sizes['one']
    # param_types['basin_area'] = 2
    # value_fields = (hru.id_field, hru.type_field, hru.area_field)
    # with arcpy.da.SearchCursor(hru.polygon_path, value_fields) as s_cursor:
    #     param_values['basin_area'][0] = sum(
    #         [float(row[2]) for row in s_cursor if int(row[1]) >= 1])
    # logging.info('  basin_area = {} acres'.format(
    #     param_values['basin_area'][0]))

    # Convert DEM_ADJ units (if necessary)
    if elev_unit_scalar != 1.0:
        logging.info('\nScaling DEM_ADJ units')
        logging.info('  DEM Units:  {}'.format(dem_units))
        logging.info('  Elev Units: {}'.format(elev_units))
        logging.info('  Multiplier: {}'.format(elev_unit_scalar))
        param_values['hru_elev'] = {
            k: v * elev_unit_scalar
            for k, v in param_values['hru_elev'].items()}

    # Calculate mean monthly maximum temperature for all active cells
    logging.info('\nCalculating tmax_index')
    logging.info('  Converting Celsius to Farenheit')
    param_names['tmax_index'] = 'tmax_index'
    param_dimen_counts['tmax_index'] = 1
    param_dimen_names['tmax_index'] = ['nmonths']
    param_value_counts['tmax_index'] = int(dimen_sizes['nmonths'])
    param_types['tmax_index'] = 2
    tmax_field_list = ['TMAX_{:02d}'.format(m) for m in range(1, 13)]
    for i, tmax_field in enumerate(tmax_field_list):
        tmax_values = [row[1] for row in arcpy.da.SearchCursor(
            hru.polygon_path, (hru.type_field, tmax_field),
            where_clause='"{}" >= 1'.format(hru.type_field))]
        tmax_c = sum(tmax_values) / len(tmax_values)
        tmax_f = 1.8 * tmax_c + 32
        param_values['tmax_index'][i] = tmax_f
        logging.info('  {} = {}'.format(
            tmax_field, param_values['tmax_index'][i]))
        del tmax_values

    #
    logging.info('\nCalculating rain_adj/snow_adj')
    ratio_field_list = ['PPT_RT_{:02d}'.format(m) for m in range(1, 13)]
    param_names['rain_adj'] = 'rain_adj'
    param_dimen_counts['rain_adj'] = 2
    param_dimen_names['rain_adj'] = ['nhru', 'nmonths']
    param_value_counts['rain_adj'] = 12 * fishnet_count
    param_types['rain_adj'] = 2

    param_names['snow_adj'] = 'snow_adj'
    param_dimen_counts['snow_adj'] = 2
    param_dimen_names['snow_adj'] = ['nhru', 'nmonths']
    param_value_counts['snow_adj'] = 12 * fishnet_count
    param_types['snow_adj'] = 2

    ratio_values = []
    for i, ratio_field in enumerate(ratio_field_list):
        ratio_values.extend([
            float(row[1]) for row in sorted(arcpy.da.SearchCursor(
                hru.polygon_path, (hru.id_field, ratio_field)))])
    for i, value in enumerate(ratio_values):
        param_values['rain_adj'][i] = value
        param_values['snow_adj'][i] = value
    del ratio_values

    #
    logging.info('\nCalculating subbasin_down')
    param_names['subbasin_down'] = 'subbasin_down'
    param_dimen_counts['subbasin_down'] = 1
    param_dimen_names['subbasin_down'] = ['nsub']
    param_value_counts['subbasin_down'] = dimen_sizes['nsub']
    param_types['subbasin_down'] = 1
    # Get list of subbasins and downstream cell for each stream/lake cell
    # Downstream is calulated from flow direction
    # logging.info('Cell out-flow dictionary')
    cell_dict = dict()
    fields = [
        hru.type_field, hru.krch_field, hru.lake_id_field,
        hru.subbasin_field, hru.flow_dir_field,
        hru.col_field, hru.row_field, hru.id_field]
    for row in arcpy.da.SearchCursor(hru.polygon_path, fields):
        # Skip inactive cells
        if int(row[0]) == 0:
            continue
        # Skip non-lake and non-stream cells
        elif (int(row[1]) == 0 and int(row[2]) == 0):
            continue
        # Read in parameters
        cell = (int(row[5]), int(row[6]))
        # support.next_row_col(FLOW_DIR, CELL)
        # HRU_ID, SUBBASIN, NEXT_CELL
        cell_dict[cell] = [
            int(row[7]), int(row[3]), support.next_row_col(int(row[4]), cell)]
        del cell

    # Get subset of cells if subbasin != next_subbasin
    subbasin_list = []
    # CELL, (HRU_ID, SUBBASIN, NEXT_CELL)
    # for cell, row in cell_dict.items():
    for cell, (hru_id, subbasin, next_cell) in cell_dict.items():
        # Skip cells that are already subbasin 0 (inactive?)
        # If next cell isn't in list, assume next cell is out of the model
        #   and set exit gauge subbasin to 0
        # If the subbasin of the current cell doesn't match the subbasin
        #   of the next cell, save the down subbasin
        if subbasin == 0:
            continue
        elif next_cell not in cell_dict.keys():
            if [subbasin, 0] not in subbasin_list:
                subbasin_list.append([subbasin, 0])
        elif subbasin != cell_dict[next_cell][1]:
            subbasin_list.append([subbasin, cell_dict[next_cell][1]])
    for i, (subbasin, subbasin_down) in enumerate(sorted(subbasin_list)):
        param_values['subbasin_down'][i] = subbasin_down
        logging.debug('  {}'.format(
            param_values['subbasin_down'][i]))
    del subbasin_list


    # # DEADBEEF - lake_hru is not used in PRMS 3.0.X or gsflow
    # #   It is used in PRMS 4.0 though
    # # lake_hru parameter
    # logging.info('\nCalculating LAKE_HRU from HRU_ID for all lake HRU\'s')
    # param_names['lake_hru'] = 'lake_hru'
    # param_dimen_counts['lake_hru'] = 1
    # param_dimen_names['lake_hru'] = ['nlake']
    # param_value_counts['lake_hru'] = dimen_sizes['nlake']
    # param_types['lake_hru'] = 1
    # lake_hru_id_list = [
    #    row[1] for row in arcpy.da.SearchCursor(
    #        hru.polygon_path, (hru.type_field, hru.id_field))
    #    if int(row[0]) == 2]
    # for i,lake_hru_id in enumerate(sorted(lake_hru_id_list)):
    #    # logging.debug('  {} {}'.format(i, lake_hru_id))
    #    param_values['lake_hru'][i] = lake_hru_id


    # Read in CRT parameters
    logging.info('\nReading CRT parameters')
    with open(crt_parameter_path, 'r') as input_f:
        crt_param_lines = [line.strip() for line in input_f.readlines()]
    input_f.close()
    # Using enumerate iterator to get .next method
    crt_param_enumerate = enumerate(crt_param_lines)
    for crt_param_line in crt_param_enumerate:
        if crt_param_line[1] == break_str:
            # Skip break string
            crt_param_line = crt_param_enumerate.next()
            # Read parameter name and get next line
            param_name = crt_param_line[1]
            param_names[param_name] = param_name
            crt_param_line = crt_param_enumerate.next()
            # Read dimension count and get next line
            param_dimen_counts[param_name] = int(crt_param_line[1])
            crt_param_line = crt_param_enumerate.next()
            # For each dimen (based on count) read in dimension name
            param_dimen_names[param_name] = []
            for dimen_i in range(param_dimen_counts[param_name]):
                param_dimen_names[param_name].append(crt_param_line[1])
                crt_param_line = crt_param_enumerate.next()
            # Read in number of parameter values
            param_value_counts[param_name] = int(crt_param_line[1])
            crt_param_line = crt_param_enumerate.next()
            # Read in parameter type
            param_types[param_name] = int(crt_param_line[1])
            # Read in parameter values
            # Get next in loop is place intentionally
            # Placing  after getting the value causes it to skip next break
            for i in range(param_value_counts[param_name]):
                crt_param_line = crt_param_enumerate.next()
                if param_types[param_name] == 1:
                    param_values[param_name][i] = int(crt_param_line[1])
                if param_types[param_name] in [2, 3]:
                    param_values[param_name][i] = float(crt_param_line[1])
                if param_types[param_name] == 4:
                    param_values[param_name][i] = crt_param_line[1]

    # Read in CRT groundwater parameters
    logging.info('Reading CRT groundwater parameters')
    with open(crt_gw_parameter_path, 'r') as input_f:
        crt_param_lines = [line.strip() for line in input_f.readlines()]
    input_f.close()
    # Using enumerate iterator to get .next method
    crt_param_enumerate = enumerate(crt_param_lines)
    for crt_param_line in crt_param_enumerate:
        if crt_param_line[1] == break_str:
            # Skip break string
            crt_param_line = crt_param_enumerate.next()
            # Read parameter name and get next line
            param_name = crt_param_line[1]
            param_names[param_name] = param_name
            crt_param_line = crt_param_enumerate.next()
            # Read dimension count and get next line
            param_dimen_counts[param_name] = int(crt_param_line[1])
            crt_param_line = crt_param_enumerate.next()
            # For each dimen (based on count) read in dimension name
            param_dimen_names[param_name] = []
            for dimen_i in range(param_dimen_counts[param_name]):
                param_dimen_names[param_name].append(crt_param_line[1])
                crt_param_line = crt_param_enumerate.next()
            # Read in number of parameter values
            param_value_counts[param_name] = int(crt_param_line[1])
            crt_param_line = crt_param_enumerate.next()
            # Read in parameter type
            param_types[param_name] = int(crt_param_line[1])
            # Read in parameter values
            # Get next in loop is place intentionally
            # Placing  after getting the value causes it to skip next break
            for i in range(param_value_counts[param_name]):
                crt_param_line = crt_param_enumerate.next()
                if param_types[param_name] == 1:
                    param_values[param_name][i] = int(crt_param_line[1])
                if param_types[param_name] in [2, 3]:
                    param_values[param_name][i] = float(crt_param_line[1])
                if param_types[param_name] == 4:
                    param_values[param_name][i] = crt_param_line[1]
    del crt_param_enumerate, crt_param_lines

    # Add lake HRU's to groundwater cascades
    logging.info('Modifying CRT groundwater parameters for all lake HRU\'s')
    logging.info('  gw_up_id = HRU_ID (lake)')
    logging.info('  gw_down_id = 0')
    # logging.info('  gw_strmseg_down_id = OUTSEG')
    logging.info('  gw_strmseg_down_id = 2')
    logging.info('  gw_pct_up = 1')
    field_list = [hru.type_field, hru.id_field, hru.outseg_field,
                 hru.outflow_field]
    lake_hru_id_dict = dict([
       (row[1], row[2])
       for row in arcpy.da.SearchCursor(hru.polygon_path, field_list)
       if int(row[0]) == 2 and int(row[3]) == 0])
    for lake_hru_id, outseg in sorted(lake_hru_id_dict.items()):
       # if lake_hru_id == 9128:
           # print lake_hru_id, outseg
       # raw_input('ENTER')
       i = dimen_sizes['ncascdgw']
       dimen_sizes['ncascdgw'] += 1
       param_values['gw_up_id'][i] = lake_hru_id
       param_values['gw_down_id'][i] = 0
       # DEADBEEF - PRMS didn't like when set to OUTSEG, but 2 worked?
       # param_values['gw_strmseg_down_id'][i] = outseg
       param_values['gw_strmseg_down_id'][i] = 2
       # DEADBEEF - Trying 0
       # param_values['gw_strmseg_down_id'][i] = 0
       param_values['gw_pct_up'][i] = 1.00
       # print param_values['gw_up_id'][i]
       # print param_values['gw_down_id'][i]
       # print param_values['gw_strmseg_down_id'][i]
       # print param_values['gw_pct_up'][i]
    param_value_counts['gw_up_id'] = int(dimen_sizes['ncascdgw'])
    param_value_counts['gw_down_id'] = int(dimen_sizes['ncascdgw'])
    param_value_counts['gw_strmseg_down_id'] = int(dimen_sizes['ncascdgw'])
    param_value_counts['gw_pct_up'] = int(dimen_sizes['ncascdgw'])
    logging.info('  ncascade = {}'.format(dimen_sizes['ncascade']))
    logging.info('  ncascdgw = {}'.format(dimen_sizes['ncascdgw']))
    # raw_input('ENTER')


    # DEADBEEF
    # Override -999 values
    # logging.info('\nChanging SOIL_MOIST_MAX nodata (-999) to 2')
    # for i,v in param_values['soil_moist_max'].items():
    #    if v == -999: param_values['soil_moist_max'][i] = 2
    # logging.info('Changing SOIL_RECHR_MAX nodata (-999) to 1')
    # for i,v in param_values['soil_rechr_max'].items():
    #    if v == -999: param_values['soil_rechr_max'][i] = 1
    # logging.info('Changing SAT_THRESHOLD nodata (-999) to 4')
    # for i,v in param_values['sat_threshold'].items():
    #    if v == -999: param_values['sat_threshold'][i] = 4

    # Override negative values
    # logging.info('Changing negative SSR2GW_RATE (< 0) to 0.1 (PRMS default)')
    # for i,v in param_values['ssr2gw_rate'].items():
    #    if v < 0: param_values['ssr2gw_rate'][i] = 0.1
    # raw_input('ENTER')


    # Write dimensions/parameters to PRMS param file
    logging.info('\nWriting parameter file(s)')
    prms_parameter_paths = sorted(list(set(
        param_files.values() + dimen_files.values())))

    for prms_parameter_path in prms_parameter_paths:
        logging.info('{}'.format(prms_parameter_path))
        if os.path.isfile(prms_parameter_path):
            logging.debug('  Removing existing file')
            os.remove(prms_parameter_path)
        # Get parameters and dimensions for each file
        param_name_list = sorted([
            p_name for p_name, p_file in param_files.items()
            if p_file == prms_parameter_path])
        dimen_name_list = sorted([
            d_name for d_name, d_file in dimen_files.items()
            if d_file == prms_parameter_path])

        with open(prms_parameter_path, 'w') as output_f:
            output_f.write(file_header_str + '\n')

            # Write dimensions
            if dimen_name_list:
                output_f.write(dimen_header_str + '\n')
                logging.debug('  Set dimensions')
            for dimen_name in dimen_name_list:
                try:
                    dimen_size = dimen_sizes[dimen_name]
                except KeyError:
                    continue
                if (type(dimen_size) is str and
                        dimen_size.lower() in ['calculated']):
                    logging.debug('    Dimension {} not calculated'.format())
                    continue
                logging.debug('    {}'.format(dimen_name))
                output_f.write(break_str + '\n')
                output_f.write(dimen_name + '\n')
                output_f.write(str(dimen_size) + '\n')

            # Then write set parameters
            if param_name_list:
                output_f.write(param_header_str + '\n')
                logging.debug('  Set parameters')
            for param_name in param_name_list:
                if param_name not in param_values.keys():
                    # logging.debug(param_name)
                    continue
                logging.debug('    {}'.format(param_name))

                output_f.write(break_str + '\n')
                output_f.write('{}\n'.format(param_name))
                output_f.write('{}\n'.format(
                    param_dimen_counts[param_name]))
                for dimen_name in param_dimen_names[param_name]:
                    output_f.write(dimen_name + '\n')
                output_f.write(str(param_value_counts[param_name]) + '\n')
                param_type = param_types[param_name]
                output_f.write(str(param_type) + '\n')

                # Get list of values sorted by parameter name
                sorted_param_values = [
                    v for i, v in sorted(param_values[param_name].items())]

                # If dimension is "nhru", write values as an array.
                # Write blocks of values for each row
                if ('nhru' in param_dimen_names[param_name] and
                        not param_column_flag):
                    n = ncol
                else:
                    n = 1

                for i in range(0, len(sorted_param_values), n):
                    values_str = ' '.join([
                        param_formats[param_type].format(v)
                        for v in sorted_param_values[i:i + n]])
                    output_f.write(values_str + '\n')

        # Close file
        output_f.close()


def prod(iterable):
    return reduce(operator.mul, iterable, 1)


def isfloat(s):
    """"""
    try:
        float(s)
        return True
    except (ValueError, TypeError):
        return False

# class dimension():
#    def __init__(self, i, data_lines):
#        self.i = i
#        self.NAME = data_lines[0]
#        self.SIZE = int(data_lines[1])
#        print self.NAME, self.SIZE

# class parameter():
#    # type_dict = dict()
#    # type_dict[1] = 'INTEGER'
#    # type_dict[2] = 'FLOAT'
#    # type_dict[3] = 'DOUBLE'
#    # type_dict[4] = 'STRING'
#    def __init__(self, i, data_lines):
#        self.i = i
#        try: self.NAME = data_lines[0]
#        except ValueError: self.NAME = data_lines[0]
#        # # There can be multiple dimensions
#        self.NO_DIMENSIONS = int(data_lines[1])
#        self.DIMENSION_NAMES = []
#        for i in range(self.NO_DIMENSIONS):
#            self.DIMENSION_NAMES.append(data_lines[2+i])
#        self.N_VALUES = int(data_lines[3+i])
#        self.TYPE = data_lines[4+i]
#        self.VALUE = data_lines[5+i:]
#        print self.NAME, self.NO_DIMENSIONS,
#        print self.DIMENSION_NAMES, self.N_VALUES, self.TYPE


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='PRMS Template Fill',
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

    # Fill PRMS Parameter Template File
    prms_template_fill(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
