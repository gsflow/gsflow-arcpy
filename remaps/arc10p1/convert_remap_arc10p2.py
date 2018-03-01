#--------------------------------
# Name:         convert_remap_arc10p2.py
# Purpose:      PRMS Remap Modfiy
# Author:       Charles Morton
# Created       2015-05-18
# Python:       2.7
#--------------------------------

import datetime as dt
import logging
import os
import re
import sys

################################################################################

def prms_remap_modify(remap_folder):
    """Calculate PRMS Remap Modify

    Args:
        remap_folder: folder path
    Returns:
        None
    """

    try:
        logging.info('\nPRMS Remap Modify')

        ## Build output folder if necessary
        input_ws = os.path.join(remap_folder, 'arc10p1')    
        output_ws = os.path.join(remap_folder, 'arc10p2')    
        if not os.path.isdir(output_ws):
            os.mkdir(output_ws)

        ## Remove comments from ASCII remap files
        for remap_name in os.listdir(input_ws):
            remap_input_path = os.path.join(input_ws, remap_name)
            remap_output_path = os.path.join(output_ws, remap_name)
            if not os.path.isfile(remap_input_path):
                continue
            if not remap_input_path.lower().endswith('.rmp'):
                continue           
            logging.debug('  Modifying: {0}'.format(remap_name))
            with open(remap_input_path, 'r') as remap_f:
                input_lines = remap_f.readlines()

            ## Separate remap values and comments
            remap_lines, comment_lines = [], []
            for i, line in enumerate(input_lines):
                ## Remove newline characters and extra leading/trailing whitespace
                line_split = [x.strip() for x in line.strip().split('/*')]
                ## Comments will result in line_split have two items
                if len(line_split) == 2:
                    comment_lines.append('# {0} - {1}'.format(
                        line_split[0].split(':')[0].strip(), line_split[1]))
                if line_split:
                    remap_lines.append(' : '.join([x.strip() for x in line_split[0].split(':')]))

            ## First write the remap values,
            ## Then write the comments (if there are any)
            with open(remap_output_path, 'w') as remap_f:
                for i, line in enumerate(remap_lines):
                    remap_f.write(line)     
                    ## Don't write newline character on last line
                    ## This causes an error in ArcGIS 10.2.2
                    if (i+1) == len(remap_lines) and not comment_lines:
                        break
                    remap_f.write('\n')

                for i, line in enumerate(comment_lines):
                    ## Replace California with abbreviation
                    ## Limit comments to 80 characters
                    remap_f.write(line.replace('California', 'CA')[:79])
                    ## Don't write newline character on last line
                    ## This causes an error in ArcGIS 10.2.2
                    if (i+1) == len(comment_lines):
                        break
                    remap_f.write('\n')

    except:
        logging.exception('Unhandled Exception Error\n\n')
        raw_input('ENTER')

################################################################################

if __name__ == '__main__':
    workspace = os.getcwd()

    ## Create Basic Logger
    logging.basicConfig(level=logging.DEBUG, format='%(message)s')

    ## Run Information
    logging.info('\n{0}'.format('#'*80))
    log_f = '{0:<20s} {1}'
    logging.info(log_f.format(
        'Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', workspace))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))

    ## Calculate PRMS Remap Modify
    prms_remap_modify(workspace)
