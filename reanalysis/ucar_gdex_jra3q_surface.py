#!/usr/bin/env python
"""
ucar_gdex_jra3q_surface.py
Written by Tyler Sutterley (05/2023)

Downloads JRA-3Q surface analysis products provided by the
    NCAR/UCAR Geoscience Data Exchange (GDEX)

JRA-3Q: Japanese Reanalysis for Three Quarters of a Century
    https://gdex.ucar.edu/datasets/d640000

CALLING SEQUENCE:
    python ucar_gdex_jra3q_surface.py

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: Working data directory
    --product X, -P X: Product to download
        * anl_surf: surface analysis fields
        * anl_surf125: Interpolated 1.25 degree product
    --variable X, -V X: Variable to download
        * pres-sfc: surface pressure
    -Y X, --year X: Years to download
    -I, --invariant: Retrieve the model invariant parameters
    -t X, --timeout X: Timeout in seconds for blocking operations
    -l, --log: Output log of files downloaded
    -M X, --mode=X: Permission mode of directories and files downloaded

PROGRAM DEPENDENCIES:
    spatial.py: spatial data class for reading, writing and processing data
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Written 07/2025
"""

from __future__ import print_function

import sys
import os
import io
import re
import copy
import time
import logging
import pathlib
import argparse
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc


# PURPOSE: sync local JRA-3Q files with UCAR/NCAR GDEX server
def ucar_gdex_download(
    base_dir,
    PRODUCT,
    VARIABLE,
    YEARS=None,
    INVARIANT=False,
    TIMEOUT=None,
    LOG=False,
    MODE=None,
):
    # directory setup
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    DIRECTORY = base_dir.joinpath('JRA-3Q')
    # check if log directory exists and recursively create if not
    DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # create log file with list of synchronized files (or print to terminal)
    if LOG:
        # format: UCAR_GDEX_JRA-3Q_2002-04-01.log
        today = time.strftime('%Y-%m-%d', time.localtime())
        output_logfile = f'UCAR_GDEX_JRA-3Q_{today}.log'
        LOGFILE = DIRECTORY.joinpath(output_logfile)
        fid = LOGFILE.open(mode='w', encoding='utf8')
        logging.basicConfig(stream=fid, level=logging.INFO)
        logging.info(f'UCAR JRA-3Q Sync Log ({today})')
    else:
        # standard output (terminal output)
        fid = sys.stdout
        logging.basicConfig(stream=fid, level=logging.INFO)

    # UCAR GDEX host
    HOST = 'https://gdex.ucar.edu/'
    # field mapping for different products
    field_mapping = {}
    field_mapping['lon'] = 'lon'
    field_mapping['lat'] = 'lat'
    field_mapping['time'] = 'time'
    # parameters for surface pressure products
    if PRODUCT == 'anl_surf':
        product_id = '7'
        invariant_id = '19'
        field_mapping['data'] = 'pres-sfc-an-gauss'
        field_mapping['points'] = (
            'original_number_of_grid_points_per_latitude_circle'
        )
        field_mapping['weight'] = 'weight'
    elif PRODUCT == 'anl_surf125':
        product_id = '25'
        invariant_id = '41'
        field_mapping['data'] = 'pres-sfc-an-ll125'
    # netCDF4 variable name for the requested product
    varname = field_mapping['data']

    # find data directories for year
    dirs, _ = mdlhmc.utilities.ucar_list(
        [HOST, 'datasets', 'd640000', 'filelist', product_id],
        tdclass='Group Name',
        timeout=TIMEOUT,
    )
    # model years in each directory
    years = [re.findall(rf'(\d+){product_id.zfill(2)}', d).pop() for d in dirs]
    # reduce list of directories to only those for the requested years
    if YEARS is not None:
        directories = copy.copy(dirs)
        dirs = [directories[years.index(y)] for y in YEARS if y in years]
    # for each directory in the (reduced) list of directories
    for d in dirs:
        # find data files for year
        dirparts = mdlhmc.utilities.url_split(d)
        cols, mods = mdlhmc.utilities.ucar_list(
            [HOST, *dirparts],
            timeout=TIMEOUT,
            pattern=VARIABLE,
            sort=True,
        )
        # for each file in the list of files
        for colname, collastmod in zip(cols, mods):
            # download file from UCAR GDEX server
            logging.debug(colname)
            fileparts = mdlhmc.utilities.url_split(colname)
            # create and submit request for file
            buffer = mdlhmc.utilities.from_http(
                colname,
                timeout=TIMEOUT,
                local=None,
                verbose=True,
                fid=fid,
                mode=MODE,
            )
            # open remote file with netCDF4
            dinput = gravtk.spatial().from_netCDF4(
                buffer,
                compression='bytes',
                field_mapping=field_mapping,
                verbose=True,
            )
            # calculate monthly mean
            mean = dinput.transpose(axes=(1, 2, 0)).mean()
            # copy additional fields
            mean.weight = dinput.weight.copy()
            mean.points = dinput.points.copy()
            # copy global attributes from input file
            mean.attributes.update(dinput.attributes)
            # copy attributes for each output field
            for field, key in field_mapping.items():
                mean.attributes[key] = dinput.attributes.get(field)
            # convert time to strings
            (datetime,) = gravtk.time.to_string(
                mean.time,
                mean.attributes['time']['units'],
                strftime='%Y%m',
            )
            # output filename
            filename = f'jra3q.{PRODUCT}.{varname}.{datetime}.nc'
            output = DIRECTORY.joinpath(filename)
            logging.info(f'\t{output}')
            # write mean to netCDF4 file
            mean.to_netCDF4(
                output,
                field_mapping=field_mapping,
                attributes=mean.attributes,
                clobber=True,
            )
            # keep remote modification time of file
            os.utime(output, (output.stat().st_atime, collastmod))

    # if retrieving the model invariant parameters
    if INVARIANT:
        # find invariant data files
        cols, mods = mdlhmc.utilities.ucar_list(
            [HOST, 'datasets', 'd640000', 'filelist', invariant_id],
            timeout=TIMEOUT,
        )
        # for each file in the list of files
        for colname, collastmod in zip(cols, mods):
            # download file from UCAR GDEX server
            logging.debug(colname)
            fileparts = mdlhmc.utilities.url_split(colname)
            # output filename for local storage
            local = DIRECTORY.joinpath(fileparts[-1])
            # Create and submit request for file
            mdlhmc.utilities.from_http(
                colname,
                timeout=TIMEOUT,
                local=local,
                verbose=True,
                fid=fid,
                mode=MODE,
            )
            # keep remote modification time of file
            os.utime(local, (local.stat().st_atime, collastmod))

    # close log file and set permissions level to MODE
    if LOG:
        fid.close()
        LOGFILE.chmod(mode=MODE)


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Downloads JRA-3Q surface analysis products
            provided by the NCAR/UCAR Geoscience Data Exchange (GDEX)
            """
    )
    # command line parameters
    # working data directory
    parser.add_argument(
        '--directory',
        '-D',
        type=pathlib.Path,
        default=pathlib.Path.cwd(),
        help='Working data directory',
    )
    # JRA-3Q surface analysis product
    parser.add_argument(
        '--product',
        '-p',
        type=str,
        default='anl_surf',
        choices=['anl_surf', 'anl_surf125'],
        help='JRA-3Q surface analysis product',
    )
    # JRA-3Q surface analysis variable
    parser.add_argument(
        '--variable',
        '-v',
        type=str,
        default='pres-sfc',
        help='JRA-3Q surface analysis variable',
    )
    # model years to download
    parser.add_argument(
        '--year',
        '-Y',
        type=str,
        nargs='+',
        help='Years to download',
    )
    # retrieve the model invariant parameters
    parser.add_argument(
        '--invariant',
        '-I',
        default=False,
        action='store_true',
        help='Retrieve model invariant parameters',
    )
    # connection timeout
    parser.add_argument(
        '--timeout',
        '-t',
        type=int,
        default=360,
        help='Timeout in seconds for blocking operations',
    )
    # Output log file in form
    # UCAR_GDEX_JRA-3Q_2002-04-01.log
    parser.add_argument(
        '--log',
        '-l',
        default=False,
        action='store_true',
        help='Output log file',
    )
    # permissions mode of the directories and files synced (number in octal)
    parser.add_argument(
        '--mode',
        '-M',
        type=lambda x: int(x, base=8),
        default=0o775,
        help='Permission mode of directories and files synced',
    )
    # return the parser
    return parser


# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args, _ = parser.parse_known_args()

    # download JRA-3Q products
    ucar_gdex_download(
        args.directory,
        args.product,
        args.variable,
        YEARS=args.year,
        INVARIANT=args.invariant,
        TIMEOUT=args.timeout,
        LOG=args.log,
        MODE=args.mode,
    )


# run main program
if __name__ == '__main__':
    main()
