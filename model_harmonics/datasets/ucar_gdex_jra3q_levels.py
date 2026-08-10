#!/usr/bin/env python
"""
ucar_gdex_jra3q_levels.py
Written by Tyler Sutterley (08/2026)

Downloads JRA-3Q monthly model-level analysis products provided by the
    NCAR/UCAR Geoscience Data Exchange (GDEX)

JRA-3Q: Japanese Reanalysis for Three Quarters of a Century Monthly Statistics
    https://gdex.ucar.edu/datasets/d640002

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: Working data directory
    --variable X, -V X: Variable to download
        * hgt-hyb: geopotential height
        * spfh-hyb: specific humidity
        * tmp-hyb: temperature
        * ugrd-hyb: U-component of wind
        * vgrd-hyb: V-component of wind
        * vvel-hyb: vertical velocity (pressure)
    -Y X, --year X: Years to download
    -t X, --timeout X: Timeout in seconds for blocking operations
    -l, --log: Output log of files downloaded
    -M X, --mode=X: Permission mode of directories and files downloaded

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 08/2026: change formatting and output of file logs
    Written 07/2026
"""

from __future__ import print_function

import sys
import os
import re
import copy
import time
import logging
import pathlib
import argparse
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc


# PURPOSE: sync local JRA-3Q files with UCAR/NCAR GDEX server
def ucar_gdex_download(
    base_dir,
    VARIABLE,
    YEARS=None,
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
        # output to log file
        # format: UCAR_GDEX_JRA-3Q_2002-04-01.log
        today = time.strftime('%Y-%m-%d', time.localtime())
        output_logfile = f'UCAR_GDEX_JRA-3Q_{today}.log'
        LOGFILE = gravtk.utilities.get_cache_path(output_logfile)
        # create a unique log and open the log file
        fid1 = gravtk.utilities.create_unique_file(LOGFILE, mode='x')
        # build logger for outputting to log file
        logger = gravtk.utilities.build_logger(
            __name__,
            level=logging.INFO,
            stream=fid1,
            format='%(message)s (%(levelname)s)',
        )
        logger.info(f'UCAR JRA-3Q Sync Log ({today})')
        logger.info(f'Filename: {pathlib.Path(sys.argv[0]).name}')
    else:
        # build logger for standard output (terminal output)
        fid1 = sys.stdout
        logger = gravtk.utilities.build_logger(
            __name__, level=logging.INFO, stream=fid1
        )

    # UCAR GDEX host
    HOST = 'https://gdex.ucar.edu/'
    dataset_id = 'd640002'
    product_id = '3'
    PRODUCT = 'anl_mdl'
    TIMENAME = 'time'
    INTERFACE = 'hybrid_level'
    LEVELNAME = 'hybrid_half_level'
    LATNAME = 'lat'
    LONNAME = 'lon'
    # netCDF4 variable name for the requested product (monthly means)
    VARNAME = f'{VARIABLE}-an-gauss-mn'
    POINTS = 'original_number_of_grid_points_per_latitude_circle'
    WEIGHT = 'weight'
    ANAME, BNAME = ('a_hybrid_half_level', 'b_hybrid_half_level')
    AINTERFACE, BINTERFACE = ('a_hybrid_level', 'b_hybrid_level')
    # build a file filter for the product and variable
    query_params = {}
    query_params['filter_wfile'] = VARNAME
    kwargs = '?' + mdlhmc.utilities.urlencode(query_params)
    # regular expression pattern for extracting datetimes
    rx = re.compile(rf'{VARNAME}.(\d{{4}})(\d{{2}})(\d{{2}})(\d{{2}})', re.I)

    # dictionary defining output structure
    struct = dict(
        dimensions=(TIMENAME, INTERFACE, LEVELNAME, LATNAME, LONNAME),
        variables={
            VARNAME: (TIMENAME, LEVELNAME, LATNAME, LONNAME),
            POINTS: (LATNAME,),
            WEIGHT: (LATNAME,),
            ANAME: (LEVELNAME,),
            BNAME: (LEVELNAME,),
            AINTERFACE: (INTERFACE,),
            BINTERFACE: (INTERFACE,),
        },
    )

    # find data directories for year
    dirs, _ = mdlhmc.utilities.ucar_list(
        [HOST, 'datasets', dataset_id, 'filelist', product_id],
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
            [HOST, *dirparts, kwargs],
            timeout=TIMEOUT,
            pattern=VARNAME,
            sort=True,
        )
        # skip to next year if no files were found
        if not any(cols):
            continue
        # for each file in the list of files
        for colname, collastmod in zip(cols, mods):
            # output filename
            YY, MM, DD, HH = rx.findall(colname).pop()
            filename = f'jra3q.{PRODUCT}.{VARNAME}.{YY}{MM}.nc'
            output_file = DIRECTORY.joinpath(filename)
            # check if local file exists and is valid
            if mdlhmc.spatial.validate_netCDF4(output_file, struct):
                continue
            # download file from UCAR GDEX server
            logger.info(colname)
            # create and submit request for file
            mdlhmc.utilities.from_http(
                colname,
                timeout=TIMEOUT,
                local=output_file,
                verbose=True,
                fid=fid1,
                mode=MODE,
            )
            # keep remote modification time of file and local access time
            os.utime(output_file, (output_file.stat().st_atime, collastmod))
            # set permissions mode to MODE
            output_file.chmod(mode=MODE)

    # close log file and set permissions level to MODE
    if LOG:
        fid1.close()
        os.chmod(fid1.name, mode=MODE)


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Downloads JRA-3Q model-level analysis products
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
    # JRA-3Q model-level analysis variable
    choices = [
        'hgt-hyb',
        'spfh-hyb',
        'tmp-hyb',
        'ugrd-hyb',
        'vgrd-hyb',
        'vvel-hyb',
    ]
    parser.add_argument(
        '--variable',
        '-v',
        type=str,
        nargs='+',
        metavar='VARIABLE',
        choices=choices,
        default=['spfh-hyb', 'tmp-hyb'],
        help='JRA-3Q model-level analysis variable',
    )
    # model years to download
    parser.add_argument(
        '--year',
        '-Y',
        type=str,
        nargs='+',
        help='Years to download',
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
    for VARIABLE in args.variable:
        ucar_gdex_download(
            args.directory,
            VARIABLE,
            YEARS=args.year,
            TIMEOUT=args.timeout,
            LOG=args.log,
            MODE=args.mode,
        )


# run main program
if __name__ == '__main__':
    main()
