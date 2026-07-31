#!/usr/bin/env python
"""
ucar_gdex_jra55_surface.py
Written by Tyler Sutterley (07/2026)

Downloads 6-Hourly 1.25 Degree Surface Analysis Fields provided by the
    NCAR/UCAR Geoscience Data Exchange (GDEX)

JRA-55: Japanese 55-year Reanalysis, Daily 3-Hourly and 6-Hourly Data
    https://gdex.ucar.edu/datasets/d628000

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: Working data directory
    -Y X, --year X: Years to download
    -I, --invariant: Retrieve the model invariant parameters
    -t X, --timeout X: Timeout in seconds for blocking operations
    -l, --log: Output log of files downloaded
    -M X, --mode=X: Permission mode of directories and files downloaded

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    pygrib: Python interface for reading and writing GRIB data
        https://pypi.python.org/pypi/pygrib

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    utilities.py: download and management utilities for files
    spatial.py: spatial data class for reading, writing and processing data

UPDATE HISTORY:
    Updated 07/2026: converted to use the NCAR/UCAR Geoscience Data Exchange
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: debug-level logging for lines within download file
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: adjust ordering of file line search in csh file
    Updated 05/2021: added option for connection timeout (in seconds)
        use try/except for retrieving netrc credentials
        define int/float precision to prevent deprecation warning
    Updated 04/2021: set a default netrc file and check access
        default credentials from environmental variables
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 12/2020: using time, spatial and utilities modules
    Updated 08/2019: option GZIP to specify if data in links list is compressed
    Updated 06/2018: using python3 compatible octal, input and urllib
    Written 03/2018
"""

from __future__ import print_function

import sys
import os
import re
import time
import logging
import pathlib
import argparse
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc


# PURPOSE: sync local JRA-55 files with UCAR/NCAR GDEX server
def ucar_gdex_download(
    base_dir,
    YEARS=None,
    INVARIANT=False,
    TIMEOUT=None,
    LOG=False,
    MODE=None,
):
    # directory setup
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    DIRECTORY = base_dir.joinpath('JRA-55')
    # check if log directory exists and recursively create if not
    DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # create log file with list of synchronized files (or print to terminal)
    if LOG:
        # format: UCAR_GDEX_JRA-55_2002-04-01.log
        today = time.strftime('%Y-%m-%d', time.localtime())
        output_logfile = f'UCAR_GDEX_JRA-55_{today}.log'
        LOGFILE = DIRECTORY.joinpath(output_logfile)
        fid = LOGFILE.open(mode='w', encoding='utf8')
        logging.basicConfig(stream=fid, level=logging.INFO)
        logging.info(f'UCAR JRA-55 Sync Log ({today})')
    else:
        # standard output (terminal output)
        fid = sys.stdout
        logging.basicConfig(stream=fid, level=logging.INFO)

    # UCAR GDEX host
    HOST = 'https://gdex.ucar.edu/'
    # interpolated surface pressure products
    PRODUCT = 'anl_surf125'
    VARIABLE = '001_pres'
    dataset_id = 'd628000'
    product_id = '13'
    invariant_id = '32'

    # build a file filter for the product and variable
    query_params = {}
    query_params['filter_wfile'] = VARIABLE
    kwargs = product_id + '?' + gravtk.utilities.urlencode(query_params)
    # find data files for product and variable
    cols, mods = mdlhmc.utilities.ucar_list(
        [HOST, 'datasets', dataset_id, 'filelist', kwargs],
        timeout=TIMEOUT,
        pattern=VARIABLE,
        sort=True,
    )

    # reduce to files for the requested years
    if YEARS:
        # compile regular expression operators
        pattern = rf'{VARIABLE}\.(' + r'|'.join(YEARS) + r')'
        i = [i for i, f in enumerate(cols) if re.search(pattern, f)]
        # reduce list of column names and last modified times
        cols = [cols[indice] for indice in i]
        mods = [mods[indice] for indice in i]

    # for each file in the list of files
    for colname, collastmod in zip(cols, mods):
        # download file from UCAR GDEX server
        logging.debug(colname)
        # create and submit request for file
        buffer = mdlhmc.utilities.from_http(
            colname,
            timeout=TIMEOUT,
            local=None,
            verbose=True,
            fid=fid,
            mode=MODE,
        )
        # open remote file with pygrib
        dinput = mdlhmc.spatial.grib().from_file(buffer)
        varname = dinput.attributes['data']['cfVarName']
        # convert GRACE months to calendar month indices
        iMON = np.mod(dinput.month - 1, 12)
        # for each month
        for mm1 in range(12):
            # find indices for the month
            tt = np.flatnonzero(iMON == mm1)
            if not np.any(tt):
                continue
            # calculate monthly mean
            mean = dinput.index(tt).mean()
            # copy attributes from input file
            mean.attributes.update(dinput.attributes)
            mean.attributes[varname] = dinput.attributes['data']
            # convert time to strings
            (datetime,) = gravtk.time.to_string(
                mean.time,
                mean.attributes['time']['units'],
                strftime='%Y%m',
            )
            # output filename
            filename = f'{PRODUCT}.{VARIABLE}.{datetime}.nc'
            local_file = DIRECTORY.joinpath(filename)
            logging.info(f'\t{local_file}')
            # write mean to netCDF4 file
            mean.to_netCDF4(
                local_file,
                varname=varname,
                attributes=mean.attributes,
                clobber=True,
            )
            # keep remote modification time of file
            os.utime(local_file, (local_file.stat().st_atime, collastmod))

    # if retrieving the model invariant parameters
    if INVARIANT:
        # build a file filter for the product and variable
        query_params['filter_wfile'] = 'll125.*.2000'
        kwargs = invariant_id + '?' + gravtk.utilities.urlencode(query_params)
        # find invariant data files
        cols, mods = mdlhmc.utilities.ucar_list(
            [HOST, 'datasets', dataset_id, 'filelist', kwargs],
            timeout=TIMEOUT,
        )
        # for each file in the list of files
        for colname, collastmod in zip(cols, mods):
            # download file from UCAR GDEX server
            logging.debug(colname)
            # create and submit request for file
            buffer = mdlhmc.utilities.from_http(
                colname,
                timeout=TIMEOUT,
                local=None,
                verbose=True,
                fid=fid,
                mode=MODE,
            )
            # open remote file with pygrib
            dinput = mdlhmc.spatial.grib().from_file(buffer)
            varname = dinput.attributes['data']['cfVarName']
            dinput.attributes[varname] = dinput.attributes['data']
            # output filename
            fileparts = mdlhmc.utilities.url_split(colname)
            output = DIRECTORY.joinpath(f'{fileparts[-1]}.nc')
            logging.info(f'\t{output}')
            # write invariant data to netCDF4 file
            dinput.to_netCDF4(
                output,
                varname=varname,
                attributes=dinput.attributes,
                clobber=True,
            )
            # keep remote modification time of file
            os.utime(output, (output.stat().st_atime, collastmod))

    # close log file and set permissions level to MODE
    if LOG:
        fid.close()
        LOGFILE.chmod(mode=MODE)


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Downloads JRA-55 surface analysis products
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
    # UCAR_GDEX_JRA-55_2002-04-01.log
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

    # download JRA-55 1.25 degree products
    ucar_gdex_download(
        args.directory,
        YEARS=args.year,
        INVARIANT=args.invariant,
        TIMEOUT=args.timeout,
        LOG=args.log,
        MODE=args.mode,
    )


# run main program
if __name__ == '__main__':
    main()
