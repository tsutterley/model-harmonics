#!/usr/bin/env python
"""
ucar_gdex_cfsr_surface.py
Written by Tyler Sutterley (07/2026)

Downloads monthly NCEP-CFSR products provided by the
    NCAR/UCAR Geoscience Data Exchange (GDEX)

NCEP Climate Forecast System Reanalysis (CFSR) Monthly Products
    https://gdex.ucar.edu/datasets/d093002/
NCEP Climate Forecast System Version 2 (CFSv2) Monthly Products
    https://gdex.ucar.edu/datasets/d094002/

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: Working data directory
    -S X, --spacing X: Data spatial resolution to download
        - 'h': high resolution
        - 'l': low resolution
    -Y X, --year X: Years to download
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
    Updated 05/2021: added option for connection timeout (in seconds)
        use try/except for retrieving netrc credentials
        define int/float precision to prevent deprecation warning
    Updated 04/2021: set a default netrc file and check access
        default credentials from environmental variables
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 12/2020: using time, spatial and utilities modules
    Written 08/2019
"""

from __future__ import print_function

import sys
import os
import io
import re
import time
import logging
import pathlib
import tarfile
import argparse
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc


# PURPOSE: sync local NCEP-CFSR files with UCAR/NCAR GDEX server
def ucar_gdex_download(
    base_dir,
    resolution,
    YEARS=None,
    INVARIANT=False,
    TIMEOUT=None,
    LOG=False,
    MODE=None,
):
    # directory setup
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    DIRECTORY = base_dir.joinpath('NCEP-CFSR')
    # check if log directory exists and recursively create if not
    DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # create log file with list of synchronized files (or print to terminal)
    if LOG:
        # format: UCAR_GDEX_NCEP-CFSR_2002-04-01.log
        today = time.strftime('%Y-%m-%d', time.localtime())
        output_logfile = f'UCAR_GDEX_NCEP-CFSR_{today}.log'
        LOGFILE = DIRECTORY.joinpath(output_logfile)
        fid = LOGFILE.open(mode='w', encoding='utf8')
        logging.basicConfig(stream=fid, level=logging.INFO)
        logging.info(f'UCAR NCEP-CFSR Sync Log ({today})')
        logging.info('PRODUCT: NCEP-DOE-2')
    else:
        # standard output (terminal output)
        fid = sys.stdout
        logging.basicConfig(stream=fid, level=logging.INFO)

    # UCAR GDEX host
    HOST = 'https://gdex.ucar.edu/'

    # NCEP-CFSR (v1)
    dataset_id = 'd093002'
    product_id = '2'
    # surface level products for NCEP-CFSR
    surface_resolution = 'a' if resolution == 'h' else 'l'
    VARIABLE = f'spl{surface_resolution}nl'
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
        pattern = r'|'.join(YEARS)
        i = [i for i, f in enumerate(cols) if re.search(pattern, f)]
        # reduce list of column names and last modified times
        cols = [cols[indice] for indice in i]
        mods = [mods[indice] for indice in i]

    # for each file in the list of files
    for colname, collastmod in zip(cols, mods):
        # download file from UCAR GDEX server
        logging.debug(colname)
        # create and submit request for file
        response = mdlhmc.utilities.from_http(
            colname,
            timeout=TIMEOUT,
            local=None,
            verbose=True,
            fid=fid,
            mode=MODE,
        )
        # open remote file with tarfile
        tar = tarfile.open(fileobj=response, mode='r')
        for member in tar.getmembers():
            # file information for member
            fileinfo = member.get_info()
            # open buffer file with pygrib
            tar_contents = tar.extractfile(member).read()
            buffer = io.BytesIO(tar_contents)
            dinput = mdlhmc.spatial.grib().from_file(
                buffer, varname='Time-mean surface pressure'
            )
            varname = dinput.attributes['data']['cfVarName']
            dinput.attributes[varname] = dinput.attributes['data']
            # convert time to strings
            (datetime,) = gravtk.time.to_string(
                dinput.time,
                dinput.attributes['time']['units'],
                strftime='%Y%m',
            )
            # output filename
            filename = f'{VARIABLE}.gdas.{datetime}.nc'
            output = DIRECTORY.joinpath(filename)
            logging.info(f'\t{output}')
            # write mean to netCDF4 file
            dinput.to_netCDF4(
                output,
                varname=varname,
                attributes=dinput.attributes,
                clobber=True,
            )
            # keep tar modification time of file
            os.utime(output, (output.stat().st_atime, fileinfo['mtime']))

    # NCEP-CFSv2 (v2)
    dataset_id = 'd094002'
    product_id = '2'
    # isentropic surface pressure products for NCEP-CFSv2
    VARIABLE = f'ipv{resolution}'
    # output product after reducing to surface pressure
    PRODUCT = f'spl{surface_resolution}nl'
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
        pattern = r'|'.join(YEARS)
        i = [i for i, f in enumerate(cols) if re.search(pattern, f)]
        # reduce list of column names and last modified times
        cols = [cols[indice] for indice in i]
        mods = [mods[indice] for indice in i]

    # for each file in the list of files
    for colname, collastmod in zip(cols, mods):
        # download file from UCAR GDEX server
        logging.debug(colname)
        # create and submit request for file
        response = mdlhmc.utilities.from_http(
            colname,
            timeout=TIMEOUT,
            local=None,
            verbose=True,
            fid=fid,
            mode=MODE,
        )
        # open remote file with tarfile
        tar = tarfile.open(fileobj=response, mode='r')
        # find analysis file in tarfile
        (member,) = [
            m for m in tar.getmembers() if re.search(rf'{VARIABLE}(nl)', m.name)
        ]
        # file information for member
        fileinfo = member.get_info()
        # open buffer file with pygrib
        tar_contents = tar.extractfile(member).read()
        buffer = io.BytesIO(tar_contents)
        dinput = mdlhmc.spatial.grib().from_file(
            buffer, varname='Time-mean surface pressure'
        )
        varname = dinput.attributes['data']['cfVarName']
        dinput.attributes[varname] = dinput.attributes['data']
        # convert time to strings
        (datetime,) = gravtk.time.to_string(
            dinput.time,
            dinput.attributes['time']['units'],
            strftime='%Y%m',
        )
        # output filename
        filename = f'{PRODUCT}.gdas.{datetime}.nc'
        output = DIRECTORY.joinpath(filename)
        logging.info(f'\t{output}')
        # write mean to netCDF4 file
        dinput.to_netCDF4(
            output,
            varname=varname,
            attributes=dinput.attributes,
            clobber=True,
        )
        # keep tar modification time of file
        os.utime(output, (output.stat().st_atime, fileinfo['mtime']))

    # if retrieving the model invariant parameters
    if INVARIANT:
        # NCEP-CFSR (v1)
        dataset_id = 'd093002'
        product_id = '3'
        VARIABLE = f'pgb{resolution}01.*.SFC'
        pattern = r'(HGT|LAND)'
        # build a file filter for the product and variable
        query_params['filter_wfile'] = VARIABLE
        kwargs = product_id + '?' + gravtk.utilities.urlencode(query_params)
        # find invariant data files
        cols, mods = mdlhmc.utilities.ucar_list(
            [HOST, 'datasets', dataset_id, 'filelist', kwargs],
            timeout=TIMEOUT,
            pattern=pattern,
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
            # calculate long-term average
            mean = dinput.mean()
            # copy attributes from input file
            mean.attributes.update(dinput.attributes)
            mean.attributes[varname] = dinput.attributes['data']
            # output filename
            (PRODUCT,) = re.findall(pattern, colname, re.I)
            filename = f'{PRODUCT.lower()}.gdas.nc'
            output = DIRECTORY.joinpath(filename)
            logging.info(f'\t{output}')
            # write mean to netCDF4 file
            mean.to_netCDF4(
                output,
                varname=varname,
                attributes=mean.attributes,
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
        description="""Downloads NCEP-CFSR surface analysis products
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
    # data spatial resolution
    parser.add_argument(
        '--spacing',
        '-S',
        type=str,
        choices=['h', 'l'],
        default='h',
        help='Data spatial resolution to download',
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

    # download NCEP-CFSR products
    ucar_gdex_download(
        args.directory,
        args.spacing,
        YEARS=args.year,
        INVARIANT=args.invariant,
        TIMEOUT=args.timeout,
        LOG=args.log,
        MODE=args.mode,
    )


# run main program
if __name__ == '__main__':
    main()
