#!/usr/bin/env python
u"""
gesdisc_gldas_sync.py
Written by Tyler Sutterley (05/2023)

Syncs GLDAS monthly datafiles from the Goddard Earth Sciences Data and
    Information Server Center (GES DISC)
    http://ldas.gsfc.nasa.gov/gldas/
    http://disc.sci.gsfc.nasa.gov/hydrology/documentation
    https://hydro1.gesdisc.eosdis.nasa.gov/data/GLDAS_V1/README.GLDAS.pdf
    https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python

Register with NASA Earthdata Login system:
    https://urs.earthdata.nasa.gov

Add "NASA GESDISC DATA ARCHIVE" to Earthdata Applications:
    https://urs.earthdata.nasa.gov/approve_app?client_id=e2WVk8Pw6weeLUKZYOxvTQ

CALLING SEQUENCE:
    python gesdisc_gldas_sync.py --user <username> CLSM NOAH VIC
    where <username> is your NASA Earthdata username

INPUTS:
    CLM: GLDAS Common Land Model (CLM)
    CLSM: GLDAS Catchment Land Surface Model (CLSM)
    MOS: GLDAS Mosaic model
    NOAH: GLDAS Noah model
    VIC: GLDAS Variable Infiltration Capacity (VIC) model

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user X: username for NASA Earthdata Login
    -W X, --password X: password for NASA Earthdata Login
    -N X, --netrc X: path to .netrc file for authentication
    -D X, --directory X: working data directory
    -Y X, --year X: years to sync
    -S X, --spacing X: spatial resolution of models to sync
        10: 1.0 degrees latitude/longitude
        025: 0.25 degrees latitude/longitude
    -T X, --temporal X: temporal resolution of models to sync
        M: Monthly
        3H: 3-hourly
    -v X, --version X: GLDAS model version to sync
    -e, --early: Sync GLDAS early products
    -t X, --timeout X: Timeout in seconds for blocking operations
    --log: output log of files downloaded
    --list: print files to be transferred, but do not execute transfer
    --clobber: Overwrite existing data in transfer
    -M X, --mode X: permissions mode of the directories and files synced

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    lxml: Pythonic XML and HTML processing library using libxml2/libxslt
        https://lxml.de/
        https://github.com/lxml/lxml
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 06/2022: use CMR queries to find model granules
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 10/2021: using python logging for handling verbose output
    Updated 06/2021: fixed timeout for monthly GLDAS download
    Updated 05/2021: added option for connection timeout (in seconds)
        use try/except for retrieving netrc credentials
    Updated 04/2021: set a default netrc file and check access
        default credentials from environmental variables
    Updated 03/2021: automatically update years to run based on current time
    Updated 01/2021: moved gesdisc_list to utilities module
    Updated 10/2020: use argparse to set command line parameters
    Updated 09/2020: using utilities program to build opener
        use gesdisc_list command to list remote directories
    Updated 08/2020: flake8 compatible regular expression strings
    Updated 06/2020: added netrc option for alternative authentication method
    Updated 04/2020: added option to download the GLDAS Early Products (EP)
    Updated 09/2019: added ssl context to urlopen headers
    Updated 06/2018: using python3 compatible octal, input and urllib
    Updated 03/2018: added option --version to set the GLDAS model version
    Updated 08/2017: use raw_input() to enter NASA Earthdata credentials rather
        than exiting with error
    Updated 05/2017: added options to change the spatial and temporal resolution
        added exception if NASA Earthdata credentials weren't entered
        using os.makedirs to recursively create directories
        using getpass to enter server password securely (remove --password)
    Updated 04/2017: using lxml to parse HTML for files and modification dates
        minor changes to check_connection function to parallel other programs
    Written 02/2017
"""
from __future__ import print_function

import sys
import os
import re
import time
import netrc
import shutil
import getpass
import logging
import pathlib
import argparse
import builtins
import posixpath
import model_harmonics as mdlhmc

# GLDAS models
gldas_products = {}
gldas_products['CLM'] = 'GLDAS Common Land Model (CLM)'
gldas_products['CLSM'] = 'GLDAS Catchment Land Surface Model (CLSM)'
gldas_products['MOS'] = 'GLDAS Mosaic model'
gldas_products['NOAH'] = 'GLDAS Noah model'
gldas_products['VIC'] = 'GLDAS Variable Infiltration Capacity (VIC) model'

# PURPOSE: sync local GLDAS files with GESDISC server
def gesdisc_gldas_sync(DIRECTORY, MODEL, YEARS, SPATIAL='', TEMPORAL='',
    VERSION='', EARLY=False, TIMEOUT=0, LOG=False, LIST=False, MODE=0o775,
    CLOBBER=False):

    # check if directory exists and recursively create if not
    DIRECTORY = pathlib.Path(DIRECTORY).expanduser().absolute()
    DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # create log file with list of synchronized files (or print to terminal)
    if LOG:
        # output to log file
        # format: NASA_GESDISC_GLDAS_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        output_logfile = f'NASA_GESDISC_GLDAS_{MODEL}_sync_{today}.log'
        LOGFILE = DIRECTORY.joinpath(output_logfile)
        logging.basicConfig(filename=LOGFILE, level=logging.INFO)
        logging.info(f'NASA GLDAS Sync Log ({today})')
    else:
        # standard output (terminal output)
        logging.basicConfig(level=logging.INFO)

    # Version and product flags
    V1,V2 = (f'_V{VERSION}','') if (VERSION == '1') else ('',f'.{VERSION}')
    PRODUCT = f'GLDAS_{MODEL}{SPATIAL}_{TEMPORAL}{V2}'
    EP = '_EP' if EARLY else ''
    # shortname for model on GESDISC server
    SHORTNAME = f'GLDAS_{MODEL}{SPATIAL}_{TEMPORAL}{EP}'

    # print header text to log/standard output
    logging.info(f'GLDAS MODEL={gldas_products[MODEL]}')
    logging.info(f'RESOLUTION={TEMPORAL},{SPATIAL}')
    logging.info(f'VERSION={VERSION}{EP}')

    # for each year to sync
    for Y in map(str,YEARS):
        # start and end date for query
        start_date,end_date = (f'{Y}-01-01',  f'{Y}-12-31')
        # query CMR for model granules
        ids,urls,mtimes = mdlhmc.utilities.cmr(SHORTNAME,
            version=VERSION, start_date=start_date, end_date=end_date,
            provider='GES_DISC', verbose=True)
        # recursively create local directory for data
        local_dir = DIRECTORY.joinpath(PRODUCT, Y)
        local_dir.mkdir(mode=MODE, parents=True, exist_ok=True)
        # sync model granules
        for id,url,mtime in zip(ids,urls,mtimes):
            # local version of the granule
            local_file = local_dir.joinpath(id)
            # copy file from remote directory comparing modified dates
            http_pull_file(url, mtime, local_file,
                TIMEOUT=TIMEOUT, LIST=LIST, CLOBBER=CLOBBER,
                MODE=MODE)

    # close log file and set permissions level to MODE
    if LOG:
        LOGFILE.chmod(mode=MODE)

# PURPOSE: pull file from a remote host checking if file exists locally
# and if the remote file is newer than the local file
def http_pull_file(remote_file,remote_mtime,local_file,
    TIMEOUT=0,LIST=False,CLOBBER=False,MODE=0o775):
    # if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    # check if local version of file exists
    if local_file.exists():
        # check last modification time of local file
        local_mtime = local_file.stat().st_mtime
        # if remote file is newer: overwrite the local file
        if (mdlhmc.utilities.even(remote_mtime) >
            mdlhmc.utilities.even(local_mtime)):
            TEST = True
            OVERWRITE = ' (overwrite)'
    else:
        TEST = True
        OVERWRITE = ' (new)'
    # if file does not exist locally, is to be overwritten, or CLOBBER is set
    if TEST or CLOBBER:
        # Printing files transferred
        logging.info(f'{remote_file} --> ')
        logging.info(f'\t{str(local_file)}{OVERWRITE}\n')
        # if executing copy command (not only printing the files)
        if not LIST:
            # Create and submit request. There are a wide range of exceptions
            # that can be thrown here, including HTTPError and URLError.
            request = mdlhmc.utilities.urllib2.Request(remote_file)
            response = mdlhmc.utilities.urllib2.urlopen(request,
                timeout=TIMEOUT)
            # chunked transfer encoding size
            CHUNK = 16 * 1024
            # copy contents to local file using chunked transfer encoding
            # transfer should work properly with ascii and binary data formats
            with local_file.open(mode='wb') as f:
                shutil.copyfileobj(response, f, CHUNK)
            # keep remote modification time of file and local access time
            os.utime(local_file, (local_file.stat().st_atime, remote_mtime))
            local_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Syncs GLDAS monthly datafiles from the Goddard Earth
            Sciences Data and Information Server Center (GES DISC)
            """
    )
    # command line parameters
    parser.add_argument('model',
        type=str, nargs='+', choices=gldas_products.keys(),
        help='GLDAS land surface model')
    # NASA Earthdata credentials
    parser.add_argument('--user','-U',
        type=str, default=os.environ.get('EARTHDATA_USERNAME'),
        help='Username for NASA Earthdata Login')
    parser.add_argument('--password','-W',
        type=str, default=os.environ.get('EARTHDATA_PASSWORD'),
        help='Password for NASA Earthdata Login')
    parser.add_argument('--netrc','-N',
        type=pathlib.Path, default=pathlib.Path.home().joinpath('.netrc'),
        help='Path to .netrc file for authentication')
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # years to download
    now = time.gmtime()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.tm_year+1),
        help='Years of model outputs to sync')
    # GLDAS model version
    parser.add_argument('--version','-v',
        type=str, default='2.1',
        help='GLDAS model version')
    # model spatial resolution
    # 10: 1.0 degrees latitude/longitude
    # 025: 0.25 degrees latitude/longitude
    parser.add_argument('--spacing','-S',
        type=str, default='10', choices=['10','025'],
        help='Spatial resolution of models to sync')
    # model temporal resolution
    # M: Monthly data products
    # 3H: 3-hourly data products
    parser.add_argument('--temporal','-T',
        type=str, default='M', choices=['M','3H'],
        help='Temporal resolution of models to sync')
    # GLDAS early products
    parser.add_argument('--early','-e',
        default=False, action='store_true',
        help='Sync GLDAS early products')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    # Output log file in form
    # NASA_GESDISC_GLDAS_sync_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    # sync options
    parser.add_argument('--list','-L',
        default=False, action='store_true',
        help='Only print files that could be transferred')
    parser.add_argument('--clobber','-C',
        default=False, action='store_true',
        help='Overwrite existing data in transfer')
    # permissions mode of the directories and files synced (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files synced')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # NASA Earthdata hostname
    URS = 'urs.earthdata.nasa.gov'
    # get NASA Earthdata credentials
    try:
        args.user,_,args.password = netrc.netrc(args.netrc).authenticators(URS)
    except:
        # check that NASA Earthdata credentials were entered
        if not args.user:
            prompt = f'Username for {URS}: '
            args.user = builtins.input(prompt)
        # enter password securely from command-line
        if not args.password:
            prompt = f'Password for {args.user}@{URS}: '
            args.password = getpass.getpass(prompt)

    # build a urllib opener for NASA GESDISC
    # Add the username and password for NASA Earthdata Login system
    mdlhmc.utilities.build_opener(args.user, args.password,
        password_manager=True, authorization_header=False)

    # check internet connection before attempting to run program
    HOST = posixpath.join('https://hydro1.gesdisc.eosdis.nasa.gov','data')
    if mdlhmc.utilities.check_credentials(HOST):
        # for each GLDAS model
        for MODEL in args.model:
            gesdisc_gldas_sync(args.directory, MODEL, args.year,
                VERSION=args.version, EARLY=args.early,
                SPATIAL=args.spacing, TEMPORAL=args.temporal,
                TIMEOUT=args.timeout, LOG=args.log, LIST=args.list,
                CLOBBER=args.clobber, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
