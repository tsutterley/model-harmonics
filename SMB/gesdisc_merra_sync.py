#!/usr/bin/env python
u"""
gesdisc_merra_sync.py
Written by Tyler Sutterley (12/2022)

Syncs MERRA-2 surface mass balance (SMB) related products from the Goddard
    Earth Sciences Data and Information Server Center (GES DISC)
    https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/
    https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python

Register with NASA Earthdata Login system:
    https://urs.earthdata.nasa.gov

Add "NASA GESDISC DATA ARCHIVE" to Earthdata Applications:
    https://urs.earthdata.nasa.gov/approve_app?client_id=e2WVk8Pw6weeLUKZYOxvTQ

tavgM_2d_int (Vertically Integrated Diagnostics) collection:
    PRECCU (convective rain)
    PRECLS (large-scale rain)
    PRECSN (snow)
    and EVAP (evaporation)
tavgM_2d_glc (Land Ice Surface Diagnostics) collection:
    RUNOFF (runoff over glaciated land)

CALLING SEQUENCE:
    python gesdisc_merra_sync.py --user <username>
    where <username> is your NASA Earthdata username

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user X: username for NASA Earthdata Login
    -W X, --password X: password for NASA Earthdata Login
    -N X, --netrc X: path to .netrc file for authentication
    -D X, --directory X: working data directory
    -v X, --version X: MERRA-2 version
    -Y X, --year X: years to sync
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
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 06/2022: use CMR queries to find reanalysis granules
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: lower case keyword arguments to output spatial
    Updated 10/2021: using python logging for handling verbose output
    Updated 06/2021: new last modified date format on GESDISC servers
    Updated 05/2021: added option for connection timeout (in seconds)
        use try/except for retrieving netrc credentials
    Updated 04/2021: set a default netrc file and check access
        default credentials from environmental variables
    Updated 02/2021: add back MERRA-2 invariant parameters sync
    Updated 01/2021: use argparse to set command line parameters
        using utilities program to build opener and list remote files
    Updated 09/2019: added ssl context to urlopen headers
    Updated 06/2018: using python3 compatible octal, input and urllib
    Updated 03/2018: --directory sets base directory similar to other programs
    Updated 08/2017: use raw_input() to enter NASA Earthdata credentials rather
        than exiting with error
    Updated 05/2017: exception if NASA Earthdata credentials weren't entered
        using os.makedirs to recursively create directories
        using getpass to enter server password securely (remove --password)
    Updated 04/2017: using lxml to parse HTML for files and modification dates
        minor changes to check_connection function to parallel other programs
    Written 11/2016
"""
from __future__ import print_function

import sys
import os
import time
import netrc
import shutil
import getpass
import logging
import argparse
import builtins
import posixpath
import model_harmonics as mdlhmc

# PURPOSE: sync local MERRA-2 files with GESDISC server
def gesdisc_merra_sync(DIRECTORY, YEARS=None, VERSION=None, TIMEOUT=None,
    LOG=False, LIST=False, MODE=None, CLOBBER=False):

    # check if directory exists and recursively create if not
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None
    # create log file with list of synchronized files (or print to terminal)
    if LOG:
        # output to log file
        # format: NASA_GESDISC_MERRA2_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = f'NASA_GESDISC_MERRA2_sync_{today}.log'
        logging.basicConfig(filename=os.path.join(DIRECTORY,LOGFILE),
            level=logging.INFO)
        logging.info(f'NASA MERRA-2 Sync Log ({today})')
    else:
        # standard output (terminal output)
        logging.basicConfig(level=logging.INFO)

    # query CMR for model MERRA-2 invariant products
    ids,urls,mtimes = mdlhmc.utilities.cmr('M2C0NXASM',
        version=VERSION, provider='GES_DISC', verbose=True)
    # sync model granules
    for id,url,mtime in zip(ids,urls,mtimes):
        # copy file from remote directory comparing modified dates
        http_pull_file(url, mtime, os.path.join(DIRECTORY,id),
            TIMEOUT=TIMEOUT, LIST=LIST, CLOBBER=CLOBBER,
            MODE=MODE)

    # for each MERRA-2 product to sync
    for SHORTNAME in ['M2TMNXINT','M2TMNXGLC']:
        PRODUCT = f'{SHORTNAME}.{VERSION}'
        logging.info(f'PRODUCT={PRODUCT}')
        # for each year to sync
        for Y in map(str,YEARS):
            # start and end date for query
            start_date = f'{Y}-01-01'
            end_date = f'{Y}-12-31'
            ids,urls,mtimes = mdlhmc.utilities.cmr(SHORTNAME,
                version=VERSION, start_date=start_date, end_date=end_date,
                provider='GES_DISC', verbose=True)
            # recursively create local directory for data
            if (not os.access(os.path.join(DIRECTORY,PRODUCT,Y), os.F_OK)):
                os.makedirs(os.path.join(DIRECTORY,PRODUCT,Y), MODE)
            # sync model granules
            for id,url,mtime in zip(ids,urls,mtimes):
                # copy file from remote directory comparing modified dates
                local_file = os.path.join(DIRECTORY,PRODUCT,Y,id)
                http_pull_file(url, mtime, local_file,
                    TIMEOUT=TIMEOUT, LIST=LIST, CLOBBER=CLOBBER,
                    MODE=MODE)

    # close log file and set permissions level to MODE
    if LOG:
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

# PURPOSE: pull file from a remote host checking if file exists locally
# and if the remote file is newer than the local file
def http_pull_file(remote_file, remote_mtime, local_file,
    TIMEOUT=None, LIST=False, CLOBBER=False, MODE=0o775):
    # if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    # check if local version of file exists
    if os.access(local_file, os.F_OK):
        # check last modification time of local file
        local_mtime = os.stat(local_file).st_mtime
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
        logging.info(f'\t{local_file}{OVERWRITE}\n')
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
            with open(local_file, 'wb') as f:
                shutil.copyfileobj(response, f, CHUNK)
            # keep remote modification time of file and local access time
            os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
            os.chmod(local_file, MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Syncs MERRA-2 surface mass balance (SMB) related
            products from the Goddard Earth Sciences Data and Information
            Server Center (GES DISC)
            """
    )
    # command line parameters
    # NASA Earthdata credentials
    parser.add_argument('--user','-U',
        type=str, default=os.environ.get('EARTHDATA_USERNAME'),
        help='Username for NASA Earthdata Login')
    parser.add_argument('--password','-W',
        type=str, default=os.environ.get('EARTHDATA_PASSWORD'),
        help='Password for NASA Earthdata Login')
    parser.add_argument('--netrc','-N',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.path.join(os.path.expanduser('~'),'.netrc'),
        help='Path to .netrc file for authentication')
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # MERRA-2 version
    parser.add_argument('--version','-v',
        type=str, default='5.12.4',
        help='MERRA-2 version')
    # years to download
    now = time.gmtime()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(1980,now.tm_year+1),
        help='Years of model outputs to sync')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    # Output log file in form
    # NASA_GESDISC_MERRA2_sync_2002-04-01.log
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
    HOST = posixpath.join('http://goldsmr4.gesdisc.eosdis.nasa.gov','data')
    if mdlhmc.utilities.check_connection(HOST):
        gesdisc_merra_sync(args.directory,
            YEARS=args.year,
            VERSION=args.version,
            TIMEOUT=args.timeout,
            LOG=args.log,
            LIST=args.list,
            CLOBBER=args.clobber,
            MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
