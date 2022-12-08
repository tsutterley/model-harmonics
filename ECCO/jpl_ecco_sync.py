#!/usr/bin/env python
u"""
jpl_ecco_sync.py
Written by Tyler Sutterley (12/2022)

Syncs ECCO Near Real-Time model outputs from the NASA JPL ECCO Drive server:
    https://ecco.jpl.nasa.gov/drive/files/NearRealTime/Readme
    https://ecco.jpl.nasa.gov/drive/files/NearRealTime/KalmanFilter/
    https://ecco.jpl.nasa.gov/drive/files/NearRealTime/Smoother/

https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
https://nsidc.org/support/faq/what-options-are-available-bulk-downloading-data-
    https-earthdata-login-enabled
http://www.voidspace.org.uk/python/articles/authentication.shtml#base64

Register with NASA Earthdata Login system:
https://urs.earthdata.nasa.gov

DATA DESCRIPTION:
    ECCO-JPL employs the MITgcm in a near-global domain (78S--78N).  Model
        resolution is 1 degree horizontally except within the tropics where
        meridional resolution gradually decreases to 0.3 degree within 10 degree
        of the Equator.  There are 46 levels in the vertical with 10m-resolution
        within 150m of the surface.
    Total grid dimension is 360*224*46 = 4*10^6.
    The GM isentropic mixing scheme (Gent and McWilliams, 1990) and the KPP
        mixed-layer formulation (Large et al., 1994) are employed. The model is
        forced by NCEP reanalysis products (12-hourly wind stress, daily
        diabatic air-sea fluxes) with time-means replaced by those of COADS.
        Temperature and salinity at the model sea surface are further relaxed
        towards observed values.
    Model fields are available at 10-day intervals (10-day averages).
    Sea level and ocean bottom pressure are also available at 12-hour intervals.

CALLING SEQUENCE:
    python jpl_ecco_sync.py --year 2015 2016 --user <username> kf080i
    where <username> is your NASA Earthdata username

INPUTS:
    ECCO Near Real-Time models
        kf080i: Kalman filter analysis
        dr080i: RTS smoother analysis

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user X: username for NASA Earthdata Login
    -W X, --webdav X: WebDAV password for JPL ECCO Drive Login
    -N X, --netrc X: path to .netrc file for authentication
    -D X, --directory X: working data directory
    -Y X, --year X: Years to sync
    -P X, --product X: Product to sync
    -t X, --timeout X: Timeout in seconds for blocking operations
    -l, --log: output log of files downloaded
    -L, --list: print files to be transferred, but do not execute transfer
    -C, --clobber: Overwrite existing data in transfer
    --checksum: compare hashes to check if overwriting existing data
    -M X, --mode X: Local permissions mode of the directories and files synced

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
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 10/2021: using python logging for handling verbose output
    Updated 05/2021: added option for connection timeout (in seconds)
        use try/except for retrieving netrc credentials
    Updated 04/2021: set a default netrc file and check access
        default credentials from environmental variables
    Updated 01/2021: added option to generalize for different products
    Updated 12/2020 for public release.
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: flake8 compatible regular expression strings
        moved urllib opener to utilities. add credential check
        moved urllib directory listing to utilities
    Updated 06/2020: increased timeout to 2 minutes
    Updated 05/2020: simplified JPL ECCO Drive login
        added netrc option for alternative authentication method
    Updated 12/2019: convert last modified time in a function for format
    Updated 09/2019: replacing ftp with https JPL ECCO Drive
        added ssl context to urlopen headers
        added checksum option to not overwrite existing data files
    Updated 06/2019: recommending kf080i for the Kalman filtered solution
    Updated 03/2018: updated with new JPL ECCO model ftp servers
    Updated 01/2018: new directory structure under NearRealTime on ftp server
    Updated 05/2017: using os.makedirs to recursively create directories
    Updated 04/2017: minor changes to check_connection function to use ftplib
    Updated 01/2017: using ftplib for syncing, input model as argument
        added --directory, --year, --mode and --clobber command-line options
    Updated 09/2016: print the model synchronized even if not using LOG
    Updated 07/2016: can run for multiple ECCO models (separate by commas)
    Updated 06/2016: added --list option for a dry-run (do not transfer files)
        added --model option for different ECCO models (kf080g vs. dr080i)
    Updated 05-06/2016: using __future__ print function. format log line
    Written 03/2016
"""
from __future__ import print_function

import sys
import os
import re
import io
import time
import netrc
import shutil
import getpass
import logging
import argparse
import builtins
import posixpath
import lxml.etree
import gravity_toolkit as gravtk

# PURPOSE: sync ECCO Near Real-Time model data from JPL ECCO drive server
def jpl_ecco_sync(DIRECTORY, MODEL, YEAR=None, PRODUCT=None, TIMEOUT=None,
    LOG=False, LIST=False, CLOBBER=False, CHECKSUM=False, MODE=None):

    # check if directory exists and recursively create if not
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None

    # remote https server for ECCO data
    HOST = 'https://ecco.jpl.nasa.gov'
    # compile HTML parser for lxml
    parser = lxml.etree.HTMLParser()

    # create log file with list of synchronized files (or print to terminal)
    if LOG:
        # format: JPL_ECCO_kf080i_OBP_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        args = (MODEL,PRODUCT,today)
        LOGFILE = f'JPL_ECCO_{MODEL}_{PRODUCT}_{today}.log'
        logging.basicConfig(filename=os.path.join(DIRECTORY,LOGFILE),
            level=logging.INFO)
        logging.info(f'ECCO Near Real-Time {PRODUCT} Sync Log ({today})')
    else:
        # standard output (terminal output)
        logging.basicConfig(level=logging.INFO)

    # print the model synchronized
    logging.info(f'MODEL: {MODEL}\n')

    # path to model files
    model_path = {}
    model_path['kf080i'] = ['NearRealTime','KalmanFilter']
    model_path['dr080i'] = ['NearRealTime','Smoother']
    # compile regular expression operator for years to sync
    # will not include the PREV, misc, forcing directories or aux files
    if YEAR is None:
        regex_years = r'\d+'
    else:
        regex_years = r'|'.join([rf'{y:d}' for y in YEAR])
    R1 = re.compile(rf'{MODEL}_({regex_years})', re.VERBOSE)
    # compile regular expression operator for subdirectories
    R2 = re.compile(r'n10day_(\d+)_(\d+)', re.VERBOSE)
    # compile regular expression operator for model product
    R3 = re.compile(rf'{PRODUCT}_(.*?).cdf$', re.VERBOSE)

    # remote subdirectory for MODEL on JPL ECCO data server
    PATH = [HOST,'drive','files',*model_path[MODEL]]
    # open connection with ECCO drive server at remote directory
    # find remote yearly directories for MODEL
    years,mtimes = gravtk.utilities.drive_list(PATH,
        timeout=TIMEOUT,build=False,parser=parser,pattern=R1,sort=True)
    for yr in years:
        # print string for year
        logging.info(yr)
        # add the year directory to the path
        PATH.append(yr)
        # open connection with ECCO drive server at remote directory
        # read and parse request for remote subdirectories
        subdirs,mtimes = gravtk.utilities.drive_list(PATH,
            timeout=TIMEOUT, build=False, parser=parser,
            pattern=R2, sort=True)
        # for each remote subdirectory
        for sd in subdirs:
            # add the subdirecotry directory to the path
            PATH.append(sd)
            # full path to remote directory
            remote_dir = posixpath.join(*PATH)
            # local directory for exact data product
            local_dir = os.path.join(DIRECTORY, yr, sd)
            # check if directory exists and recursively create if not
            if not os.path.exists(local_dir):
                os.makedirs(local_dir,MODE)
            # read and parse request for files (find names and modified dates)
            colnames,mtimes = gravtk.utilities.drive_list(PATH,
                timeout=TIMEOUT, build=False, parser=parser,
                pattern=R3, sort=True)
            # for each file on the remote server
            for colname,remote_mtime in zip(colnames,mtimes):
                # remote and local versions of the file
                remote_file = posixpath.join(remote_dir,colname)
                local_file = os.path.join(local_dir,colname)
                http_pull_file(remote_file, remote_mtime,
                    local_file, TIMEOUT=TIMEOUT, LIST=LIST,
                    CLOBBER=CLOBBER, CHECKSUM=CHECKSUM, MODE=MODE)
            # remove the directory from the path
            PATH.remove(sd)
        # remove the year directory from the path
        PATH.remove(yr)

    # close log file and set permissions level to MODE
    if LOG:
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

# PURPOSE: pull file from a remote host checking if file exists locally
# and if the remote file is newer than the local file
def http_pull_file(remote_file, remote_mtime, local_file,
    TIMEOUT=None, LIST=False, CLOBBER=False, CHECKSUM=False, MODE=0o775):
    # if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    # check if local version of file exists
    if CHECKSUM and os.access(local_file, os.F_OK):
        # generate checksum hash for local file
        # open the local_file in binary read mode
        local_hash = gravtk.utilities.get_hash(local_file)
        # Create and submit request.
        # There are a wide range of exceptions that can be thrown here
        # including HTTPError and URLError.
        request = gravtk.utilities.urllib2.Request(remote_file)
        response = gravtk.utilities.urllib2.urlopen(request,
            timeout=TIMEOUT)
        # copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO(response.read())
        remote_buffer.seek(0)
        # generate checksum hash for remote file
        remote_hash = gravtk.utilities.get_hash(remote_buffer)
        # compare checksums
        if (local_hash != remote_hash):
            TEST = True
            OVERWRITE = f' (checksums: {local_hash} {remote_hash})'
    elif os.access(local_file, os.F_OK):
        # check last modification time of local file
        local_mtime = os.stat(local_file).st_mtime
        # if remote file is newer: overwrite the local file
        if (remote_mtime > local_mtime):
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
            # chunked transfer encoding size
            CHUNK = 16 * 1024
            # copy bytes or transfer file
            if CHECKSUM and os.access(local_file, os.F_OK):
                # store bytes to file using chunked transfer encoding
                remote_buffer.seek(0)
                with open(local_file, 'wb') as f:
                    shutil.copyfileobj(remote_buffer, f, CHUNK)
            else:
                # Create and submit request.
                # There are a wide range of exceptions that can be thrown here
                # including HTTPError and URLError.
                request = gravtk.utilities.urllib2.Request(remote_file)
                response = gravtk.utilities.urllib2.urlopen(request,
                    timeout=TIMEOUT)
                # copy contents to local file using chunked transfer encoding
                # transfer should work properly with ascii and binary formats
                with open(local_file, 'wb') as f:
                    shutil.copyfileobj(response, f, CHUNK)
            # keep remote modification time of file and local access time
            os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
            os.chmod(local_file, MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Syncs ECCO Near Real-Time model outputs from the
        NASA JPL ECCO Drive server
        """
    )
    # command line parameters
    parser.add_argument('model',
        type=str, nargs='+',
        default=['kf080i','dr080i'], choices=['kf080i','dr080i'],
        help='ECCO Near Real-Time Model')
    # NASA Earthdata credentials
    parser.add_argument('--user','-U',
        type=str, default=os.environ.get('EARTHDATA_USERNAME'),
        help='Username for NASA Earthdata Login')
    parser.add_argument('--webdav','-W',
        type=str, default=os.environ.get('ECCO_PASSWORD'),
        help='WebDAV password for JPL ECCO Drive Login')
    parser.add_argument('--netrc','-N',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.path.join(os.path.expanduser('~'),'.netrc'),
        help='Path to .netrc file for authentication')
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # ECCO model years to sync
    parser.add_argument('--year','-Y',
        type=int, nargs='+',
        help='Years to sync')
    # ECCO model product to sync
    parser.add_argument('--product', '-P',
        type=str, default='OBP',
        help='Product to sync')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    # Output log file in form
    # JPL_ECCO_kf080i_OBP_sync_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    # sync options
    parser.add_argument('--list','-L',
        default=False, action='store_true',
        help='Only print files that could be transferred')
    parser.add_argument('--checksum',
        default=False, action='store_true',
        help='Compare hashes to check for overwriting existing data')
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

    # JPL ECCO drive hostname
    HOST = 'ecco.jpl.nasa.gov'
    # get NASA Earthdata and JPL ECCO drive credentials
    try:
        args.user,_,args.webdav = netrc.netrc(args.netrc).authenticators(HOST)
    except:
        # check that NASA Earthdata credentials were entered
        if not args.user:
            prompt = f'Username for {HOST}: '
            args.user = builtins.input(prompt)
        # enter WebDAV password securely from command-line
        if not args.webdav:
            prompt = f'Password for {args.user}@{HOST}: '
            args.webdav = getpass.getpass(prompt)

    # build a urllib opener for JPL ECCO Drive
    # Add the username and password for NASA Earthdata Login system
    gravtk.utilities.build_opener(args.user,args.webdav)

    # check internet connection before attempting to run program
    # check JPL ECCO Drive credentials before attempting to run program
    DRIVE = f'https://{HOST}/drive/files'
    if gravtk.utilities.check_credentials(DRIVE):
        for MODEL in args.model:
            jpl_ecco_sync(args.directory, MODEL, YEAR=args.year,
                PRODUCT=args.product, TIMEOUT=args.timeout, LOG=args.log,
                LIST=args.list, CLOBBER=args.clobber, CHECKSUM=args.checksum,
                MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
