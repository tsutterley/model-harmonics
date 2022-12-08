#!/usr/bin/env python
u"""
jpl_ecco_llc_sync.py
Written by Tyler Sutterley (12/2022)

Syncs ECCO LLC tile model outputs from the NASA JPL ECCO Drive server:
https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/nctiles_monthly
https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/nctiles_monthly

https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
https://nsidc.org/support/faq/what-options-are-available-bulk-downloading-data-
    https-earthdata-login-enabled
http://www.voidspace.org.uk/python/articles/authentication.shtml#base64

Register with NASA Earthdata Login system:
https://urs.earthdata.nasa.gov

DATA DESCRIPTION:
    ECCO employs the MITgcm in a global domain incorporating most of
        the available satellite and in-situ data to produce a physically
        consistent ocean estimate
    Model fields are available at monthly-day intervals

CALLING SEQUENCE:
    python jpl_ecco_llc_sync.py --year 2015 2016 --user <username> V4r4
    where <username> is your NASA Earthdata username

INPUTS:
    ECCO version 4 or 5 models
        V4r4: Version 4, Revision 4
        V5alpha: Version 5, Alpha release

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
    Updated 07/2021: add warning for Version 4, Revision 4
    Updated 05/2021: added option for connection timeout (in seconds)
        use try/except for retrieving netrc credentials
    Updated 04/2021: set a default netrc file and check access
        default credentials from environmental variables
    Written 02/2021
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
import warnings
import posixpath
import lxml.etree
import gravity_toolkit as gravtk

# PURPOSE: sync ECCO LLC tile data from JPL ECCO drive server
def jpl_ecco_llc_sync(ddir, MODEL, YEAR=None, PRODUCT=None, TIMEOUT=None,
    LOG=False, LIST=False, CLOBBER=False, CHECKSUM=False, MODE=None):

    # check if directory exists and recursively create if not
    DIRECTORY = os.path.join(ddir,f'ECCO-{MODEL}','nctiles_monthly')
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None

    # remote https server for ECCO data
    HOST = 'https://ecco.jpl.nasa.gov'
    # compile HTML parser for lxml
    parser = lxml.etree.HTMLParser()

    # create log file with list of synchronized files (or print to terminal)
    if LOG:
        # format: JPL_ECCO_V5alpha_PHIBOT_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = f'JPL_ECCO_{MODEL}_{PRODUCT}_{today}.log'
        logging.basicConfig(filename=os.path.join(DIRECTORY,LOGFILE),
            level=logging.INFO)
        logging.info(f'ECCO LLC {MODEL} {PRODUCT} Sync Log ({today})')
    else:
        # standard output (terminal output)
        logging.basicConfig(level=logging.INFO)

    # print the model synchronized
    logging.info(f'MODEL: {MODEL}\n')

    # print warning for Version 4, Revision 4
    # https://ecco-group.org/docs/ECCO_V4r4_errata.pdf
    if MODEL in ('V4r4',):
        warnings.filterwarnings("always")
        warnings.warn("See Errata for V4r4 Atmospheric Pressure Forcing")

    # download the ECCO llc grid file
    grid_path = {}
    grid_path['V4r4'] = ['Version4','Release4','nctiles_grid']
    grid_path['V5alpha'] = ['Version5','Alpha','nctiles_grid']
    # remote subdirectory for MODEL on JPL ECCO data server
    PATH = [HOST,'drive','files',*grid_path[MODEL]]
    colnames,mtimes = gravtk.utilities.drive_list(PATH,
        timeout=TIMEOUT, build=False, parser=parser,
        pattern=r'ECCO-GRID.nc$')
    # full path to remote directory
    remote_dir = posixpath.join(*PATH)
    # remote and local versions of the file
    for colname,remote_mtime in zip(colnames,mtimes):
        remote_file = posixpath.join(remote_dir,colname)
        local_file = os.path.join(DIRECTORY,colname)
        http_pull_file(remote_file, remote_mtime, local_file,
            TIMEOUT=TIMEOUT, LIST=LIST, CLOBBER=CLOBBER,
            CHECKSUM=CHECKSUM, MODE=MODE)

    # path to model llc tile files
    model_path = {}
    model_path['V4r4'] = ['Version4','Release4','nctiles_monthly',PRODUCT]
    model_path['V5alpha'] = ['Version5','Alpha','nctiles_monthly',PRODUCT]
    # compile regular expression operator for years to sync
    if YEAR is None:
        regex_years = r'\d+'
    else:
        regex_years = r'|'.join([rf'{y:d}' for y in YEAR])
    # compile regular expression operator finding years
    if MODEL in ('V5alpha'):
        R1 = re.compile(rf'{PRODUCT}([\.\_])({regex_years})(_\d+)?.nc$')
    elif MODEL in ('V4r4',):
        R1 = re.compile(regex_years, re.VERBOSE)

    # open connection with ECCO drive server at remote directory
    # find remote yearly directories for MODEL
    years,mtimes = gravtk.utilities.drive_list(PATH,
        timeout=TIMEOUT, build=False, parser=parser,
        pattern=R1, sort=True)
    for yr in years:
        # extract year from file
        if MODEL in ('V5alpha'):
            _,YY,_ = R1.findall(yr).pop()
        elif MODEL in ('V4r4',):
            # print string for year
            logging.info(yr)
            # add the year directory to the path
            YY, = R1.findall(yr)
            PATH.append(yr)
        # compile regular expression operator for model product files
        R3 = re.compile(rf'{PRODUCT}([\.\_])({YY})(_\d+)?.nc$', re.VERBOSE)
        # full path to remote directory
        remote_dir = posixpath.join(*PATH)
        # read and parse request for files (find names and modified dates)
        colnames,mtimes = gravtk.utilities.drive_list(PATH,
            timeout=TIMEOUT, build=False, parser=parser,
            pattern=R3, sort=True)
        # for each file on the remote server
        for colname,remote_mtime in zip(colnames,mtimes):
            # remote and local versions of the file
            remote_file = posixpath.join(remote_dir,colname)
            local_file = os.path.join(DIRECTORY,colname)
            http_pull_file(remote_file, remote_mtime, local_file,
                TIMEOUT=TIMEOUT, LIST=LIST, CLOBBER=CLOBBER,
                CHECKSUM=CHECKSUM, MODE=MODE)
        # remove the year directory from the path
        if MODEL in ('V4r4',):
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
        description="""Syncs ECCO Version 4 and 5 LLC tile outputs
        from the NASA JPL ECCO Drive server
        """
    )
    # command line parameters
    parser.add_argument('model',
        type=str, nargs='+',
        default=['V5alpha'], choices=['V4r4','V5alpha'],
        help='ECCO Version 4 or 5 Model')
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
    parser.add_argument('--product','-P',
        type=str, default='PHIBOT',
        help='Product to sync')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    # Output log file in form
    # JPL_ECCO_V5alpha_PHIBOT_sync_2002-04-01.log
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
            jpl_ecco_llc_sync(args.directory, MODEL, YEAR=args.year,
                PRODUCT=args.product, TIMEOUT=args.timeout, LOG=args.log,
                LIST=args.list, CLOBBER=args.clobber, CHECKSUM=args.checksum,
                MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
