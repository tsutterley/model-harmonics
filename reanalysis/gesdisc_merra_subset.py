#!/usr/bin/env python
u"""
gesdisc_merra_subset.py
Written by Tyler Sutterley (12/2022)

Subsets monthly MERRA-2 products for specific variables from the
    Goddard Earth Sciences Data and Information Server Center (GES DISC)
    https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/
    https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python

Register with NASA Earthdata Login system:
    https://urs.earthdata.nasa.gov

Add "NASA GESDISC DATA ARCHIVE" to Earthdata Applications:
    https://urs.earthdata.nasa.gov/approve_app?client_id=e2WVk8Pw6weeLUKZYOxvTQ

CALLING SEQUENCE:
    python gesdisc_merra_subset.py --user <username> --variables PS
    where <username> is your NASA Earthdata username

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user X: username for NASA Earthdata Login
    -W X, --password X: password for NASA Earthdata Login
    -N X, --netrc X: path to .netrc file for authentication
    -D X, --directory X: Working data directory
    -s X, --shortname X: MERRA-2 product shortname
    -v X, --version X: MERRA-2 version
    -Y X, --year X: years to sync
    -V X, --variables X: variables to subset
    -t X, --timeout X: Timeout in seconds for blocking operations
    -F, --flatten: Do not create subdirectories
    -l, --log: output log of files downloaded
    -C, --clobber: Overwrite existing files
    -M X, --mode X: Local permissions mode of the files created

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
    time.py: Utilities for calculating time operations
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 10/2022: added flatten option to not create output subdirectories
    Written 06/2022
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
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc

# PURPOSE: subsets MERRA-2 files for specific variables
def gesdisc_merra_subset(base_dir, SHORTNAME, VERSION=None, YEARS=None,
    VARIABLES=None, TIMEOUT=None, FLATTEN=False, LOG=False, CLOBBER=False,
    MODE=None):
    # set up data directory
    if FLATTEN:
        DIRECTORY = os.path.expanduser(base_dir)
    else:
        DIRECTORY = os.path.join(base_dir, f'{SHORTNAME}.{VERSION}')
    # check if DIRECTORY exists and recursively create if not
    if not os.access(os.path.join(DIRECTORY), os.F_OK):
        os.makedirs(os.path.join(DIRECTORY), mode=MODE, exist_ok=True)

    # create log file with list of synchronized files (or print to terminal)
    if LOG:
        # format: NASA_GESDISC_MERRA2_subset_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = f'NASA_GESDISC_MERRA2_subset_{today}.log'
        fid = open(os.path.join(DIRECTORY,LOGFILE), mode='w', encoding='utf8')
        logging.basicConfig(stream=fid, level=logging.INFO)
        logging.info(f'NASA MERRA-2 Sync Log ({today})')
        logging.info(f'PRODUCT: {SHORTNAME}.{VERSION}')
    else:
        # standard output (terminal output)
        fid = sys.stdout
        logging.basicConfig(stream=fid, level=logging.INFO)

    # for each unique date
    for YEAR in YEARS:
        dpm = gravtk.time.calendar_days(YEAR)
        # for each month of the year
        for i,days_per_month in enumerate(dpm):
            # year and month as strings
            YY = f'{YEAR:4d}'
            MM = f'{i+1:02d}'
            # start and end date for query
            start_date = f'{YY}-{MM}-{1:02.0f}'
            end_date = f'{YY}-{MM}-{days_per_month:02.0f}'
            # query for data
            ids,urls,mtimes = mdlhmc.utilities.cmr(SHORTNAME,
                version=VERSION, start_date=start_date, end_date=end_date,
                provider='GES_DISC', verbose=True)
            # skip years and months without any data
            if not ids:
                continue
            # for each granule
            for id,url,mtime in zip(ids,urls,mtimes):
                # build filename for output
                fileBasename,_ = os.path.splitext(id)
                FILE = f'{fileBasename}.SUB.nc'
                local_file = os.path.join(DIRECTORY, FILE)
                # get subsetting API url for granule
                request_url = mdlhmc.utilities.build_request(
                    SHORTNAME, VERSION, url, variables=VARIABLES,
                    bbox=[-90,-180,90,180], LABEL=FILE)
                # copy subsetted file and update modified dates
                http_pull_file(request_url, mtime, local_file,
                    TIMEOUT=TIMEOUT, CLOBBER=CLOBBER, MODE=0o775)

    # close log file and set permissions level to MODE
    if LOG:
        fid.close()
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

# PURPOSE: pull file from a remote host checking if file exists locally
# and if the remote file is newer than the local file
def http_pull_file(remote_file, remote_mtime, local_file,
    TIMEOUT=None, CLOBBER=False, MODE=0o775):
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
        logging.info(f'{remote_file} -->')
        logging.info(f'\t{local_file}{OVERWRITE}\n')
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
        description="""Creates monthly MERRA-2 3D model level
            products syncing data from the Goddard Earth Sciences
            Data and Information Server Center (GES DISC)
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
    # MERRA-2 product shortname
    parser.add_argument('--shortname','-s',
        type=str, default='M2TMNXSLV',
        help='MERRA-2 product shortname')
    # MERRA-2 version
    parser.add_argument('--version','-v',
        type=str, default='5.12.4',
        help='MERRA-2 version')
    # years to download
    now = time.gmtime()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.tm_year+1),
        help='Years of model outputs to sync')
    # variable names to subset from product
    parser.add_argument('--variables','-V',
        type=str, nargs='+',
        help='Variables to subset from product')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    # output subdirectories
    parser.add_argument('--flatten','-F',
        default=False, action='store_true',
        help='Do not create subdirectories')
    # Output log file in form
    # NASA_GESDISC_MERRA2_subset_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    # sync options
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
    except Exception as e:
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
    HOST = 'https://goldsmr5.gesdisc.eosdis.nasa.gov/'
    if mdlhmc.utilities.check_credentials(HOST):
        gesdisc_merra_subset(args.directory, args.shortname,
            VERSION=args.version,
            YEARS=args.year,
            VARIABLES=args.variables,
            TIMEOUT=args.timeout,
            FLATTEN=args.flatten,
            LOG=args.log,
            CLOBBER=args.clobber,
            MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
