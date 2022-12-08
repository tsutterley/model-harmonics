#!/usr/bin/env python
u"""
noaa_cdc_ncep_ftp.py
Written by Tyler Sutterley (12/2022)

Syncs NOAA-DOE-2 surface reanalysis outputs with the NOAA CDC ftp server
    ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis2.dailyavgs/surface/

CALLING SEQUENCE:
    noaa_cdc_ncep_ftp(base_dir, YEAR=YEAR, MASK=True, INVARIANT=True)

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory=X: working data directory
    -Y X, --year X: years to download
    --mask: Download land-sea mask file (land.nc)
    --invariant: Download invariant parameters file (hgt.sfc.nc)
    -t X, --timeout X: Timeout in seconds for blocking operations
    -l, --log: output log of files downloaded
    -M X, --mode X: Permissions mode of the directories and files downloaded

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
    Updated 03/2021: automatically update years to run based on current time
    Updated 12/2020: using utilities module to list and retrieve from ftp
    Written 09/2019
"""
from __future__ import print_function

import sys
import os
import re
import time
import logging
import argparse
import posixpath
import gravity_toolkit as gravtk

# PURPOSE: sync local NCEP-DOE-2 reanalysis files with NOAA CDC server
def noaa_cdc_ncep_ftp(base_dir, YEAR=None, MASK=False, INVARIANT=False,
    TIMEOUT=None, LOG=False, MODE=None):

    # directory setup
    DIRECTORY = os.path.join(base_dir,'NCEP-DOE-2')
    # check if log directory exists and recursively create if not
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None
    # create log file with list of synchronized files (or print to terminal)
    if LOG:
        # output to log file
        # format: NOAA_CDC_NCEP-DOE-2_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = f'NOAA_CDC_NCEP-DOE-2_sync_{today}.log'
        fid1 = open(os.path.join(DIRECTORY,LOGFILE), mode='w', encoding='utf8')
        logging.basicConfig(stream=fid1, level=logging.INFO)
        logging.info(f'NOAA CDC Sync Log ({today})')
        logging.info('PRODUCT: NCEP-DOE-2')
    else:
        # standard output (terminal output)
        fid1 = sys.stdout
        logging.basicConfig(stream=fid1, level=logging.INFO)

    # remote directory for data
    HOST = ['ftp.cdc.noaa.gov','Datasets','ncep.reanalysis2.dailyavgs','surface']

    # compile the regular expression operator to find files
    regex_years = r'|'.join([rf'{y:4d}' for y in YEAR])
    R1 = re.compile(rf'pres.sfc.({regex_years}).nc', re.VERBOSE)
    # list filenames and modification times from remote directory
    remote_files,remote_mtimes = gravtk.utilities.ftp_list(HOST,
        timeout=TIMEOUT, basename=True, pattern=R1, sort=True)
    for fi,mtime in zip(remote_files,remote_mtimes):
        # extract filename from regex object
        local_file = os.path.join(DIRECTORY,fi)
        MD5 = gravtk.utilities.get_hash(local_file)
        gravtk.utilities.from_ftp(posixpath.join(*HOST,fi),
            local=local_file, timeout=TIMEOUT, hash=MD5,
            verbose=True, fid=fid1, mode=MODE)

    # get mask file
    if MASK:
        # extract filename from regex object
        local_file = os.path.join(DIRECTORY,'land.nc')
        MD5 = gravtk.utilities.get_hash(local_file)
        gravtk.utilities.from_ftp(posixpath.join(*HOST,'land.nc'),
            local=local_file, timeout=TIMEOUT, hash=MD5,
            verbose=True, fid=fid1, mode=MODE)

    # get invariant file
    if INVARIANT:
        # extract filename from regex object
        local_file = os.path.join(DIRECTORY,'hgt.sfc.nc')
        MD5 = gravtk.utilities.get_hash(local_file)
        gravtk.utilities.from_ftp(posixpath.join(*HOST,'hgt.sfc.nc'),
            local=local_file, timeout=TIMEOUT, hash=MD5,
            verbose=True, fid=fid1, mode=MODE)

    # close log file and set permissions level to MODE
    if LOG:
        fid1.close()
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Syncs NOAA-DOE-2 surface reanalysis
            outputs from NOAA CDC ftp server
            """
    )
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # years to retrieve
    now = time.gmtime()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.tm_year+1),
        help='Model years to retrieve')
    # retrieve the model land surface mask (land.nc)
    parser.add_argument('--mask','-m',
        default=False, action='store_true',
        help='Retrieve model land surface mask')
    # retrieve the model invariant parameters (hgt.sfc.nc)
    parser.add_argument('--invariant','-I',
        default=False, action='store_true',
        help='Retrieve model invariant parameters')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    # Output log file in form
    # NOAA_CDC_NCEP-DOE-2_sync_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    # permissions mode of the directories and files retrieved
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files retrieved')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # run program for model
    if gravtk.utilities.check_ftp_connection('ftp.cdc.noaa.gov'):
        noaa_cdc_ncep_ftp(args.directory, YEAR=args.year, MASK=args.mask,
            INVARIANT=args.invariant, TIMEOUT=args.timeout, LOG=args.log,
            MODE=args.mode)
    else:
        raise RuntimeError('Check internet connection')

# run main program
if __name__ == '__main__':
    main()
