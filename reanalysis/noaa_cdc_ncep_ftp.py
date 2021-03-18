#!/usr/bin/env python
u"""
noaa_cdc_ncep_ftp.py
Written by Tyler Sutterley (03/2021)
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
    Updated 03/2021: automatically update years to run based on current time
    Updated 12/2020: using utilities module to list and retrieve from ftp
    Written 09/2019
"""
from __future__ import print_function

import sys
import os
import re
import time
import argparse
import gravity_toolkit.utilities

#-- PURPOSE: sync local NCEP-DOE-2 reanalysis files with NOAA CDC server
def noaa_cdc_ncep_ftp(base_dir, YEAR=None, MASK=False, INVARIANT=False,
    LOG=False, MODE=None):

    #-- directory setup
    DIRECTORY = os.path.join(base_dir,'NCEP-DOE-2')
    #-- check if log directory exists and recursively create if not
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None
    #-- create log file with list of synchronized files (or print to terminal)
    if LOG:
        #-- output to log file
        #-- format: NOAA_CDC_NCEP-DOE-2_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'NOAA_CDC_NCEP-DOE-2_sync_{0}.log'.format(today)
        fid1 = open(os.path.join(DIRECTORY,LOGFILE),'w')
        print('NOAA CDC Sync Log ({0})'.format(today), file=fid1)
        print('PRODUCT: NCEP-DOE-2',file=fid1)
    else:
        #-- standard output (terminal output)
        fid1 = sys.stdout

    #-- remote directory for data
    HOST = ['ftp.cdc.noaa.gov','Datasets','ncep.reanalysis2.dailyavgs','surface']

    #-- compile the regular expression operator to find files
    regex_years = '|'.join(['{0:4d}'.format(y) for y in YEAR])
    R1 = re.compile('pres.sfc.({0}).nc'.format(regex_years))
    #-- list filenames and modification times from remote directory
    remote_files,remote_mtimes = gravity_toolkit.utilities.ftp_list(HOST,
        basename=True, pattern=R1, sort=True)
    for fi,mtime in zip(remote_files,remote_mtimes):
        #-- extract filename from regex object
        local_file = os.path.join(DIRECTORY,fi)
        MD5 = gravity_toolkit.utilities.get_hash(local_file)
        gravity_toolkit.utilities.from_ftp(HOST + [fi],
            local=local_file, hash=MD5, verbose=True, fid=fid1, mode=MODE)

    #-- get mask file
    if MASK:
        #-- extract filename from regex object
        local_file = os.path.join(DIRECTORY,'land.nc')
        MD5 = gravity_toolkit.utilities.get_hash(local_file)
        gravity_toolkit.utilities.from_ftp(HOST + ['land.nc'],
            local=local_file, hash=MD5, verbose=True, fid=fid1, mode=MODE)

    #-- get invariant file
    if INVARIANT:
        #-- extract filename from regex object
        local_file = os.path.join(DIRECTORY,'hgt.sfc.nc')
        MD5 = gravity_toolkit.utilities.get_hash(local_file)
        gravity_toolkit.utilities.from_ftp(HOST + ['hgt.sfc.nc'],
            local=local_file, hash=MD5, verbose=True, fid=fid1, mode=MODE)

    #-- close log file and set permissions level to MODE
    if LOG:
        fid1.close()
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

#-- Main program that calls ecmwf_reanalysis_retrieve()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Syncs NOAA-DOE-2 surface reanalysis
            outputs from NOAA CDC ftp server
            """
    )
    #-- command line parameters
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- years to retrieve
    now = time.gmtime()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.tm_year+1),
        help='Model years to retrieve')
    #-- retrieve the model land surface mask (land.nc)
    parser.add_argument('--mask','-m',
        default=False, action='store_true',
        help='Retrieve model land surface mask')
    #-- retrieve the model invariant parameters (hgt.sfc.nc)
    parser.add_argument('--invariant','-I',
        default=False, action='store_true',
        help='Retrieve model invariant parameters')
    #-- Output log file in form
    #-- NOAA_CDC_NCEP-DOE-2_sync_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    #-- permissions mode of the directories and files retrieved
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files retrieved')
    args = parser.parse_args()

    #-- run program for model
    if gravity_toolkit.utilities.check_ftp_connection('ftp.cdc.noaa.gov'):
        noaa_cdc_ncep_ftp(args.directory, YEAR=args.year, MASK=args.mask,
            INVARIANT=args.invariant, LOG=args.log, MODE=args.mode)
    else:
        raise RuntimeError('Check internet connection')

#-- run main program
if __name__ == '__main__':
    main()
