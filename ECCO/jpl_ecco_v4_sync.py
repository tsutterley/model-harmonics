#!/usr/bin/env python
u"""
jpl_ecco_v4_sync.py
Written by Tyler Sutterley (12/2020)

Syncs ECCO Ocean Bottom Pressure outputs from the NASA JPL ECCO Drive server:
https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/interp_monthly/README
https://ecco-group.org/products-ECCO-V4r4.htm
https://ecco-group.org/user-guide-v4r4.htm

https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
https://nsidc.org/support/faq/what-options-are-available-bulk-downloading-data-
    https-earthdata-login-enabled
http://www.voidspace.org.uk/python/articles/authentication.shtml#base64

Register with NASA Earthdata Login system:
https://urs.earthdata.nasa.gov

DATA DESCRIPTION:
    ECCO Version4 employs the MITgcm in a global domain incorporating most of
        the available satellite and in-situ data to produce a physically
        consistent ocean estimate. Model resolution is 0.5 degrees horizontal.
    Model fields are available at monthly-day intervals

CALLING SEQUENCE:
    python jpl_ecco_v4_sync.py --year=2015,2016 --user <username> V4r4
    where <username> is your NASA Earthdata username

INPUTS:
    ECCO version 4 models
        V4r3: Version 4, Revision 3
        V4r4: Version 4, Revision 4

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user X: username for NASA Earthdata Login
    -N X, --netrc X: path to .netrc file for authentication
    -D X, --directory X: working data directory
    -Y X, --year X: Years to sync
    -L, --list: print files to be transferred, but do not execute transfer
    -l, --log: output log of files downloaded
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
    utilities: download and management utilities for syncing files

UPDATE HISTORY:
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
    Written 06/2018
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
import argparse
import builtins
import posixpath
import lxml.etree
import gravity_toolkit.utilities

#-- PURPOSE: sync ECCO Ocean Bottom Pressure data from JPL ECCO drive server
def jpl_ecco_v4_sync(ddir, MODEL, YEAR=None, LOG=False, LIST=False,
    CLOBBER=False, CHECKSUM=False, MODE=None):

    #-- check if directory exists and recursively create if not
    DIRECTORY = os.path.join(ddir, 'ECCO-{0}'.format(MODEL))
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None

    #-- remote https server for ECCO data
    HOST = 'https://ecco.jpl.nasa.gov'
    #-- compile HTML parser for lxml
    parser = lxml.etree.HTMLParser()

    #-- create log file with list of synchronized files (or print to terminal)
    if LOG:
        #-- format: JPL_ECCO_V4r4_OBP_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'JPL_ECCO_{0}_OBP_{1}.log'.format(MODEL,today)
        fid1 = open(os.path.join(DIRECTORY,LOGFILE),'w')
        print('ECCO Version 4 OBP Sync Log ({0})'.format(today), file=fid1)
    else:
        #-- standard output (terminal output)
        fid1 = sys.stdout

    #-- print the model synchronized
    print('MODEL: {0}\n'.format(MODEL), file=fid1)

    #-- path to model files
    model_path = {}
    model_path['V4r3'] = ['Version4','Release3','interp_monthly','PHIBOT']
    model_path['V4r4'] = ['Version4','Release4','interp_monthly','PHIBOT']
    #-- compile regular expression operator for years to sync
    if YEAR is None:
        regex_years = r'\d+'
    else:
        regex_years = r'|'.join('{0:d}'.format(y) for y in YEAR)
    #-- compile regular expression operator finding years
    if MODEL in ('V4r3',):
        R1 = re.compile(r'PHIBOT([\.\_])({0})(_\d+)?.nc$'.format(regex_years))
    elif MODEL in ('V4r4',):
        R1 = re.compile(regex_years)

    #-- remote subdirectory for MODEL on JPL ECCO data server
    PATH = [HOST,'drive','files',*model_path[MODEL]]
    #-- open connection with ECCO drive server at remote directory
    #-- find remote yearly directories for MODEL
    years,mtimes = gravity_toolkit.utilities.drive_list(PATH,
        timeout=120,build=False,parser=parser,pattern=R1,sort=True)
    for yr in years:
        #-- extract year from file
        if MODEL in ('V4r3',):
            _,YY,_ = R1.findall(yr).pop()
        elif MODEL in ('V4r4',):
            #-- print string for year
            print(yr, file=fid1)
            #-- add the year directory to the path
            YY, = R1.findall(yr)
            PATH.append(yr)
        #-- compile regular expression operator for ocean bottom pressure files
        R3 = re.compile(r'PHIBOT([\.\_])({0})(_\d+)?.nc$'.format(YY))
        #-- full path to remote directory
        remote_dir = posixpath.join(*PATH)
        #-- read and parse request for files (find names and modified dates)
        colnames,mtimes=gravity_toolkit.utilities.drive_list(PATH,
            timeout=120,build=False,parser=parser,pattern=R3,sort=True)
        #-- for each file on the remote server
        for colname,remote_mtime in zip(colnames,mtimes):
            #-- remote and local versions of the file
            remote_file = posixpath.join(remote_dir,colname)
            local_file = os.path.join(DIRECTORY,colname)
            http_pull_file(fid1, remote_file, remote_mtime,
                local_file, LIST, CLOBBER, CHECKSUM, MODE)
        #-- remove the year directory from the path
        if MODEL in ('V4r4',):
            PATH.remove(yr)

    #-- close log file and set permissions level to MODE
    if LOG:
        fid1.close()
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

#-- PURPOSE: pull file from a remote host checking if file exists locally
#-- and if the remote file is newer than the local file
def http_pull_file(fid, remote_file, remote_mtime, local_file, LIST, CLOBBER,
    CHECKSUM, MODE):
    #-- if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    #-- check if local version of file exists
    if CHECKSUM and os.access(local_file, os.F_OK):
        #-- generate checksum hash for local file
        #-- open the local_file in binary read mode
        local_hash = gravity_toolkit.utilities.get_hash(local_file)
        #-- Create and submit request.
        #-- There are a wide range of exceptions that can be thrown here
        #-- including HTTPError and URLError.
        req=gravity_toolkit.utilities.urllib2.Request(remote_file)
        resp=gravity_toolkit.utilities.urllib2.urlopen(req,timeout=120)
        #-- copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO(resp.read())
        remote_buffer.seek(0)
        #-- generate checksum hash for remote file
        remote_hash = gravity_toolkit.utilities.get_hash(remote_buffer)
        #-- compare checksums
        if (local_hash != remote_hash):
            TEST = True
            OVERWRITE = ' (checksums: {0} {1})'.format(local_hash,remote_hash)
    elif os.access(local_file, os.F_OK):
        #-- check last modification time of local file
        local_mtime = os.stat(local_file).st_mtime
        #-- if remote file is newer: overwrite the local file
        if (remote_mtime > local_mtime):
            TEST = True
            OVERWRITE = ' (overwrite)'
    else:
        TEST = True
        OVERWRITE = ' (new)'
    #-- if file does not exist locally, is to be overwritten, or CLOBBER is set
    if TEST or CLOBBER:
        #-- Printing files transferred
        print('{0} --> '.format(remote_file), file=fid)
        print('\t{0}{1}\n'.format(local_file,OVERWRITE), file=fid)
        #-- if executing copy command (not only printing the files)
        if not LIST:
            #-- chunked transfer encoding size
            CHUNK = 16 * 1024
            #-- copy bytes or transfer file
            if CHECKSUM and os.access(local_file, os.F_OK):
                #-- store bytes to file using chunked transfer encoding
                remote_buffer.seek(0)
                with open(local_file, 'wb') as f:
                    shutil.copyfileobj(remote_buffer, f, CHUNK)
            else:
                #-- Create and submit request.
                #-- There are a wide range of exceptions that can be thrown here
                #-- including HTTPError and URLError.
                req=gravity_toolkit.utilities.urllib2.Request(remote_file)
                resp=gravity_toolkit.utilities.urllib2.urlopen(req,timeout=120)
                #-- copy contents to local file using chunked transfer encoding
                #-- transfer should work properly with ascii and binary formats
                with open(local_file, 'wb') as f:
                    shutil.copyfileobj(resp, f, CHUNK)
            #-- keep remote modification time of file and local access time
            os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
            os.chmod(local_file, MODE)

#-- Main program that calls jpl_ecco_sync()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Syncs ECCO Ocean Bottom Pressure outputs from the
        NASA JPL ECCO Drive server
        """
    )
    #-- command line parameters
    parser.add_argument('model',
        metavar='MODEL', type=str, nargs='+',
        default=['V4r3','V4r4'], choices=['V4r3','V4r4'],
        help='ECCO Model')
    #-- NASA Earthdata credentials
    parser.add_argument('--user','-U',
        type=str, default='',
        help='Username for NASA Earthdata Login')
    parser.add_argument('--netrc','-N',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Path to .netrc file for authentication')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- ECCO model years to sync
    parser.add_argument('--year','-Y',
        type=int, nargs='+',
        help='Years to sync')
    #-- Output log file in form
    #-- JPL_ECCO_kf080i_OBP_sync_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    #-- sync options
    parser.add_argument('--list','-L',
        default=False, action='store_true',
        help='Only print files that could be transferred')
    parser.add_argument('--checksum',
        default=False, action='store_true',
        help='Compare hashes to check for overwriting existing data')
    parser.add_argument('--clobber','-C',
        default=False, action='store_true',
        help='Overwrite existing data in transfer')
    #-- permissions mode of the directories and files synced (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files synced')
    args = parser.parse_args()

    #-- JPL ECCO drive hostname
    HOST = 'ecco.jpl.nasa.gov'
    #-- get NASA Earthdata and JPL ECCO drive credentials
    if not args.user and not args.netrc:
        #-- check that NASA Earthdata credentials were entered
        args.user=builtins.input('Username for {0}: '.format(HOST))
        #-- enter password securely from command-line
        PASSWORD=getpass.getpass('Password for {0}@{1}: '.format(args.user,HOST))
    elif args.netrc:
        args.user,LOGIN,PASSWORD=netrc.netrc(args.netrc).authenticators(HOST)
    else:
        #-- enter password securely from command-line
        PASSWORD=getpass.getpass('Password for {0}@{1}: '.format(args.user,HOST))

    #-- build a urllib opener for JPL ECCO Drive
    #-- Add the username and password for NASA Earthdata Login system
    gravity_toolkit.utilities.build_opener(args.user,PASSWORD)

    #-- check internet connection before attempting to run program
    #-- check JPL ECCO Drive credentials before attempting to run program
    if gravity_toolkit.utilities.check_credentials('https://{0}'.format(HOST)):
        for MODEL in args.model:
            jpl_ecco_v4_sync(args.directory, MODEL, YEAR=args.year,
                LIST=args.list, LOG=args.log, CLOBBER=args.clobber,
                CHECKSUM=args.checksum, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
