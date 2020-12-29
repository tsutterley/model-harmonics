#!/usr/bin/env python
u"""
gesdisc_gldas_sync.py
Written by Tyler Sutterley (10/2020)

This program syncs GLDAS monthly datafiles from the Goddard Earth Sciences Data
    and Information Server Center
    http://ldas.gsfc.nasa.gov/gldas/
    http://disc.sci.gsfc.nasa.gov/hydrology/documentation
    https://hydro1.gesdisc.eosdis.nasa.gov/data/GLDAS_V1/README.GLDAS.pdf
    https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python

Register with NASA Earthdata Login system:
    https://urs.earthdata.nasa.gov

Add "NASA GESDISC DATA ARCHIVE" to Earthdata Applications:
    https://urs.earthdata.nasa.gov/approve_app?client_id=e2WVk8Pw6weeLUKZYOxvTQ

CALLING SEQUENCE:
    python gesdisc_gldas_sync.py --user=<username>
    where <username> is your NASA Earthdata username

INPUTS:
    CLM: GLDAS Common Land Model (CLM)
    MOS: GLDAS Mosaic model
    NOAH: GLDAS Noah model
    VIC: GLDAS Variable Infiltration Capacity (VIC) model

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory: working data directory
    -Y X, --year X: years to sync separated by commas
    -S X, --spacing X: spatial resolution of models to sync
        10: 1.0 degrees latitude/longitude
        025: 0.25 degrees latitude/longitude
    -T X, --temporal X: temporal resolution of models to sync
        M: Monthly
        3H: 3-hourly
    -v X, --version X: GLDAS model version
    -e, --early: GLDAS early products
    -U X, --user X: username for NASA Earthdata Login
    -N X, --netrc X: path to .netrc file for authentication
    -M X, --mode X: permissions mode of the directories and files synced
    --log: output log of files downloaded
    --list: print files to be transferred, but do not execute transfer
    --clobber: Overwrite existing data in transfer

PYTHON DEPENDENCIES:
    lxml: Pythonic XML and HTML processing library using libxml2/libxslt
        http://lxml.de/
        https://github.com/lxml/lxml
    future: Compatibility layer between Python 2 and Python 3
        http://python-future.org/

UPDATE HISTORY:
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
import argparse
import builtins
import posixpath
import lxml.etree
import gravity_toolkit.utilities

#-- GLDAS models
gldas_products = {}
gldas_products['CLM'] = 'GLDAS Common Land Model (CLM)'
gldas_products['CLSM'] = 'GLDAS Catchment Land Surface Model (CLSM)'
gldas_products['MOS'] = 'GLDAS Mosaic model'
gldas_products['NOAH'] = 'GLDAS Noah model'
gldas_products['VIC'] = 'GLDAS Variable Infiltration Capacity (VIC) model'

#-- PURPOSE: sync local GLDAS files with GESDISC server
def gesdisc_gldas_sync(DIRECTORY, MODEL, YEARS, SPATIAL='', TEMPORAL='',
    VERSION='', EARLY=False, LOG=False, LIST=False, MODE=0o775, CLOBBER=False):

    #-- check if directory exists and recursively create if not
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None
    #-- create log file with list of synchronized files (or print to terminal)
    if LOG:
        #-- output to log file
        #-- format: NASA_GESDISC_GLDAS_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'NASA_GESDISC_GLDAS_{0}_sync_{1}.log'.format(MODEL,today)
        fid = open(os.path.join(DIRECTORY,LOGFILE),'w')
        print('NASA GLDAS Sync Log ({0})'.format(today), file=fid)
    else:
        #-- standard output (terminal output)
        fid = sys.stdout

    #-- Version flags
    V1,V2 = ('_V1','') if (VERSION == '1') else ('','.{0}'.format(VERSION))
    EP = '_EP' if EARLY else ''
    #-- GLDAS model remote base directory
    HOST = ['https://hydro1.gesdisc.eosdis.nasa.gov','data','GLDAS{0}'.format(V1)]

    #-- compile regular expression operator for years to sync
    regex_pattern = '|'.join('{0:d}'.format(y) for y in YEARS)
    R1 = re.compile(r'({0})'.format(regex_pattern), re.VERBOSE)

    #-- print header text to log/standard output
    print('GLDAS MODEL={0}'.format(gldas_products[MODEL]), file=fid)
    print('RESOLUTION={0},{1}'.format(TEMPORAL,SPATIAL), file=fid)
    print('VERSION={0}{1}'.format(VERSION,EP), file=fid)
    #-- subdirectory for model on GESDISC server
    REMOTE = "GLDAS_{0}{1}_{2}{3}{4}".format(MODEL,SPATIAL,TEMPORAL,EP,V2)
    PRODUCT = "GLDAS_{0}{1}_{2}{3}".format(MODEL,SPATIAL,TEMPORAL,V2)

    #-- open connection with GESDISC server at remote directory
    #-- find remote yearly directories for MODEL
    remote_years,mtimes = gesdisc_list([*HOST,REMOTE],pattern=R1,sort=True)
    #-- compile regular expression operator for model on GESDISC server
    args = (MODEL,SPATIAL)
    R2 = re.compile(r'GLDAS_{0}{1}_(.*?).(nc4|grb)(.xml)?$'.format(*args))
    #-- if running monthly data
    if (TEMPORAL == 'M'):
        #-- for each yearly subdirectory
        for Y in remote_years:
            #-- check if local directory exists and recursively create if not
            if (not os.access(os.path.join(DIRECTORY,PRODUCT,Y), os.F_OK)):
                os.makedirs(os.path.join(DIRECTORY,PRODUCT,Y), MODE)
            #-- open connection with GESDISC server at remote directory
            #-- read and parse request for files (names and modified dates)
            #-- find remote files for MODEL and YEAR
            files,mtimes = gesdisc_list([*HOST,REMOTE,Y],pattern=R2,sort=True)
            for colname,remote_mtime in zip(files,mtimes):
                #-- local and remote versions of the file
                local_file = os.path.join(DIRECTORY,PRODUCT,Y,colname)
                remote_file = posixpath.join(*HOST,REMOTE,Y,colname)
                #-- copy file from remote directory comparing modified dates
                http_pull_file(fid, remote_file, remote_mtime, local_file,
                    LIST, CLOBBER, MODE)
    #-- if running daily data
    elif (TEMPORAL == '3H'):
        #-- for each yearly subdirectory
        for Y in remote_years:
            #-- compile regular expression operator for days
            R3 = re.compile(r'\d+',re.VERBOSE)
            #-- open connection with GESDISC server at remote directory
            #-- find remote daily directories for MODEL
            remote_days,mtimes = gesdisc_list([*HOST,REMOTE,Y],pattern=R3,sort=True)
            #-- for each daily subdirectory
            for D in remote_days:
                #-- check if local directory exists and recursively create if not
                if (not os.access(os.path.join(DIRECTORY,PRODUCT,Y,D),os.F_OK)):
                    os.makedirs(os.path.join(DIRECTORY,PRODUCT,Y,D), MODE)
                #-- open connection with GESDISC server at remote directory
                #-- read and parse request for files (names and modified dates)
                #-- find remote files for MODEL and YEAR
                files,mtimes = gesdisc_list([*HOST,REMOTE,Y,D],pattern=R2,sort=True)
                for colname,remote_mtime in zip(files,mtimes):
                    #-- local and remote versions of the file
                    local_file = os.path.join(DIRECTORY,PRODUCT,Y,D,colname)
                    remote_file = posixpath.join(*HOST,REMOTE,Y,D,colname)
                    #-- copy file from remote directory comparing modified dates
                    http_pull_file(fid, remote_file, remote_mtime, local_file,
                        LIST, CLOBBER, MODE)

    #-- close log file and set permissions level to MODE
    if LOG:
        fid.close()
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

#-- PURPOSE: list a directory on NASA GESDISC https server
def gesdisc_list(HOST,timeout=None,parser=lxml.etree.HTMLParser(),
    pattern='',sort=False):
    """
    List a directory on NASA GESDISC

    Arguments
    ---------
    HOST: remote https host path split as list

    Keyword arguments
    -----------------
    timeout: timeout in seconds for blocking operations
    parser: HTML parser for lxml
    pattern: regular expression pattern for reducing list
    sort: sort output list

    Returns
    -------
    colnames: list of column names in a directory
    collastmod: list of last modification times for items in the directory
    """
    #-- try listing from https
    try:
        #-- Create and submit request.
        request=gravity_toolkit.utilities.urllib2.Request(posixpath.join(*HOST))
        response=gravity_toolkit.utilities.urllib2.urlopen(request,timeout=timeout)
    except (gravity_toolkit.utilities.urllib2.HTTPError,
        gravity_toolkit.utilities.urllib2.URLError):
        raise Exception('List error from {0}'.format(posixpath.join(*HOST)))
    else:
        #-- read and parse request for files (column names and modified times)
        tree = lxml.etree.parse(response,parser)
        colnames = tree.xpath('//tr/td[not(@*)]//a/@href')
        #-- get the Unix timestamp value for a modification time
        lastmod=[gravity_toolkit.utilities.get_unix_time(i,format='%Y-%m-%d %H:%M')
            for i in tree.xpath('//tr/td[@align="right"][1]/text()')]
        #-- reduce using regular expression pattern
        if pattern:
            i = [i for i,f in enumerate(colnames) if re.search(pattern,f)]
            #-- reduce list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            lastmod = [lastmod[indice] for indice in i]
        #-- sort the list
        if sort:
            i = [i for i,j in sorted(enumerate(colnames), key=lambda i: i[1])]
            #-- sort list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            lastmod = [lastmod[indice] for indice in i]
        #-- return the list of column names and last modified times
        return (colnames,lastmod)

#-- PURPOSE: pull file from a remote host checking if file exists locally
#-- and if the remote file is newer than the local file
def http_pull_file(fid,remote_file,remote_mtime,local_file,LIST,CLOBBER,MODE):
    #-- if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    #-- check if local version of file exists
    if os.access(local_file, os.F_OK):
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
            #-- Create and submit request. There are a wide range of exceptions
            #-- that can be thrown here, including HTTPError and URLError.
            request = gravity_toolkit.utilities.urllib2.Request(remote_file)
            response = gravity_toolkit.utilities.urllib2.urlopen(request)
            #-- chunked transfer encoding size
            CHUNK = 16 * 1024
            #-- copy contents to local file using chunked transfer encoding
            #-- transfer should work properly with ascii and binary data formats
            with open(local_file, 'wb') as f:
                shutil.copyfileobj(response, f, CHUNK)
            #-- keep remote modification time of file and local access time
            os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
            os.chmod(local_file, MODE)

#-- Main program that calls gesdisc_gldas_sync()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Syncs GLDAS monthly datafiles from the Goddard Earth
            Sciences Data and Information Server Center (GES DISC)
            """
    )
    #-- command line parameters
    parser.add_argument('model',
        type=str, nargs='+', choices=gldas_products.keys(),
        help='GLDAS land surface model')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- years to download
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,2021),
        help='Years of model outputs to sync')
    #-- GLDAS model version
    parser.add_argument('--version','-v',
        type=str, default='2.1',
        help='GLDAS model version')
    #-- model spatial resolution
    #-- 10: 1.0 degrees latitude/longitude
    #-- 025: 0.25 degrees latitude/longitude
    parser.add_argument('--spacing','-S',
        type=str, default='10', choices=['10','025'],
        help='Spatial resolution of models to sync')
    #-- model temporal resolution
    #-- M: Monthly data products
    #-- 3H: 3-hourly data products
    parser.add_argument('--temporal','-T',
        type=str, default='M', choices=['M','3H'],
        help='Temporal resolution of models to sync')
    #-- GLDAS early products
    parser.add_argument('--early','-e',
        default=False, action='store_true',
        help='Sync GLDAS early products')
    #-- NASA Earthdata credentials
    parser.add_argument('--user','-U',
        type=str, default='',
        help='Username for NASA Earthdata Login')
    parser.add_argument('--netrc','-N',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Path to .netrc file for authentication')
    #-- Output log file in form
    #-- NASA_GESDISC_GLDAS_sync_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    #-- sync options
    parser.add_argument('--list','-L',
        default=False, action='store_true',
        help='Only print files that could be transferred')
    parser.add_argument('--clobber','-C',
        default=False, action='store_true',
        help='Overwrite existing data in transfer')
    #-- permissions mode of the directories and files synced (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files synced')
    args = parser.parse_args()

    #-- NASA Earthdata hostname
    URS = 'urs.earthdata.nasa.gov'
    #-- get NASA Earthdata credentials
    if not args.user and not args.netrc:
        #-- check that NASA Earthdata credentials were entered
        args.user = builtins.input('Username for {0}: '.format(URS))
        #-- enter password securely from command-line
        PASSWORD=getpass.getpass('Password for {0}@{1}: '.format(args.user,URS))
    elif args.netrc:
        args.user,_,PASSWORD = netrc.netrc(args.netrc).authenticators(URS)
    else:
        #-- enter password securely from command-line
        PASSWORD=getpass.getpass('Password for {0}@{1}: '.format(args.user,URS))

    #-- build a urllib opener for NASA GESDISC
    #-- Add the username and password for NASA Earthdata Login system
    gravity_toolkit.utilities.build_opener(args.user, PASSWORD,
        password_manager=True, authorization_header=False)

    #-- check internet connection before attempting to run program
    HOST = posixpath.join('https://hydro1.gesdisc.eosdis.nasa.gov','data')
    if gravity_toolkit.utilities.check_connection(HOST):
        #-- for each GLDAS model
        for MODEL in args.model:
            gesdisc_gldas_sync(args.directory, MODEL, args.year,
                VERSION=args.version, EARLY=args.early,
                SPATIAL=args.spacing, TEMPORAL=args.temporal,
                LOG=args.log, LIST=args.list, CLOBBER=args.clobber,
                MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
