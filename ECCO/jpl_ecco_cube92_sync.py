#!/usr/bin/env python
u"""
jpl_ecco_cube92_sync.py
Written by Tyler Sutterley (04/2021)

Converts ECCO2 Cube92 daily model outputs from the NASA JPL ECCO2 server
    into monthly averages
https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/readme.txt
https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/PHIBOT.nc/

DATA DESCRIPTION:
    ECCO2 Cube92 employs the MITgcm in a global domain incorporating a
    Green's function approach.  Makes an optimal adjustment of model
    parameters, initial conditions, and boundary conditions through a
    series of sensitivity numerical simulations and minimization of
    model/observation misfit.

CALLING SEQUENCE:
    python jpl_ecco_cube92_sync.py --year 2015 2016 --user <username>
    where <username> is your NASA Earthdata username

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user X: username for NASA Earthdata Login
    -W X, --webdav X: WebDAV password for JPL ECCO Drive Login
    -N X, --netrc X: path to .netrc file for authentication
    -D X, --directory X: working data directory
    -Y X, --year X: Years to sync
    -l, --log: output log of files downloaded
    -V, --verbose: Output information for each output file
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
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://h5py.org
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    spatial.py: spatial data class for reading, writing and processing data
        ncdf_read.py: reads input spatial data from netCDF4 files
        hdf5_read.py: reads input spatial data from HDF5 files
        ncdf_write.py: writes output spatial data to netCDF4
        hdf5_write.py: writes output spatial data to HDF5
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 04/2021: set a default netrc file and check access
        default credentials from environmental variables
    Updated 03/2021: automatically update years to run based on current time
    Updated 01/2021 for public release.
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
import time
import netrc
import getpass
import netCDF4
import argparse
import builtins
import lxml.etree
import numpy as np
import gravity_toolkit.time
import gravity_toolkit.spatial
import gravity_toolkit.utilities

#-- PURPOSE: sync ECCO2 Cube92 model outputs from JPL ECCO drive server
#-- combines daily files to calculate monthly averages
def jpl_ecco_cube92_sync(ddir, YEAR=None, PRODUCT=None, LOG=False,
    VERBOSE=False, MODE=None):

    #-- check if directory exists and recursively create if not
    DIRECTORY = os.path.join(ddir, 'cube92_latlon_quart_90S90N')
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None

    #-- remote subdirectory for Cube92 data on JPL ECCO data server
    PATH = ['https://ecco.jpl.nasa.gov','drive','files','ECCO2',
        'cube92_latlon_quart_90S90N','{0}.nc'.format(PRODUCT)]
    #-- compile HTML parser for lxml
    parser = lxml.etree.HTMLParser()

    #-- create log file with list of synchronized files (or print to terminal)
    if LOG:
        #-- format: JPL_ECCO2_Cube92_PHIBOT_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        args = (PRODUCT, today)
        LOGFILE = 'JPL_ECCO2_Cube92_{0}_sync_{1}.log'.format(*args)
        fid1 = open(os.path.join(DIRECTORY,LOGFILE),'w')
        print('ECCO2 Cube92 {0} Sync Log ({1})'.format(today), file=fid1)
    else:
        #-- standard output (terminal output)
        fid1 = sys.stdout

    #-- regular expression for grouping months from daily data
    regex_pattern = '{0}\.(\d+)x(\d+)\.({1:4})({2:02d})(\d{{2}}).nc$'
    #-- input and output variable names
    LONNAME = 'LONGITUDE_T'
    LATNAME = 'LATITUDE_T'
    TIMENAME = 'TIME'
    #-- output netCDF4 keyword arguments
    kwargs = dict(varname=PRODUCT, latname=LATNAME, lonname=LONNAME,
        timename=TIMENAME, title="ECCO2 cube92 monthly average",
        units='m^2/s^2', time_units='days since 1992-01-01 00:00:00',
        time_longname='center time of averaging period')

    #-- for each year
    for YY in YEAR:
        #-- days per month in the year
        dpm = gravity_toolkit.time.calendar_days(YY)
        #-- for each month
        for MM in range(12):
            #-- compile regular expression pattern for year and month
            R1 = re.compile(regex_pattern.format(YY,MM+1), re.VERBOSE)
            #-- read and parse request for files (find names and modified dates)
            colnames,mtimes=gravity_toolkit.utilities.drive_list(PATH,
                timeout=120,build=False,parser=parser,pattern=R1,sort=True)
            #-- check if all files are available for the month
            if (len(colnames) != dpm[MM]):
                continue
            #-- python list with daily data
            daily = []
            #-- for each file in the month
            for remote_file,remote_mtime in zip(colnames,mtimes):
                #-- extract dimension variables
                dim1,dim2,Y,M,D = R1.findall(remote_file).pop()
                #-- Create and submit request to retrieve bytes
                response = gravity_toolkit.utilities.from_drive(
                    [*PATH,remote_file], build=False, verbose=VERBOSE,
                    fid=fid1, mode=MODE)
                #-- open remote file with netCDF4
                #-- remove singleton dimensions
                dinput = gravity_toolkit.spatial().from_netCDF4(response,
                    compression='bytes', latname=LATNAME, lonname=LONNAME,
                    varname=PRODUCT, timename=TIMENAME).squeeze()
                #-- replace fill value with missing value attribute
                dinput.fill_value = dinput.attributes['data']['missing_value']
                dinput.update_mask()
                #-- get variable longname attribute from daily file
                kwargs['longname'] = dinput.attributes['data']['long_name']
                #-- append to daily list
                daily.append(dinput)
            #-- calculate monthly mean from list of daily files
            monthly = gravity_toolkit.spatial().from_list(daily).mean()
            #-- output to netCDF4 file
            args = (PRODUCT,dim1,dim2,YY,MM+1)
            FILE = '{0}.{1}x{2}.{3}{4:02d}.nc'.format(*args)
            monthly.to_netCDF4(os.path.join(DIRECTORY,FILE),
                date=True, verbose=VERBOSE, **kwargs)
            #-- set permissions mode to MODE
            os.chmod(os.path.join(DIRECTORY,FILE), MODE)

    #-- close log file and set permissions level to MODE
    if LOG:
        fid1.close()
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

#-- Main program that calls jpl_ecco_cube92_sync()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Converts ECCO2 Cube92 daily model outputs
            from the NASA JPL ECCO2 server into monthly averages
            """
    )
    #-- command line parameters
    #-- NASA Earthdata credentials
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
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- ECCO model years to sync
    now = gravity_toolkit.time.datetime.datetime.now()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.year+1),
        help='Years to sync')
    #-- ECCO model product to sync
    parser.add_argument('--product', '-P',
        type=str, default='PHIBOT',
        help='Product to sync')
    #-- Output log file in form
    #-- JPL_ECCO2_Cube92_OBP_sync_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    #-- print information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the directories and files synced (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files synced')
    args = parser.parse_args()

    #-- JPL ECCO drive hostname
    HOST = 'ecco.jpl.nasa.gov'
    #-- get NASA Earthdata and JPL ECCO drive credentials
    if not args.user and not os.access(args.netrc,os.F_OK):
        #-- check that NASA Earthdata credentials were entered
        args.user=builtins.input('Username for {0}: '.format(HOST))
        #-- enter password securely from command-line
        args.webdav=getpass.getpass('Password for {0}@{1}: '.format(args.user,HOST))
    elif not args.user and os.access(args.netrc,os.F_OK):
        args.user,_,args.webdav=netrc.netrc(args.netrc).authenticators(HOST)
    else:
        #-- enter password securely from command-line
        args.webdav=getpass.getpass('Password for {0}@{1}: '.format(args.user,HOST))

    #-- build a urllib opener for JPL ECCO Drive
    #-- Add the username and password for NASA Earthdata Login system
    gravity_toolkit.utilities.build_opener(args.user,args.webdav)

    #-- check internet connection before attempting to run program
    #-- check JPL ECCO Drive credentials before attempting to run program
    if gravity_toolkit.utilities.check_credentials('https://{0}'.format(HOST)):
        jpl_ecco_cube92_sync(args.directory, YEAR=args.year,
            PRODUCT=args.product, LOG=args.log, VERBOSE=args.verbose,
            MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
