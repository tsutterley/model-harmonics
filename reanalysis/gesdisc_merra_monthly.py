#!/usr/bin/env python
u"""
gesdisc_merra_monthly.py
Written by Tyler Sutterley (12/2020)

Downloads MERRA-2 products using a links list provided by the Goddard Earth
    Sciences Data and Information Server Center (GES DISC)
    https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/
    https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
Combines daily model level outputs into monthly averages

Register with NASA Earthdata Login system:
    https://urs.earthdata.nasa.gov

Add "NASA GESDISC DATA ARCHIVE" to Earthdata Applications:
    https://urs.earthdata.nasa.gov/approve_app?client_id=e2WVk8Pw6weeLUKZYOxvTQ

CALLING SEQUENCE:
    python gesdisc_merra_monthly.py --user=<username> links_list_file
    where <username> is your NASA Earthdata username

INPUTS:
    links_list_file: GES DISC generated file listing files to download

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user X: username for NASA Earthdata Login
    -N X, --netrc X: path to .netrc file for authentication
    -D X, --directory X: Working data directory
    -l, --log: output log of files downloaded
    -V, --verbose: Output information for each output file
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
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 12/2020: use argparse to set command line parameters
        using utilities program to build opener
    Updated 09/2019: added ssl context to urlopen headers
        modified regular expression to check if suffix of file is nc4
    Updated 08/2019: new GESDISC server and links list file format
    Updated 06/2018: using python3 compatible octal, input and urllib
    Written 03/2018
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
import numpy as np
import gravity_toolkit.utilities

#-- PURPOSE: sync local MERRA-2 files with GESDISC server
def gesdisc_merra_monthly(base_dir, links_list_file, LOG=False,
    VERBOSE=False, MODE=None):
    #-- full path to MERRA-2 directory
    DIRECTORY = os.path.join(base_dir,'MERRA-2')
    #-- check if DIRECTORY exists and recursively create if not
    if (not os.access(os.path.join(DIRECTORY), os.F_OK)):
        os.makedirs(os.path.join(DIRECTORY), MODE)

    #-- create log file with list of synchronized files (or print to terminal)
    if LOG:
        #-- format: NASA_GESDISC_MERRA2_monthly_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'NASA_GESDISC_MERRA2_monthly_{0}.log'.format(today)
        fid = open(os.path.join(DIRECTORY,LOGFILE),'w')
        print('NASA MERRA-2 Sync Log ({0})'.format(today), file=fid)
    else:
        #-- standard output (terminal output)
        fid = sys.stdout

    #-- read the links list file
    with open(links_list_file,'rb') as fileID:
        lines = fileID.read().decode("utf-8-sig").encode("utf-8").splitlines()

    #-- regular expression for grouping months from daily data
    regex_pattern = rb'MERRA2_(\d+).(.*?).(\d{4})(\d{2})(\d{2})(.*?).nc[4]?'
    rx1 = re.compile(regex_pattern, re.VERBOSE)

    #-- output variable names
    VARNAME = 'PS'
    TNAME = 'T'
    QNAME = 'QV'
    LONNAME = 'lon'
    LATNAME = 'lat'
    LEVELNAME = 'lev'
    TIMENAME = 'time'
    TIME_LONGNAME = 'time'
    #-- output dimensions
    nlevels,nlat,nlon = (72,361,576)

    #-- output arrays with year and month for each file
    valid_lines = []
    valid_files = []
    year = []
    month = []
    #-- for each line in the links_list_file
    for i,remote_file in enumerate(lines):
        #-- extract filename from url
        if re.search(rb'LABEL\=(.*?)\&SHORTNAME',remote_file):
            match_object = re.match(rb'LABEL\=(.*?)\&SHORTNAME', remote_file)
            valid_lines.append(remote_file.strip().decode('utf-8'))
            valid_files.append(match_object.group(0).strip().decode('utf-8'))
            year.append(match_object.group(3))
            month.append(match_object.group(4))
        elif rx1.search(remote_file):
            search_object = rx1.search(remote_file)
            valid_lines.append(remote_file.strip().decode('utf-8'))
            valid_files.append(search_object.group(0).strip().decode('utf-8'))
            year.append(search_object.group(3).decode('utf-8'))
            month.append(search_object.group(4).decode('utf-8'))

    #-- for each unique date
    DATES = sorted(set(zip(year,month)))
    for YY,MM in DATES:
        #-- regular expression for finding data in year and month
        regex_pattern = r'MERRA2_(\d+).(.*?).({0})({1})(\d{{2}})(.*?).nc[4]?$'
        rx2 = re.compile(regex_pattern.format(YY,MM),re.VERBOSE)
        remote_lines = [i for i,fi in enumerate(valid_files) if rx2.search(fi)]
        #-- python dictionary with output data
        dinput = {}
        dinput[TIMENAME] = np.zeros((1))
        dinput[VARNAME] = np.zeros((1,nlat,nlon))
        dinput[TNAME] = np.zeros((1,nlevels,nlat,nlon))
        dinput[QNAME] = np.zeros((1,nlevels,nlat,nlon))
        #-- python dictionary with count for converting totals to means
        count = {}
        count[TIMENAME] = np.zeros((1))
        count[VARNAME] = np.zeros((1,nlat,nlon))
        count[TNAME] = np.zeros((1,nlevels,nlat,nlon))
        count[QNAME] = np.zeros((1,nlevels,nlat,nlon))
        #-- for each file in the month
        for i in remote_lines:
            #-- Create and submit request. There are a wide range of exceptions
            #-- that can be thrown here, including HTTPError and URLError.
            URL = gravity_toolkit.utilities.url_split(valid_lines[i])
            response = gravity_toolkit.utilities.from_http(URL, context=None,
                verbose=VERBOSE, fid=fid)
            response.seek(0)
            #-- open remote file with netCDF4
            fileID = netCDF4.Dataset(valid_files[i],'r',memory=response.read())
            MOD,DATASET,Y,M,D,AUX = rx2.findall(valid_files[i]).pop()
            #-- extract dimension variables
            dinput[LEVELNAME] = fileID.variables[LEVELNAME][:].copy()
            dinput[LATNAME] = fileID.variables[LATNAME][:].copy()
            dinput[LONNAME] = fileID.variables[LONNAME][:].copy()
            nt, = fileID.variables[TIMENAME].shape
            #-- bad value
            fill_value = fileID.variables[VARNAME]._FillValue
            #-- add over time slices products to monthly output
            for t in range(nt):
                dinput[TIMENAME][0] += fileID.variables[TIMENAME][t].astype('f')
                count[TIMENAME][0] += 1.0
                #-- surface pressure
                PS = fileID.variables[VARNAME][t,:,:].copy()
                ii,jj = np.nonzero(PS != fill_value)
                dinput[VARNAME][0,ii,jj] += PS[ii,jj]
                count[VARNAME][0,ii,jj] += 1.0
                #-- air temperature
                T = fileID.variables[TNAME][t,:,:,:].copy()
                ii,jj,kk = np.nonzero(T != fill_value)
                dinput[TNAME][0,ii,jj,kk] += T[ii,jj,kk]
                count[TNAME][0,ii,jj,kk] += 1.0
                #-- specific humidity
                QV = fileID.variables[QNAME][t,:,:,:].copy()
                ii,jj,kk = np.nonzero(QV != fill_value)
                dinput[QNAME][0,ii,jj,kk] += QV[ii,jj,kk]
                count[QNAME][0,ii,jj,kk] += 1.0
            #-- close the input file from remote url
            fileID.close()
        #-- calculate mean from totals
        dinput[TIMENAME] /= count[TIMENAME]
        for key in [VARNAME,TNAME,QNAME]:
            #-- find valid values
            valid_indices = np.nonzero(count[key] > 0)
            dinput[key][valid_indices] /= count[key][valid_indices]
            #-- replace points where no values with fill_value
            complementary_indices = np.nonzero(count[key] == 0)
            dinput[key][complementary_indices] = fill_value
        #-- output to netCDF4 file (replace hour variable with monthly)
        DATASET = re.sub(r'inst\d+_3d',r'instM_3d',DATASET)
        TIME_UNITS = 'minutes since {0}-{1}-01 00:00:00'.format(YY,MM)
        local_file = 'MERRA2_{0}.{1}.{2}{3}.SUB.nc'.format(MOD,DATASET,YY,MM)
        ncdf_model_write(dinput, fill_value, VARNAME=VARNAME, TNAME=TNAME,
            QNAME=QNAME, LONNAME=LONNAME, LATNAME=LATNAME, LEVELNAME=LEVELNAME,
            TIMENAME=TIMENAME,TIME_UNITS=TIME_UNITS,TIME_LONGNAME=TIME_LONGNAME,
            FILENAME=os.path.join(DIRECTORY,local_file), VERBOSE=VERBOSE)
        #-- set permissions mode to MODE
        os.chmod(os.path.join(DIRECTORY,local_file), MODE)

    #-- close log file and set permissions level to MODE
    if LOG:
        fid.close()
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

#-- PURPOSE: write output model layer fields data to file
def ncdf_model_write(dinput, fill_value, VARNAME=None, TNAME=None, QNAME=None,
    LONNAME=None, LATNAME=None, LEVELNAME=None, TIMENAME=None, TIME_UNITS=None,
    TIME_LONGNAME=None, FILENAME=None, VERBOSE=False):
    #-- opening NetCDF4 file for writing
    fileID = netCDF4.Dataset(FILENAME, 'w', format="NETCDF4")

    #-- Defining the NetCDF4 dimensions and creating dimension variables
    nc = {}
    for key in [LONNAME,LATNAME,TIMENAME,LEVELNAME]:
        fileID.createDimension(key, len(dinput[key]))
        nc[key] = fileID.createVariable(key,dinput[key].dtype,(key,))
    #-- creating the layered NetCDF4 variables
    for key in [TNAME,QNAME]:
        nc[key] = fileID.createVariable(key, dinput[key].dtype,
            (TIMENAME,LEVELNAME,LATNAME,LONNAME,), fill_value=fill_value, zlib=True)
    #-- creating the surface NetCDF4 variables
    for key in [VARNAME]:
        nc[key] = fileID.createVariable(key, dinput[key].dtype,
            (TIMENAME,LATNAME,LONNAME,), fill_value=fill_value, zlib=True)

    #-- filling NetCDF4 variables
    for key,val in dinput.items():
        nc[key][:] = val.copy()

    #-- Defining attributes for longitude and latitude
    nc[LONNAME].long_name = 'Longitude'
    nc[LONNAME].units = 'degrees_east'
    nc[LATNAME].long_name = 'Latitude'
    nc[LATNAME].units = 'degrees_north'
    #-- Defining attributes for time
    nc[TIMENAME].units = TIME_UNITS
    nc[TIMENAME].long_name = TIME_LONGNAME
    #-- Definining attributes for model levels
    nc[LEVELNAME].long_name = 'vertical_level'
    nc[LEVELNAME].units = 'layer'
    nc[LEVELNAME].positive = 'down'
    #-- Defining attributes for air temperature
    nc[TNAME].long_name = 'Air_Temperature'
    nc[TNAME].units = 'K'
    #-- Defining attributes for specific humidity
    nc[QNAME].long_name = 'Specific_Humidity'
    nc[QNAME].units = 'kg/kg'
    #-- Defining attributes for surface pressure
    nc[VARNAME].long_name = 'Surface_Air_Pressure'
    nc[VARNAME].units = 'Pa'
    #-- date created
    fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())

    #-- Output NetCDF structure information
    if VERBOSE:
        print(os.path.basename(FILENAME))
        print(list(fileID.variables.keys()))

    #-- Closing the NetCDF file
    fileID.close()

#-- Main program that calls gesdisc_merra_monthly()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Downloads MERRA-2 products using a links list
            provided by the Goddard Earth Sciences Data and Information
            Server Center (GES DISC)
            """
    )
    #-- command line parameters
    parser.add_argument('file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='GESDISC links list file')
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
    #-- Output log file in form
    #-- NASA_GESDISC_MERRA2_monthly_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    #-- sync options
    parser.add_argument('--list','-L',
        default=False, action='store_true',
        help='Only print files that could be transferred')
    #-- print information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
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
        args.user,_,PASSWORD=netrc.netrc(args.netrc).authenticators(URS)
    else:
        #-- enter password securely from command-line
        PASSWORD=getpass.getpass('Password for {0}@{1}: '.format(args.user,URS))

    #-- build a urllib opener for NASA GESDISC
    #-- Add the username and password for NASA Earthdata Login system
    gravity_toolkit.utilities.build_opener(args.user, PASSWORD,
        password_manager=True, authorization_header=False)

    #-- check internet connection before attempting to run program
    HOST = 'http://disc.sci.gsfc.nasa.gov/'
    if gravity_toolkit.utilities.check_connection(HOST):
        #-- for each links list file from GESDISC
        for FILE in args.file:
            gesdisc_merra_monthly(args.directory, FILE, LOG=args.log,
                VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
