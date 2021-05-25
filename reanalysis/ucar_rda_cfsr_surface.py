#!/usr/bin/env python
u"""
ucar_rda_cfsr_surface.py
Written by Tyler Sutterley (05/2021)

Downloads NCEP-CFSR products using a links list provided by the
    NCAR/UCAR Research Data Archive (RDA): https://rda.ucar.edu/

NCEP Climate Forecast System Reanalysis (CFSR) Monthly Products
    https://rda.ucar.edu/datasets/ds093.2/
NCEP Climate Forecast System Version 2 (CFSv2) Monthly Products
    https://rda.ucar.edu/datasets/ds094.2/

Parameters for links list file:
    Pressure
    Subset regular monthly means converted to netCDF4
    Monthly Mean (4 per day) of Analyses
    Ground or Water Surface

CALLING SEQUENCE:
    python ucar_rda_cfsr_surface.py --user <username> links_list_file

INPUTS:
    links_list_file: UCAR links list file (csh)

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user X: Username for UCAR/NCAR RDA login
    -P X, --password X: Password for UCAR/NCAR RDA login
    -N X, --netrc X: Path to .netrc file for authentication
    -D X, --directory X: Full path to output directory
    -Y X, --year X: Years to download from input links file
    -I, --isentropic: Input data is over isentropic levels
    -G, --gzip: Input data is compressed
    -t X, --timeout X: Timeout in seconds for blocking operations
    -l, --log: Output log of files downloaded
    -M X, --mode=X: Permission mode of directories and files downloaded

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        https://www.h5py.org/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    utilities.py: download and management utilities for files
    spatial.py: spatial data class for reading, writing and processing data
        ncdf_read.py: reads input spatial data from netCDF4 files
        hdf5_read.py: reads input spatial data from HDF5 files
        ncdf_write.py: writes output spatial data to netCDF4
        hdf5_write.py: writes output spatial data to HDF5

UPDATE HISTORY:
    Updated 05/2021: added option for connection timeout (in seconds)
    Updated 04/2021: set a default netrc file and check access
        default credentials from environmental variables
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 12/2020: using time, spatial and utilities modules
    Written 08/2019
"""
from __future__ import print_function

import sys
import os
import re
import io
import time
import gzip
import netrc
import netCDF4
import getpass
import builtins
import argparse
import numpy as np
import gravity_toolkit.time
import gravity_toolkit.utilities
from gravity_toolkit.spatial import spatial

#-- PURPOSE: sync local NCEP-CFSR files with UCAR/NCAR RDA server
def ucar_rda_download(links_list_file, DIRECTORY=None, YEARS=None,
    ISENTROPIC=False, GZIP=False, TIMEOUT=None, LOG=False, MODE=None):
    #-- check if directory exists and recursively create if not
    if (not os.access(os.path.join(DIRECTORY), os.F_OK)):
        os.makedirs(os.path.join(DIRECTORY), MODE)

    #-- create log file with list of synchronized files (or print to terminal)
    if LOG:
        #-- format: UCAR_RDA_NCEP-CFSR_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'UCAR_RDA_NCEP-CFSR_{0}.log'.format(today)
        fid = open(os.path.join(DIRECTORY,LOGFILE),'w')
        print('UCAR NCEP-CFSR Sync Log ({0})'.format(today), file=fid)
    else:
        #-- standard output (terminal output)
        fid = sys.stdout

    #-- read the links list file
    with open(links_list_file,'rb') as fileID:
        lines = fileID.read().decode("utf-8-sig").encode("utf-8").splitlines()

    #-- regular expression pattern for finding netCDF4 files
    prefix = r'ipvh' if ISENTROPIC else r'pgbh'
    suffix = r'\.gz' if GZIP else r''
    regex_pattern = r'{0}(nl)?\.(.*?)\.({1})(\d{{2}})?\.(.*?)\.nc{2}'
    #-- compile regular expression operators
    rx1=re.compile(regex_pattern.format(prefix,r'\d{4}',suffix).encode('utf-8'))
    rx2=re.compile(rb'(http(s)?\:\/\/rda.ucar.edu\/dsrqst\/(.*?))\"', re.VERBOSE)
    #-- output arrays with year for each file
    valid_lines = []
    year = []
    #-- for each line in the links_list_file
    for i,input_lines in enumerate(lines):
        #-- extract filename from url
        if rx1.search(input_lines):
            match_object = rx1.search(input_lines)
            valid_lines.append(match_object.group(0).decode('utf-8'))
            year.append(match_object.group(3).decode('utf-8'))
        if rx2.search(input_lines):
            match_object = rx2.search(input_lines)
            HOST=match_object.group(1).decode('utf-8').replace('http:','https:')

    #-- output variable names
    VARNAME = 'PRES_L1_Avg'
    LONNAME = 'lon'
    LATNAME = 'lat'
    TIMENAME = 'time'
    #-- output a regular grid
    nlat,nlon = (361,720)

    #-- for each unique date
    YEARS = np.unique(year).astype(np.float) if (YEARS is None) else YEARS
    for YY in YEARS:
        #-- compile regular expression operators for finding files in year
        rx = re.compile(regex_pattern.format(prefix,str(YY),suffix), re.VERBOSE)
        remote_lines = [i for i,fi in enumerate(valid_lines) if rx.search(fi)]
        #-- create a list object for each month
        p_month = [[] for m in range(12)]
        #-- for each file in the year
        for i in sorted(remote_lines):
            #-- print input file to be read into memory
            print(valid_lines[i])
            #-- Create and submit request. There are a wide range of exceptions
            # local = os.path.join(DIRECTORY,valid_lines[i])
            response=gravity_toolkit.utilities.from_http([HOST,valid_lines[i]],
                timeout=TIMEOUT, context=None, verbose=True, fid=fid, local=None)
            #-- extract information from monthly files
            NL,OP,YEAR,MONTH,AUX = rx.findall(valid_lines[i]).pop()
            #-- decompress file or convert stream to BytesIO object
            if GZIP:
                fd = gzip.GzipFile(fileobj=io.BytesIO(response.read()))
            else:
                fd = io.BytesIO(response.read()).seek(0)
            #-- open remote file with netCDF4
            dinput = spatial().from_netCDF4(fd, compression='bytes',
                varname=VARNAME, latname=LATNAME, lonname=LONNAME,
                timename=TIMENAME, verbose=False).transpose(axes=(1,2,0))
            #-- read variable for hours since start of file
            epoch,to_secs = gravity_toolkit.time.parse_date_string(
                dinput.attributes['time']['units'])
            #-- calculate the date in Julian Days
            JD = 2400000.5 + gravity_toolkit.time.convert_delta_time(
                dinput.time*to_secs, epoch1=epoch,
                epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)
            #-- convert from Julian days to calendar dates
            Y,M,D,h,m,s = gravity_toolkit.time.convert_julian(JD,
                FORMAT='tuple')
            #-- output delta time in hours since start of year
            dinput.time = gravity_toolkit.time.convert_calendar_dates(Y,
                M, D, hour=h, minute=m, second=s,
                epoch=(YY,1,1,0,0,0), scale=24.0)
            #-- for each month
            for tt,mn in enumerate(M):
                #-- convert month number to variable indice
                m1 = np.int(mn - 1)
                #-- extract data for the month
                p_month[m1].append(dinput.index(tt))

        #-- reduce to valid months
        valid_months = [m for m in range(12) if p_month[m]]
        nmon = len(valid_months)
        #-- python dictionary with output data
        output = {}
        output[VARNAME] = np.ma.zeros((nmon,nlat,nlon))
        output[VARNAME].mask = np.ones((nmon,nlat,nlon),dtype=bool)
        output[TIMENAME] = np.zeros((len(valid_months)))
        #-- copy dimension variables
        output[LATNAME] = dinput.lat.copy()
        output[LONNAME] = dinput.lon.copy()
        #-- for each valid month
        for m in valid_months:
            p_mean = spatial().from_list(p_month[m]).mean()
            output[VARNAME].data[m,:,:] = p_mean.data.copy()
            output[VARNAME].mask[m,:,:] = p_mean.mask.copy()
            output[TIMENAME][m] = p_mean.time.copy()
        #-- output netCDF4 filename
        FILE = '{0}.{1}.{2:4.0f}.nc'.format(prefix,OP,YY)
        #-- output time format: hours since start of file
        TIME_UNITS = 'hours since {0:4.0f}-01-01 00:00:0.0'.format(YY)
        TIME_LONGNAME = 'time'
        #-- output to netCDF4 file
        ncdf_model_write(output, dinput.fill_value, VARNAME=VARNAME,
            LONNAME=LONNAME, LATNAME=LATNAME, TIMENAME=TIMENAME,
            TIME_UNITS=TIME_UNITS, TIME_LONGNAME=TIME_LONGNAME,
            FILENAME=os.path.join(DIRECTORY,FILE), VERBOSE=True)
        #-- set permissions mode to MODE
        os.chmod(os.path.join(DIRECTORY,FILE), MODE)

    #-- close log file and set permissions level to MODE
    if LOG:
        fid.close()
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

#-- PURPOSE: write output model layer fields data to file
def ncdf_model_write(dinput, fill_value, VARNAME=None, LONNAME=None,
    LATNAME=None, TIMENAME=None, TIME_UNITS=None, TIME_LONGNAME=None,
    FILENAME=None, VERBOSE=False):
    #-- opening NetCDF file for writing
    fileID = netCDF4.Dataset(FILENAME, 'w', format="NETCDF4")

    #-- Defining the NetCDF dimensions and creating dimension variables
    nc = {}
    for key in [LONNAME,LATNAME,TIMENAME]:
        fileID.createDimension(key, len(dinput[key]))
        nc[key] = fileID.createVariable(key,dinput[key].dtype,(key,))
    #-- creating the surface NetCDF variables
    nc[VARNAME] = fileID.createVariable(VARNAME, dinput[VARNAME].dtype,
        (TIMENAME,LATNAME,LONNAME,), fill_value=fill_value, zlib=True)

    #-- filling NetCDF variables
    for key,val in dinput.items():
        nc[key][:] = np.copy(val)

    #-- Defining attributes for longitude and latitude
    nc[LONNAME].long_name = 'Longitude'
    nc[LONNAME].units = 'degrees_east'
    nc[LATNAME].long_name = 'Latitude'
    nc[LATNAME].units = 'degrees_north'
    #-- Defining attributes for time
    nc[TIMENAME].units = TIME_UNITS
    nc[TIMENAME].long_name = TIME_LONGNAME
    #-- Defining attributes for surface pressure
    nc[VARNAME].long_name = 'Surface_Air_Pressure'
    nc[VARNAME].standard_name = 'air_pressure'
    nc[VARNAME].product_description = 'Monthly Mean (4 per day) of Analyses'
    nc[VARNAME].level = 'Ground or water surface'
    nc[VARNAME].units = 'Pa'
    #-- date created
    fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())

    #-- Output NetCDF structure information
    if VERBOSE:
        print(os.path.basename(FILENAME))
        print(list(fileID.variables.keys()))

    #-- Closing the NetCDF file
    fileID.close()

#-- Main program that calls ucar_rda_download()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Downloads NCEP-CFSR products using a links list
            provided by the NCAR/UCAR Research Data Archive (RDA)
            """
    )
    #-- command line parameters
    parser.add_argument('file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='UCAR links list file')
    #-- UCAR/NCAR RDA credentials
    parser.add_argument('--user','-U',
        type=str, default=os.environ.get('UCAR_RDA_USERNAME'),
        help='Username for UCAR/NCAR RDA Login')
    parser.add_argument('--password','-P',
        type=str, default=os.environ.get('UCAR_RDA_PASSWORD'),
        help='Password for UCAR/NCAR RDA Login')
    parser.add_argument('--netrc','-N',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.path.join(os.path.expanduser('~'),'.netrc'),
        help='Path to .netrc file for authentication')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- model years to download
    parser.add_argument('--year','-Y',
        type=int, nargs='+',
        help='Years to download from input links file')
    #-- Download isentropic level data
    parser.add_argument('--isentropic','-I',
        default=False, action='store_true',
        help='Input data is over isentropic levels')
    #-- Download compressed data
    parser.add_argument('--gzip','-G',
        default=False, action='store_true',
        help='Input data is compressed')
    #-- connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    #-- Output log file in form
    #-- UCAR_RDA_NCEP-CFSR_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    #-- permissions mode of the directories and files synced (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files synced')
    args = parser.parse_args()

    #-- UCAR/NCAR RDA hostname
    HOST = 'rda.ucar.edu'
    #-- get UCAR/NCAR RDA credentials
    if not args.user and not os.access(args.netrc,os.F_OK):
        #-- check that UCAR/NCAR RDA credentials were entered
        args.user=builtins.input('Username for {0}: '.format(HOST))
        #-- enter password securely from command-line
        args.password=getpass.getpass('Password for {0}@{1}: '.format(args.user,HOST))
    elif not args.user and os.access(args.netrc,os.F_OK):
        args.user,_,args.password=netrc.netrc(args.netrc).authenticators(HOST)
    elif not args.password:
        #-- enter password securely from command-line
        args.password=getpass.getpass('Password for {0}@{1}: '.format(args.user,HOST))

    #-- Build opener with cookie jar for storing cookies
    #-- This is used to store and return the session cookie given
    #-- to use by the data server
    gravity_toolkit.utilities.build_opener(args.user, args.password,
        password_manager=False, get_ca_certs=False, redirect=False,
        authorization_header=False, urs=HOST)
    #-- post authorization header to retrieve cookies
    data = gravity_toolkit.utilities.urlencode(
        dict(email=args.user, passwd=args.password, action='login')).encode("utf-8")
    request = gravity_toolkit.utilities.urllib2.Request(
        'https://rda.ucar.edu/cgi-bin/login',data=data)
    response = gravity_toolkit.utilities.urllib2.urlopen(request,timeout=20)
    #-- All calls to urllib2.urlopen will now use handler

    #-- check that connection to UCAR RDA was successful
    if (response.read().startswith(b'Authentication successful')):
        #-- for each links list file from UCAR/NCAR RDA
        for FILE in args.file:
            ucar_rda_download(FILE, DIRECTORY=args.directory, YEARS=args.year,
                ISENTROPIC=args.isentropic, GZIP=args.gzip, TIMEOUT=args.timeout,
                LOG=args.log, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
