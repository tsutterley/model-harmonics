#!/usr/bin/env python
u"""
racmo_downscaled_cumulative.py
Written by Tyler Sutterley (12/2022)
Calculates cumulative anomalies of RACMO surface mass balance products

COMMAND LINE OPTIONS:
    --help: list the command line options
    --directory=X: set the base data directory
    --version=X: Downscaled RACMO Version
        1.0: RACMO2.3/XGRN11
        2.0: RACMO2.3p2/XGRN11
        3.0: RACMO2.3p2/FGRN055
        4.0: RACMO2.3p2/FGRN055
    --product: RACMO product to calculate
        SMB: Surface Mass Balance
        PRECIP: Precipitation
        RUNOFF: Melt Water Runoff
        SNOWMELT: Snowmelt
        REFREEZE: Melt Water Refreeze
    --mean: Start and end year of mean (separated by commas)
    -G, --gzip: Input netCDF data files are compressed
    -M X, --mode=X: Permission mode of directories and files created
    -V, --verbose: Verbose output of netCDF4 variables

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 10/2022: added version 4.0 (RACMO2.3p2 for 1958-2022 from FGRN055)
    Updated 08/2022: updated docstrings to numpy documentation format
    Updated 02/2021: using argparse to set parameters
    Forked 09/2019 from downscaled_cumulative_netcdf.py
    Updated 07/2019: added version 3.0 (RACMO2.3p2 for 1958-2018 from FGRN055)
    Updated 06/2018: using python3 compatible octal and input
    Updated 11/2017: added version 2.0 (RACMO2.3p2 for 1958-2016)
    Updated 02/2017: using getopt to set base directory
    Written 10/2016
"""
from __future__ import print_function

import sys
import os
import re
import uuid
import gzip
import netCDF4
import argparse
import numpy as np
from datetime import date
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc

# data product longnames
longname = {}
longname['SMB'] = 'Cumulative Surface Mass Balance Anomalies'
longname['PRECIP'] = 'Cumulative Precipitation Anomalies'
longname['RUNOFF'] = 'Cumulative Runoff Anomalies'
longname['SNOWMELT'] = 'Cumulative Snowmelt Anomalies'
longname['REFREEZE'] = 'Cumulative Melt Water Refreeze Anomalies'
# netcdf variable names
input_products = {}
input_products['SMB'] = 'SMB_rec'
input_products['PRECIP'] = 'precip'
input_products['RUNOFF'] = 'runoff'
input_products['SNOWMELT'] = 'snowmelt'
input_products['REFREEZE'] = 'refreeze'

# PURPOSE: get the dimensions for the input data matrices
def get_dimensions(input_dir, VERSION, PRODUCT, GZIP=False):
    """Get the total dimensions of the input data

    Parameters
    ----------
    input_dir: str
        Working data directory
    VERSION: str
        Downscaled RACMO Version
    VARIABLE: str
        RACMO product to run

            - ``SMB``: Surface Mass Balance
            - ``PRECIP``: Precipitation
            - ``RUNOFF``: Melt Water Runoff
            - ``SNOWMELT``: Snowmelt
            - ``REFREEZE``: Melt Water Refreeze
    GZIP: bool, default False
        netCDF data files are compressed
    """
    # names within netCDF4 files
    VARIABLE = input_products[PRODUCT]
    # variable of interest
    if PRODUCT in ('SMB','PRECIP') and (VERSION == '2.0'):
        VARNAME = VARIABLE
    else:
        VARNAME = f'{VARIABLE}corr'
    # if reading yearly files or compressed files
    if VERSION in ('1.0','4.0'):
        # find input files
        pattern = rf'{VARIABLE}.(\d+).BN_(.*?).MM.nc(\.gz)?'
        rx = re.compile(pattern, re.VERBOSE | re.IGNORECASE)
        infiles = sorted([f for f in os.listdir(input_dir) if rx.match(f)])
        nt = 12*len(infiles)
        # read netCDF file for dataset (could also set memory=None)
        if GZIP:
            # read bytes from compressed file
            fd = gzip.open(os.path.join(input_dir,infiles[0]),'rb')
            # read netCDF file for dataset from bytes
            fileID = netCDF4.Dataset(uuid.uuid4().hex,mode='r',memory=fd.read())
        else:
            fileID = netCDF4.Dataset(os.path.join(input_dir,infiles[0]), mode='r')
        # shape of the input data matrix
        nm,ny,nx = fileID.variables[VARIABLE].shape
        fileID.close()
    elif VERSION in ('2.0','3.0'):
        # if reading bytes from compressed file or netcdf file directly
        gz = '.gz' if GZIP else ''
        # input dataset for variable
        file_format = {}
        file_format['2.0'] = '{0}.1958-2016.BN_RACMO2.3p2_FGRN11_GrIS.MM.nc{1}'
        file_format['3.0'] = '{0}.1958-2016.BN_RACMO2.3p2_FGRN055_GrIS.MM.nc{1}'
        f = file_format[VERSION].format(VARIABLE.lower(),gz)
        if GZIP:
            # read bytes from compressed file
            fd = gzip.open(os.path.join(input_dir,f),'rb')
            # read netCDF file for dataset from bytes
            fileID = netCDF4.Dataset(uuid.uuid4().hex,mode='r',memory=fd.read())
        else:
            # read netCDF file for dataset (could also set memory=None)
            fileID = netCDF4.Dataset(os.path.join(input_dir,f), mode='r')
        # shape of the input data matrix
        nt,ny,nx = fileID.variables[VARNAME].shape
        fd.close() if GZIP else fileID.close()
    # return the data dimensions
    return (nt,ny,nx)

# PURPOSE: read individual yearly netcdf files and calculate anomalies
def yearly_file_cumulative(input_dir, VERSION, PRODUCT, MEAN, GZIP=False):
    """Read individual yearly netcdf files and calculate
    cumulative anomalies

    Parameters
    ----------
    input_dir: str
        Working data directory
    VERSION: str
        Downscaled RACMO Version
    PRODUCT: str
        RACMO product to run

            - ``SMB``: Surface Mass Balance
            - ``PRECIP``: Precipitation
            - ``RUNOFF``: Melt Water Runoff
            - ``SNOWMELT``: Snowmelt
            - ``REFREEZE``: Melt Water Refreeze
    MEAN: list
        Start and end year for mean
    GZIP: bool, default False
        netCDF data files are compressed
    """
    # names within netCDF4 files
    VARIABLE = input_products[PRODUCT]
    # regular expression operator for finding variables
    regex = re.compile(VARIABLE, re.VERBOSE | re.IGNORECASE)
    # find input files for years of interest
    pattern = rf'{VARIABLE}.(\d+).BN_(.*?).MM.nc(\.gz)?'
    rx = re.compile(pattern, re.VERBOSE | re.IGNORECASE)
    input_files = sorted([fi for fi in os.listdir(input_dir) if rx.match(fi)])
    # number of input files
    n_files = len(input_files)
    # input dimensions and counter variable
    # get dimensions for input VERSION
    nt,ny,nx = get_dimensions(input_dir,VERSION,PRODUCT,GZIP=GZIP)
    # create counter variable
    c = 0
    # allocate for all data
    dinput = {}
    dinput['LON'] = np.zeros((ny,nx))
    dinput['LAT'] = np.zeros((ny,nx))
    dinput['x'] = np.zeros((nx))
    dinput['y'] = np.zeros((ny))
    dinput['TIME'] = np.zeros((nt))
    dinput['MASK'] = np.zeros((ny,nx),dtype=np.int8)
    dinput[VARIABLE] = np.zeros((nt,ny,nx))
    CUMULATIVE = np.zeros((ny,nx))

    # if reading bytes from compressed file or netcdf file directly
    gz = '.gz' if GZIP else ''
    # input area file with ice mask and model topography
    if (VERSION == '4.0'):
        f1 = f'Icemask_Topo_Iceclasses_lon_lat_average_1km_GrIS.nc{gz}'
        if GZIP:
            # read bytes from compressed file
            fd = gzip.open(os.path.join(input_dir,f1),'rb')
            # read netCDF file for topography and ice classes from bytes
            fileID = netCDF4.Dataset(uuid.uuid4().hex, mode='r', memory=fd.read())
        else:
            # read netCDF file for topography and ice classes
            fileID = netCDF4.Dataset(os.path.join(input_dir,f1), mode='r')
        # Getting the data from each netCDF variable
        dinput['LON'] = np.array(fileID.variables['LON'][:,:])
        dinput['LAT'] = np.array(fileID.variables['LAT'][:,:])
        dinput['x'] = np.array(fileID.variables['x'][:])
        dinput['y'] = np.array(fileID.variables['y'][:])
        promicemask = np.array(fileID.variables['Promicemask'][:,:])
        topography = np.array(fileID.variables['Topography'][:,:])
        # close the compressed file objects
        fd.close() if GZIP else fileID.close()
        # find ice sheet points from promicemask that valid
        ii,jj = np.nonzero((promicemask >= 1) & (promicemask <= 3))
        dinput['MASK'] = np.zeros((ny,nx),dtype=np.int8)
        dinput['MASK'][ii,jj] = 1

    # for each file of interest
    for t in range(n_files):
        # Open the NetCDF file for reading
        if GZIP:
            # read bytes from compressed file
            fd = gzip.open(os.path.join(input_dir,input_files[t]),'rb')
            # read netCDF file for dataset from bytes
            fileID = netCDF4.Dataset(uuid.uuid4().hex, mode='r', memory=fd.read())
        else:
            # read netCDF file for dataset (could also set memory=None)
            fileID = netCDF4.Dataset(os.path.join(input_dir,input_files[t]), mode='r')
        # Getting the data from each netCDF variable
        if (VERSION == '1.0'):
            dinput['LON'][:,:] = fileID.variables['LON'][:,:].copy()
            dinput['LAT'][:,:] = fileID.variables['LAT'][:,:].copy()
            dinput['x'][:] = fileID.variables['x'][:].copy()
            dinput['y'][:] = fileID.variables['y'][:].copy()
            dinput['MASK'][:,:] = fileID.variables['icemask'][:,:].astype(np.int8)
        # calculate dates from delta times
        delta_time = fileID.variables['time'][:].copy()
        date_string = fileID.variables['time'].units
        epoch,to_secs = gravtk.time.parse_date_string(date_string)
        # calculate time array in Julian days
        JD = gravtk.time.convert_delta_time(delta_time*to_secs, epoch1=epoch,
            epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0) + 2400000.5
        # for each month
        for m in range(12):
            # convert from Julian days to calendar dates
            YY,MM,DD,hh,mm,ss = gravtk.time.convert_julian(JD[m],
                format='tuple')
            # calculate time in year-decimal
            dinput['TIME'][c] = gravtk.time.convert_calendar_decimal(YY, MM,
                day=DD, hour=hh, minute=mm, second=ss)
            # find variable of interest
            ncvar, = [v for v in fileID.variables.keys() if regex.match(v)]
            # extract data and add to total cumulative matrix
            CUMULATIVE += (fileID.variables[ncvar][m,:,:].copy() - MEAN)
            dinput[VARIABLE][c,:,:] = CUMULATIVE.copy()
            # add to counter
            c += 1
        # close the NetCDF file
        fileID.close()

    # return the cumulative anomalies
    return dinput

# PURPOSE: read compressed netCDF4 files and calculate cumulative anomalies
def compressed_file_cumulative(input_dir, VERSION, PRODUCT, MEAN, GZIP=False):
    """Read compressed netCDF4 files and calculate cumulative anomalies

    Parameters
    ----------
    input_dir: str
        Working data directory
    VERSION: str
        Downscaled RACMO Version
    PRODUCT: str
        RACMO product to run

            - ``SMB``: Surface Mass Balance
            - ``PRECIP``: Precipitation
            - ``RUNOFF``: Melt Water Runoff
            - ``SNOWMELT``: Snowmelt
            - ``REFREEZE``: Melt Water Refreeze
    MEAN: list
        Start and end year for mean
    GZIP: bool, default False
        netCDF data files are compressed
    """
    # names within netCDF4 files
    VARIABLE = input_products[PRODUCT]
    # variable of interest
    if (PRODUCT == 'SMB') or ((PRODUCT == 'PRECIP') and (VERSION == '2.0')):
        VARNAME = VARIABLE
    else:
        VARNAME = f'{VARIABLE}corr'

    # if reading bytes from compressed file or netcdf file directly
    gz = '.gz' if GZIP else ''
    # allocate for all data
    dinput = {}

    # input area file with ice mask and model topography
    f1 = f'Icemask_Topo_Iceclasses_lon_lat_average_1km_GrIS.nc{gz}'
    if GZIP:
        # read bytes from compressed file
        fd = gzip.open(os.path.join(input_dir,f1),'rb')
        # read netCDF file for topography and ice classes from bytes
        fileID = netCDF4.Dataset(uuid.uuid4().hex, mode='r', memory=fd.read())
    else:
        # read netCDF file for topography and ice classes
        fileID = netCDF4.Dataset(os.path.join(input_dir,f1), mode='r')
    # Getting the data from each netCDF variable
    dinput['LON'] = np.array(fileID.variables['LON'][:,:])
    dinput['LAT'] = np.array(fileID.variables['LAT'][:,:])
    dinput['x'] = np.array(fileID.variables['x'][:])
    dinput['y'] = np.array(fileID.variables['y'][:])
    promicemask = np.array(fileID.variables['Promicemask'][:,:])
    topography = np.array(fileID.variables['Topography'][:,:])
    # close the compressed file objects
    fd.close() if GZIP else fileID.close()

    # file format for each version
    file_format = {}
    file_format['2.0'] = '{0}.1958-2016.BN_RACMO2.3p2_FGRN11_GrIS.MM.nc{1}'
    file_format['3.0'] = '{0}.1958-2018.BN_RACMO2.3p2_FGRN055_GrIS.MM.nc{1}'

    # input dataset for variable
    f2 = file_format[VERSION].format(VARIABLE.lower(),gz)
    if GZIP:
        # read bytes from compressed file
        fd = gzip.open(os.path.join(input_dir,f2),'rb')
        # read netCDF file for dataset from bytes
        fileID = netCDF4.Dataset(uuid.uuid4().hex, mode='r', memory=fd.read())
    else:
        # read netCDF file for dataset (could also set memory=None)
        fileID = netCDF4.Dataset(os.path.join(input_dir,f2), mode='r')
    # shape of the input data matrix
    nt,ny,nx = fileID.variables[VARNAME].shape

    # find ice sheet points from promicemask that valid
    ii,jj = np.nonzero((promicemask >= 1) & (promicemask <= 3))
    dinput['MASK'] = np.zeros((ny,nx),dtype=np.int8)
    dinput['MASK'][ii,jj] = 1

    # calculate dates from delta times
    # Months since 1958-01-15 at 00:00:00
    delta_time = fileID.variables['time'][:].copy()
    date_string = fileID.variables['time'].units
    epoch,to_secs = gravtk.time.parse_date_string(date_string)
    # calculate time array in Julian days
    JD = gravtk.time.convert_delta_time(delta_time*to_secs, epoch1=epoch,
        epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0) + 2400000.5
    # convert from Julian days to calendar dates
    YY,MM,DD,hh,mm,ss = gravtk.time.convert_julian(JD, format='tuple')
    # calculate time in year-decimal
    dinput['TIME'] = gravtk.time.convert_calendar_decimal(YY, MM,
        day=DD, hour=hh, minute=mm, second=ss)

    # calculate cumulative
    CUMULATIVE = np.zeros((ny,nx))
    dinput[VARNAME] = np.zeros((nt,ny,nx))
    for t in range(nt):
        # extract data and add to total cumulative matrix
        CUMULATIVE += (fileID.variables[VARNAME][t,:,:].copy() - MEAN)
        dinput[VARNAME][t,:,:] = CUMULATIVE.copy()

    # close the compressed file objects
    fd.close() if GZIP else fileID.close()

    # return the cumulative anomalies
    return dinput

# PURPOSE: write RACMO downscaled data to netCDF4
def ncdf_racmo(dinput, FILENAME=None, UNITS=None, LONGNAME=None, VARNAME=None,
    LONNAME=None, LATNAME=None, XNAME=None, YNAME=None, TIMENAME=None,
    MASKNAME=None, TIME_UNITS='years', TIME_LONGNAME='Date_in_Decimal_Years',
    TITLE = None, CLOBBER = False, VERBOSE=False):

    # setting NetCDF clobber attribute
    if CLOBBER:
        clobber = 'w'
    else:
        clobber = 'a'

    # opening NetCDF file for writing
    # Create the NetCDF file
    fileID = netCDF4.Dataset(FILENAME, clobber, format="NETCDF4")

    # Dimensions of parameters
    n_time = 1 if (np.ndim(dinput[TIMENAME]) == 0) else len(dinput[TIMENAME])
    # Defining the NetCDF dimensions
    fileID.createDimension(XNAME, len(dinput[XNAME]))
    fileID.createDimension(YNAME, len(dinput[YNAME]))
    fileID.createDimension(TIMENAME, n_time)

    # python dictionary with netCDF4 variables
    nc = {}

    # defining the NetCDF variables
    nc[XNAME] = fileID.createVariable(XNAME, dinput[XNAME].dtype, (XNAME,))
    nc[YNAME] = fileID.createVariable(YNAME, dinput[YNAME].dtype, (YNAME,))
    nc[TIMENAME] = fileID.createVariable(TIMENAME, dinput[TIMENAME].dtype,
        (TIMENAME,))
    nc[LONNAME] = fileID.createVariable(LONNAME, dinput[LONNAME].dtype,
        (YNAME,XNAME,))
    nc[LATNAME] = fileID.createVariable(LATNAME, dinput[LATNAME].dtype,
        (YNAME,XNAME,))
    nc[MASKNAME] = fileID.createVariable(MASKNAME, dinput[MASKNAME].dtype,
        (YNAME,XNAME,), fill_value=0, zlib=True)
    if (n_time > 1):
        nc[VARNAME] = fileID.createVariable(VARNAME, dinput[VARNAME].dtype,
            (TIMENAME,YNAME,XNAME,), zlib=True)
    else:
        nc[VARNAME] = fileID.createVariable(VARNAME, dinput[VARNAME].dtype,
            (YNAME,XNAME,), zlib=True)

    # filling NetCDF variables
    for key,val in dinput.items():
        nc[key][:] = val.copy()

    # Defining attributes for longitude and latitude
    nc[LONNAME].long_name = 'longitude'
    nc[LONNAME].units = 'degrees_east'
    nc[LATNAME].long_name = 'latitude'
    nc[LATNAME].units = 'degrees_north'
    # Defining attributes for x and y coordinates
    nc[XNAME].long_name = 'easting'
    nc[XNAME].units = 'meters'
    nc[YNAME].long_name = 'northing'
    nc[YNAME].units = 'meters'
    # Defining attributes for dataset
    nc[VARNAME].long_name = LONGNAME
    nc[VARNAME].units = UNITS
    nc[MASKNAME].long_name = 'mask'
    # Defining attributes for date
    nc[TIMENAME].long_name = TIME_LONGNAME
    nc[TIMENAME].units = TIME_UNITS

    # add software information
    fileID.software_reference = mdlhmc.version.project_name
    fileID.software_version = mdlhmc.version.full_version
    fileID.reference = f'Output from {os.path.basename(sys.argv[0])}'
    # global variable of NetCDF file
    fileID.TITLE = TITLE
    fileID.date_created = date.isoformat(date.today())

    # Output NetCDF structure information
    if VERBOSE:
        print(FILENAME)
        print(list(fileID.variables.keys()))

    # Closing the NetCDF file
    fileID.close()

# PURPOSE: calculate RACMO cumulative anomalies with respect to a mean field
def racmo_downscaled_cumulative(base_dir, VERSION, PRODUCT,
    RANGE=[1961,1990], GZIP=False, VERBOSE=False, MODE=0o775):
    """
    Calculate RACMO cumulative anomalies with respect to a mean field

    Parameters
    ----------
    base_dir: str
        Working data directory
    VERSION: str
        Downscaled RACMO Version

            - ``1.0``: RACMO2.3/XGRN11
            - ``2.0``: RACMO2.3p2/XGRN11
            - ``3.0``: RACMO2.3p2/FGRN055
    PRODUCT: str
        RACMO product to calculate

            - ``SMB``: Surface Mass Balance
            - ``PRECIP``: Precipitation
            - ``RUNOFF``: Melt Water Runoff
            - ``SNOWMELT``: Snowmelt
            - ``REFREEZE``: Melt Water Refreeze
    RANGE: list, default [1961,1990]
        Start and end year of mean
    GZIP: bool, default False
        netCDF data files are compressed
    VERBOSE: bool, default False
        Verbose output of netCDF4 variables
    MODE: oct, default 0o775
        Permission mode of directories and files created
    """

    # Full Directory Setup
    DIRECTORY = f'SMB1km_v{VERSION}'

    # versions 1 and 4 are in separate files for each year
    if (VERSION == '1.0'):
        RACMO_MODEL = ['XGRN11','2.3']
        VARNAME = input_products[PRODUCT]
        SUBDIRECTORY = f'{VARNAME}_v{VERSION}'
        input_dir = os.path.join(base_dir, DIRECTORY, SUBDIRECTORY)
    elif (VERSION == '2.0'):
        RACMO_MODEL = ['XGRN11','2.3p2']
        var = input_products[PRODUCT]
        VARNAME = var if PRODUCT in ('SMB','PRECIP') else f'{var}corr'
        input_dir = os.path.join(base_dir, DIRECTORY)
    elif (VERSION == '3.0'):
        RACMO_MODEL = ['FGRN055','2.3p2']
        var = input_products[PRODUCT]
        VARNAME = var if (PRODUCT == 'SMB') else f'{var}corr'
        input_dir = os.path.join(base_dir, DIRECTORY)
    elif (VERSION == '4.0'):
        RACMO_MODEL = ['FGRN055','2.3p2']
        var = input_products[PRODUCT]
        VARNAME = var if (PRODUCT == 'SMB') else f'{var}corr'
        input_dir = os.path.join(base_dir, DIRECTORY)

    # read mean from netCDF4 file
    arg = (RACMO_MODEL[0],RACMO_MODEL[1],VERSION,PRODUCT,RANGE[0],RANGE[1])
    mean_file = '{0}_RACMO{1}_DS1km_v{2}_{3}_Mean_{4:4d}-{5:4d}.nc'.format(*arg)
    with netCDF4.Dataset(os.path.join(input_dir,mean_file),'r') as fileID:
        MEAN = fileID[VARNAME][:,:].copy()

    # calculate cumulative
    if VERSION in ('1.0','4.0'):
        dinput = yearly_file_cumulative(input_dir, VERSION, PRODUCT, MEAN,
            GZIP=GZIP)
    elif VERSION in ('2.0','3.0'):
        dinput = compressed_file_cumulative(input_dir, VERSION, PRODUCT,
            MEAN, GZIP=GZIP)

    # output cumulative as netCDF4 file
    args = (RACMO_MODEL[0],RACMO_MODEL[1],VERSION,PRODUCT)
    output_file = '{0}_RACMO{1}_DS1km_v{2}_{3}_cumul.nc'.format(*args)
    TITLE = 'Downscaled_cumulative_anomalies_relative_to_{0:4d}-{1:4d}_Mean'
    ncdf_racmo(dinput, FILENAME=os.path.join(input_dir,output_file),
        UNITS='mmWE', LONGNAME=longname[PRODUCT], VARNAME=VARNAME,
        LONNAME='LON', LATNAME='LAT', XNAME='x', YNAME='y', TIMENAME='TIME',
        MASKNAME='MASK', TITLE=TITLE.format(RANGE[0],RANGE[1]),
        CLOBBER=True, VERBOSE=VERBOSE)
    # change the permissions mode
    os.chmod(os.path.join(input_dir,output_file),MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates cumulative anomalies of
            RACMO downscaled surface mass balance products
            with respect to a mean field
            """
    )
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # Downscaled version
    # 1.0: RACMO2.3/XGRN11
    # 2.0: RACMO2.3p2/XGRN11
    # 3.0: RACMO2.3p2/FGRN055
    # 4.0: RACMO2.3p2/FGRN055
    parser.add_argument('--version','-v',
        type=str, default='4.0', choices=['1.0','2.0','3.0','4.0'],
        help='Downscaled RACMO Version')
    # Products to calculate cumulative
    parser.add_argument('--product','-p',
        metavar='PRODUCT', type=str, nargs='+',
        default=['SMB'], choices=input_products.keys(),
        help='RACMO product to calculate')
    # start and end years to run for mean
    parser.add_argument('--mean','-m',
        metavar=('START','END'), type=int, nargs=2,
        default=[1961,1990],
        help='Start and end year range for mean')
    # netCDF4 files are gzip compressed
    parser.add_argument('--gzip','-G',
        default=False, action='store_true',
        help='netCDF4 file is locally gzip compressed')
    # verbose output of processing run
    # print information about each input and output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files')
    # return the parser
    return parser

# Main program that calls racmo_downscaled_cumulative()
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args = parser.parse_args()

    # run program for each input product
    for PRODUCT in args.product:
        # run downscaled cumulative program with parameters
        racmo_downscaled_cumulative(args.directory, args.version, PRODUCT,
            RANGE=args.mean, GZIP=args.gzip, VERBOSE=args.verbose,
            MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()

