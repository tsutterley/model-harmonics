#!/usr/bin/env python
u"""
racmo_downscaled_mean.py
Written by Tyler Sutterley (06/2024)
Calculates the temporal mean of downscaled RACMO
surface mass balance products

COMMAND LINE OPTIONS:
    --help: list the command line options
    --directory=X: set the base data directory
    -R X, --region X: Region to calculate (gris, ais)
    --version=X: Downscaled RACMO Version
        gris:
            1.0: RACMO2.3/XGRN11/DS1km
            2.0: RACMO2.3p2/XGRN11/DS1km
            3.0: RACMO2.3p2/FGRN055/DS1km
            4.0: RACMO2.3p2/FGRN055/DS1km
            5.0: RACMO2.3p2/FGRN055/DS1km
            6.0: RACMO2.3p2/FGRN055/DS1km
        ais:
            6.0: RACMO2.3/XANT27/DS2km
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
    Updated 06/2024: added version 6.0 for Greenland and Antarctica
        RACMO2.3p2/FGRN055/DS1km for 1958-2023
        RACMO2.3p2/XANT27/DS2km for 1979-2023
    Updated 06/2023: added version 5.0 (RACMO2.3p2 for 1958-2023 from FGRN055)
    Updated 10/2022: added version 5.0 (RACMO2.3p2 for 1958-2023 from FGRN055)
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 10/2022: added version 4.0 (RACMO2.3p2 for 1958-2022 from FGRN055)
    Updated 08/2022: updated docstrings to numpy documentation format
    Updated 02/2021: using argparse to set parameters
    Forked 09/2019 from downscaled_mean_netcdf.py
    Updated 07/2019: added version 3.0 (RACMO2.3p2 for 1958-2018 from FGRN055)
    Updated 06/2018: using python3 compatible octal and input
    Updated 11/2017: added version 2.0 (RACMO2.3p2 for 1958-2016)
    Updated 02/2017: using getopt to set base directory
    Written 11/2016
"""
from __future__ import print_function

import sys
import re
import uuid
import gzip
import logging
import netCDF4
import pathlib
import argparse
import numpy as np
from datetime import date
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc

# data product longnames
longname = {}
longname['SMB'] = 'Surface Mass Balance'
longname['PRECIP'] = 'Precipitation'
longname['RUNOFF'] = 'Runoff'
longname['SNOWMELT'] = 'Snowmelt'
longname['REFREEZE'] = 'Melt Water Refreeze'
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
    # verify input directory
    input_dir = pathlib.Path(input_dir).expanduser().absolute()
    # names within netCDF4 files
    VARIABLE = input_products[PRODUCT]
    # regular expression operator for finding variables
    regex = re.compile(VARIABLE, re.VERBOSE | re.IGNORECASE)
    # if reading yearly files or compressed files
    if VERSION in ('1.0','4.0','5.0','6.0'):
        # find input files
        pattern = r'{0}.(\d+).BN_(.*?).MM.nc(\.gz)?'.format(VARIABLE)
        rx = re.compile(pattern, re.VERBOSE | re.IGNORECASE)
        infiles = sorted([f for f in input_dir.iterdir() if rx.match(f.name)])
        nt = 12*len(infiles)
        # read netCDF file for dataset (could also set memory=None)
        if GZIP:
            # read bytes from compressed file
            fd = gzip.open(str(infiles[0]), 'rb')
            # read netCDF file for dataset from bytes
            fileID = netCDF4.Dataset(uuid.uuid4().hex,mode='r',memory=fd.read())
        else:
            fileID = netCDF4.Dataset(infiles[0], mode='r')
        # shape of the input data matrix
        ncvar, = [v for v in fileID.variables.keys() if regex.match(v)]
        nm,ny,nx = fileID.variables[ncvar].shape
        fileID.close()
    elif VERSION in ('2.0','3.0'):
        # if reading bytes from compressed file or netcdf file directly
        gz = '.gz' if GZIP else ''
        # input dataset for variable
        file_format = {}
        file_format['2.0'] = '{0}.1958-2016.BN_RACMO2.3p2_FGRN11_GrIS.MM.nc{1}'
        file_format['3.0'] = '{0}.1958-2016.BN_RACMO2.3p2_FGRN055_GrIS.MM.nc{1}'
        f = input_dir.joinpath(file_format[VERSION].format(VARIABLE.lower(),gz))
        if GZIP:
            # read bytes from compressed file
            fd = gzip.open(str(f), 'rb')
            # read netCDF file for dataset from bytes
            fileID = netCDF4.Dataset(uuid.uuid4().hex,mode='r',memory=fd.read())
        else:
            # read netCDF file for dataset (could also set memory=None)
            fileID = netCDF4.Dataset(f, mode='r')
        # shape of the input data matrix
        ncvar, = [v for v in fileID.variables.keys() if regex.match(v)]
        nt,ny,nx = fileID.variables[ncvar].shape
        fd.close() if GZIP else fileID.close()
    # return the data dimensions
    return (nt,ny,nx)

# PURPOSE: read individual yearly netcdf files and calculate mean over period
def yearly_file_mean(
        input_dir,
        REGION,
        VERSION,
        PRODUCT,
        START,
        END,
        GZIP=False
    ):
    """Read individual yearly netcdf files and calculate mean

    Parameters
    ----------
    input_dir: str
        Working data directory
    REGION: str
        Region of RACMO model to calculate
    VERSION: str
        Downscaled RACMO Version
    PRODUCT: str
        RACMO product to run

            - ``SMB``: Surface Mass Balance
            - ``PRECIP``: Precipitation
            - ``RUNOFF``: Melt Water Runoff
            - ``SNOWMELT``: Snowmelt
            - ``REFREEZE``: Melt Water Refreeze
    START: int
        Starting year to run
    END: int
        Ending year to run
    GZIP: bool, default False
        netCDF data files are compressed
    """
    # verify input directory
    input_dir = pathlib.Path(input_dir).expanduser().absolute()
    # names within netCDF4 files
    VARIABLE = input_products[PRODUCT]
    # regular expression operator for finding variables
    regex = re.compile(VARIABLE, re.VERBOSE | re.IGNORECASE)
    # find input files for years of interest
    regex_years = r'|'.join(rf'{Y:4d}' for Y in range(START,END+1))
    pattern = rf'{VARIABLE}.({regex_years}).BN_(.*?).MM.nc(\.gz)?'
    rx = re.compile(pattern, re.VERBOSE | re.IGNORECASE)
    input_files = sorted([f for f in input_dir.iterdir() if rx.match(f.name)])
    # number of input files
    n_files = len(input_files)
    # input dimensions and counter variable
    # get dimensions for input VERSION
    nt,ny,nx = get_dimensions(input_dir, VERSION, PRODUCT, GZIP=GZIP)
    # create counter variable
    c = 0
    # allocate for all data
    dinput = {}
    dinput['LON'] = np.zeros((ny,nx))
    dinput['LAT'] = np.zeros((ny,nx))
    dinput['x'] = np.zeros((nx))
    dinput['y'] = np.zeros((ny))
    dinput['MASK'] = np.zeros((ny,nx),dtype=np.int8)
    # calculate total
    dinput[VARIABLE] = np.zeros((ny,nx))
    # date in year-decimal form
    tdec = np.zeros((nt))

    # if reading bytes from compressed file or netcdf file directly
    gz = '.gz' if GZIP else ''
    # input area file with ice mask and model topography
    if VERSION in ('4.0','5.0'):
        f1 = f'Icemask_Topo_Iceclasses_lon_lat_average_1km_GrIS.nc{gz}'
        input_mask_file = input_dir.joinpath(f1)
        if GZIP:
            # read bytes from compressed file
            fd = gzip.open(str(input_mask_file), 'rb')
            # read netCDF file for topography and ice classes from bytes
            fileID = netCDF4.Dataset(uuid.uuid4().hex, mode='r', memory=fd.read())
        else:
            # read netCDF file for topography and ice classes
            fileID = netCDF4.Dataset(input_mask_file, mode='r')
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
            fd = gzip.open(str(input_files[t]), 'rb')
            # read netCDF file for dataset from bytes
            fileID = netCDF4.Dataset(uuid.uuid4().hex, mode='r', memory=fd.read())
        else:
            # read netCDF file for dataset (could also set memory=None)
            fileID = netCDF4.Dataset(input_files[t], 'r')
        # check if ERA5 3-hourly data
        ERA5_3h = re.search(r'ERA5_3h', input_files[t].name)
        # Getting the data from each netCDF variable
        if (VERSION == '1.0'):
            dinput['LON'][:,:] = fileID.variables['LON'][:,:].copy()
            dinput['LAT'][:,:] = fileID.variables['LAT'][:,:].copy()
            dinput['x'][:] = fileID.variables['x'][:].copy()
            dinput['y'][:] = fileID.variables['y'][:].copy()
            dinput['MASK'][:,:] = fileID.variables['icemask'][:,:].astype(np.int8)
        elif (VERSION == '6.0'):
            # extract coordinates
            if ERA5_3h:
                dinput['x'][:] = fileID.variables['x'][:].copy()
                dinput['y'][:] = fileID.variables['y'][:].copy()
            else:
                dinput['x'][:] = fileID.variables['lon'][:].copy()
                dinput['y'][:] = fileID.variables['lat'][:].copy()
            EPSG = 3413 if REGION == 'gris' else 3031
            # convert from edges to centers
            dx = np.diff(dinput['x'])[0]
            dy = np.diff(dinput['y'])[0]
            # calculate latitude and longitude of grid cells
            dinput['LON'][:,:], dinput['LAT'][:,:] = mdlhmc.spatial.get_latlon(
                dinput['x'] - dx/2.0, dinput['y'] - dy/2.0,
                srs_epsg=EPSG)
            # find variable of interest
            ncvar, = [v for v in fileID.variables.keys() if regex.match(v)]
            valid = np.any(fileID.variables[ncvar], axis=0)
            dinput['MASK'] = valid.astype(np.int8)
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
            tdec[c] = gravtk.time.convert_calendar_decimal(YY, MM,
                day=DD, hour=hh, minute=mm, second=ss)
            # find variable of interest
            ncvar, = [v for v in fileID.variables.keys() if regex.match(v)]
            # read product of interest and add to total
            dinput[VARIABLE] += fileID.variables[ncvar][m,:,:].copy()
            # add to counter
            c += 1
        # close the NetCDF file
        fileID.close()

    # calculate mean time over period
    dinput['TIME'] = np.mean(tdec)
    # convert from total to mean
    dinput[VARIABLE] /= np.float64(c)

    # return the mean variables
    return dinput

# PURPOSE: read compressed netCDF4 files and calculate mean over period
def compressed_file_mean(
        input_dir,
        REGION,
        VERSION,
        PRODUCT,
        START,
        END,
        GZIP=False
    ):
    """Read compressed netCDF4 files and calculate mean

    Parameters
    ----------
    input_dir: str
        Working data directory
    REGION: str
        Region of RACMO model to calculate
    VERSION: str
        Downscaled RACMO Version
    PRODUCT: str
        RACMO product to run

            - ``SMB``: Surface Mass Balance
            - ``PRECIP``: Precipitation
            - ``RUNOFF``: Melt Water Runoff
            - ``SNOWMELT``: Snowmelt
            - ``REFREEZE``: Melt Water Refreeze
    START: int
        Starting year to run
    END: int
        Ending year to run
    GZIP: bool, default False
        netCDF data files are compressed
    """
    # names within netCDF4 files
    VARIABLE = input_products[PRODUCT]
    # variable of interest
    if (PRODUCT == 'SMB') or ((PRODUCT == 'PRECIP') and (VERSION == '2.0')):
        VARNAME = VARIABLE
    else:
        VARNAME = '{0}corr'.format(VARIABLE)

    # if reading bytes from compressed file or netcdf file directly
    gz = '.gz' if GZIP else ''
    # allocate for all data
    dinput = {}

    # input area file with ice mask and model topography
    f1 = f'Icemask_Topo_Iceclasses_lon_lat_average_1km_GrIS.nc{gz}'
    input_mask_file = input_dir.joinpath(f1)
    if GZIP:
        # read bytes from compressed file
        fd = gzip.open(str(input_mask_file), 'rb')
        # read netCDF file for topography and ice classes from bytes
        fileID = netCDF4.Dataset(uuid.uuid4().hex, mode='r', memory=fd.read())
    else:
        # read netCDF file for topography and ice classes
        fileID = netCDF4.Dataset(input_mask_file, mode='r')
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
    input_file = input_dir.joinpath(f2)
    if GZIP:
        # read bytes from compressed file
        fd = gzip.open(str(input_file), 'rb')
        # read netCDF file for dataset from bytes
        fileID = netCDF4.Dataset(uuid.uuid4().hex, mode='r', memory=fd.read())
    else:
        # read netCDF file for dataset (could also set memory=None)
        fileID = netCDF4.Dataset(input_file, mode='r')
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
    tdec = gravtk.time.convert_calendar_decimal(YY, MM,
        day=DD, hour=hh, minute=mm, second=ss)

    # calculate total
    dinput[VARNAME] = np.zeros((ny,nx))
    # find indices for dates of interest
    indices, = np.nonzero((tdec >= START) & (tdec < END+1))
    c = np.count_nonzero((tdec >= START) & (tdec < END+1))
    for t in indices:
        # read product of interest and add to total
        dinput[VARNAME] += fileID.variables[VARNAME][t,:,:].copy()
    # convert from total to mean
    dinput[VARNAME] /= np.float64(c),
    # calculate mean time over period
    dinput['TIME'] = np.mean(tdec[indices])

    # close the compressed file objects
    fd.close() if GZIP else fileID.close()
    # return the mean variables
    return dinput

# PURPOSE: write RACMO downscaled data to netCDF4
def ncdf_racmo(dinput, FILENAME=None, REGION=None,
    UNITS=None, LONGNAME=None, VARNAME=None,
    LONNAME=None, LATNAME=None, XNAME=None, YNAME=None,
    TIMENAME=None, MASKNAME=None,
    TIME_UNITS='years', TIME_LONGNAME='Date_in_Decimal_Years',
    TITLE = None, CLOBBER = False):

    # setting NetCDF clobber attribute
    if CLOBBER:
        clobber = 'w'
    else:
        clobber = 'a'

    # opening NetCDF file for writing
    FILENAME = pathlib.Path(FILENAME).expanduser().absolute()
    fileID = netCDF4.Dataset(FILENAME, clobber, format="NETCDF4")

    # Dimensions of parameters
    n_time = 1 if (np.ndim(dinput[TIMENAME]) == 0) else len(dinput[TIMENAME])
    # Defining the NetCDF dimensions
    fileID.createDimension(XNAME, len(dinput[XNAME]))
    fileID.createDimension(YNAME, len(dinput[YNAME]))
    fileID.createDimension(TIMENAME, n_time)

    # python dictionary with netCDF4 variables
    nc = {}

    # create variable and attributes for projection
    if REGION in ('gris',):
        crs = fileID.createVariable('Polar_Stereographic',np.byte,())
        crs.standard_name = 'Polar_Stereographic'
        crs.grid_mapping_name = 'polar_stereographic'
        crs.straight_vertical_longitude_from_pole = -45.0
        crs.latitude_of_projection_origin = 90.0
        crs.standard_parallel = 70.0
        crs.scale_factor_at_projection_origin = 1.
        crs.false_easting = 0.0
        crs.false_northing = 0.0
        crs.semi_major_axis = 6378.137
        crs.semi_minor_axis = 6356.752
        crs.inverse_flattening = 298.257223563
        crs.spatial_epsg = '3413'
    elif REGION in ('ais',):
        crs = fileID.createVariable('Polar_Stereographic',np.byte,())
        crs.standard_name = 'Polar_Stereographic'
        crs.grid_mapping_name = 'polar_stereographic'
        crs.straight_vertical_longitude_from_pole = 0.0
        crs.latitude_of_projection_origin = -90.0
        crs.standard_parallel = -71.0
        crs.scale_factor_at_projection_origin = 1.
        crs.false_easting = 0.0
        crs.false_northing = 0.0
        crs.semi_major_axis = 6378.137
        crs.semi_minor_axis = 6356.752
        crs.inverse_flattening = 298.257223563
        crs.spatial_epsg = '3031'

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
    nc[XNAME].grid_mapping = 'Polar_Stereographic'
    nc[YNAME].long_name = 'northing'
    nc[YNAME].units = 'meters'
    nc[YNAME].grid_mapping = 'Polar_Stereographic'
    # Defining attributes for dataset
    nc[VARNAME].long_name = LONGNAME
    nc[VARNAME].units = UNITS
    nc[VARNAME].grid_mapping = 'Polar_Stereographic'
    nc[MASKNAME].long_name = 'mask'
    nc[MASKNAME].grid_mapping = 'Polar_Stereographic'
    # Defining attributes for date
    nc[TIMENAME].long_name = TIME_LONGNAME
    nc[TIMENAME].units = TIME_UNITS
    # add software information
    fileID.software_reference = mdlhmc.version.project_name
    fileID.software_version = mdlhmc.version.full_version
    fileID.reference = f'Output from {pathlib.Path(sys.argv[0]).name}'
    # global variable of NetCDF file
    fileID.TITLE = TITLE
    fileID.date_created = date.isoformat(date.today())

    # Output NetCDF structure information
    logging.info(str(FILENAME))
    logging.info(list(fileID.variables.keys()))

    # Closing the NetCDF file
    fileID.close()

# PURPOSE: calculate RACMO mean data over a polar stereographic grid
def racmo_downscaled_mean(base_dir, REGION, VERSION, PRODUCT,
    RANGE=[1961,1990], GZIP=False, MODE=0o775):
    """
    Calculates the temporal mean of downscaled RACMO
    surface mass balance products

    Parameters
    ----------
    base_dir: str
        Working data directory
    REGION: str
        RACMO model region

            - ``gris``: Greenland Ice Sheet
            - ``ais``: Antarctic Ice Sheet
    VERSION: str
        Downscaled RACMO Version

            - ``'gris'``:
                - ``1.0``: RACMO2.3/XGRN11/DS1km
                - ``2.0``: RACMO2.3p2/XGRN11/DS1km
                - ``3.0``: RACMO2.3p2/FGRN055/DS1km
                - ``4.0``: RACMO2.3p2/FGRN055/DS1km
                - ``5.0``: RACMO2.3p2/FGRN055/DS1km
                - ``6.0``: RACMO2.3p2/FGRN055/DS1km
            - ``'ais'``:
                - ``6.0``: RACMO2.3/XANT27/DS2km
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
    MODE: oct, default 0o775
        Permission mode of directories and files created
    """

    # Full Directory Setup
    base_dir = pathlib.Path(base_dir).expanduser().absolute()

    # versions 1 and 4 are in separate files for each year
    if (VERSION == '1.0') and (REGION == 'gris'):
        RACMO_MODEL = ['XGRN11','2.3','DS1km']
        VARNAME = input_products[PRODUCT]
        SUBDIRECTORY = f'{VARNAME}_v{VERSION}'
        input_dir = base_dir.joinpath(f'SMB1km_v{VERSION}', SUBDIRECTORY)
        file_type = 'yearly'
    elif (VERSION == '2.0') and (REGION == 'gris'):
        RACMO_MODEL = ['XGRN11','2.3p2','DS1km']
        var = input_products[PRODUCT]
        VARNAME = var if PRODUCT in ('SMB','PRECIP') else f'{var}corr'
        input_dir = base_dir.joinpath(f'SMB1km_v{VERSION}')
        file_type = 'compressed'
    elif (VERSION == '3.0') and (REGION == 'gris'):
        RACMO_MODEL = ['FGRN055','2.3p2','DS1km']
        var = input_products[PRODUCT]
        VARNAME = var if (PRODUCT == 'SMB') else f'{var}corr'
        input_dir = base_dir.joinpath(f'SMB1km_v{VERSION}')
        file_type = 'compressed'
    elif (VERSION == '4.0') and (REGION == 'gris'):
        RACMO_MODEL = ['FGRN055','2.3p2','DS1km']
        var = input_products[PRODUCT]
        VARNAME = var if (PRODUCT == 'SMB') else f'{var}corr'
        SUBDIRECTORY = PRODUCT.lower()
        input_dir = base_dir.joinpath(f'SMB1km_v{VERSION}', SUBDIRECTORY)
        file_type = 'yearly'
    elif (VERSION == '5.0') and (REGION == 'gris'):
        RACMO_MODEL = ['FGRN055','2.3p2','DS1km']
        var = input_products[PRODUCT]
        VARNAME = var if (PRODUCT == 'SMB') else f'{var}corr'
        SUBDIRECTORY = PRODUCT.lower()
        input_dir = base_dir.joinpath(f'SMB1km_v{VERSION}', SUBDIRECTORY)
        file_type = 'yearly'
    elif (VERSION == '6.0') and (REGION == 'gris'):
        RACMO_MODEL = ['FGRN055','2.3p2','DS1km']
        var = input_products[PRODUCT]
        VARNAME = var if (PRODUCT == 'SMB') else f'{var}corr'
        input_dir = base_dir.joinpath('GrIS-1km')
        file_type = 'yearly'
    elif (VERSION == '6.0') and (REGION == 'ais'):
        RACMO_MODEL = ['XANT27','2.3p2','DS2km']
        var = input_products[PRODUCT]
        VARNAME = var if (PRODUCT == 'SMB') else f'{var}corr'
        input_dir = base_dir.joinpath('AIS-2km')
        file_type = 'yearly'

    # read data and calculate mean
    if (file_type == 'yearly'):
        dinput = yearly_file_mean(input_dir, REGION, VERSION, PRODUCT,
            RANGE[0], RANGE[1], GZIP=GZIP)
    elif (file_type == 'compressed'):
        dinput = compressed_file_mean(input_dir, REGION, VERSION, PRODUCT,
            RANGE[0], RANGE[1], GZIP=GZIP)

    # output mean as netCDF4 file
    arg = (*RACMO_MODEL, VERSION, PRODUCT, *RANGE)
    FILE = '{0}_RACMO{1}_{2}_v{3}_{4}_Mean_{5:4d}-{6:4d}.nc'.format(*arg)
    output_file = input_dir.joinpath(FILE)
    ncdf_racmo(dinput, FILENAME=output_file, REGION=REGION,
        UNITS='mmWE', LONGNAME=longname[PRODUCT], VARNAME=VARNAME,
        LONNAME='LON', LATNAME='LAT', XNAME='x', YNAME='y', TIMENAME='TIME',
        MASKNAME='MASK', TITLE='Mean_downscaled_field', CLOBBER=True)
    # change the permission mode
    output_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates the temporal mean of downscaled
            RACMO surface mass balance products
            """
    )
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # region of RACMO model
    parser.add_argument('--region','-R',
        type=str, default='gris', choices=['gris','ais'],
        help='Region of RACMO model to calculate')
    # Downscaled version
    # 1.0: RACMO2.3/XGRN11/DS1km
    # 2.0: RACMO2.3p2/XGRN11/DS1km
    # 3.0: RACMO2.3p2/FGRN055/DS1km
    # 4.0: RACMO2.3p2/FGRN055/DS1km
    # 5.0: RACMO2.3p2/FGRN055/DS1km
    # 6.0: RACMO2.3p2/FGRN055/DS1km
    # 6.0: RACMO2.3p2/XANT27/DS2km
    choices = ['1.0','2.0','3.0','4.0','5.0','6.0']
    parser.add_argument('--version','-v',
        type=str, default='6.0', choices=choices,
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
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files')
    # return the parser
    return parser

# Main program that calls racmo_downscaled_mean()
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args = parser.parse_args()

    # create logger for verbosity level
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # run program for each input product
    for PRODUCT in args.product:
        # run downscaled cumulative program with parameters
        racmo_downscaled_mean(
            args.directory,
            args.region,
            args.version,
            PRODUCT,
            RANGE=args.mean,
            GZIP=args.gzip,
            MODE=args.mode
        )

# run main program
if __name__ == '__main__':
    main()
