#!/usr/bin/env python
"""
racmo_downscaled_cumulative.py
Written by Tyler Sutterley (06/2024)
Calculates cumulative anomalies of RACMO surface mass balance products

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
            6.1: RACMO2.3p2/FGRN055/DS1km
        ais:
            6.0: RACMO2.3/XANT27/DS2km
            6.1: RACMO2.3/ANT27/DS2km
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
    Updated 07/2026: added version 6.1 for Greenland and Antarctica
        output using structured dictionary with netCDF4 parameters
    Updated 06/2024: added version 6.0 for Greenland and Antarctica
        RACMO2.3p2/FGRN055/DS1km for 1958-2023
        RACMO2.3p2/XANT27/DS2km for 1979-2023
    Updated 06/2023: added version 5.0 (RACMO2.3p2 for 1958-2023 from FGRN055)
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

# dictionary defining the spatial reference system for each region
crs = dict(gris={}, ais={})
# Greenland (EPSG:3413)
crs['gris']['standard_name'] = 'Polar_Stereographic'
crs['gris']['grid_mapping_name'] = 'polar_stereographic'
crs['gris']['straight_vertical_longitude_from_pole'] = -45.0
crs['gris']['latitude_of_projection_origin'] = 90.0
crs['gris']['standard_parallel'] = 70.0
crs['gris']['scale_factor_at_projection_origin'] = 1.0
crs['gris']['false_easting'] = 0.0
crs['gris']['false_northing'] = 0.0
crs['gris']['semi_major_axis'] = 6378.137
crs['gris']['semi_minor_axis'] = 6356.752
crs['gris']['inverse_flattening'] = 298.257223563
crs['gris']['spatial_epsg'] = '3413'
# Antarctica (EPSG:3031)
crs['ais']['standard_name'] = 'Polar_Stereographic'
crs['ais']['grid_mapping_name'] = 'polar_stereographic'
crs['ais']['straight_vertical_longitude_from_pole'] = 0.0
crs['ais']['latitude_of_projection_origin'] = -90.0
crs['ais']['standard_parallel'] = -71.0
crs['ais']['scale_factor_at_projection_origin'] = 1.0
crs['ais']['false_easting'] = 0.0
crs['ais']['false_northing'] = 0.0
crs['ais']['semi_major_axis'] = 6378.137
crs['ais']['semi_minor_axis'] = 6356.752
crs['ais']['inverse_flattening'] = 298.257223563
crs['ais']['spatial_epsg'] = '3031'


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
    if VERSION in ('1.0', '4.0', '5.0', '6.0', '6.1'):
        # find input files
        pattern = rf'{VARIABLE}.(\d+).BN_(.*?).MM.nc(\.gz)?'
        rx = re.compile(pattern, re.VERBOSE | re.IGNORECASE)
        infiles = sorted([f for f in input_dir.iterdir() if rx.match(f.name)])
        nt = 12 * len(infiles)
        # read netCDF file for dataset (could also set memory=None)
        if GZIP:
            # read bytes from compressed file
            fd = gzip.open(str(infiles[0]), 'rb')
            # read netCDF file for dataset from bytes
            fileID = netCDF4.Dataset(
                uuid.uuid4().hex, mode='r', memory=fd.read()
            )
        else:
            fileID = netCDF4.Dataset(infiles[0], mode='r')
        # shape of the input data matrix
        (ncvar,) = [v for v in fileID.variables.keys() if regex.match(v)]
        _, ny, nx = fileID.variables[ncvar].shape
        # close the (compressed) file objects
        fileID.close()
        if GZIP:
            fd.close()
    elif VERSION in ('2.0', '3.0'):
        # if reading bytes from compressed file or netcdf file directly
        gz = '.gz' if GZIP else ''
        # input dataset for variable
        file_format = {}
        file_format['2.0'] = '{0}.1958-2016.BN_RACMO2.3p2_FGRN11_GrIS.MM.nc{1}'
        file_format['3.0'] = '{0}.1958-2016.BN_RACMO2.3p2_FGRN055_GrIS.MM.nc{1}'
        f = input_dir.joinpath(
            file_format[VERSION].format(VARIABLE.lower(), gz)
        )
        if GZIP:
            # read bytes from compressed file
            fd = gzip.open(str(f), 'rb')
            # read netCDF file for dataset from bytes
            fileID = netCDF4.Dataset(
                uuid.uuid4().hex, mode='r', memory=fd.read()
            )
        else:
            # read netCDF file for dataset (could also set memory=None)
            fileID = netCDF4.Dataset(f, mode='r')
        # shape of the input data matrix
        (ncvar,) = [v for v in fileID.variables.keys() if regex.match(v)]
        nt, ny, nx = fileID.variables[ncvar].shape
        # close the (compressed) file objects
        fileID.close()
        if GZIP:
            fd.close()
    # return the data dimensions
    return (nt, ny, nx)


# PURPOSE: read individual yearly netcdf files and calculate anomalies
def yearly_file_cumulative(
    input_dir, REGION, VERSION, PRODUCT, MEAN, GZIP=False
):
    """Read individual yearly netcdf files and calculate
    cumulative anomalies

    Parameters
    ----------
    input_dir: str
        Working data directory
    REGION: str
        RACMO model region
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
    # get logger
    logger = logging.getLogger(__name__)
    # verify input directory
    input_dir = pathlib.Path(input_dir).expanduser().absolute()
    # names within netCDF4 files
    VARIABLE = input_products[PRODUCT]
    # regular expression operator for finding variables
    regex = re.compile(VARIABLE, re.VERBOSE | re.IGNORECASE)
    # find input files for years of interest
    pattern = rf'{VARIABLE}.(\d+).BN_(.*?).MM.nc(\.gz)?'
    rx = re.compile(pattern, re.VERBOSE | re.IGNORECASE)
    input_files = sorted([f for f in input_dir.iterdir() if rx.match(f.name)])
    # number of input files
    n_files = len(input_files)
    # input dimensions and counter variable
    # get dimensions for input VERSION
    nt, ny, nx = get_dimensions(input_dir, VERSION, PRODUCT, GZIP=GZIP)
    # create counter variable
    c = 0
    # allocate for all data
    dinput = {}
    # dictionary defining file-level and variable attributes
    attributes = dict(ROOT={})
    attributes_list = ['long_name', 'units', 'standard_name', 'axis']

    # create variable and attributes for projection
    dinput['Polar_Stereographic'] = np.byte()
    # add projection attributes to dictionary
    attributes['Polar_Stereographic'] = crs[REGION]
    # create variables for coordinates and attributes
    dinput['LON'] = np.zeros((ny, nx))
    dinput['LAT'] = np.zeros((ny, nx))
    dinput['x'] = np.zeros((nx))
    dinput['y'] = np.zeros((ny))
    dinput['TIME'] = np.zeros((nt))
    dinput['MASK'] = np.zeros((ny, nx), dtype=np.int8)
    dinput[VARIABLE] = np.zeros((nt, ny, nx))
    CUMULATIVE = np.zeros((ny, nx))

    # if reading bytes from compressed file or netcdf file directly
    gz = '.gz' if GZIP else ''
    # input area file with ice mask and model topography
    use_icemask = False
    use_icemask |= VERSION in ('4.0', '5.0')
    use_icemask |= (VERSION == '6.1') and (REGION.lower() == 'gris')
    if use_icemask:
        if VERSION in ('4.0', '5.0'):
            # input area file with ice mask and model topography
            f1 = f'Icemask_Topo_Iceclasses_lon_lat_average_1km_GrIS.nc{gz}'
            # full path to input mask file
            input_mask_file = input_dir.joinpath(f1)
        elif (VERSION == '6.1') and (REGION.lower() == 'gris'):
            # input area file with ice mask and model topography
            f1 = f'GrIS_topo_icemask_lsm_lon_lat_1km.nc{gz}'
            # full path to input mask file
            input_mask_file = input_dir.parents[1].joinpath(f1)
        # log the input mask file
        logger.debug(str(input_mask_file))
        # read netCDF file for topography and ice classes
        if GZIP:
            # read bytes from compressed file
            fd = gzip.open(input_mask_file, 'rb')
            # read netCDF file for topography and ice classes from bytes
            fileID = netCDF4.Dataset(
                uuid.uuid4().hex, mode='r', memory=fd.read()
            )
        else:
            # read netCDF file for topography and ice classes
            fileID = netCDF4.Dataset(input_mask_file, mode='r')
        # get coordinate data from each netCDF variable
        for var in ['LON', 'LAT', 'x', 'y']:
            dinput[var] = np.array(fileID.variables[var][:])
            attributes[var] = ncdf_attributes(
                fileID.variables[var], attributes_list
            )
            attributes[var]['grid_mapping'] = 'Polar_Stereographic'
        # find ice sheet points from promicemask that valid
        promicemask = fileID.variables['Promicemask'][:, :].copy()
        # close the (compressed) file objects
        fileID.close()
        if GZIP:
            fd.close()
        # find ice sheet points from promicemask that valid
        ii, jj = np.nonzero((promicemask >= 1) & (promicemask <= 3))
        # create ice mask for output
        dinput['MASK'] = np.zeros((ny, nx), dtype=np.int8)
        dinput['MASK'][ii, jj] = 1

    # for each file of interest
    for t in range(n_files):
        # log the input file
        logger.debug(str(input_files[t]))
        # Open the NetCDF file for reading
        if GZIP:
            # read bytes from compressed file
            fd = gzip.open(str(input_files[t]), 'rb')
            # read netCDF file for dataset from bytes
            fileID = netCDF4.Dataset(
                uuid.uuid4().hex, mode='r', memory=fd.read()
            )
        else:
            # read netCDF file for dataset (could also set memory=None)
            fileID = netCDF4.Dataset(input_files[t], mode='r')
        # check if ERA5 3-hourly data
        ERA5_3h = re.search(r'ERA5_3h', input_files[t].name)
        # Getting the data from each netCDF variable
        if VERSION == '1.0':
            # get coordinate data from each netCDF variable
            for var in ['LON', 'LAT', 'x', 'y']:
                dinput[var] = np.array(fileID.variables[var][:])
                attributes[var] = ncdf_attributes(
                    fileID.variables[var], attributes_list
                )
                attributes[var]['grid_mapping'] = 'Polar_Stereographic'
            icemask = fileID.variables['icemask'][:, :]
            dinput['MASK'][:, :] = icemask.astype(np.int8)
        elif (VERSION == '6.0') or (VERSION == '6.1' and REGION == 'ais'):
            # extract coordinates
            if ERA5_3h or (VERSION == '6.1' and REGION == 'ais'):
                for dim in ['x', 'y']:
                    dinput[dim][:] = fileID.variables[dim][:].copy()
                    attributes[dim] = ncdf_attributes(
                        fileID.variables[dim], attributes_list
                    )
            else:
                for dim, var in zip(['x', 'y'], ['lon', 'lat']):
                    dinput[dim][:] = fileID.variables[var][:].copy()
                    attributes[dim] = ncdf_attributes(
                        fileID.variables[var], attributes_list
                    )
            # coordianate reference system
            EPSG = 3413 if REGION == 'gris' else 3031
            # convert from edges to centers
            dx = np.diff(dinput['x'])[0]
            dy = np.diff(dinput['y'])[0]
            # calculate latitude and longitude of grid cells
            dinput['LON'], dinput['LAT'] = mdlhmc.spatial.get_latlon(
                dinput['x'] - dx / 2.0, dinput['y'] - dy / 2.0, srs_epsg=EPSG
            )
            # find variable of interest
            (ncvar,) = [v for v in fileID.variables.keys() if regex.match(v)]
            valid = np.any(fileID.variables[ncvar], axis=0)
            dinput['MASK'] = valid.astype(np.int8)

        # calculate dates from delta times
        delta_time = fileID.variables['time'][:].copy()
        date_string = fileID.variables['time'].units
        epoch, to_secs = gravtk.time.parse_date_string(date_string)
        # calculate time array in Julian days
        JD = 2400000.5 + gravtk.time.convert_delta_time(
            to_secs * delta_time,
            epoch1=epoch,
            epoch2=(1858, 11, 17, 0, 0, 0),
            scale=1.0 / 86400.0,
        )
        # for each month
        for m in range(12):
            # convert from Julian days to calendar dates
            YY, MM, DD, hh, mm, ss = gravtk.time.convert_julian(
                JD[m], format='tuple'
            )
            # calculate time in year-decimal
            dinput['TIME'][c] = gravtk.time.convert_calendar_decimal(
                YY, MM, day=DD, hour=hh, minute=mm, second=ss
            )
            # find variable of interest
            (ncvar,) = [v for v in fileID.variables.keys() if regex.match(v)]
            # extract data and add to total cumulative matrix
            CUMULATIVE += fileID.variables[ncvar][m, :, :].copy() - MEAN
            dinput[VARIABLE][c, :, :] = CUMULATIVE.copy()
            # add to counter
            c += 1
        # close the (compressed) file objects
        fileID.close()
        if GZIP:
            fd.close()

    # return the cumulative anomalies
    return dinput, attributes


# PURPOSE: read compressed netCDF4 files and calculate cumulative anomalies
def compressed_file_cumulative(
    input_dir, REGION, VERSION, PRODUCT, MEAN, GZIP=False
):
    """Read compressed netCDF4 files and calculate cumulative anomalies

    Parameters
    ----------
    input_dir: str
        Working data directory
    REGION: str
        RACMO model region
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
    # get logger
    logger = logging.getLogger(__name__)
    # verify input directory
    input_dir = pathlib.Path(input_dir).expanduser().absolute()
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

    # dictionary defining file-level and variable attributes
    attributes = dict(ROOT={})
    attributes_list = ['long_name', 'units', 'standard_name', 'axis']

    # create variable and attributes for projection
    dinput['Polar_Stereographic'] = np.byte()
    # add projection attributes to dictionary
    attributes['Polar_Stereographic'] = crs[REGION]

    # input area file with ice mask and model topography
    f1 = f'Icemask_Topo_Iceclasses_lon_lat_average_1km_GrIS.nc{gz}'
    input_mask_file = input_dir.joinpath(f1)
    logger.debug(str(input_mask_file))
    if GZIP:
        # read bytes from compressed file
        fd = gzip.open(input_mask_file, 'rb')
        # read netCDF file for topography and ice classes from bytes
        fileID = netCDF4.Dataset(uuid.uuid4().hex, mode='r', memory=fd.read())
    else:
        # read netCDF file for topography and ice classes
        fileID = netCDF4.Dataset(input_mask_file, mode='r')
    # Getting the data from each netCDF variable
    for var in ['LON', 'LAT', 'x', 'y']:
        dinput[var] = np.array(fileID.variables[var][:])
        attributes[var] = ncdf_attributes(
            fileID.variables[var], attributes_list
        )
        attributes[var]['grid_mapping'] = 'Polar_Stereographic'
    # find ice sheet points from promicemask that valid
    promicemask = fileID.variables['Promicemask'][:, :].copy()

    # close the (compressed) file objects
    fileID.close()
    if GZIP:
        fd.close()

    # file format for each version
    file_format = {}
    file_format['2.0'] = '{0}.1958-2016.BN_RACMO2.3p2_FGRN11_GrIS.MM.nc{1}'
    file_format['3.0'] = '{0}.1958-2018.BN_RACMO2.3p2_FGRN055_GrIS.MM.nc{1}'

    # input dataset for variable
    f2 = file_format[VERSION].format(VARIABLE.lower(), gz)
    input_file = input_dir.joinpath(f2)
    logger.debug(str(input_file))
    if GZIP:
        # read bytes from compressed file
        fd = gzip.open(input_file, 'rb')
        # read netCDF file for dataset from bytes
        fileID = netCDF4.Dataset(uuid.uuid4().hex, mode='r', memory=fd.read())
    else:
        # read netCDF file for dataset (could also set memory=None)
        fileID = netCDF4.Dataset(input_file, mode='r')
    # shape of the input data matrix
    nt, ny, nx = fileID.variables[VARNAME].shape

    # find ice sheet points from promicemask that valid
    ii, jj = np.nonzero((promicemask >= 1) & (promicemask <= 3))
    dinput['MASK'] = np.zeros((ny, nx), dtype=np.int8)
    dinput['MASK'][ii, jj] = 1
    attributes['MASK']['grid_mapping'] = 'Polar_Stereographic'

    # calculate dates from delta times
    # Months since 1958-01-15 at 00:00:00
    delta_time = fileID.variables['time'][:].copy()
    date_string = fileID.variables['time'].units
    epoch, to_secs = gravtk.time.parse_date_string(date_string)
    # calculate time array in Julian days
    JD = 2400000.5 + gravtk.time.convert_delta_time(
        to_secs * delta_time,
        epoch1=epoch,
        epoch2=(1858, 11, 17, 0, 0, 0),
        scale=1.0 / 86400.0,
    )
    # convert from Julian days to calendar dates
    YY, MM, DD, hh, mm, ss = gravtk.time.convert_julian(JD, format='tuple')
    # calculate time in year-decimal
    dinput['TIME'] = gravtk.time.convert_calendar_decimal(
        YY, MM, day=DD, hour=hh, minute=mm, second=ss
    )

    # calculate cumulative
    CUMULATIVE = np.zeros((ny, nx))
    dinput[VARNAME] = np.zeros((nt, ny, nx))
    # get attributes for variable of interest
    attributes[VARNAME] = ncdf_attributes(
        fileID.variables[VARNAME], attributes_list
    )
    attributes[VARNAME]['grid_mapping'] = 'Polar_Stereographic'
    for t in range(nt):
        # extract data and add to total cumulative matrix
        CUMULATIVE += fileID.variables[VARNAME][t, :, :].copy() - MEAN
        dinput[VARNAME][t, :, :] = CUMULATIVE.copy()

    # close the (compressed) file objects
    fileID.close()
    if GZIP:
        fd.close()

    # return the cumulative anomalies
    return dinput, attributes


# PURPOSE: get attributes for netCDF4 files and variable
def ncdf_attributes(nc, attributes_list):
    # get logger
    logger = logging.getLogger(__name__)
    # output dictionary of attributes
    attributes = {}
    # for each attribute to try to get
    for att_name in attributes_list:
        rx = re.compile(att_name, re.IGNORECASE)
        try:
            # use case-insensitive regex to find attribute name
            (ncattr,) = [s for s in nc.ncattrs() if rx.match(s)]
            att_val = nc.getncattr(ncattr)
        except Exception as exc:
            ncvar = getattr(nc, 'name', 'ROOT')
            logger.debug(f'Attribute {att_name} not found in {ncvar}')
            continue
        # strip whitespace for string attributes
        if isinstance(att_val, str):
            att_val = att_val.strip()
        # add attribute to dictionary
        attributes[att_name] = att_val
    # return the dictionary of attributes
    return attributes


# PURPOSE: calculate RACMO cumulative anomalies with respect to a mean field
def racmo_downscaled_cumulative(
    base_dir,
    REGION,
    VERSION,
    PRODUCT,
    RANGE=[1961, 1990],
    GZIP=False,
    MODE=0o775,
):
    """
    Calculate RACMO cumulative anomalies with respect to a mean field

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
                - ``6.1``: RACMO2.3p2/FGRN055/DS1km
            - ``'ais'``:
                - ``6.0``: RACMO2.3/XANT27/DS2km
                - ``6.1``: RACMO2.3/ANT27/DS2km
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
        RACMO_MODEL = ['XGRN11', '2.3', 'DS1km']
        VARNAME = input_products[PRODUCT]
        SUBDIRECTORY = f'{VARNAME}_v{VERSION}'
        input_dir = base_dir.joinpath(f'SMB1km_v{VERSION}', SUBDIRECTORY)
        file_type = 'yearly'
    elif (VERSION == '2.0') and (REGION == 'gris'):
        RACMO_MODEL = ['XGRN11', '2.3p2', 'DS1km']
        var = input_products[PRODUCT]
        VARNAME = var if PRODUCT in ('SMB', 'PRECIP') else f'{var}corr'
        input_dir = base_dir.joinpath(f'SMB1km_v{VERSION}')
        file_type = 'compressed'
    elif (VERSION == '3.0') and (REGION == 'gris'):
        RACMO_MODEL = ['FGRN055', '2.3p2', 'DS1km']
        var = input_products[PRODUCT]
        VARNAME = var if (PRODUCT == 'SMB') else f'{var}corr'
        input_dir = base_dir.joinpath(f'SMB1km_v{VERSION}')
        file_type = 'compressed'
    elif (VERSION == '4.0') and (REGION == 'gris'):
        RACMO_MODEL = ['FGRN055', '2.3p2', 'DS1km']
        var = input_products[PRODUCT]
        VARNAME = var if (PRODUCT == 'SMB') else f'{var}corr'
        SUBDIRECTORY = PRODUCT.lower()
        input_dir = base_dir.joinpath(f'SMB1km_v{VERSION}', SUBDIRECTORY)
        file_type = 'yearly'
    elif (VERSION == '5.0') and (REGION == 'gris'):
        RACMO_MODEL = ['FGRN055', '2.3p2', 'DS1km']
        var = input_products[PRODUCT]
        VARNAME = var if (PRODUCT == 'SMB') else f'{var}corr'
        SUBDIRECTORY = PRODUCT.lower()
        input_dir = base_dir.joinpath(f'SMB1km_v{VERSION}', SUBDIRECTORY)
        file_type = 'yearly'
    elif (VERSION == '6.0') and (REGION == 'gris'):
        RACMO_MODEL = ['FGRN055', '2.3p2', 'DS1km']
        var = input_products[PRODUCT]
        VARNAME = var if (PRODUCT == 'SMB') else f'{var}corr'
        input_dir = base_dir.joinpath('GrIS-1km')
        file_type = 'yearly'
    elif (VERSION == '6.0') and (REGION == 'ais'):
        RACMO_MODEL = ['XANT27', '2.3p2', 'DS2km']
        var = input_products[PRODUCT]
        VARNAME = var if (PRODUCT == 'SMB') else f'{var}corr'
        input_dir = base_dir.joinpath('AIS-2km')
        file_type = 'yearly'
    elif (VERSION == '6.1') and (REGION == 'gris'):
        RACMO_MODEL = ['FGRN055', '2.3p2', 'DS1km']
        var = input_products[PRODUCT]
        VARNAME = var if (PRODUCT == 'SMB') else f'{var}corr'
        input_dir = base_dir.joinpath('GrIS-1km', 'Monthly-1km', var.lower())
        file_type = 'yearly'
    elif (VERSION == '6.1') and (REGION == 'ais'):
        RACMO_MODEL = ['ANT27', '2.3p2', 'DS2km']
        var = input_products[PRODUCT]
        VARNAME = var if (PRODUCT == 'SMB') else f'{var}corr'
        input_dir = base_dir.joinpath('AIS-2km', 'Monthly-2km', var.lower())
        file_type = 'yearly'

    # read mean from netCDF4 file
    arg = (*RACMO_MODEL, VERSION, PRODUCT, *RANGE)
    mean_file = '{0}_RACMO{1}_{2}_v{3}_{4}_Mean_{5:4d}-{6:4d}.nc'.format(*arg)
    with netCDF4.Dataset(input_dir.joinpath(mean_file), mode='r') as fileID:
        MEAN = fileID[VARNAME][:, :].copy()

    # dictionary defining output structure
    struct = dict(
        dimensions=('TIME', 'y', 'x'),
        variables={
            VARNAME: ('TIME', 'y', 'x'),
            'MASK': ('y', 'x'),
            'LAT': ('y', 'x'),
            'LON': ('y', 'x'),
            'Polar_Stereographic': (),
        },
    )

    # read data and calculate cumulative
    if file_type == 'yearly':
        dinput, attributes = yearly_file_cumulative(
            input_dir, REGION, VERSION, PRODUCT, MEAN, GZIP=GZIP
        )
    elif file_type == 'compressed':
        dinput, attributes = compressed_file_cumulative(
            input_dir, REGION, VERSION, PRODUCT, MEAN, GZIP=GZIP
        )
    # edit file-level and variable attributes
    TITLE = 'Downscaled_cumulative_anomalies_relative_to_{0:4d}-{1:4d}_mean'
    REFERENCE = f'Output from {pathlib.Path(sys.argv[0]).name}'
    attributes['ROOT']['title'] = TITLE.format(*RANGE)
    attributes['ROOT']['model'] = ''.join(RACMO_MODEL)
    attributes['ROOT']['version'] = VERSION
    attributes['ROOT']['region'] = REGION
    attributes['ROOT']['reference'] = REFERENCE
    attributes[VARNAME]['long_name'] = longname[PRODUCT]
    attributes[VARNAME]['units'] = 'mmWE'
    # Defining attributes for date
    attributes['TIME'] = dict(
        long_name='time',
        units='Date_in_Decimal_Years',
        calendar='standard',
    )

    # output cumulative as netCDF4 file
    args = (*RACMO_MODEL, VERSION, PRODUCT)
    FILE = '{0}_RACMO{1}_{2}_v{3}_{4}_cumul.nc'.format(*args)
    output_file = input_dir.joinpath(FILE)
    mdlhmc.spatial.to_netCDF4(output_file, dinput, attributes, struct)
    # change the permission mode
    output_file.chmod(mode=MODE)


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
    parser.add_argument(
        '--directory',
        '-D',
        type=pathlib.Path,
        default=pathlib.Path.cwd(),
        help='Working data directory',
    )
    # region of RACMO model
    parser.add_argument(
        '--region',
        '-R',
        type=str,
        default='gris',
        choices=['gris', 'ais'],
        help='Region of RACMO model to calculate',
    )
    # Downscaled version
    # 1.0: RACMO2.3/XGRN11/DS1km
    # 2.0: RACMO2.3p2/XGRN11/DS1km
    # 3.0: RACMO2.3p2/FGRN055/DS1km
    # 4.0: RACMO2.3p2/FGRN055/DS1km
    # 5.0: RACMO2.3p2/FGRN055/DS1km
    # 6.0: RACMO2.3p2/FGRN055/DS1km
    # 6.0: RACMO2.3p2/XANT27/DS2km
    # 6.1: RACMO2.3p2/FGRN055/DS1km
    # 6.1: RACMO2.3p2/ANT27/DS2km
    parser.add_argument(
        '--version',
        '-v',
        type=str,
        default='6.0',
        help='Downscaled RACMO Version',
    )
    # Products to calculate cumulative
    parser.add_argument(
        '--product',
        '-p',
        metavar='PRODUCT',
        type=str,
        nargs='+',
        default=['SMB'],
        choices=input_products.keys(),
        help='RACMO product to calculate',
    )
    # start and end years to run for mean
    parser.add_argument(
        '--mean',
        '-m',
        metavar=('START', 'END'),
        type=int,
        nargs=2,
        default=[1961, 1990],
        help='Start and end year range for mean',
    )
    # netCDF4 files are gzip compressed
    parser.add_argument(
        '--gzip',
        '-G',
        default=False,
        action='store_true',
        help='netCDF4 file is locally gzip compressed',
    )
    # verbose output of processing run
    parser.add_argument(
        '--verbose',
        '-V',
        action='count',
        default=0,
        help='Verbose output of processing run',
    )
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument(
        '--mode',
        '-M',
        type=lambda x: int(x, base=8),
        default=0o775,
        help='Permission mode of directories and files',
    )
    # return the parser
    return parser


# Main program that calls racmo_downscaled_cumulative()
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args = parser.parse_args()

    # create logger for verbosity level
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logger = gravtk.utilities.build_logger(
        __name__, level=loglevels[args.verbose]
    )

    # run program for each input product
    for PRODUCT in args.product:
        # run downscaled cumulative program with parameters
        racmo_downscaled_cumulative(
            args.directory,
            args.region,
            args.version,
            PRODUCT,
            RANGE=args.mean,
            GZIP=args.gzip,
            MODE=args.mode,
        )


# run main program
if __name__ == '__main__':
    main()
