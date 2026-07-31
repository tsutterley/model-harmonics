#!/usr/bin/env python
"""
reanalysis_geopotential_heights.py
Written by Tyler Sutterley (07/2026)
Reads temperature and specific humidity data to calculate geopotential height
    and pressure difference fields at half levels from reanalysis

Reanalysis models:
    ERA-Interim:
        http://apps.ecmwf.int/datasets/data/interim-full-moda
    ERA5:
        http://apps.ecmwf.int/data-catalogues/era5/?class=ea
    MERRA-2:
        https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -Y X, --year X: years to run
    -V, --verbose: Output information for each output file
    -M X, --mode X: Permission mode of directories and files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 07/2026: use struct dictionary to define netCDF4 parameters
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: use full path to output file in verbose logging
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: can use input files to define command line arguments
        added check for ERA5 expver dimension (denotes mix of ERA5 and ERA5T)
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 03/2021: automatically update years to run based on current time
    Updated 01/2021: read from netCDF4 file in slices to reduce memory load
    Updated 12/2020: using argparse to set command line options
    Updated 01/2020: outputs variables as 32-bit floats instead of 64-bit floats
        clear variables and iterate years to reduce required memory
        iterate over time variable to calculate heights using incomplete files
    Updated 08/2019: adjust time scale variable for MERRA-2
    Updated 07/2018: added parameters for ERA5
    Written 03/2018
"""

from __future__ import print_function

import sys
import re
import time
import logging
import netCDF4
import pathlib
import argparse
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc


# PURPOSE: reads temperature and specific humidity data to calculate
# geopotential height fields at half levels from reanalysis
def reanalysis_geopotential_heights(base_dir, MODEL, YEAR=None, MODE=0o775):
    # directory setup
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    ddir = base_dir.joinpath(MODEL)

    # set model specific parameters
    WEIGHT = None
    input_coordinate_file = None
    if MODEL == 'ERA-Interim':
        # invariant parameters file
        input_invariant_file = 'ERA-Interim-Invariant-Parameters.nc'
        # coordinate parameters file
        input_coordinate_file = 'ERA-Interim_coordvars.nc'
        # surface pressure file format
        input_pressure_file = 'ERA-Interim-Monthly-SP-{0}.nc'
        # regular expression pattern for finding files
        regex_pattern = r'ERA\-Interim\-Monthly\-Levels\-({0})\.nc$'
        # output file format
        output_file_format = 'ERA-Interim-GPH-Levels-{0}.nc'
        SURFNAME = 'z'
        ZNAME = 'z'
        VARNAME = 'sp'
        TNAME = 't'
        QNAME = 'q'
        DIFFNAME = 'dp'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'time'
        LEVELNAME = 'lvl'
        ANAME, BNAME = ('a_model_alt', 'b_model_alt')
        AINTERFACE, BINTERFACE = ('a_interface', 'b_interface')
        # hours since 1900-01-01 00:00:0.0
        TIME_LONGNAME = 'Time'
        UNITS = 'm**2 s**-2'
        GRAVITY = 1.0
    elif MODEL == 'ERA5':
        # invariant parameters file
        input_invariant_file = 'ERA5-Invariant-Parameters.nc'
        # coordinate parameters file
        input_coordinate_file = 'ERA5_coordvars.nc'
        # surface pressure file format
        input_pressure_file = 'ERA5-Monthly-SP-{0}.nc'
        # regular expression pattern for finding files
        regex_pattern = r'ERA5\-Monthly\-Levels\-({0})\.nc$'
        # output file format
        output_file_format = 'ERA5-GPH-Levels-{0}.nc'
        SURFNAME = 'z'
        ZNAME = 'z'
        VARNAME = 'sp'
        TNAME = 't'
        QNAME = 'q'
        DIFFNAME = 'dp'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'valid_time'
        LEVELNAME = 'lvl'
        ANAME, BNAME = ('a_half', 'b_half')
        AINTERFACE, BINTERFACE = ('a_interface', 'b_interface')
        # hours since 1900-01-01 00:00:0.0
        TIME_LONGNAME = 'Time'
        UNITS = 'm**2 s**-2'
        GRAVITY = 1.0
    elif MODEL == 'MERRA-2':
        # invariant parameters file
        input_invariant_file = 'MERRA2_101.const_2d_asm_Nx.00000000.nc4'
        # coordinate parameters file
        input_coordinate_file = 'MERRA2_101.const_3d_coords_Nx.00000000.nc4'
        # regular expression pattern for finding files
        regex_pattern = r'MERRA2_(\d+).tavgM_3d_ana_Nv.({0})(\d{{2}}).SUB.nc$'
        # output file format
        output_file_format = 'MERRA2_{0}.tavgM_3d_PHIS.{1}{2}.SUB.nc'
        SURFNAME = 'PHIS'
        ZNAME = 'PHIS'
        VARNAME = 'PS'
        TNAME = 'T'
        QNAME = 'QV'
        DIFFNAME = 'dP'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'
        LEVELNAME = 'lev'
        ANAME, BNAME = ('a_half', 'b_half')
        AINTERFACE, BINTERFACE = ('a_interface', 'b_interface')
        # minutes since start of file
        TIME_LONGNAME = 'Time'
        UNITS = 'm+2 s-2'
        GRAVITY = 1.0
    elif MODEL == 'JRA-3Q':
        # invariant parameters file
        input_invariant_file = (
            'jra3q.tl479_surf.0_3_4.gp-sfc-cn-gauss.1947090100_1947090100.nc'
        )
        # surface pressure and specific humidity file formats
        input_pressure_file = 'jra3q.anl_surf.pres-sfc-an-gauss.{0}{1}.nc'
        input_humidity_file = 'jra3q.anl_mdl.spfh-hyb-an-gauss.{0}{1}.nc'
        # regular expression pattern for finding files
        regex_pattern = r'jra3q\.anl_mdl\.tmp-hyb-an-gauss\.({0})(\d+).nc$'
        # output file format
        output_file_format = 'jra3q.anl_mdl.hgt-hyb-an-gauss.{0}{1}.nc'
        SURFNAME = 'gp-sfc-cn-gauss'
        ZNAME = 'hgt-hyb-an-gauss'
        VARNAME = 'pres-sfc-an-gauss'
        TNAME = 'tmp-hyb-an-gauss'
        QNAME = 'spfh-hyb-an-gauss'
        DIFFNAME = 'dpres-hyb-an-gauss'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'
        LEVELNAME = 'hybrid_level'
        WEIGHT = 'weight'
        ANAME, BNAME = ('a_hybrid_half_level', 'b_hybrid_half_level')
        AINTERFACE, BINTERFACE = ('a_hybrid_level', 'b_hybrid_level')
        # hours since 1900-01-01 00:00:00
        TIME_LONGNAME = 'time'
        UNITS = 'm2 s-2'
        GRAVITY = 1.0

    # dictionary defining output structure
    struct = dict(
        dimensions=(TIMENAME, LEVELNAME, LATNAME, LONNAME),
        variables={
            ZNAME: (TIMENAME, LEVELNAME, LATNAME, LONNAME),
            DIFFNAME: (TIMENAME, LEVELNAME, LATNAME, LONNAME),
        },
    )

    # dictionary defining file-level and variable attributes
    attributes = dict(ROOT={})
    # reference attribute
    REFERENCE = f'Output from {pathlib.Path(sys.argv[0]).name}'
    attributes['ROOT']['reference'] = REFERENCE
    # variable attributes
    attributes[LONNAME] = dict(
        long_name='Longitude',
        units='degrees_east',
    )
    attributes[LATNAME] = dict(
        long_name='Latitude',
        units='degrees_north',
    )
    attributes[TIMENAME] = dict(
        long_name=TIME_LONGNAME,
        standard_name='time',
        calendar='standard',
    )
    # Defining attributes for model levels
    attributes[LEVELNAME] = dict(
        long_name='Model_Level_Number',
    )
    # Defining attributes for geopotential height
    attributes[ZNAME] = dict(
        long_name='Geopotential_Heights_on_Model_Levels',
        units=UNITS,
    )
    # Defining attributes for pressure differences
    attributes[DIFFNAME] = dict(
        long_name='Pressure_Differences_between_Levels',
        units='Pa',
    )
    # append additional weight variable for gaussian grid models
    if MODEL in ('JRA-3Q',):
        struct['variables'][WEIGHT] = (LATNAME,)
        attributes[WEIGHT] = dict(
            long_name='gaussian weight', short_name='wgt', units='1'
        )

    # read model orography for dimensions
    geopotential, lon, lat = ncdf_invariant(
        ddir.joinpath(input_invariant_file), LONNAME, LATNAME, SURFNAME
    )
    # read parameters for calculating pressures at levels
    # flips the order of the levels so that bottom=0
    if input_coordinate_file is not None:
        lev, A, B, AI, BI = ncdf_coordinates(
            ddir.joinpath(input_coordinate_file),
            LEVELNAME,
            ANAME,
            BNAME,
            AINTERFACE,
            BINTERFACE,
        )
    # Gas constant for dry air
    R_dry = 287.06
    # minimum allowable pressure at the top of the atmosphere
    Pmin = 0.1

    # read each reanalysis pressure field for each year
    regex_years = r'\d{4}' if (YEAR is None) else '|'.join(map(str, YEAR))
    rx = re.compile(regex_pattern.format(regex_years), re.VERBOSE)
    input_files = sorted([f for f in ddir.iterdir() if rx.match(f.name)])
    # for each reanalysis file
    for i, temperature_file in enumerate(input_files):
        # read input temperature and specific humidity data
        logging.debug(str(temperature_file))
        fid = [None] * 3
        fid[0] = netCDF4.Dataset(temperature_file, mode='r')
        # extract shape from temperature variable
        ntime, nlevels, nlat, nlon = fid[0].variables[TNAME].shape
        # invalid value
        fill_value = fid[0].variables[TNAME]._FillValue
        # dictionary with output variables and dimensions
        dinput = {}
        # allocate for output variables
        for var in (ZNAME, DIFFNAME):
            dinput[var] = np.ma.zeros((ntime, nlevels, nlat, nlon), dtype='f')
            dinput[var].set_fill_value(fill_value)
        # extract time and time units
        dinput[TIMENAME] = np.copy(fid[0].variables[TIMENAME][:])
        attributes[TIMENAME]['units'] = fid[0].variables[TIMENAME].units
        # copy latitude and longitude
        dinput[LONNAME] = lon.copy()
        dinput[LATNAME] = lat.copy()

        if MODEL in ('MERRA-2'):
            # extract date from temperature files
            MOD, YEAR, MONTH = rx.findall(temperature_file.name).pop()
            # output monthly filename
            FILENAME = output_file_format.format(MOD, YEAR, MONTH)
            output_file = ddir.joinpath(FILENAME)
            # model levels
            dinput[LEVELNAME] = lev[:].copy()
            # specific humidity from same file as temperature
            fid[1] = fid[0]
            # read surface pressure
            pressure = fid[0].variables[VARNAME][:]
        elif MODEL in ('ERA-Interim', 'ERA5'):
            # extract year from temperature files
            (YEAR,) = rx.findall(temperature_file.name)
            # output yearly filename
            FILENAME = output_file_format.format(YEAR)
            output_file = ddir.joinpath(FILENAME)
            # model levels
            dinput[LEVELNAME] = lev[:].copy()
            # specific humidity from same file as temperature
            fid[1] = fid[0]
            # read input surface pressure data
            pressure_file = ddir.joinpath(input_pressure_file.format(YEAR))
            logging.debug(str(pressure_file))
            with netCDF4.Dataset(pressure_file, 'r') as fid[2]:
                pressure = np.copy(fid[2].variables[VARNAME][:])
        elif MODEL in ('JRA-3Q'):
            # extract date from temperature files
            YEAR, MONTH = rx.findall(temperature_file.name).pop()
            # output monthly filename
            FILENAME = output_file_format.format(YEAR, MONTH)
            output_file = ddir.joinpath(FILENAME)
            # extract the A and B coefficients for the hybrid levels
            dinput[LEVELNAME] = fid[0].variables[LEVELNAME][:]
            A = fid[0].variables[ANAME][:]
            B = fid[0].variables[BNAME][:]
            AI = fid[0].variables[AINTERFACE][:]
            BI = fid[0].variables[BINTERFACE][:]
            # extract the gaussian grid variables
            dinput[WEIGHT] = fid[0].variables[WEIGHT][:]
            # specific humidity data file
            humidity_file = ddir.joinpath(
                input_humidity_file.format(YEAR, MONTH)
            )
            logging.debug(str(humidity_file))
            fid[1] = netCDF4.Dataset(humidity_file, 'r')
            # read input surface pressure data
            pressure_file = ddir.joinpath(
                input_pressure_file.format(YEAR, MONTH)
            )
            logging.debug(str(pressure_file))
            with netCDF4.Dataset(pressure_file, 'r') as fid[2]:
                # reorder dimensions to match the required order
                dims = fid[2].variables[VARNAME].dimensions
                order = [dims.index(d) for d in (TIMENAME, LATNAME, LONNAME)]
                pressure = fid[2].variables[VARNAME][:].transpose(order)

        # iterate over dates
        for t in range(ntime):
            # extract temperature and specific humidity for time t
            if fid[0].variables[TNAME].ndim == 5:
                # check dimensions for expver slice
                t_time, expver = ncdf_expver(fid[0], t, TNAME)
                q_time, _ = ncdf_expver(fid[1], t, QNAME, index=expver)
            else:
                # temperature and specific humidity
                # reverse layers so bottom=0
                t_time = np.flipud(fid[0].variables[TNAME][t, :, :, :])
                q_time = np.flipud(fid[1].variables[QNAME][t, :, :, :])
            # convert to masked arrays
            t_time = np.ma.masked_equal(t_time, fill_value)
            q_time = np.ma.masked_equal(q_time, fill_value)
            # calculate geopotential over model levels
            geopotential_height = np.empty((nlat, nlon), dtype=np.float32)
            # start with surface geopotential (m^2/s^2)
            geopotential_height[:, :] = geopotential * GRAVITY
            # extract surface pressure for time t
            SP = pressure[t, :, :]
            # Integrate the model layers in the atmosphere
            for k in range(nlevels):
                # specific humidity and temperature for level k and time t
                T = t_time[k, :, :]
                QV = q_time[k, :, :]
                # calculate virtual temperature
                Tv = (1.0 + 0.609133 * QV) * T
                # calculate numerator and denominator for pressure ratio
                Pnum = A[k] + B[k] * SP
                if (k + 1) == len(A):
                    Pdom = Pmin
                else:
                    # add a threshold to avoid dividing by zero
                    Pdom = np.maximum(A[k + 1] + B[k + 1] * SP, Pmin)
                # calculate geopotential difference using hypsometric equation
                # add level to geopotential_levels
                geopotential_height[:, :] += R_dry * Tv * np.log(Pnum / Pdom)
                # save level to output variable
                dinput[ZNAME][t, k, :, :] = geopotential_height
                # calculate pressure difference between levels (at interfaces)
                Plower = AI[k] + BI[k] * SP
                # check if there is an upper bound
                if (k + 1) == len(AI):
                    Pupper = Pmin
                else:
                    # add a threshold limit to the atmospheric pressure
                    Pupper = np.maximum(AI[k + 1] + BI[k + 1] * SP, Pmin)
                dinput[DIFFNAME][t, k, :, :] = Pupper - Plower
        # replace invalid values with fill_value
        for var in (ZNAME, DIFFNAME):
            dinput[var].data[dinput[var].mask] = dinput[var].fill_value
        # write structured data to netCDF4 file
        mdlhmc.spatial.to_netCDF4(output_file, dinput, attributes, struct)
        # set the permissions level of the output file to MODE
        output_file.chmod(mode=MODE)
        # clear dinput dictionary variable
        dinput = None
        # close the input netCDF4 files
        ncdf_close_all(fid)


# PURPOSE: Compute the Specific Humidity from parameters (Bolton 1980)
# http://cires1.colorado.edu/~voemel/vp.html
# https://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
# https://github.com/NCAR/ncl/blob/master/ni/src/lib/nfpfort/mixhum_ptrh.f
def calculate_specific_humidity(P, T, RH):
    # ratio of the molecular weights of water vapor to dry air
    epsilon = 0.622
    # calibration pressure and temperature
    pc = 6.112
    tc = 243.5
    # calculate the saturation vapor pressure in mb
    Es = pc * np.exp((17.67 * T) / (T + tc))
    # calculate the vapor pressure in mb
    Ev = Es * (RH / 100.0)
    # calculate the dew point temperature
    Td = np.log(Ev / pc) * tc / (17.67 - np.log(Ev / pc))
    # calculate the specific humidity
    Q = (epsilon * Ev) / (P / 100.0 - (0.378 * Ev))
    return (Q, Td)


# PURPOSE: extract temperature and specific humidity variables
# from a 5d netCDF4 dataset
# ERA5 expver dimension (denotes mix of ERA5 and ERA5T)
def ncdf_expver(fileID, slice, VARNAME, index=None):
    ntime, nexp, nlevel, nlat, nlon = fileID.variables[VARNAME].shape
    # if expver is specified, check if data is valid for expver
    if index is not None and np.any(fileID.variables[VARNAME][slice, index, :]):
        # reverse layers so bottom=0
        variable = np.flipud(fileID.variables[VARNAME][slice, index, :])
        # return the reduced variables
        return variable, index
    # reduced variable (temperature or specific humidity) for time
    variable = np.zeros((nlevel, nlat, nlon))
    # iterate over expver slices to find valid outputs
    for j in range(nexp):
        # check if any are valid for expver
        if np.any(fileID.variables[VARNAME][slice, j, :]):
            # reverse layers so bottom=0
            variable = np.flipud(fileID.variables[VARNAME][slice, j, :])
            index = j
            break
    # return the reduced variables and the valid expver slice
    return variable, index


# PURPOSE: read reanalysis invariant parameters (geopotential,lat,lon)
def ncdf_invariant(FILENAME, LONNAME, LATNAME, ZNAME):
    logging.debug(str(FILENAME))
    with netCDF4.Dataset(FILENAME, mode='r') as fileID:
        geopotential = fileID.variables[ZNAME][:].squeeze()
        longitude = fileID.variables[LONNAME][:].copy()
        latitude = fileID.variables[LATNAME][:].copy()
    return (geopotential, longitude, latitude)


# PURPOSE: read reanalysis coordinate parameters
# reverse order to go from surface to top-of-atmosphere
def ncdf_coordinates(FILENAME, LEVELNAME, ANAME, BNAME, AINTERFACE, BINTERFACE):
    logging.debug(str(FILENAME))
    with netCDF4.Dataset(FILENAME, mode='r') as fileID:
        # reverse layers so bottom=0
        levels = np.flip(fileID.variables[LEVELNAME][:])
        A = np.flip(fileID.variables[ANAME][:])
        B = np.flip(fileID.variables[BNAME][:])
        AI = np.flip(fileID.variables[AINTERFACE][:])
        BI = np.flip(fileID.variables[BINTERFACE][:])
    return (levels, A, B, AI, BI)


# PURPOSE: attempt to close all open netCDF4 files
def ncdf_close_all(fileID: list[__loader__]):
    [fid.close() for fid in fileID if fid._isopen]


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads temperature and specific humidity data
            to calculate geopotential height and pressure difference
            fields at half levels from reanalysis
            """,
        fromfile_prefix_chars='@',
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    choices = ['ERA-Interim', 'ERA5', 'JRA-3Q', 'MERRA-2']
    parser.add_argument(
        'model',
        type=str,
        nargs='+',
        default=['ERA5', 'MERRA-2'],
        choices=choices,
        help='Reanalysis Model',
    )
    # working data directory
    parser.add_argument(
        '--directory',
        '-D',
        type=pathlib.Path,
        default=pathlib.Path.cwd(),
        help='Working data directory',
    )
    # years to run
    now = time.gmtime()
    parser.add_argument(
        '--year',
        '-Y',
        type=int,
        nargs='+',
        default=range(2000, now.tm_year + 1),
        help='Years of model outputs to run',
    )
    # print information about each input and output file
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


# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args, _ = parser.parse_known_args()

    # create logger
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # for each reanalysis model
    for MODEL in args.model:
        # run program
        reanalysis_geopotential_heights(
            args.directory, MODEL, YEAR=args.year, MODE=args.mode
        )


# run main program
if __name__ == '__main__':
    main()
