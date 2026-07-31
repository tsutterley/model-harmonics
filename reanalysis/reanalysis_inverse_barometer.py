#!/usr/bin/env python
"""
reanalysis_inverse_barometer.py
Written by Tyler Sutterley (07/2026)
Reads hourly mean sea level pressure fields from reanalysis and
    calculates the inverse-barometer response

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
    -m X, --mean X: Start and end year range for mean
    -d X, --density X: Density of seawater in kg/m^3
    -V, --verbose: Output information about each created file
    -M X, --mode X: Permission mode of directories and files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for files

REFERENCES:
    Wunsch and Stammer. Atmospheric loading and the oceanic "inverted
        barometer" effect. Reviews of Geophysics, 35(1), 79-107, (1997).
        https://doi.org/10.1029/96RG03037

    Hofmann-Wellenhof and Moritz. Physical Geodesy, (2005).
        https://doi.org/10.1007/978-3-211-33545-1

UPDATE HISTORY:
    Updated 07/2026: use struct dictionary to define netCDF4 parameters
        use authalic area for the grid cell areas
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: use full path to output file in verbose logging
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: can use input files to define command line arguments
    Written 03/2021
"""

from __future__ import print_function

import sys
import re
import logging
import netCDF4
import pathlib
import argparse
import datetime
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc


# PURPOSE: read land sea mask to get indices of oceanic values
def ncdf_landmask(FILENAME, MASKNAME, OCEAN):
    logging.debug(str(FILENAME))
    with netCDF4.Dataset(FILENAME, mode='r') as fileID:
        landsea = np.squeeze(fileID.variables[MASKNAME][:].copy())
    return landsea == OCEAN


# PURPOSE: read reanalysis mean sea level pressure
def ncdf_mean_pressure(FILENAME, VARNAME, LONNAME, LATNAME):
    logging.debug(str(FILENAME))
    with netCDF4.Dataset(FILENAME, mode='r') as fileID:
        # extract pressure and remove singleton dimensions
        mean_pressure = np.array(fileID.variables[VARNAME][:].squeeze())
        longitude = fileID.variables[LONNAME][:].squeeze()
        latitude = fileID.variables[LATNAME][:].squeeze()
    return (mean_pressure, longitude, latitude)


# PURPOSE:  calculate the instantaneous inverse barometer response
def reanalysis_inverse_barometer(
    base_dir, MODEL, YEAR=None, RANGE=None, DENSITY=None, MODE=0o775
):
    # directory setup
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    ddir = base_dir.joinpath(MODEL)

    # set model specific parameters
    if MODEL == 'ERA-Interim':
        # regular expression pattern for finding files
        regex_pattern = (
            r'ERA\-Interim\-Hourly\-MSL\-'
            r'({0})-(\d{{2}})-(\d{{2}})\.nc$'
        )
        # mean sea level pressure file
        input_mean_file = 'ERA-Interim-Mean-MSL-{0:4d}-{1:4d}.nc'
        # input land-sea mask for ocean redistribution
        input_mask_file = 'ERA-Interim-Invariant-Parameters.nc'
        # output file format
        output_file_format = 'ERA-Interim-Hourly-IB-{0}-{1}-{2}.nc'
        VARNAME = 'msl'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'valid_time'
        IBNAME = 'ib'
        UNITS = 'm'
        # hours since 1900-01-01 00:00:0.0
        TIME_LONGNAME = 'Time'
        # land-sea mask variable name and value of oceanic points
        MASKNAME = 'lsm'
        OCEAN = 0
    elif MODEL == 'ERA5':
        # regular expression pattern for finding files
        regex_pattern = r'ERA5\-Hourly\-MSL\-({0})-(\d{{2}})-(\d{{2}})\.nc$'
        # mean sea level pressure file
        input_mean_file = 'ERA5-Mean-MSL-{0:4d}-{1:4d}.nc'
        # input land-sea mask for ocean redistribution
        input_mask_file = 'ERA5-Invariant-Parameters.nc'
        # output file format
        output_file_format = 'ERA5-Hourly-IB-{0}-{1}-{2}.nc'
        VARNAME = 'msl'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'time'
        IBNAME = 'ib'
        UNITS = 'm'
        # hours since 1900-01-01 00:00:0.0
        TIME_LONGNAME = 'Time'
        # land-sea mask variable name and value of oceanic points
        MASKNAME = 'lsm'
        OCEAN = 0
    elif MODEL == 'MERRA-2':
        # regular expression pattern for finding files
        regex_pattern = (
            r'MERRA2_(\d{{3}}).tavg1_2d_slv_Nx.'
            r'({0})(\d{{2}})(\d{{2}}).(.*?).nc$'
        )
        # mean sea level pressure file
        input_mean_file = 'MERRA2.Mean_SLP.{0:4d}-{1:4d}.nc'
        # input land-sea mask for ocean redistribution
        input_mask_file = 'MERRA2_101.const_2d_asm_Nx.00000000.nc4'
        # output file format
        output_file_format = 'MERRA2_{0}.tavg1_2d_IB.{1}{2}{3}.{4}.nc'
        VARNAME = 'SLP'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'
        IBNAME = 'IB'
        UNITS = 'm'
        # minutes since start of file
        TIME_LONGNAME = 'Time'
        # land-sea mask variable name and value of oceanic points
        MASKNAME = 'FROCEAN'
        OCEAN = 1

    # dictionary defining output structure
    struct = dict(
        dimensions=(TIMENAME, LATNAME, LONNAME),
        variables={
            IBNAME: (TIMENAME, LATNAME, LONNAME),
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
    )
    attributes[IBNAME] = dict(
        long_name='Instantaneous_inverse_barometer_(IB)_response',
        units=UNITS,
        density=DENSITY,
    )

    # read mean pressure field
    mean_file = ddir.joinpath(input_mean_file.format(*RANGE))
    mean_pressure, lon, lat = ncdf_mean_pressure(
        mean_file, VARNAME, LONNAME, LATNAME
    )
    # shape of mean pressure field
    ny, nx = np.shape(mean_pressure)

    # grid step size in radians
    dphi = np.radians(np.abs(lon[1] - lon[0]))
    dth = np.radians(np.abs(lat[1] - lat[0]))
    # calculate meshgrid from latitude and longitude
    gridlon, gridlat = np.meshgrid(lon, lat)
    gridphi = np.radians(gridlon)
    # calculate colatitude
    gridtheta = np.radians(90.0 - gridlat)

    # get reference parameters for ellipsoid
    ellipsoid_params = mdlhmc.datum(ellipsoid='WGS84')
    # semimajor and semiminor axes of the ellipsoid [m]
    a_axis = ellipsoid_params.a_axis
    # first numerical eccentricity
    ecc1 = ellipsoid_params.ecc1
    e12 = ecc1**2.0
    # convert from geodetic latitude to geocentric latitude
    # radius of curvature in prime vertical direction (east-west)
    N = a_axis / np.sqrt(1.0 - e12 * np.cos(gridtheta) ** 2.0)
    # radius of curvature in meridional direction (north-south)
    M = a_axis * (1.0 - e12) / (1.0 - e12 * np.cos(gridtheta) ** 2) ** 1.5
    # calculate area of each grid cell
    AREA = (M * dth) * (N * np.sin(gridtheta) * dphi)

    # read land-sea mask to find ocean values
    # ocean pressure points will be based on reanalysis mask
    MASK = ncdf_landmask(ddir.joinpath(input_mask_file), MASKNAME, OCEAN)
    # calculate total area of reanalysis ocean
    # ocean pressure points will be based on reanalysis mask
    ii, jj = np.nonzero(MASK)
    ocean_area = np.sum(AREA[ii, jj])

    # gravitational acceleration at mean sea level at the equator
    ge = 9.780356
    # gravitational acceleration at mean sea level over colatitudes
    # from Heiskanen and Moritz, Physical Geodesy, (1967)
    gs = ge * (
        1.0
        + 5.2885e-3 * np.cos(gridtheta) ** 2
        - 5.9e-6 * np.cos(2.0 * gridtheta) ** 2
    )

    # read each reanalysis pressure field for each year
    regex_years = r'\d{4}' if (YEAR is None) else '|'.join(map(str, YEAR))
    rx = re.compile(regex_pattern.format(regex_years), re.VERBOSE)
    input_files = sorted([f for f in ddir.iterdir() if rx.match(f.name)])
    # for each reanalysis file
    for i, input_file in enumerate(input_files):
        # extract parameters from filename
        if MODEL in ('MERRA-2'):
            # extract date from hourly files
            MOD, YEAR, MONTH, DAY, AUX = rx.findall(input_file.name).pop()
            # output inverse barometer filename
            FILENAME = output_file_format.format(MOD, YEAR, MONTH, DAY, AUX)
            output_file = ddir.joinpath(FILENAME)
        elif MODEL in ('ERA-Interim', 'ERA5'):
            # extract date from hourly files
            YEAR, MONTH, DAY = rx.findall(input_file.name).pop()
            # output inverse barometer filename
            FILENAME = output_file_format.format(YEAR, MONTH, DAY)
            output_file = ddir.joinpath(FILENAME)
        # read netCDF4 mean sea level file
        with netCDF4.Dataset(input_file, mode='r') as fileID:
            # number of time points in file
            (nt,) = fileID.variables[TIMENAME].shape
            # extract time and time units
            dinput = {}
            dinput[TIMENAME] = np.copy(fileID.variables[TIMENAME][:])
            TIME_UNITS = fileID.variables[TIMENAME].units
            attributes[TIMENAME]['units'] = TIME_UNITS
            # copy latitude and longitude
            dinput[LONNAME] = lon.copy()
            dinput[LATNAME] = lat.copy()
            # invalid value
            fill_value = fileID.variables[VARNAME]._FillValue
            # reduced sea level pressure field
            SLP = np.ma.zeros((nt, ny, nx), fill_value=fill_value)
            # calculate reduced sea level pressure for each time
            for t, dt in enumerate(dinput[TIMENAME]):
                # check dimensions for expver slice
                if fileID.variables[VARNAME].ndim == 4:
                    _, nexp, _, _ = fileID.variables[VARNAME].shape
                    # sea level pressure for time
                    pressure = fileID.variables[VARNAME][t, :, :, :].copy()
                    # iterate over expver slices to find valid outputs
                    for j in range(nexp):
                        # check if any are valid for expver
                        if np.any(pressure[j, :, :]):
                            # remove average with respect to time
                            AveRmvd = pressure[j, :, :] - mean_pressure
                            break
                else:
                    # sea level pressure for time
                    pressure = fileID.variables[VARNAME][t, :, :].copy()
                    # remove average with respect to time
                    AveRmvd = pressure - mean_pressure
                # calculate average oceanic pressure values
                AVERAGE = np.sum(AveRmvd[ii, jj] * AREA[ii, jj]) / ocean_area
                # calculate sea level pressure anomalies
                SLP[t, :, :] = AveRmvd - AVERAGE
                # clear temp variables for iteration to free up memory
                pressure, AveRmvd = (None, None)
            # calculate inverse barometer response
            dinput[IBNAME] = -SLP * (DENSITY * gs) ** -1
            # replace masks
            dinput[IBNAME].data[dinput[IBNAME].mask] = fill_value
            # write structured data to netCDF4 file
            mdlhmc.spatial.to_netCDF4(output_file, dinput, attributes, struct)
            # change permissions mode
            output_file.chmod(mode=MODE)


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads hourly mean sea level pressure
            fields from reanalysis and calculates the
            inverse-barometer response
            """,
        fromfile_prefix_chars='@',
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    choices = ['ERA-Interim', 'ERA5', 'MERRA-2']
    parser.add_argument(
        'model',
        type=str,
        nargs='+',
        default=['ERA5', 'MERRA-2'],
        choices=choices,
        help='Reanalysis Model',
    )
    # directory with reanalysis data
    parser.add_argument(
        '--directory',
        '-D',
        type=pathlib.Path,
        default=pathlib.Path.cwd(),
        help='Working data directory',
    )
    # years to run
    now = datetime.datetime.now()
    parser.add_argument(
        '--year',
        '-Y',
        type=int,
        nargs='+',
        default=range(2000, now.year + 1),
        help='Years of model outputs to run',
    )
    # start and end years to run for mean
    parser.add_argument(
        '--mean',
        '-m',
        metavar=('START', 'END'),
        type=int,
        nargs=2,
        default=[2000, 2020],
        help='Start and end year range for mean',
    )
    # ocean fluidic density [kg/m^3]
    parser.add_argument(
        '--density',
        '-d',
        metavar='RHO',
        type=float,
        default=1030.0,
        help='Density of seawater in kg/m^3',
    )
    # verbosity settings
    # verbose will output information about each output file
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
        reanalysis_inverse_barometer(
            args.directory,
            MODEL,
            YEAR=args.year,
            RANGE=args.mean,
            DENSITY=args.density,
            MODE=args.mode,
        )


# run main program
if __name__ == '__main__':
    main()
