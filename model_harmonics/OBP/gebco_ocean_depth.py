#!/usr/bin/env python
"""
gebco_ocean_depth.py
Written by Tyler Sutterley (05/2023)

Interpolates GEBCO bathymetry to model grids

GEBCO 2014/2020 Gridded bathymetry data:
https://www.bodc.ac.uk/data/hosted_data_systems/gebco_gridded_bathymetry_data/

INPUTS:
    model_file: GEBCO zipped bathymetry file

COMMAND LINE OPTIONS:
    -D X, --directory X: working data directory
    -S X, --spacing X: spatial resolution of output data (dlon,dlat)
    -M X, --mode X: Permission mode of directories and files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html

UPDATE HISTORY:
    Updated 07/2026: generalized program for any grid spacing
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: lower case keyword arguments to output spatial
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 12/2020: use argparse to set command line parameters
        using spatial module for read and write routines
    Written 10/2018
"""

from __future__ import print_function

import sys
import re
import pathlib
import argparse
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc


# PURPOSE: interpolate GEBCO bathymetry to model grids
def gebco_ocean_depth(ddir, FILE, SPACING=[0.5, 0.5], MODE=0o775):
    # verify output directory
    ddir = pathlib.Path(ddir).expanduser().absolute()
    ddir.mkdir(parents=True, exist_ok=True, mode=MODE)
    # verify input file path
    FILE = pathlib.Path(FILE).expanduser().absolute()

    # input bathymetry model parameters
    rx = re.compile(r'gebco_(\d+)', re.IGNORECASE | re.VERBOSE)
    (VERSION,) = rx.findall(FILE.name)
    # read zipped file and extract file into in-memory file object
    bathymetry = gravtk.spatial().from_netCDF4(
        FILE, date=False, varname='elevation', compression='zip'
    )
    bathymetry.data = extend_matrix(bathymetry.data, 1)
    bathymetry.lon = extend_array(bathymetry.lon, 1)
    # verify validity of bathymetry data and mask invalid values
    bathymetry.mask |= np.isnan(bathymetry.data)
    bathymetry.mask |= bathymetry.data >= 0.0
    bathymetry.update_mask()
    # convert to meshgrids
    gridlon, gridlat = np.meshgrid(bathymetry.lon, bathymetry.lat)
    # find valid bathymetry points
    valid = np.logical_not(bathymetry.mask)
    # calculate grid spacing in radians
    dphi, dth = np.radians(bathymetry.spacing)

    # get reference parameters for WGS84 ellipsoid
    ellipsoid_params = mdlhmc.datum(ellipsoid='WGS84')
    # semimajor axis of the ellipsoid [m]
    a_axis = ellipsoid_params.a_axis
    # square of first numerical eccentricity
    e12 = ellipsoid_params.ecc1**2.0
    # colatitude in radians
    theta = np.radians(90.0 - bathymetry.lat)
    # radius of curvature in prime vertical direction (east-west)
    N = a_axis / np.sqrt(1.0 - e12 * np.cos(theta) ** 2.0)
    # radius of curvature in meridional direction (north-south)
    M = (a_axis * (1.0 - e12)) / np.power(1.0 - e12 * np.cos(theta) ** 2, 1.5)
    # calculate area of grid cells
    area = (M * dth) * (N * np.sin(theta) * dphi)

    # grid spacing
    dx, dy = np.broadcast_to(np.atleast_1d(SPACING), (2,))
    # output dimensions and extents
    xmin = -180 + dx / 2.0
    xmax = 180 - dx / 2.0
    ymin = -90 + dy / 2.0
    ymax = 90 - dy / 2.0
    nlon = int((xmax - xmin) / dx) + 1
    nlat = int((ymax - ymin) / dy) + 1
    extent = [xmin, xmax, ymin, ymax]
    # bad value
    fill_value = 99999.0

    # create output data
    interp = gravtk.spatial(fill_value=fill_value)
    # calculate dimension variables
    interp.lon = np.arange(extent[0], extent[1] + dx, dx)
    interp.lat = np.arange(extent[2], extent[3] + dy, dy)
    interp.data = np.zeros((nlat, nlon))
    interp.mask = np.ones((nlat, nlon), dtype=bool)
    for j, lat in enumerate(interp.lat):
        for i, lon in enumerate(interp.lon):
            # find valid bathymetry points within grid cell
            ilon = np.abs(gridlon - lon) <= (dx / 2.0)
            ilat = np.abs(gridlat - lat) <= (dy / 2.0)
            if not np.any(ilon & ilat & valid):
                continue
            # calculate the area-weighted mean of grid cell bathymetries
            jj, ii = np.nonzero(ilon & ilat & valid)
            total_weighted = np.sum(bathymetry.data[jj, ii] * area[jj, ii])
            total_area = np.sum(area[jj, ii])
            # convert from bathymetry to depth and assign to grid
            interp.data[j, i] = -(total_weighted / total_area)
            interp.mask[j, i] = np.all(bathymetry.mask[jj, ii])
    # update the mask
    interp.update_mask()

    # output file attributes
    attributes = {}
    attributes['varname'] = 'depth'
    attributes['units'] = 'm'
    attributes['longname'] = 'Depth'
    attributes['title'] = (
        'General Depth Chart of the Oceans, produced by the'
        'International Hydrographic Organization (IHO) and the United Nations '
        '(UNESCO) Intergovernmental Oceanographic Commission (IOC)'
    )
    attributes['source'] = (
        'https://www.gebco.net/data_and_products/gridded_bathymetry_data/'
    )
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'
    # output netCDF4 dataset
    bathymetry_file = ddir.joinpath(f'DEPTH.{VERSION}.{nlon:d}x{nlat:d}.nc')
    interp.to_netCDF4(bathymetry_file, date=False, **attributes)
    # change the permissions mode to MODE
    bathymetry_file.chmod(mode=MODE)


# wrapper function to extend a matrix
def extend_matrix(input_matrix, count):
    ny, nx = np.shape(input_matrix)
    temp = np.zeros((ny, nx + 2 * count), dtype=input_matrix.dtype)
    temp[:, 0:count] = input_matrix[:, -count:]
    temp[:, count:-count] = input_matrix[:, :]
    temp[:, -count:] = input_matrix[:, 0:count]
    return temp


# wrapper function to linearly extend an array
def extend_array(input_array, count):
    n = len(input_array)
    step_size = np.abs(input_array[1] - input_array[0])
    temp = np.zeros((n + 2 * count), dtype=input_array.dtype)
    temp[0:count] = input_array[0] - step_size * np.arange(count, 0, -1)
    temp[count:-count] = input_array[:]
    temp[-count:] = input_array[-1] + step_size * np.arange(1, count + 1)
    return temp


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Interpolates GEBCO bathymetry to model grids
            """
    )
    # command line parameters
    parser.add_argument(
        'file',
        type=pathlib.Path,
        help='GEBCO bathymetry file',
    )
    # working data directory
    parser.add_argument(
        '--directory',
        '-D',
        type=pathlib.Path,
        default=mdlhmc.utilities.get_data_path('data'),
        help='Working data directory',
    )
    # output grid parameters
    parser.add_argument(
        '--spacing',
        '-S',
        type=float,
        nargs=2,
        default=[0.5, 0.5],
        metavar=('dlon', 'dlat'),
        help='Spatial resolution of output data',
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

    # run program
    gebco_ocean_depth(
        args.directory,
        args.file,
        SPACING=args.spacing,
        MODE=args.mode,
    )


# run main program
if __name__ == '__main__':
    main()
