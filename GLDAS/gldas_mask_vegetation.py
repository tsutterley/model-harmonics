#!/usr/bin/env python
u"""
gldas_mask_vegetation.py
Written by Tyler Sutterley (12/2022)

Creates a mask for GLDAS data using the GLDAS vegetation type binary files
    https://ldas.gsfc.nasa.gov/gldas/GLDASvegetation.php
    https://ldas.gsfc.nasa.gov/gldas/vegetation-class-mask

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -S X, --spacing X: spatial resolution of models to run
        10: 1.0 degrees latitude/longitude
        025: 0.25 degrees latitude/longitude
    -V, --verbose: Output information for each output file
    -M X, --mode X: Permission mode of directories and files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 01/2021: using argparse to set parameters
    Updated 06/2018: using python3 compatible octal and input
    Written 03/2018
"""
from __future__ import print_function

import sys
import os
import time
import logging
import netCDF4
import argparse
import numpy as np
import model_harmonics as mdlhmc

# Read the GLDAS vegetation index and create a mask defining each type
def gldas_mask_vegetation(ddir, SPACING=None, MODE=0o775):

    # parameters for each grid spacing
    if (SPACING == '025'):
        dx, dy = (0.25, 0.25)
        nx, ny = (1440, 600)
        latlimit_south = -59.875
        longlimit_west = -179.875
        input_file = f'modmodis_domveg20_{dx:4.2f}.bin'
        output_file = f'modmodis_domveg20_{SPACING}.nc'
    elif (SPACING == '10'):
        dx, dy = (1.0, 1.0)
        nx, ny = (360, 150)
        latlimit_south = -59.5
        longlimit_west = -179.5
        input_file = f'modmodis_domveg20_{dx:3.1f}.bin'
        output_file = f'modmodis_domveg20_{SPACING}.nc'

    # python dictionary with input data
    dinput = {}
    # latitude and longitude
    dinput['longitude'] = longlimit_west + np.arange(nx)*dx
    dinput['latitude'] = latlimit_south + np.arange(ny)*dy
    # read MODIS vegetation index binary file
    mask_input = np.fromfile(os.path.join(ddir,input_file),'>f4')
    dinput['index'] = np.zeros((ny,nx),dtype=np.uint16)
    dinput['index'][:,:] = mask_input.reshape(ny,nx)
    # write to output netCDF4 (.nc)
    ncdf_index_write(dinput, FILENAME=os.path.join(ddir,output_file))
    # change the permission level to MODE
    os.chmod(os.path.join(ddir,output_file),MODE)

# PURPOSE: write vegetation index data to netCDF4 file
def ncdf_index_write(dinput, FILENAME=None):
    # opening NetCDF file for writing
    fileID = netCDF4.Dataset(FILENAME, 'w', format="NETCDF4")

    # Defining the NetCDF dimensions
    LATNAME,LONNAME = ('latitude','longitude')
    for key in [LONNAME,LATNAME]:
        fileID.createDimension(key, len(dinput[key]))

    # defining the NetCDF variables
    nc = {}
    nc[LATNAME]=fileID.createVariable(LATNAME,dinput[LATNAME].dtype,(LATNAME,))
    nc[LONNAME]=fileID.createVariable(LONNAME,dinput[LONNAME].dtype,(LONNAME,))
    nc['index'] = fileID.createVariable('index', dinput['index'].dtype,
        (LATNAME,LONNAME,), fill_value=0, zlib=True)
    # filling NetCDF variables
    for key,val in dinput.items():
        nc[key][:] = np.copy(val)

    # Defining attributes for longitude and latitude
    nc[LONNAME].long_name = 'longitude'
    nc[LONNAME].units = 'degrees_east'
    nc[LATNAME].long_name = 'latitude'
    nc[LATNAME].units = 'degrees_north'
    nc['index'].long_name = 'vegetation_index'
    description = []
    description.append('1: Evergreen Needleleaf Forest')
    description.append('2: Evergreen Broadleaf Forest')
    description.append('3: Deciduous Needleleaf Forest')
    description.append('4: Deciduous Broadleaf Forest')
    description.append('5: Mixed Forest')
    description.append('6: Closed Shrublands')
    description.append('7: Open Shrublands')
    description.append('8: Woody Savannas')
    description.append('9: Savannas')
    description.append('10: Grassland')
    description.append('11: Permanent Wetland')
    description.append('12: Cropland')
    description.append('13: Urban and Built-Up')
    description.append('14: Cropland/Natural Vegetation Mosaic')
    description.append('15: Snow and Ice')
    description.append('16: Barren or Sparsely Vegetated')
    description.append('17: Ocean')
    description.append('18: Wooded Tundra')
    description.append('19: Mixed Tundra')
    description.append('20: Bare Ground Tundra')
    nc['index'].description = ', '.join(description)

    # add software information
    fileID.software_reference = mdlhmc.version.project_name
    fileID.software_version = mdlhmc.version.full_version
    fileID.reference = f'Output from {os.path.basename(sys.argv[0])}'
    # date created
    fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())

    # Output NetCDF structure information
    logging.info(os.path.basename(FILENAME))
    logging.info(list(fileID.variables.keys()))

    # Closing the NetCDF file
    fileID.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Creates a mask for GLDAS data using
            the GLDAS vegetation type binary files
            """
    )
    # command line parameters
    # working data directory for location of GLDAS data
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # model spatial resolution
    # 10: 1.0 degrees latitude/longitude
    # 025: 0.25 degrees latitude/longitude
    parser.add_argument('--spacing','-S',
        type=str, default='10', choices=['10','025'],
        help='Spatial resolution of models to run')
    # verbosity settings
    # verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output files')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # create logger
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # run program
    gldas_mask_vegetation(args.directory, SPACING=args.spacing,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
