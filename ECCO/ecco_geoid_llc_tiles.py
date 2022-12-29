#!/usr/bin/env python
u"""
ecco_geoid_llc_tiles.py
Written by Tyler Sutterley (12/2022)

Calculates geoid heights for ECCO ocean model LLC tiles using model
    coefficients from the GFZ International Centre for Global Earth
    Models (ICGEM)

ECCO Version 4, Revision 4
    https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/nctiles_monthly

ECCO Version 5, Alpha
    https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/nctiles_monthly

INPUTS:
    path to ECCO LLC tile grid file
    path to output geoid LLC tile file

COMMAND LINE OPTIONS:
    -G X, --geoid X: gfc file from the GFZ ICGEM
    -l X, --lmax X: maximum spherical harmonic degree
    -M X, --mode X: Permission mode of directories and files
    -V, --verbose: Output information for each output file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html

PROGRAM DEPENDENCIES:
    geoid_undulation.py: geoidal undulation at a given latitude and longitude
    read_ICGEM_harmonics.py: reads the coefficients for a given gravity model file
    calculate_tidal_offset.py: calculates the C20 offset for a tidal system
    real_potential.py: real potential at a latitude and height for gravity model
    norm_potential.py: normal potential of an ellipsoid at a latitude and height
    norm_gravity.py: normal gravity of an ellipsoid at a latitude and height
    ref_ellipsoid.py: Computes parameters for a reference ellipsoid
    gauss_weights.py: Computes Gaussian weights as a function of degree

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Written 02/2021
"""
from __future__ import print_function

import sys
import os
import logging
import netCDF4
import datetime
import argparse
import numpy as np
import geoid_toolkit as geoidtk
import model_harmonics as mdlhmc

# PURPOSE: read ECCO tiled ocean bottom pressure data and calculate mean
def ecco_geoid_llc_tiles(input_file, output_file, GEOID=None,
    LMAX=None, MODE=0o775):

    # input variable names for each model
    LONNAME = 'XC'
    LATNAME = 'YC'
    ZNAME = 'Depth'
    MASKNAME = 'maskC'
    # read ECCO tile grid file
    invariant = ncdf_invariant(os.path.expanduser(input_file),
        i='i',j='j',tile='tile',lon=LONNAME,lat=LATNAME,
        depth=ZNAME,mask=MASKNAME)
    Nt,Nj,Ni = np.shape(invariant['depth'])
    invalid_mask = np.logical_not(invariant['mask'][0,:,:,:]) | \
        (invariant['depth'] == 0.0)
    # bad value
    fill_value = -1e+10

    # read gravity model spherical harmonics
    Ylms = geoidtk.read_ICGEM_harmonics(os.path.expanduser(GEOID), LMAX=LMAX)
    # extract parameters
    R = np.float64(Ylms['radius'])
    GM = np.float64(Ylms['earth_gravity_constant'])
    LMAX = np.int64(Ylms['max_degree']) if not LMAX else LMAX
    # dictionary with output variables
    output = {}
    # copy from invariant
    for key in ('i','j','tile','lon','lat'):
        output[key] = np.copy(invariant[key])
    # allocate for output geoid height
    output['geoid'] = np.ma.zeros((Nt,Nj,Ni),fill_value=fill_value)
    output['geoid'].mask = np.ones((Nt,Nj,Ni),dtype=bool)
    # calculate geoid for each tile
    for k in range(Nt):
        # find valid points for tile
        indj,indi = np.nonzero(np.logical_not(invalid_mask[k,:,:]))
        # latitude and longitude for tile
        lat = invariant['lat'][k,indj,indi]
        lon = invariant['lon'][k,indj,indi]
        # geoid height for valid points in tile
        output['geoid'].data[k,indj,indi] = geoidtk.geoid_undulation(
            lat, lon, 'WGS84', Ylms['clm'], Ylms['slm'], LMAX, R, GM,
            GAUSS=0)
        output['geoid'].mask[k,indj,indi] = False
    # replace invalid data with fill value
    output['geoid'].data[output['geoid'].mask] = output['geoid'].fill_value

    # Defining output attributes
    attributes = {}
    attributes['title'] = os.path.basename(GEOID)
    attributes['radius'] = Ylms['radius']
    attributes['earth_gravity_constant'] = Ylms['earth_gravity_constant']
    attributes['max_degree'] = str(LMAX)
    # dimension attributes
    attributes['i'] = {}
    attributes['i']['long_name'] = 'x-dimension of the t grid'
    attributes['i']['axis'] = 'X'
    attributes['j'] = {}
    attributes['j']['long_name'] = 'y-dimension of the t grid'
    attributes['j']['axis'] = 'Y'
    attributes['tile'] = {}
    attributes['tile']['long_name'] = 'index of llc grid tile'
    # longitude and latitude
    attributes['lon'] = {}
    attributes['lon']['long_name'] = 'longitude'
    attributes['lon']['units'] = 'degrees_east'
    attributes['lat'] = {}
    attributes['lat']['long_name'] = 'latitude'
    attributes['lat']['units'] = 'degrees_north'
    # output geoid height
    attributes['geoid'] = {}
    attributes['geoid']['long_name'] = 'geoidal_undulation'
    attributes['geoid']['units'] = 'meter'
    attributes['geoid']['reference'] = \
        'https://doi.org/10.5194/essd-11-647-2019'

    # netcdf (.nc)
    ncdf_tile_write(output, attributes, FILENAME=output_file,
        LONNAME='lon', LATNAME='lat', VARNAME='geoid')
    # change the permissions mode of the output file to MODE
    os.chmod(output_file,MODE)

# PURPOSE: read ECCO invariant grid file
def ncdf_invariant(invariant_file,**kwargs):
    # output dictionary with invariant parameters
    invariant = {}
    # open netCDF4 file for reading
    with netCDF4.Dataset(os.path.expanduser(invariant_file),'r') as fileID:
        # extract latitude, longitude, depth, area and valid mask
        for key,val in kwargs.items():
            invariant[key] = fileID.variables[val][:].copy()
    # return the invariant parameters
    return invariant

# PURPOSE: write tiled data to a netCDF4 flie
def ncdf_tile_write(output, attributes, FILENAME=None, LONNAME=None,
    LATNAME=None, VARNAME=None):

    # opening NetCDF file for writing
    fileID = netCDF4.Dataset(os.path.expanduser(FILENAME),'w')

    # python dictionary with NetCDF variables
    nc = {}
    # Defining the NetCDF dimensions and variables
    for key in ('i','j','tile'):
        fileID.createDimension(key, len(np.atleast_1d(output[key])))
        nc[key] = fileID.createVariable(key, output[key].dtype, (key,))
        # filling NetCDF variables
        nc[key][:] = np.copy(output[key])
        # Defining attributes for variable
        for att_name,att_val in attributes[key].items():
            setattr(nc[key],att_name,att_val)

    # Defining the NetCDF variables
    for key in (LONNAME,LATNAME,VARNAME):
        if hasattr(output[key],'fill_value'):
            nc[key] = fileID.createVariable(key, output[key].dtype,
                ('tile','j','i'), fill_value=output[key].fill_value,
                zlib=True)
        else:
            nc[key] = fileID.createVariable(key, output[key].dtype,
                ('tile','j','i'))
        # filling NetCDF variables
        nc[key][:] = np.copy(output[key])
        # Defining attributes for variable
        for att_name,att_val in attributes[key].items():
            setattr(nc[key],att_name,att_val)
    # add attribute for date created
    fileID.date_created = datetime.datetime.now().isoformat()
    fileID.title = attributes['title']
    fileID.radius = attributes['radius']
    fileID.earth_gravity_constant = attributes['earth_gravity_constant']
    fileID.max_degree = attributes['max_degree']
    # add software information
    fileID.software_reference = mdlhmc.version.project_name
    fileID.software_version = mdlhmc.version.full_version
    fileID.reference = f'Output from {os.path.basename(sys.argv[0])}'
    # Output NetCDF structure information
    logging.info(FILENAME)
    logging.info(list(fileID.variables.keys()))
    # Closing the NetCDF file
    fileID.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates geoid heights for ECCO ocean
            model LLC tiles
            """
    )
    # command line options
    # input and output file
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='Input file')
    parser.add_argument('outfile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='Output file')
    # path to static gravity harmonics file
    parser.add_argument('--geoid','-G',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='gfc file from the GFZ ICGEM')
    # maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=None,
        help='Maximum spherical harmonic degree')
    # print information about each output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files')
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
    ecco_geoid_llc_tiles(args.infile, args.outfile, GEOID=args.geoid,
        LMAX=args.lmax, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
