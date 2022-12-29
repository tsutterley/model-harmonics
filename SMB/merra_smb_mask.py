#!/usr/bin/env python
u"""
merra_smb_mask.py
Written by Tyler Sutterley (12/2022)

Creates a mask for MERRA-2 land ice data using a set of shapefiles
https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/
    M2C0NXASM.5.12.4/1980/MERRA2_101.const_2d_asm_Nx.00000000.nc4

INPUTS:
    path to MERRA-2 invariant file
    path to output mask file

COMMAND LINE OPTIONS:
    -v X, --variable X: Variable from input netCDF4 file to extract
    -F X, --shapefile X: Shapefiles to run
    -A X, --area X: Minimum area threshold for polygons
    -B X, --buffer X: Distance to buffer polygons
    -V, --verbose: Output information for each output file
    -M X, --mode X: Permission mode of directories and files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    fiona: Python wrapper for vector data access functions from the OGR library
        https://fiona.readthedocs.io/en/latest/manual.html
    shapely: PostGIS-ish operations outside a database context for Python
        http://toblerity.org/shapely/index.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        https://www.h5py.org/

PROGRAM DEPENDENCIES:
    spatial.py: spatial data class for reading, writing and processing data

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 07/2022: place some imports behind try/except statements
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 02/2021: use spatial class to read input mask file
        use fiona to read from shapefiles. convert to projection of shapefile
    Updated 01/2020: Using pyproj for coordinate conversion
    Updated 02/2019: shapely updates for python3 compatibility
    Updated 06/2018: using python3 compatible octal and input
    Updated 03/2019: using bivariate splines instead of basemap interp
    Written 09/2018
"""
from __future__ import print_function

import sys
import os
import time
import pyproj
import logging
import netCDF4
import argparse
import warnings
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc

# attempt imports
try:
    import fiona
except ModuleNotFoundError:
    warnings.filterwarnings("always")
    warnings.warn("fiona not available")
    warnings.warn("Some functions will throw an exception if called")
try:
    import shapely.ops
    import shapely.geometry
except ModuleNotFoundError:
    warnings.filterwarnings("always")
    warnings.warn("shapely not available")
    warnings.warn("Some functions will throw an exception if called")
# ignore warnings
warnings.filterwarnings("ignore")

# PURPOSE: read shapefile to find points within a specified region
def read_shapefile(input_shapefile, AREA=None, BUFFER=None):
    # reading shapefile
    shape = fiona.open(input_shapefile)
    # create projection object from shapefile
    crs = pyproj.CRS.from_string(shape.crs['init'])
    # list of polygons
    poly_list = []
    # for each entity
    for ent in shape.values():
        # extract coordinates for entity
        for coords in ent['geometry']['coordinates']:
            # extract shapefile points
            x,y = np.transpose(coords)
            # create shapely polygon
            poly_obj = shapely.geometry.Polygon(np.c_[x,y])
            # cannot have overlapping exterior or interior rings
            poly_obj = poly_obj.buffer(BUFFER)
            # add to list if area is above threshold value
            # and resultant polygon is valid
            if poly_obj.is_valid and (poly_obj.area > AREA):
                poly_list.append(poly_obj)
    # return the shapely multipolygon object and the projection
    mpoly_obj = shapely.geometry.MultiPolygon(poly_list)
    return (shapely.ops.unary_union(mpoly_obj), crs)

# PURPOSE: create a mask for MERRA-2 surface mass balance
def merra_smb_mask(input_file, output_file, VARNAME=None,
    SHAPEFILES=None, AREA=None, BUFFER=None, VERBOSE=False,
    MODE=0o775):

    # create logger for verbosity level
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[VERBOSE])

    # create output directory if non-existent
    ddir = os.path.dirname(output_file)
    os.makedirs(ddir) if not os.access(ddir, os.F_OK) else None

    # read input mask file
    dinput = gravtk.spatial().from_netCDF4(input_file,
        lonname='lon', latname='lat', varname=VARNAME,
        verbose=VERBOSE)
    # remove singleton dimensions
    dinput.squeeze()
    # update mask and replace fill value
    dinput.mask = np.where(dinput.data > 0.9, True, False)
    dinput.replace_invalid(1.0)
    # grid spacing
    dlon,dlat = dinput.spacing
    # sign to convert from center to patch
    lon_sign = np.array([-0.5,0.5,0.5,-0.5,-0.5])
    lat_sign = np.array([-0.5,-0.5,0.5,0.5,-0.5])

    # find valid points from mask input
    ii,jj = np.nonzero(~dinput.mask)
    valid_count = np.count_nonzero(~dinput.mask)
    intersection_mask = np.zeros((valid_count),dtype=np.uint8)
    # create meshgrid of lat and long
    gridlon,gridlat = np.meshgrid(dinput.lon, dinput.lat)
    # projection object for converting from latitude/longitude
    crs1 = pyproj.CRS.from_epsg(4326)

    # dictionary with output variables
    output = {}
    # copy geolocation variables from input file
    for key in ['lat','lon']:
        output[key] = getattr(dinput,key)
    # reshape intersection mask to output
    output['mask'] = np.zeros_like(dinput.mask,dtype=np.uint8)
    # iterate over shapefiles
    for i,SHAPEFILE in enumerate(SHAPEFILES):
        # read shapefile to find points within region
        poly_obj,crs2 = read_shapefile(SHAPEFILE, AREA=AREA, BUFFER=BUFFER)
        logging.info(f'Polygon Count: {len(poly_obj):d}')
        # pyproj transformer for converting from latitude/longitude
        # to projection of input shapefile
        transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
        # convert from points to polygon patch
        for indy,indx in zip(ii,jj):
            patchlon = gridlon[indy,indx] + dlon*lon_sign
            patchlat = gridlat[indy,indx] + dlat*lat_sign
            np.clip(patchlon, -180.0, 180.0, out=patchlon)
            np.clip(patchlat, -90.0, 90.0, out=patchlat)
            # convert projection from latitude/longitude to output
            X,Y = transformer.transform(patchlon, patchlat)
            # shapely polygon object for patch
            # cannot have overlapping exterior or interior rings
            patch = shapely.geometry.Polygon(np.c_[X,Y]).buffer(0.0)
            # testing for intersection of points and polygon
            if poly_obj.intersects(patch):
                # extract intersected points
                # if area of intersection is greater than 15%
                int_area = poly_obj.intersection(patch).area
                output['mask'][indy,indx] = (int_area/patch.area) > 0.15
    # write to output netCDF4 (.nc)
    ncdf_mask_write(output, FILENAME=output_file)
    # change the permission level to MODE
    os.chmod(output_file, MODE)

# PURPOSE: write land sea mask to netCDF4 file
def ncdf_mask_write(dinput, FILENAME=None):
    # opening NetCDF file for writing
    fileID = netCDF4.Dataset(FILENAME, 'w', format="NETCDF4")

    # Defining the NetCDF dimensions
    LATNAME,LONNAME = ('lat','lon')
    for key in [LONNAME,LATNAME]:
        fileID.createDimension(key, len(dinput[key]))

    # defining the NetCDF variables
    nc = {}
    nc[LATNAME]=fileID.createVariable(LATNAME,dinput[LATNAME].dtype,(LATNAME,))
    nc[LONNAME]=fileID.createVariable(LONNAME,dinput[LONNAME].dtype,(LONNAME,))
    nc['mask'] = fileID.createVariable('mask', dinput['mask'].dtype,
        (LATNAME,LONNAME,), fill_value=0, zlib=True)
    # filling NetCDF variables
    for key,val in dinput.items():
        nc[key][:] = np.copy(val)

    # Defining attributes for longitude and latitude
    nc[LONNAME].long_name = 'longitude'
    nc[LONNAME].units = 'degrees_east'
    nc[LATNAME].long_name = 'latitude'
    nc[LATNAME].units = 'degrees_north'
    nc['mask'].long_name = 'land_sea_mask'

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
        description="""Creates a mask for MERRA-2 land ice data
            using a set of shapefiles
            """
    )
    # command line parameters
    # input and output file
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='Input file')
    parser.add_argument('outfile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='Output file')
    # variable from input file to extract as mask
    parser.add_argument('--variable','-v',
        type=str, default='FROCEAN',
        help='Variable from input netCDF4 file to extract')
    # input shapefiles to run
    parser.add_argument('--shapefile','-F',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        nargs='+', help='Shapefiles to run')
    # minimum area threshold for polygons within shapefiles
    parser.add_argument('--area','-A',
        type=float, default=0.0,
        help='Minimum area threshold for polygons')
    # distance to buffer polygons within shapefiles
    parser.add_argument('--buffer','-B',
        type=float, default=0.0,
        help='Distance to buffer polygons')
    # verbosity settings
    # verbose will output information about each output file
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

    # run program
    merra_smb_mask(args.infile, args.outfile, VARNAME=args.variable,
        SHAPEFILES=args.shapefile, AREA=args.area, BUFFER=args.buffer,
        VERBOSE=args.verbose, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
