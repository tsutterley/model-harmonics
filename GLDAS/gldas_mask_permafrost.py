#!/usr/bin/env python
u"""
gldas_mask_permafrost.py
Written by Tyler Sutterley (12/2022)

Creates a mask for GLDAS data based on the permafrost/surface classification
from the NSIDC Circum-Arctic Map of Permafrost and Ground-Ice Conditions
    1: Continuous Permafrost
    2: Discontinuous Permafrost
    3: Isolated Permafrost
    4: Sporadic Permafrost
    5: Glaciated Area

Projection: Modification of Lambert Azimuthal projection, EASE-Grid
    http://nsidc.org/data/docs/fgdc/ggd318_map_circumarctic/index.html
    http://nsidc.org/data/ggd318.html

GLDAS land mask:
    https://ldas.gsfc.nasa.gov/gldas/vegetation-class-mask

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -S X, --spacing X: spatial resolution of models to run
        10: 1.0 degrees latitude/longitude
        025: 0.25 degrees latitude/longitude
    -F X, --shapefile X: Shapefile to run
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

REFERENCES:
    Brown, J., O.J. Ferrians, Jr., J.A. Heginbottom, and E.S. Melnikov.
        1998, revised February 2001. Circum-arctic map of permafrost and
        ground ice conditions. Boulder, CO: National Snow and Ice Data Center.

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 07/2022: place some imports behind try/except statements
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 01/2021: fiona for shapefile read. pyproj for coordinate conversion
    Updated 02/2019: shapely updates for python3 compatibility
    Updated 08/2018: use getopt to set parameters. output a spatial grid
    Written 07/2013
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
import model_harmonics as mdlhmc

# attempt imports
try:
    import fiona
except ModuleNotFoundError:
    warnings.filterwarnings("always")
    warnings.warn("fiona not available")
    warnings.warn("Some functions will throw an exception if called")
try:
    import shapely.geometry
except ModuleNotFoundError:
    warnings.filterwarnings("always")
    warnings.warn("shapely not available")
    warnings.warn("Some functions will throw an exception if called")
# ignore warnings
warnings.filterwarnings("ignore")

# Read the NSIDC Circum-Arctic Map of Permafrost and Ground-Ice Conditions
# and create a mask for continuous/discontinuous permafrost
def gldas_mask_permafrost(ddir, SPACING=None, SHAPEFILE=None, MODE=0o775):

    # parameters for each grid spacing
    if (SPACING == '025'):
        nx,ny = (1440,600)
        dx,dy = (0.25,0.25)
        latlimit_south = -59.875
        longlimit_west = -179.875
    elif (SPACING == '10'):
        nx,ny = (360,150)
        dx,dy = (1.0,1.0)
        latlimit_south = -59.5
        longlimit_west = -179.5
    # input binary land mask and output netCDF4 mask
    input_file = f'landmask_mod44w_{SPACING}.1gd4r'
    output_file = f'permafrost_mod44w_{SPACING}.nc'

    # python dictionary with input data
    dinput = {}
    # latitude and longitude
    dinput['longitude'] = longlimit_west + np.arange(nx)*dx
    dinput['latitude'] = latlimit_south + np.arange(ny)*dy
    # read GLDAS mask binary file
    binary_input = np.fromfile(os.path.join(ddir,input_file),'>f4')
    mask_input = binary_input.reshape(ny,nx).astype(bool)
    # create meshgrid of lat and long
    gridlon,gridlat = np.meshgrid(dinput['longitude'],dinput['latitude'])
    # find valid northern hemisphere points from mask input
    ii,jj = np.nonzero(mask_input & (gridlat >= 26) & (gridlat <= 86))

    # reading shapefile
    shape = fiona.open(SHAPEFILE)
    # pyproj transformer for converting from latitude/longitude
    # into NSIDC EASE-Grid North
    crs1 = pyproj.CRS.from_epsg(4326)
    crs2 = pyproj.CRS.from_dict(shape.crs)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    # convert projection from latitude/longitude to polar stereographic
    X,Y = transformer.transform(gridlon[ii,jj], gridlat[ii,jj])
    # shapely multipoint object for points
    xy_point = shapely.geometry.MultiPoint(np.c_[X,Y])

    # sparse intersection array
    valid_count=np.count_nonzero(mask_input & (gridlat >= 26) & (gridlat <= 86))
    intersection_mask=np.zeros((valid_count),dtype=np.uint8)
    # iterate over shapefile entities of interest
    attribute_keys = ['NUM_CODE','EXTENT','EXTENT','EXTENT','EXTENT']
    # C: continuous, D: discontinuous, I: isolated, S: sporadic, 21: glaciers
    for j,val in enumerate(['21','S','I','D','C']):
        # find entities
        key = attribute_keys[j]
        entities = [v for v in shape.values() if (v['properties'][key] == val)]
        # for each entity of interest
        for ent in entities:
            # extract coordinates for entity
            poly_list = []
            for coords in ent['geometry']['coordinates']:
                # convert points to latitude/longitude
                x,y = np.transpose(coords)
                poly_list.append(np.c_[x,y])
            # try intersecting polygon with input points
            try:
                # create shapely polygon
                poly_obj=shapely.geometry.Polygon(poly_list[0],poly_list[1:])
            except:
                continue
            else:
                # testing for intersection of points and polygon
                int_test = poly_obj.intersects(xy_point)
                # if there is an intersection
                if int_test:
                    # extract intersected points
                    int_map = list(map(poly_obj.intersects,xy_point))
                    int_indices, = np.nonzero(int_map)
                    intersection_mask[int_indices] = (5-j)
    # fill larger data mask
    dinput['mask'] = np.zeros((ny,nx),dtype=np.uint8)
    dinput['mask'][ii,jj] = intersection_mask[:]
    # write to output netCDF4 (.nc)
    ncdf_mask_write(dinput, FILENAME=os.path.join(ddir,output_file))
    # change the permission level to MODE
    os.chmod(os.path.join(ddir,output_file), MODE)

# PURPOSE: write permafrost mask to netCDF4 file
def ncdf_mask_write(output_data, FILENAME=None):
    # opening NetCDF file for writing
    fileID = netCDF4.Dataset(FILENAME, 'w', format="NETCDF4")

    # python dictionary with the NetCDF4 data variables
    nc = {}
    # Defining the NetCDF4 dimensions
    LATNAME,LONNAME = ('latitude','longitude')
    for key in [LONNAME,LATNAME]:
        fileID.createDimension(key, len(output_data[key]))
        nc[key] = fileID.createVariable(key,output_data[key].dtype,(key,))
    # create the NetCDF4 data variables
    nc['mask'] = fileID.createVariable('mask', output_data['mask'].dtype,
        (LATNAME,LONNAME,), fill_value=0, zlib=True)
    # filling NetCDF variables
    for key,val in output_data.items():
        nc[key][:] = np.copy(val)

    # Defining attributes
    nc[LONNAME].long_name = 'longitude'
    nc[LONNAME].units = 'degrees_east'
    nc[LATNAME].long_name = 'latitude'
    nc[LATNAME].units = 'degrees_north'
    nc['mask'].long_name = 'permafrost_mask'
    nc['mask'].reference = 'https://nsidc.org/data/ggd318'
    description = []
    description.append('1: Continuous Permafrost')
    description.append('2: Discontinuous Permafrost')
    description.append('3: Isolated Permafrost')
    description.append('4: Sporadic Permafrost')
    description.append('5: Glaciated Area')
    nc['mask'].description = ', '.join(description)

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
        description="""Creates a mask for GLDAS data based on the
            permafrost/surface classification from the NSIDC
            Circum-Arctic Map of Permafrost and Ground-Ice Conditions
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
    # input shapefile to run
    parser.add_argument('--shapefile','-F',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Shapefile to run')
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
    gldas_mask_permafrost(args.directory, SPACING=args.spacing,
        SHAPEFILE=args.shapefile, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
