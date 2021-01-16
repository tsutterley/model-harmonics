#!/usr/bin/env python
u"""
gldas_mask_arctic.py
Written by Tyler Sutterley (01/2021)

Creates a mask for GLDAS data for Greenland, Svalbard, Iceland and the
    Russian High Arctic defined by a set of shapefiles
    https://ldas.gsfc.nasa.gov/gldas/vegetation-class-mask

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -S X, --spacing X: spatial resolution of models to run
        10: 1.0 degrees latitude/longitude
        025: 0.25 degrees latitude/longitude
    -F X, --shapefile X: Shapefiles to run
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

UPDATE HISTORY:
    Updated 01/2021: fiona for shapefile read. pyproj for coordinate conversion
    Updated 02/2019: shapely updates for python3 compatibility
    Updated 06/2018: using python3 compatible octal and input
    Updated 05/2018: include outputs of Svalbard and Iceland.
    Written 03/2018
"""
from __future__ import print_function

import os
import fiona
import pyproj
import netCDF4
import argparse
import numpy as np
import shapely.geometry

#-- PURPOSE: read coarse shapefile to find points within a specified region
def read_shapefile(input_shapefile):
    #-- reading shapefile
    shape = fiona.open(input_shapefile)
    crs1 = pyproj.CRS.from_string(shape.crs['init'])
    crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(3413))
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    #-- list of polygons
    poly_list = []
    #-- for each entity
    for ent in shape.values():
        #-- extract coordinates for entity
        for coords in ent['geometry']['coordinates']:
            #-- convert points to latitude/longitude
            x,y = np.transpose(coords)
            #-- convert points to polar stereographic
            xi,yi = transformer.transform(x, y)
            #-- create shapely polygon
            poly_obj = shapely.geometry.Polygon(np.c_[xi,yi])
            #-- cannot have overlapping exterior or interior rings
            if (not poly_obj.is_valid):
                poly_obj = poly_obj.buffer(0)
            poly_list.append(poly_obj)
    #-- return the shapely multipolygon object
    return shapely.geometry.MultiPolygon(poly_list)

#-- PURPOSE: create a mask for Greenland, Svalbard and Iceland
def arctic_gldas_landmask(ddir, SPACING=None, SHAPEFILES=None, VERBOSE=False,
    MODE=0o775):
    #-- parameters for each grid spacing
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
    #-- input binary land mask and output netCDF4 mask
    input_file = 'landmask_mod44w_{0}.1gd4r'.format(SPACING)
    output_file = 'arcticmask_mod44w_{0}.nc'.format(SPACING)

    #-- python dictionary with input data
    dinput = {}
    #-- latitude and longitude
    dinput['longitude'] = longlimit_west + np.arange(nx)*dx
    dinput['latitude'] = latlimit_south + np.arange(ny)*dy
    #-- read GLDAS mask binary file
    binary_input = np.fromfile(os.path.join(ddir,input_file),'>f4')
    mask_input = binary_input.reshape(ny,nx).astype(np.bool)
    #-- find valid points from mask input
    ii,jj = np.nonzero(mask_input)
    valid_count = np.count_nonzero(mask_input)
    intersection_mask = np.zeros((valid_count),dtype=np.uint8)
    #-- create meshgrid of lat and long
    gridlon,gridlat = np.meshgrid(dinput['longitude'],dinput['latitude'])
    #-- pyproj transformer for converting from latitude/longitude
    #-- into polar stereographic north
    crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(4326))
    crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(3413))
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    #-- convert projection from latitude/longitude to polar stereographic
    X,Y = transformer.transform(gridlon[ii,jj], gridlat[ii,jj])
    #-- shapely multipoint object for points
    xy_point = shapely.geometry.MultiPoint(np.c_[X,Y])

    #-- iterate over Greenland, Svalbard and Iceland
    for i,SHAPEFILE in enumerate(SHAPEFILES):
        #-- read shapefile to find points within region
        poly_obj = read_shapefile(SHAPEFILE)
        #-- testing for intersection of points and polygon
        int_test = poly_obj.intersects(xy_point)
        #-- if there is an intersection
        if int_test:
            #-- extract intersected points
            int_map = list(map(poly_obj.intersects,xy_point))
            int_indices, = np.nonzero(int_map)
            intersection_mask[int_indices] = i+1
    #-- create larger data mask
    dinput['mask'] = np.zeros((ny,nx),dtype=np.uint8)
    dinput['mask'][ii,jj] = intersection_mask[:]
    #-- write to output netCDF4 (.nc)
    ncdf_mask_write(dinput, FILENAME=os.path.join(ddir,output_file),
        VERBOSE=VERBOSE)
    #-- change the permission level to MODE
    os.chmod(os.path.join(ddir,output_file),MODE)

#-- PURPOSE: write land sea mask to netCDF4 file
def ncdf_mask_write(dinput, FILENAME=None, VERBOSE=False):
    #-- opening NetCDF file for writing
    fileID = netCDF4.Dataset(FILENAME, 'w', format="NETCDF4")

    #-- Defining the NetCDF dimensions
    LATNAME,LONNAME = ('latitude','longitude')
    for key in [LONNAME,LATNAME]:
        fileID.createDimension(key, len(dinput[key]))

    #-- defining the NetCDF variables
    nc = {}
    nc[LATNAME]=fileID.createVariable(LATNAME,dinput[LATNAME].dtype,(LATNAME,))
    nc[LONNAME]=fileID.createVariable(LONNAME,dinput[LONNAME].dtype,(LONNAME,))
    nc['mask'] = fileID.createVariable('mask', dinput['mask'].dtype,
        (LATNAME,LONNAME,), fill_value=0, zlib=True)
    #-- filling NetCDF variables
    for key,val in dinput.items():
        nc[key][:] = np.copy(val)

    #-- Defining attributes for longitude and latitude
    nc[LONNAME].long_name = 'longitude'
    nc[LONNAME].units = 'degrees_east'
    nc[LATNAME].long_name = 'latitude'
    nc[LATNAME].units = 'degrees_north'
    nc['mask'].long_name = 'land_sea_mask'

    #-- Output NetCDF structure information
    if VERBOSE:
        print(os.path.basename(FILENAME))
        print(list(fileID.variables.keys()))

    #-- Closing the NetCDF file
    fileID.close()

#-- Main program that calls gldas_mask_arctic()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Creates a mask for GLDAS data for
            Greenland, Svalbard, Iceland and the Russian
            High Arctic defined by a set of shapefiles
            """
    )
    #-- command line parameters
    #-- working data directory for location of GLDAS data
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- model spatial resolution
    #-- 10: 1.0 degrees latitude/longitude
    #-- 025: 0.25 degrees latitude/longitude
    parser.add_argument('--spacing','-S',
        type=str, default='10', choices=['10','025'],
        help='Spatial resolution of models to run')
    #-- input shapefiles to run
    parser.add_argument('--shapefile','-F',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        nargs='+', help='Shapefiles to run')
    #-- verbosity settings
    #-- verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the local files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output files')
    args = parser.parse_args()

    #-- run program
    arctic_gldas_landmask(args.directory,SPACING=args.spacing,
        SHAPEFILES=args.shapefile,VERBOSE=args.verbose,MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
