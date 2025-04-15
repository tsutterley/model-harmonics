#!/usr/bin/env python
u"""
era5_land_mask_arctic.py
Written by Tyler Sutterley (05/2023)

Creates a mask for ERA5-Land data for Greenland, Svalbard, Iceland and
    the Russian High Arctic defined by a set of shapefiles

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
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

UPDATE HISTORY:
    Written 04/2025
"""
from __future__ import print_function

import sys
import time
import pyproj
import logging
import netCDF4
import pathlib
import argparse
import warnings
import numpy as np
import model_harmonics as mdlhmc

# attempt imports
try:
    import fiona
except ModuleNotFoundError:
    warnings.filterwarnings("module")
    warnings.warn("fiona not available")
    warnings.warn("Some functions will throw an exception if called")
try:
    import shapely.geometry
except ModuleNotFoundError:
    warnings.filterwarnings("module")
    warnings.warn("shapely not available")
    warnings.warn("Some functions will throw an exception if called")
# ignore warnings
warnings.filterwarnings("ignore")

# PURPOSE: read shapefile to find points within a specified region
def read_shapefile(input_shapefile, AREA=None, BUFFER=None):
    # reading shapefile
    input_shapefile = pathlib.Path(input_shapefile).expanduser().absolute()
    logging.debug(str(input_shapefile))
    shape = fiona.open(str(input_shapefile))
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
    return (shapely.geometry.MultiPolygon(poly_list), crs)

# PURPOSE: create a mask for Greenland, Svalbard and Iceland
def era5_land_mask_arctic(base_dir,
        SHAPEFILES=None,
        AREA=None,
        BUFFER=None,
        MODE=0o775
    ):

    # directory models
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    # ERA5-land products
    MODEL = 'ERA5-Land'
    ddir = base_dir.joinpath(MODEL)

    # input land mask and output arctic mask
    input_file = ddir.joinpath(f'lsm.nc')
    output_file = ddir.joinpath(f'carctic.nc')

    # python dictionary with input data
    dinput = {}
    with netCDF4.Dataset(input_file, 'r') as fileID:
        # find valid points from land mask
        mask_input = (fileID.variables['lsm'][0,:,:] > 0.0)
        ntime, nlat, nlon = fileID.variables['lsm'].shape
        # read the latitude and longitude
        dinput['latitude'] = fileID.variables['latitude'][:]
        dinput['longitude'] = fileID.variables['longitude'][:]
        dinput['time'] = fileID.variables['time'][:]

    # create meshgrid of lat and long
    gridlon,gridlat = np.meshgrid(dinput['longitude'], dinput['latitude'])
    gridlon[gridlon > 180.0] -= 360.0
    # latitude range for valid points
    latmin, latmax = (26.0, 86.0)
    # find valid northern hemisphere points from mask input
    ii,jj = np.nonzero(mask_input & (gridlat >= latmin) & (gridlat <= latmax))
    # projection object for converting from latitude/longitude
    crs1 = pyproj.CRS.from_epsg(4326)

    # sparse intersection array
    count = np.count_nonzero(mask_input & (gridlat >= latmin) & (gridlat <= latmax))
    intersection_mask = np.zeros((count),dtype=np.uint8)
    # iterate over shapefiles
    for i,SHAPEFILE in enumerate(SHAPEFILES):
        # read shapefile to find points within region
        poly_obj,crs2 = read_shapefile(SHAPEFILE, AREA=AREA, BUFFER=BUFFER)
        # pyproj transformer for converting from latitude/longitude
        # to projection of input shapefile
        transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
        # convert projection from latitude/longitude to output
        X,Y = transformer.transform(gridlon[ii,jj], gridlat[ii,jj])
        # shapely multipoint object for points
        xy_point = shapely.geometry.MultiPoint(np.c_[X,Y])
        # testing for intersection of points and polygon
        int_test = poly_obj.intersects(xy_point)
        # if there is an intersection
        if int_test:
            # extract intersected points
            int_map = list(map(poly_obj.intersects, xy_point.geoms))
            int_indices, = np.nonzero(int_map)
            intersection_mask[int_indices] = i+1
    # create larger data mask
    dinput['mask'] = np.zeros((ntime, nlat, nlon), dtype=np.uint8)
    dinput['mask'][0,ii,jj] = intersection_mask[:]
    # write to output netCDF4 (.nc)
    ncdf_mask_write(dinput, FILENAME=output_file)
    # change the permission level to MODE
    output_file.chmod(mode=MODE)

# PURPOSE: write land sea mask to netCDF4 file
def ncdf_mask_write(output_data, FILENAME=None):
    # opening NetCDF file for writing
    FILENAME = pathlib.Path(FILENAME).expanduser().absolute()
    fileID = netCDF4.Dataset(FILENAME, 'w', format="NETCDF4")

    # python dictionary with the NetCDF4 data variables
    nc = {}
    # Defining the NetCDF4 dimensions
    TIMENAME,LATNAME,LONNAME = ('time','latitude','longitude')
    for key in [TIMENAME,LONNAME,LATNAME]:
        fileID.createDimension(key, len(output_data[key]))
        nc[key] = fileID.createVariable(key,output_data[key].dtype,(key,))
    # create the NetCDF4 data variables
    nc['mask'] = fileID.createVariable('mask', output_data['mask'].dtype,
        (TIMENAME,LATNAME,LONNAME,), fill_value=0, zlib=True)
    # filling NetCDF variables
    for key,val in output_data.items():
        nc[key][:] = np.copy(val)

    # Defining attributes
    nc[TIMENAME].long_name = 'time'
    nc[TIMENAME].units = 'hours since 1900-01-01 00:00:00.0'
    nc[LONNAME].long_name = 'longitude'
    nc[LONNAME].units = 'degrees_east'
    nc[LATNAME].long_name = 'latitude'
    nc[LATNAME].units = 'degrees_north'
    nc['mask'].long_name = 'land_sea_mask'

    # add software information
    fileID.software_reference = mdlhmc.version.project_name
    fileID.software_version = mdlhmc.version.full_version
    fileID.reference = f'Output from {pathlib.Path(sys.argv[0]).name}'
    # date created
    fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())

    # Output NetCDF structure information
    logging.info(str(FILENAME))
    logging.info(list(fileID.variables.keys()))

    # Closing the NetCDF file
    fileID.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Creates a mask for ERA5-Land data for
            Greenland, Svalbard, Iceland and the Russian
            High Arctic defined by a set of shapefiles
            """
    )
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # input shapefiles to run
    parser.add_argument('--shapefile','-F',
        type=pathlib.Path,
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
        help='permissions mode of output files')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # create logger
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # run program
    era5_land_mask_arctic(args.directory,
        SHAPEFILES=args.shapefile, AREA=args.area,
        BUFFER=args.buffer,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
