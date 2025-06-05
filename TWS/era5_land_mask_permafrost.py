#!/usr/bin/env python
u"""
era5_land_mask_permafrost.py
Written by Tyler Sutterley (04/2025)

Creates a mask for ERA5-Land based on the permafrost/surface classification
from the NSIDC Circum-Arctic Map of Permafrost and Ground-Ice Conditions
    1: Continuous Permafrost
    2: Discontinuous Permafrost
    3: Isolated Permafrost
    4: Sporadic Permafrost
    5: Glaciated Area

Projection: Modification of Lambert Azimuthal projection, EASE-Grid
    http://nsidc.org/data/docs/fgdc/ggd318_map_circumarctic/index.html
    http://nsidc.org/data/ggd318.html

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
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

# Read the NSIDC Circum-Arctic Map of Permafrost and Ground-Ice Conditions
# and create a mask for continuous/discontinuous permafrost
def era5_land_mask_permafrost(base_dir,
        SHAPEFILE=None,
        BUFFER=None,
        MODE=0o775
    ):

    # directory models
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    # ERA5-land products
    MODEL = 'ERA5-Land'
    ddir = base_dir.joinpath(MODEL)

    # input land mask and output permafrost mask
    input_file = ddir.joinpath(f'lsm.nc')
    output_file = ddir.joinpath(f'cpfrost.nc')

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
    # latitude range for valid points
    latmin, latmax = (26.0, 86.0)
    # find valid northern hemisphere points from mask input
    ii,jj = np.nonzero(mask_input & (gridlat >= latmin) & (gridlat <= latmax))

    # reading shapefile
    SHAPEFILE = pathlib.Path(SHAPEFILE).expanduser().absolute()
    logging.debug(str(SHAPEFILE))
    shape = fiona.open(str(SHAPEFILE))
    # pyproj transformer for converting from latitude/longitude
    # into NSIDC EASE-Grid North
    crs1 = pyproj.CRS.from_epsg(4326)
    crs2 = pyproj.CRS.from_user_input(shape.crs)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    # convert projection from latitude/longitude to polar stereographic
    X,Y = transformer.transform(gridlon[ii,jj], gridlat[ii,jj])
    # shapely multipoint object for points
    xy_point = shapely.geometry.MultiPoint(np.c_[X,Y])

    # sparse intersection array
    count = np.count_nonzero(mask_input & (gridlat >= latmin) & (gridlat <= latmax))
    intersection_mask = np.zeros((count),dtype=np.uint8)
    # iterate over shapefile entities of interest
    attribute_keys = ['NUM_CODE','EXTENT','EXTENT','EXTENT','EXTENT']
    # C: continuous, D: discontinuous, I: isolated, S: sporadic, 21: glaciers
    for j,val in enumerate(['21','S','I','D','C']):
        # find entities
        key = attribute_keys[j]
        entities = [v for v in shape.values() if (v['properties'][key] == val)]
        # list of polygon objects for permafrost type
        poly_list = []
        # for each entity of interest
        for ent in entities:
            # extract coordinates for entity
            coord_list = []
            for coords in ent['geometry']['coordinates']:
                # convert points to latitude/longitude
                x,y = np.transpose(coords)
                coord_list.append(np.c_[x,y])
            # try creating a polygon object for entity
            try:
                # create shapely polygon
                poly_obj = shapely.geometry.Polygon(coord_list[0],coord_list[1:])
            except:
                continue
            else:
                # buffer polygon and add to list
                poly_list.append(poly_obj.buffer(BUFFER))
        # create shapely multipolygon object using a union
        mpoly_obj = shapely.unary_union(poly_list)
        # testing for intersection of points and multipolygon
        int_test = mpoly_obj.intersects(xy_point)
        # if there is an intersection
        if int_test:
            # extract intersected points
            int_map = list(map(mpoly_obj.intersects, xy_point.geoms))
            int_indices, = np.nonzero(int_map)
            intersection_mask[int_indices] = (5-j)
    # fill output data mask
    dinput['pf'] = np.zeros((ntime, nlat, nlon), dtype=np.uint8)
    dinput['pf'][0,ii,jj] = intersection_mask[:]
    # find valid southern hemisphere points and replace Antarctica
    ii,jj = np.nonzero(mask_input & (gridlat <= -60))
    dinput['pf'][0,ii,jj] = 5
    # write to output netCDF4 (.nc)
    ncdf_mask_write(dinput, FILENAME=output_file)
    # change the permission level to MODE
    output_file.chmod(mode=MODE)

# PURPOSE: write permafrost mask to netCDF4 file
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
    nc['pf'] = fileID.createVariable('pf', output_data['pf'].dtype,
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
    nc['pf'].long_name = 'permafrost_mask'
    nc['pf'].reference = 'https://nsidc.org/data/ggd318'
    description = []
    description.append('1: Continuous Permafrost')
    description.append('2: Discontinuous Permafrost')
    description.append('3: Isolated Permafrost')
    description.append('4: Sporadic Permafrost')
    description.append('5: Glaciated Area')
    nc['pf'].description = ', '.join(description)

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
        description="""Creates a mask for ERA5-Land data based on the
            permafrost/surface classification from the NSIDC
            Circum-Arctic Map of Permafrost and Ground-Ice Conditions
            """
    )
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # input shapefile to run
    parser.add_argument('--shapefile','-F',
        type=pathlib.Path,
        help='Shapefile to run')
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
    era5_land_mask_permafrost(args.directory, 
        SHAPEFILE=args.shapefile,
        BUFFER=args.buffer,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
