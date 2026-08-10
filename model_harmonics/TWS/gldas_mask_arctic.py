#!/usr/bin/env python
"""
gldas_mask_arctic.py
Written by Tyler Sutterley (07/2026)

Creates a mask for GLDAS data for Greenland, Svalbard, Iceland and the
    Russian High Arctic defined by a set of shapefiles
    https://ldas.gsfc.nasa.gov/gldas/vegetation-class-mask

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -S X, --spacing X: spatial resolution of models to run
        10: 1.0 degrees latitude/longitude
        025: 0.25 degrees latitude/longitude
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
    Updated 07/2026: output using structured dictionary with netCDF4 parameters
    Updated 09/2025: use importlib to attempt to import dependencies
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: use full path to output file in verbose logging
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 07/2022: place some imports behind try/except statements
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 02/2021: convert to projection of each shapefile
        replaced numpy bool to prevent deprecation warning
    Updated 01/2021: use fiona to read from shapefiles
        use pyproj for coordinate conversion
    Updated 02/2019: shapely updates for python3 compatibility
    Updated 06/2018: using python3 compatible octal and input
    Updated 05/2018: include outputs of Svalbard and Iceland.
    Written 03/2018
"""

from __future__ import print_function

import sys
import pyproj
import logging
import pathlib
import argparse
import numpy as np
import model_harmonics as mdlhmc
from gravity_toolkit.utilities import build_logger

# attempt imports
fiona = mdlhmc.utilities.import_dependency('fiona')
shapely = mdlhmc.utilities.import_dependency('shapely')
shapely.ops = mdlhmc.utilities.import_dependency('shapely.ops')
shapely.geometry = mdlhmc.utilities.import_dependency('shapely.geometry')


# PURPOSE: read shapefile to find points within a specified region
def read_shapefile(input_shapefile, AREA=None, BUFFER=None):
    # get logger
    logger = logging.getLogger(__name__)
    # reading shapefile
    logger.debug(str(input_shapefile))
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
            x, y = np.transpose(coords)
            # create shapely polygon
            poly_obj = shapely.geometry.Polygon(np.c_[x, y])
            # cannot have overlapping exterior or interior rings
            poly_obj = poly_obj.buffer(BUFFER)
            # add to list if area is above threshold value
            # and resultant polygon is valid
            if poly_obj.is_valid and (poly_obj.area > AREA):
                poly_list.append(poly_obj)
    # return the shapely multipolygon object and the projection
    return (shapely.geometry.MultiPolygon(poly_list), crs)


# PURPOSE: create a mask for Greenland, Svalbard and Iceland
def gldas_mask_arctic(
    ddir, SPACING=None, SHAPEFILES=None, AREA=None, BUFFER=None, MODE=0o775
):
    # get logger
    logger = logging.getLogger(__name__)
    # parameters for each grid spacing
    if SPACING == '025':
        nx, ny = (1440, 600)
        dx, dy = (0.25, 0.25)
        latlimit_south = -59.875
        longlimit_west = -179.875
    elif SPACING == '10':
        nx, ny = (360, 150)
        dx, dy = (1.0, 1.0)
        latlimit_south = -59.5
        longlimit_west = -179.5
    # input binary land mask and output netCDF4 mask
    ddir = pathlib.Path(ddir).expanduser().absolute()
    input_file = ddir.joinpath(f'landmask_mod44w_{SPACING}.1gd4r')
    output_file = ddir.joinpath(f'arcticmask_mod44w_{SPACING}.nc')

    # python dictionary with input data
    dinput = {}
    # latitude and longitude
    dinput['longitude'] = longlimit_west + np.arange(nx) * dx
    dinput['latitude'] = latlimit_south + np.arange(ny) * dy

    # dictionary describing the output netCDF4 structure
    struct = dict(
        dimensions=('latitude', 'longitude'),
        variables={
            'mask': ('latitude', 'longitude'),
        },
    )
    # list of input files for provenance
    lineage = []
    lineage.append(input_file.name)
    # netCDF4 attributes for output files
    attributes = dict(ROOT={})
    reference = f'Output from {pathlib.Path(sys.argv[0]).name}'
    attributes['ROOT']['reference'] = reference
    attributes['ROOT']['title'] = 'GLDAS Arctic land-sea mask'
    attributes['longitude'] = dict(long_name='longitude', units='degrees_east')
    attributes['latitude'] = dict(long_name='latitude', units='degrees_north')
    attributes['mask'] = dict(long_name='land_sea_mask')

    # read GLDAS mask binary file
    logger.info(str(input_file))
    binary_input = np.fromfile(input_file, '>f4')
    mask_input = binary_input.reshape(ny, nx).astype(bool)
    # find valid points from mask input
    ii, jj = np.nonzero(mask_input)
    valid_count = np.count_nonzero(mask_input)
    intersection_mask = np.zeros((valid_count), dtype=np.uint8)
    # create meshgrid of lat and long
    gridlon, gridlat = np.meshgrid(dinput['longitude'], dinput['latitude'])
    # projection object for converting from latitude/longitude
    crs1 = pyproj.CRS.from_epsg(4326)

    # iterate over shapefiles
    for i, SHAPEFILE in enumerate(SHAPEFILES):
        # read shapefile to find points within region
        SHAPEFILE = pathlib.Path(SHAPEFILE).expanduser().absolute()
        poly_obj, crs2 = read_shapefile(SHAPEFILE, AREA=AREA, BUFFER=BUFFER)
        lineage.append(SHAPEFILE.name)
        # pyproj transformer for converting from latitude/longitude
        # to projection of input shapefile
        transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
        # convert projection from latitude/longitude to output
        X, Y = transformer.transform(gridlon[ii, jj], gridlat[ii, jj])
        # shapely multipoint object for points
        xy_point = shapely.geometry.MultiPoint(np.c_[X, Y])
        # testing for intersection of points and polygon
        int_test = poly_obj.intersects(xy_point)
        # if there is an intersection
        if int_test:
            # extract intersected points
            int_map = list(map(poly_obj.intersects, xy_point.geoms))
            (int_indices,) = np.nonzero(int_map)
            intersection_mask[int_indices] = i + 1
    # create larger data mask
    dinput['mask'] = np.zeros((ny, nx), dtype=np.uint8)
    dinput['mask'][ii, jj] = intersection_mask[:]
    # write to output netCDF4 (.nc)
    attributes['ROOT']['lineage'] = lineage
    mdlhmc.spatial.to_netCDF4(output_file, dinput, attributes, struct)
    # change the permission level to MODE
    output_file.chmod(mode=MODE)


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Creates a mask for GLDAS data for
            Greenland, Svalbard, Iceland and the Russian
            High Arctic defined by a set of shapefiles
            """
    )
    # command line parameters
    # working data directory for location of GLDAS data
    parser.add_argument(
        '--directory',
        '-D',
        type=pathlib.Path,
        default=pathlib.Path.cwd(),
        help='Working data directory',
    )
    # model spatial resolution
    # 10: 1.0 degrees latitude/longitude
    # 025: 0.25 degrees latitude/longitude
    parser.add_argument(
        '--spacing',
        '-S',
        type=str,
        default='10',
        choices=['10', '025'],
        help='Spatial resolution of models to run',
    )
    # input shapefiles to run
    parser.add_argument(
        '--shapefile',
        '-F',
        type=pathlib.Path,
        nargs='+',
        help='Shapefiles to run',
    )
    # minimum area threshold for polygons within shapefiles
    parser.add_argument(
        '--area',
        '-A',
        type=float,
        default=0.0,
        help='Minimum area threshold for polygons',
    )
    # distance to buffer polygons within shapefiles
    parser.add_argument(
        '--buffer',
        '-B',
        type=float,
        default=0.0,
        help='Distance to buffer polygons',
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
        help='permissions mode of output files',
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
    logger = build_logger(__name__, level=loglevels[args.verbose])

    # run program
    gldas_mask_arctic(
        args.directory,
        SPACING=args.spacing,
        SHAPEFILES=args.shapefile,
        AREA=args.area,
        BUFFER=args.buffer,
        MODE=args.mode,
    )


# run main program
if __name__ == '__main__':
    main()
