#!/usr/bin/env python
u"""
GIA_uplift_ICESat2_ATL15.py
Written by Tyler Sutterley (09/2023)
Calculates GIA-induced crustal uplift over polar stereographic grids for
    correcting ICESat-2 ATL15 gridded land ice height change data
Calculated directly from GIA spherical harmonics

INPUTS:
    ATL15 file to read

COMMAND LINE OPTIONS:
    -h, --help: list the command line options
    -l X, --lmax X: maximum spherical harmonic degree
    -G X, --gia X: GIA model type to read
        IJ05-R2: Ivins R2 GIA Models
        W12a: Whitehouse GIA Models
        SM09: Simpson/Milne GIA Models
        ICE6G: ICE-6G GIA Models
        Wu10: Wu (2010) GIA Correction
        AW13-ICE6G: Geruo A ICE-6G GIA Models
        Caron: Caron JPL GIA Assimilation
        ICE6G-D: ICE-6G Version-D GIA Models
        ascii: reformatted GIA in ascii format
        netCDF4: reformatted GIA in netCDF4 format
        HDF5: reformatted GIA in HDF5 format
    --gia-file X: GIA file to read
    -I, --iterate: Iterations over mask to solve for large matrices
    -V, --verbose: verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    h5py: Pythonic interface to the HDF5 binary data format
        (http://www.h5py.org/)
    netCDF4: Python interface to the netCDF C library
         (https://unidata.github.io/netcdf4-python/netCDF4/index.html)
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/
    future: Compatibility layer between Python 2 and Python 3
        (http://python-future.org/)

PROGRAM DEPENDENCIES:
    read_GIA_model.py: reads harmonics for a glacial isostatic adjustment model
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    clenshaw_summation.py: calculate spatial field from spherical harmonics
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    harmonics.py: data class for working with spherical harmonics
    utilities.py: download and management utilities for files

REFERENCE:
    Holmes and Featherstone, "A Unified Approach to the Clenshaw Summation and
        the Recursive Computation of Very High Degree and Order Normalised
        Associated Legendre Functions", Journal of Geodesy (2002)
        http://dx.doi.org/10.1007/s00190-002-0216-2
    Tscherning and Poder, "Some Geodetic Applications of Clenshaw Summation",
        Bollettino di Geodesia e Scienze (1982)

UPDATE HISTORY:
    Updated 09/2023: simplify input arguments
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 02/2023: added iterate option to reduce memory load for 1km files
        add debug logging for ATL15 netCDF4 variables
        use love numbers class with additional attributes
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose output
    Updated 05/2022: use argparse descriptions within documentation
        use GIA reference and citation output from GIA read program
    Forked 04/2022 from calculate_GIA_uplift_ICESat2.py from Smith et al. (2020)
    Updated 04/2020: updates to reading load love numbers
    Updated 10/2019: set spatial reference using polar stereographic EPSG codes
        using pyproj for coordinate conversion
    Forked 08/2019 from gia_stereographic.py
    Updated 07/2019: added AW13 ICE6G models and option to set data directory
    Updated 06/2018: using python3 compatible octal and input
    Updated 08/2017: use clenshaw summation to calculate at each point
    Written 10/2016
"""
from __future__ import print_function

import sys
import os
import re
import time
import pyproj
import logging
import pathlib
import netCDF4
import argparse
import traceback
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: read a variable group from ICESat-2 ATL15
def read_ATL15(infile, group='delta_h'):
    # dictionary with ATL15 variables
    ATL15 = {}
    attributes = {}
    infile = pathlib.Path(infile).expanduser().absolute()
    with netCDF4.Dataset(infile, mode='r') as fileID:
        # check if reading from root group or sub-group
        ncf = fileID.groups[group] if group else fileID
        # netCDF4 structure information
        logging.debug(str(infile))
        logging.debug(list(ncf.variables.keys()))
        for key,val in ncf.variables.items():
            ATL15[key] = val[:].copy()
            attributes[key] = {}
            for att_name in val.ncattrs():
                attributes[key][att_name] = val.getncattr(att_name)
    # return the data and attributes
    return (ATL15, attributes)

# PURPOSE: calculate spatial fields of GIA crustal uplift to correct
# ICESat-2 ATL15 Gridded Land Ice Height Change data
def calculate_GIA_uplift(input_file, LMAX,
    GIA=None,
    GIA_FILE=None,
    ITERATIONS=1,
    MODE=0o775):
    """
    Calculate spatial fields of GIA crustal uplift to correct
    ICESat-2 ATL15 Gridded Land Ice Height Change data

    Parameters
    ----------
    input_file: str
        full path to ATL15 gridded land ice height change file
    LMAX: int
        Maximum degree and order of spherical harmonics
    GIA: str or NoneType, default None
        GIA model type to read and output

            - ``'IJ05-R2'``: Ivins R2 GIA Models [Ivins2013]_
            - ``'W12a'``: Whitehouse GIA Models [Whitehouse2012]_
            - ``'SM09'``: Simpson/Milne GIA Models [Simpson2009]_
            - ``'ICE6G'``: ICE-6G GIA Models [Peltier2015]_
            - ``'Wu10'``: Wu (2010) GIA Correction [Wu2010]_
            - ``'AW13-ICE6G'``: Geruo A ICE-6G GIA Models [A2013]_
            - ``'Caron'``: Caron JPL GIA Assimilation [Caron2018]_
            - ``'ICE6G-D'``: ICE-6G Version-D GIA Models [Peltier2018]_
            - ``'ascii'``: reformatted GIA in ascii format
            - ``'netCDF4'``: reformatted GIA in netCDF4 format
            - ``'HDF5'``: reformatted GIA in HDF5 format
    GIA_FILE: str or NoneType, default None
        full path to input GIA file
    ITERATE: bool, default False
        Iterate over mask to solve for large matrices
    MODE: oct, default 0o775
        Permissions mode of output files

    References
    ----------
    .. [Blewett2003] G. Blewitt, "Self-consistency in reference frames, geocenter
        definition, and surface loading of the solid Earth",
        *Journal of Geophysical Research: Solid Earth*, 108(B2), 2103, (2003).
        `doi: 10.1029/2002JB002082 <https://doi.org/10.1029/2002JB002082>`_

    .. [Dziewonski1981] A. M. Dziewonski and D. L. Anderson,
        "Preliminary reference Earth model",
        *Physics of the Earth and Planetary Interiors*, 25(4), 297--356, (1981).
        `doi: 10.1016/0031-9201(81)90046-7 <https://doi.org/10.1016/0031-9201(81)90046-7>`_

    .. [Gegout2010] P. Gegout, J. Boehm, and D. Wijaya,
        "Practical numerical computation of love numbers and applications",
        Workshop of the COST Action ES0701, (2010).
        `doi: 10.13140/RG.2.1.1866.7045 <https://doi.org/10.13140/RG.2.1.1866.7045>`_

    .. [Han1995] D. Han and J. Wahr, "The viscoelastic relaxation of a
        realistically stratified earth, and a further analysis of postglacial
        rebound", *Geophysical Journal International*, 120(2), 287--311, (1995).
        `doi: 10.1111/j.1365-246X.1995.tb01819.x <https://doi.org/10.1111/j.1365-246X.1995.tb01819.x>`_

    .. [Wahr1998] J. Wahr, M. Molenaar, and F. Bryan,
        "Time variability of the Earth's gravity field: Hydrological and
        oceanic effects and their possible detection using GRACE",
        *Journal of Geophysical Research*, 103(B12), 30205--30229, (1998).
        `doi: 10.1029/98JB02844 <https://doi.org/10.1029/98JB02844>`_

    .. [Wahr2000] J. Wahr, D. Wingham, and C. Bentley,
        "A method of combining ICESat and GRACE satellite data
        to constrain Antarctic mass balance",
        *Journal of Geophysical Research: Solid Earth*,
        105(B7), 16279--16294, (2000).
        `doi: 10.1029/2000JB900113 <https://doi.org/10.1029/2000JB900113>`_

    .. [Wang2012] H. Wang et al., "Load Love numbers and Green's
        functions for elastic Earth models PREM, iasp91, ak135, and
        modified models with refined crustal structure from Crust 2.0",
        *Computers & Geosciences*, 49, 190--199, (2012).
        `doi: 10.1016/j.cageo.2012.06.022 <https://doi.org/10.1016/j.cageo.2012.06.022>`_

    """
    # verify path to input ATL15 file
    input_file = pathlib.Path(input_file).expanduser().absolute()
    # parse ATL15 file
    pattern = r'(ATL\d{2})_(.*?)_(\d{2})(\d{2})_(.*?)_(\d{3})_(\d{2}).nc$'
    rx = re.compile(pattern, re.VERBOSE)
    PRD, RGN, SCYC, ECYC, RES, RL, VERS = rx.findall(input_file.name).pop()
    # read ATL15 data for group
    ATL15, attrib = read_ATL15(input_file, group='delta_h')
    nt, ny, nx = ATL15['delta_h'].shape
    grid_mapping_name = attrib['delta_h']['grid_mapping']
    crs_wkt = attrib[grid_mapping_name]['crs_wkt']
    fill_value = attrib['delta_h']['_FillValue']
    # coordinate reference systems for converting from projection
    crs1 = pyproj.CRS.from_wkt(crs_wkt)
    crs2 = pyproj.CRS.from_epsg(4326)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    # dictionary of coordinate reference system variables
    cs_to_cf = crs1.cs_to_cf()
    crs_to_dict = crs1.to_dict()

    # read arrays of kl, hl, and ll Love Numbers
    # these Love numbers are not used in the GIA uplift calculation
    # but are a required input for the clenshaw summation
    LOVE = gravtk.load_love_numbers(LMAX,
        LOVE_NUMBERS=0, REFERENCE='CF',
        FORMAT='class')

    # input GIA spherical harmonic datafiles
    GIA_Ylms_rate = gravtk.gia(lmax=LMAX).from_GIA(GIA_FILE, GIA=GIA)

    # convert x and y axes to grid
    gridx, gridy = np.meshgrid(ATL15['x'], ATL15['y'])
    gridlon, latitude_geodetic = transformer.transform(gridx, gridy)

    # WGS84 ellipsoid parameters
    # semimajor axis of the ellipsoid [m]
    a_axis = crs2.ellipsoid.semi_major_metre
    # flattening of the ellipsoid
    flat = crs2.ellipsoid.inverse_flattening**-1
    # first numerical eccentricity
    ecc1 = np.sqrt((2.0*flat - flat**2)*a_axis**2)/a_axis

    # convert from geodetic latitude to geocentric latitude
    # geodetic latitude in radians
    latitude_geodetic_rad = np.pi*latitude_geodetic/180.0
    # prime vertical radius of curvature
    N = a_axis/np.sqrt(1.0 - ecc1**2.*np.sin(latitude_geodetic_rad)**2.)
    # calculate X, Y and Z from geodetic latitude and longitude
    X = N * np.cos(latitude_geodetic_rad) * np.cos(np.pi*gridlon/180.0)
    Y = N * np.cos(latitude_geodetic_rad) * np.sin(np.pi*gridlon/180.0)
    Z = (N * (1.0 - ecc1**2.0)) * np.sin(latitude_geodetic_rad)
    # calculate geocentric latitude and convert to degrees
    latitude_geocentric = 180.0*np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))/np.pi

    # output data and attributes dictionaries
    output_data = {}
    attributes = dict(x={}, y={})
    # x and y
    output_data['x'] = ATL15['x'].copy()
    output_data['y'] = ATL15['y'].copy()
    for att_name in ['long_name','standard_name','units']:
        attributes['x'][att_name] = cs_to_cf[0][att_name]
        attributes['y'][att_name] = cs_to_cf[1][att_name]

    # find valid points (bare-ice, near-shore ocean or ground)
    if (ITERATIONS > 1):
        valid_mask = np.logical_not(ATL15['delta_h'].mask)
    else:
        valid_mask = np.ones_like(ATL15['delta_h'].mask, dtype=bool)
    # indices of valid points
    indy, indx = np.nonzero(np.any(valid_mask, axis=0))
    nind = len(indy)
    # allocate for output GIA uplift
    output_data['dhdt_gia'] = np.empty((ny, nx))
    output_data['dhdt_gia'][:,:] = fill_value
    for i in range(ITERATIONS):
        ind = slice(i, nind, ITERATIONS)
        iY, iX = (indy[ind], indx[ind])
        # Converting GIA rates to crustal uplift (Wahr et al., 2000)
        # convert from cm to meters uplift
        output_data['dhdt_gia'][iY,iX] = 0.01*gravtk.clenshaw_summation(
            GIA_Ylms_rate.clm, GIA_Ylms_rate.slm,
            gridlon[iY,iX], latitude_geocentric[iY,iX],
            LMAX=LMAX, RAD=0, UNITS=6, LOVE=LOVE,
            ASTYPE=np.float64, SCALE=1e-32)

    # GIA correction attributes
    attributes['dhdt_gia'] = {}
    attributes['dhdt_gia']['long_name'] = 'Viscoelastic Crustal Uplift'
    attributes['dhdt_gia']['model'] = str(GIA_Ylms_rate.citation)
    attributes['dhdt_gia']['description'] = str(GIA_Ylms_rate.title)
    attributes['dhdt_gia']['reference'] = str(GIA_Ylms_rate.reference)
    attributes['dhdt_gia']['coordinates'] = 'y x'
    attributes['dhdt_gia']['units'] = 'm/yr'
    attributes['dhdt_gia']['grid_mapping'] = 'crs'
    attributes['dhdt_gia']['_FillValue'] = fill_value

    # output file
    FILE = f'{PRD}_{RGN}_{RES}_GIA_{GIA_Ylms_rate.title}_L{LMAX:d}.nc'
    output_file = input_file.with_name(FILE)
    fileID = netCDF4.Dataset(output_file, mode='w')

    # dictionary with netCDF4 variables
    nc = {}
    # netCDF4 dimension variables
    dimensions = []
    dimensions.append('y')
    dimensions.append('x')
    dims = tuple(dimensions)

    # create projection variable
    nc['crs'] = fileID.createVariable('crs', np.byte, ())
    # add projection attributes
    nc['crs'].setncattr('standard_name', 'Polar_Stereographic')
    # nc['crs'].setncattr('spatial_epsg', crs1.to_epsg())
    nc['crs'].setncattr('spatial_ref', crs1.to_wkt())
    nc['crs'].setncattr('proj4_params', crs1.to_proj4())
    nc['crs'].setncattr('latitude_of_projection_origin', crs_to_dict['lat_0'])
    for att_name, att_val in crs1.to_cf().items():
        nc['crs'].setncattr(att_name, att_val)

    # netCDF4 dimensions
    for i,key in enumerate(dimensions):
        val = output_data[key]
        fileID.createDimension(key, len(val))
        nc[key] = fileID.createVariable(key, val.dtype, (key,))
        # filling netCDF4 dimension variables
        nc[key][:] = val
        # Defining attributes for variable
        for att_name,att_val in attributes[key].items():
            nc[key].setncattr(att_name,att_val)

    # netCDF4 spatial variables
    variables = set(output_data.keys()) - set(dimensions)
    for key in sorted(variables):
        val = output_data[key]
        if '_FillValue' in attributes[key].keys():
            nc[key] = fileID.createVariable(key, val.dtype, dims,
                fill_value=attributes[key]['_FillValue'], zlib=True)
            attributes[key].pop('_FillValue')
        elif val.shape:
            nc[key] = fileID.createVariable(key, val.dtype, dims,
                zlib=True)
        else:
            nc[key] = fileID.createVariable(key, val.dtype, ())
        # filling netCDF4 variables
        nc[key][:] = val
        # Defining attributes for variable
        for att_name,att_val in attributes[key].items():
            nc[key].setncattr(att_name, att_val)

    # add root level attributes
    fileID.setncattr('title', 'ATL15_GIA_Correction')
    fileID.setncattr('summary', 'Glacial_Isostatic_Adjustment_Corrections_'
        'for_NASA_ICESat-2_ATL15_Gridded_Land_Ice_Height_Change_data.')
    fileID.setncattr('GDAL_AREA_OR_POINT', 'Area')
    fileID.setncattr('Conventions', 'CF-1.6')
    today = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    fileID.setncattr('date_created', today)
    # add software information
    fileID.setncattr('software_reference', gravtk.version.project_name)
    fileID.setncattr('software_version', gravtk.version.full_version)
    # add geospatial and temporal attributes
    fileID.setncattr('geospatial_lat_min', latitude_geodetic.min())
    fileID.setncattr('geospatial_lat_max', latitude_geodetic.max())
    fileID.setncattr('geospatial_lon_min', gridlon.min())
    fileID.setncattr('geospatial_lon_max', gridlon.max())
    fileID.setncattr('geospatial_lat_units', "degrees_north")
    fileID.setncattr('geospatial_lon_units', "degrees_east")
    # Output NetCDF structure information
    logging.info(str(output_file))
    logging.info(list(fileID.variables.keys()))
    # Closing the netCDF4 file
    fileID.close()
    # change the permissions mode
    output_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates GIA-induced crustal uplift over polar
            stereographic grids for correcting ICESat-2 ATL15 gridded
            land ice height change data
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = \
        gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('infile',
        type=pathlib.Path,
        help='ICESat-2 ATL15 file to run')
    # maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    # GIA model type list
    models = {}
    models['IJ05-R2'] = 'Ivins R2 GIA Models'
    models['W12a'] = 'Whitehouse GIA Models'
    models['SM09'] = 'Simpson/Milne GIA Models'
    models['ICE6G'] = 'ICE-6G GIA Models'
    models['Wu10'] = 'Wu (2010) GIA Correction'
    models['AW13-ICE6G'] = 'Geruo A ICE-6G GIA Models'
    models['Caron'] = 'Caron JPL GIA Assimilation'
    models['ICE6G-D'] = 'ICE-6G Version-D GIA Models'
    models['ascii'] = 'reformatted GIA in ascii format'
    models['netCDF4'] = 'reformatted GIA in netCDF4 format'
    models['HDF5'] = 'reformatted GIA in HDF5 format'
    # GIA model type
    parser.add_argument('--gia','-G',
        type=str, metavar='GIA', choices=models.keys(),
        help='GIA model type to read')
    # full path to GIA file
    parser.add_argument('--gia-file',
        type=pathlib.Path,
        help='GIA file to read')
    # iterate over mask to solve for large matrices
    parser.add_argument('--iterate','-I',
        type=uint, default=1,
        help='Iterations over mask to solve for large matrices')
    # print information about each input and output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of run')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
    # return the parser
    return parser

def uint(x):
    """Argument parser type for positive integers
    """
    x = int(x)
    if (x <= 0):
        raise argparse.ArgumentTypeError('must be a positive integer')
    return x

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # create logger
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # try to run the analysis with listed parameters
    try:
        info(args)
        # run algorithm with parameters
        calculate_GIA_uplift(args.infile, args.lmax,
            GIA=args.gia,
            GIA_FILE=args.gia_file,
            ITERATIONS=args.iterate,
            MODE=args.mode)
    except Exception as exc:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())

# run main program
if __name__ == '__main__':
    main()
