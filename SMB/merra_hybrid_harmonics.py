#!/usr/bin/env python
u"""
merra_hybrid_harmonics.py
Written by Tyler Sutterley (04/2024)
Read MERRA-2 hybrid variables and converts to spherical harmonics
MERRA-2 Hybrid firn model outputs provided by Brooke Medley at GSFC

CALLING SEQUENCE:
    python merra_hybrid_harmonics.py --directory <path> --region gris \
        --version v1.1 --product SMB_a --lmax 60 --verbose

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -R X, --region X: Region to calculate (gris, ais)
    -v X, --version X: Version of firn model to calculate
        v0
        v1
        v1.0
        v1.1
        v1.2
        v1.2.1
    -P X, --product X: MERRA-2 hybrid product to convert to harmonics
    -Y X, --year X: Years to run
    --mask X: netCDF4 mask files for reducing to regions
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -n X, --love X: Treatment of the Love Love numbers
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
        3: Wang et al. (2012) values from PREM with hard sediment
        4: Wang et al. (2012) values from PREM with soft sediment
    --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -F X, --format X: Output data format
        ascii
        netcdf
        HDF5
    -G, --gzip: input netCDF4 file is gzip compressed
    -V, --verbose: Output information for each output file
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    utilities.py: download and management utilities for files
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    datum.py: calculate reference parameters for common ellipsoids
    gen_point_load.py: calculates spherical harmonics from point masses
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: spatial data class for reading, writing and processing data

UPDATE HISTORY:
    Updated 04/2024: changed polar stereographic area function to scale_factors
    Updated 04/2023: added option to convert firn air content variables
    Updated 03/2023: add root attributes to output netCDF4 and HDF5 files
        use spatial function for calculating geocentric latitude
    Updated 02/2023: use love numbers class with additional attributes
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 10/2022: move polar stereographic scaling function to spatial
        add Greenland and Antarctic versions v1.2.1
    Updated 06/2022: change default variables to include firn height anomaly
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: use wrapper function for reading load Love numbers
    Updated 12/2021: open MERRA-2 hybrid product command line options
        added GSFC MERRA-2 Hybrid Greenland v1.2
        can use variable loglevels for verbose output
    Updated 10/2021: add pole case in stereographic area scale calculation
        using python logging for handling verbose output
    Updated 09/2021: use original FDM file for ais products
    Written 08/2021
"""
from __future__ import print_function

import sys
import copy
import gzip
import uuid
import pyproj
import logging
import netCDF4
import pathlib
import argparse
import datetime
import warnings
import numpy as np
import scipy.interpolate
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc
# ignore pyproj and divide by zero warnings
warnings.filterwarnings("ignore")

# PURPOSE: set the projection parameters based on the region name
def set_projection(REGION):
    if (REGION == 'ais'):
        projection_flag = 'EPSG:3031'
    elif (REGION == 'gris'):
        projection_flag = 'EPSG:3413'
    return projection_flag

# PURPOSE: read MERRA-2 hybrid SMB estimates and convert to spherical harmonics
def merra_hybrid_harmonics(base_dir, REGION, VARIABLE, YEARS,
    VERSION='v1',
    MASKS=None,
    LMAX=None,
    MMAX=None,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    DATAFORM=None,
    GZIP=False,
    MODE=0o775):

    # MERRA-2 hybrid directory
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    DIRECTORY = base_dir.joinpath(VERSION)

    # suffix if compressed
    suffix = '.gz' if GZIP else ''
    # set the input netCDF4 file for the variable of interest
    if VARIABLE in ('cum_smb_anomaly',):
        FILE_VERSION = copy.copy(VERSION)
        args = (VERSION,REGION.lower(),suffix)
        hybrid_file = 'gsfc_fdm_{0}_{1}.nc{2}'.format(*args)
    elif VARIABLE in ('FAC','height','h_a') or (REGION.lower() == 'ais'):
        FILE_VERSION = VERSION.replace('.','_')
        args = (FILE_VERSION,REGION.lower(),suffix)
        hybrid_file = 'gsfc_fdm_{0}_{1}.nc{2}'.format(*args)
    elif VARIABLE in ('Me_a','Ra_a','Ru_a','Sn-Ev_a','SMB_a'):
        FILE_VERSION = VERSION.replace('.','_')
        args = (FILE_VERSION,REGION.lower(),suffix)
        hybrid_file = 'gsfc_fdm_smb_cumul_{0}_{1}.nc{2}'.format(*args)
    else:
        raise ValueError(f'Unknown variable {VARIABLE}')

    # Open the MERRA-2 Hybrid NetCDF file for reading
    input_file = DIRECTORY.joinpath(hybrid_file)
    if GZIP:
        # read as in-memory (diskless) netCDF4 dataset
        with gzip.open(str(input_file), mode='r') as f:
            fileID = netCDF4.Dataset(uuid.uuid4().hex, memory=f.read())
    else:
        # read netCDF4 dataset
        fileID = netCDF4.Dataset(input_file, mode='r')

    # Output NetCDF file information
    logging.info(str(input_file))
    logging.info(list(fileID.variables.keys()))

    # Get data from each netCDF variable and remove singleton dimensions
    fd = {}
    # time is year decimal at time step 5 days
    time_step = 5.0/365.25
    # reduce grids to time period of input buffered by time steps
    tmin = np.min(YEARS) - 2.0*time_step
    tmax = np.max(YEARS) + 1.0 + 2.0*time_step
    # find indices to times
    nt, = fileID.variables['time'].shape
    f = scipy.interpolate.interp1d(fileID.variables['time'][:],
        np.arange(nt), kind='nearest', bounds_error=False,
        fill_value=(0,nt))
    imin,imax = f((tmin,tmax)).astype(np.int64)
    # read reduced time variables
    fd['time'] = fileID.variables['time'][imin:imax+1].copy()
    # read reduced dataset and remove singleton dimensions
    fd[VARIABLE] = np.squeeze(fileID.variables[VARIABLE][imin:imax+1,:,:])
    # invalid data value
    fv = np.float64(fileID.variables[VARIABLE]._FillValue)
    # input variable units
    variable_units = fileID.variables[VARIABLE].units
    # input shape of MERRA-2 Hybrid firn data
    nt,nx,ny = np.shape(fd[VARIABLE])
    # extract x and y coordinate arrays from grids if applicable
    # else create meshgrids of coordinate arrays
    if (np.ndim(fileID.variables['x'][:]) == 2):
        xg = fileID.variables['x'][:].copy()
        yg = fileID.variables['y'][:].copy()
        fd['x'],fd['y'] = (xg[:,0],yg[0,:])
    else:
        fd['x'] = fileID.variables['x'][:].copy()
        fd['y'] = fileID.variables['y'][:].copy()
        xg,yg = np.meshgrid(fd['x'],fd['y'],indexing='ij')
    # extract area of each grid cell if applicable
    # calculate using dimensions if not possible
    fd['area'] = np.zeros((nx,ny))
    try:
        # ice covered area
        fd['area'][:,:] = fileID.variables['iArea'][:,:].copy()
    except:
        # calculate grid areas (assume fully ice covered)
        dx = np.abs(fd['x'][1] - fd['x'][0])
        dy = np.abs(fd['y'][1] - fd['y'][0])
        fd['area'][:,:] = dx*dy
    # close the NetCDF files
    fileID.close()

    # create mask object for reducing data
    if not MASKS:
        fd['mask'] = np.ones((nx,ny),dtype=bool)
    else:
        fd['mask'] = np.zeros((nx,ny),dtype=bool)
    # read masks for reducing regions before converting to harmonics
    for mask_file in MASKS:
        logging.info(mask_file)
        fileID = netCDF4.Dataset(mask_file,'r')
        fd['mask'] |= fileID.variables['mask'][:].astype(bool)
        fileID.close()
    # indices of valid MERRA hybrid data
    fd['mask'] &= (fd[VARIABLE].data[0,:,:] != fv)

    # pyproj transformer for converting to input coordinates (EPSG)
    MODEL_EPSG = set_projection(REGION)
    crs1 = pyproj.CRS.from_string('EPSG:4326')
    crs2 = pyproj.CRS.from_string(MODEL_EPSG)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    direction = pyproj.enums.TransformDirection.INVERSE
    # convert projection from model coordinates
    gridlon,gridlat = transformer.transform(xg, yg, direction=direction)
    # polar stereographic standard parallel (latitude of true scale)
    reference_latitude = crs2.to_dict().pop('lat_ts')

    # get reference parameters for ellipsoid
    ellipsoid_params = mdlhmc.datum(ellipsoid='WGS84')
    # semimajor axis of ellipsoid [m]
    a_axis = ellipsoid_params.a_axis
    # ellipsoidal flattening
    flat = ellipsoid_params.flat
    # Average Radius of the Earth with equal surface area [m]
    rad_e = ellipsoid_params.rad_e
    # calculate geocentric latitude and convert to degrees
    latitude_geocentric = mdlhmc.spatial.geocentric_latitude(gridlon, gridlat,
        a_axis=a_axis, flat=flat)

    # reduce latitude and longitude to valid and masked points
    indx,indy = np.nonzero(fd['mask'])
    lon,lat = (gridlon[indx,indy],latitude_geocentric[indx,indy])
    # scaled areas
    ps_scale = mdlhmc.spatial.scale_factors(gridlat[indx,indy], flat=flat,
        reference_latitude=reference_latitude)
    # unit parameters for each input variable type
    if VARIABLE in ('FAC','height','h_a'):
        # areas in terms of solid angle (steradians)
        scaling_factor = ps_scale*fd['area'][indx,indy]/(rad_e**2)
        # use custom UNITS to keep as inputs but use 4-pi norm
        UNITS = np.ones((LMAX+1))/(4.0*np.pi)
        # output spherical harmonic units
        harmonic_units = copy.copy(variable_units)
    else:
        # densities of meteoric ice [kg/m^3]
        rho_ice = 917.0
        # scaling factor to convert inputs from from meters ice eq to g
        scaling_factor = 1000.0*rho_ice*ps_scale*fd['area'][indx,indy]
        # use named point mass units code (grams)
        UNITS = 1
        # output spherical harmonic units
        harmonic_units = 'Geodesy_Normalization'

    # read load love numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE, FORMAT='class')
    # upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    # output string for both LMAX == MMAX and LMAX != MMAX cases
    order_str = 'M{MMAX:d}' if (MMAX != LMAX) else ''

    # allocate for output spherical harmonics
    Ylms = gravtk.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX+1,MMAX+1,nt))
    Ylms.slm = np.zeros((LMAX+1,MMAX+1,nt))
    Ylms.time = np.zeros((nt))
    Ylms.month = np.zeros((nt), dtype=np.int64)
    # for each time step
    for t in range(nt):
        # reduce data for date and scale
        scaled = scaling_factor*fd[VARIABLE][t,indx,indy]
        # convert to spherical harmonics
        YLMS = gravtk.gen_point_load(scaled, lon, lat,
            LMAX=LMAX, MMAX=MMAX, UNITS=UNITS, LOVE=LOVE)
        # copy harmonics for time step
        Ylms.clm[:,:,t] = YLMS.clm[:,:].copy()
        Ylms.slm[:,:,t] = YLMS.slm[:,:].copy()
        # copy date parameters for time step
        Ylms.time[t] = fd['time'][t].copy()
        Ylms.month[t] = gravtk.time.calendar_to_grace(Ylms.time[t])

    # output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')
    # attributes for output files
    attributes = {}
    attributes['institution'] = 'NASA Goddard Space Flight Center (GSFC)'
    attributes['project'] = 'MERRA2-hybrid'
    attributes['product_region'] = REGION
    attributes['product_version'] = VERSION
    attributes['product_name'] = VARIABLE
    attributes['product_type'] = 'gravity_field'
    # add attributes for earth parameters
    attributes['earth_model'] = LOVE.model
    attributes['earth_love_numbers'] = LOVE.citation
    attributes['reference_frame'] = LOVE.reference
    # add attributes for maximum degree and order
    attributes['max_degree'] = LMAX
    attributes['max_order'] = MMAX
    attributes['lineage'] = input_file.name
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'
    # add attributes to output harmonics
    Ylms.attributes['ROOT'] = attributes
    # output spherical harmonic data file
    args = (FILE_VERSION,REGION.lower(),VARIABLE,LMAX,order_str,suffix[DATAFORM])
    FILE = 'gsfc_fdm_{0}_{1}_{2}_CLM_L{3:d}{4}.{5}'.format(*args)
    output_file = DIRECTORY.joinpath(FILE)
    Ylms.to_file(output_file, format=DATAFORM,
        date=True, units=harmonic_units)
    # change the permissions mode of the output file to MODE
    output_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Read MERRA-2 hybrid variables and
            converts to spherical harmonics
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # region of firn model
    parser.add_argument('--region','-R',
        type=str, default='gris', choices=['gris','ais'],
        help='Region of firn model to calculate')
    # version of firn model
    versions = ['v0','v1','v1.0','v1.1','v1.2','v1.2.1']
    parser.add_argument('--version','-v',
        type=str, default='v1.2.1', choices=versions,
        help='Version of firn model to calculate')
    # products from firn model
    parser.add_argument('--product','-P',
        type=str, default='SMB_a',
        help='MERRA-2 hybrid product to calculate')
    # years to run
    now = datetime.datetime.now()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.year+1),
        help='Years of model outputs to run')
    # mask file for reducing to regions
    parser.add_argument('--mask',
        type=pathlib.Path,
        nargs='+', default=[],
        help='netCDF4 masks file for reducing to regions')
    # maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
    # different treatments of the load Love numbers
    # 0: Han and Wahr (1995) values from PREM
    # 1: Gegout (2005) values from PREM
    # 2: Wang et al. (2012) values from PREM
    # 3: Wang et al. (2012) values from PREM with hard sediment
    # 4: Wang et al. (2012) values from PREM with soft sediment
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2,3,4],
        help='Treatment of the Load Love numbers')
    # option for setting reference frame for gravitational load love number
    # reference frame options (CF, CM, CE)
    parser.add_argument('--reference',
        type=str.upper, default='CF', choices=['CF','CM','CE'],
        help='Reference frame for load Love numbers')
    # input and output data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input and output data format')
    # print information about each input and output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
    # netCDF4 files are gzip compressed
    parser.add_argument('--gzip','-G',
        default=False, action='store_true',
        help='netCDF4 file is locally gzip compressed')
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
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # run program
    merra_hybrid_harmonics(args.directory, args.region, args.product, args.year,
        VERSION=args.version,
        MASKS=args.mask,
        LMAX=args.lmax,
        MMAX=args.mmax,
        LOVE_NUMBERS=args.love,
        REFERENCE=args.reference,
        DATAFORM=args.format,
        GZIP=args.gzip,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
