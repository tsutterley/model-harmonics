#!/usr/bin/env python
u"""
gemb_smb_harmonics.py
Written by Tyler Sutterley (04/2024)
Read GEMB SMB variables and convert to spherical harmonics
Shifts dates of SMB point masses to mid-month values to correspond with GRACE

CALLING SEQUENCE:
    python gemb_smb_harmonics.py --lmax 60 --verbose <path_to_gemb_file>

COMMAND LINE OPTIONS:
    -P X, --product X: GEMB product to convert to harmonics
    --mask X: netCDF4 mask files for reducing to regions
    --area X: netCDF4 area files for calculating mass
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
        add option for adding area file for calculating mass
    Updated 02/2023: use love numbers class with additional attributes
    Updated 12/2022: single implicit import of spherical harmonic tools
        use constants class in place of geoid-toolkit ref_ellipsoid
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Written 11/2022
"""
from __future__ import print_function

import sys
import re
import copy
import pyproj
import logging
import netCDF4
import pathlib
import argparse
import warnings
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc
# ignore pyproj and divide by zero warnings
warnings.filterwarnings("ignore")

# PURPOSE: set the projection parameters based on the region name
def set_projection(REGION):
    if (REGION == 'Antarctica'):
        projection_flag = 'EPSG:3031'
    elif (REGION == 'Greenland'):
        projection_flag = 'EPSG:3413'
    return projection_flag

# PURPOSE: read GEMB SMB estimates and convert to spherical harmonics
def gemb_smb_harmonics(model_file,
    VARIABLE='accum_SMB',
    MASKS=None,
    AREA=None,
    LMAX=None,
    MMAX=None,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    DATAFORM=None,
    MODE=0o775):

    # regular expression pattern for extracting parameters
    pattern = r'GEMB_(Greenland|Antarctica)_(.*?)_(v.*?).nc$'
    model_file = pathlib.Path(model_file).expanduser().absolute()
    region, _, version = re.findall(pattern, model_file.name).pop()

    # Open the GEMB NetCDF file for reading
    fileID = netCDF4.Dataset(model_file, mode='r')

    # Output NetCDF file information
    logging.info(model_file)
    logging.info(list(fileID.variables.keys()))

    # Get data from each netCDF variable and remove singleton dimensions
    fd = {}
    # read reduced time variables
    fd['time'] = fileID.variables['time'][:].copy()
    # read surface mass balance or firn variables
    fd[VARIABLE] = fileID.variables[VARIABLE][:,:,:].copy()
    # invalid data value
    try:
        fv = np.float64(fileID.variables[VARIABLE]._FillValue)
    except (ValueError,AttributeError):
        fv = np.nan
    # input variable units
    variable_units = fileID.variables[VARIABLE].units
    # input shape of GEMB surface mass balance or firn data
    nt,ny,nx = np.shape(fd[VARIABLE])
    # extract x and y coordinate arrays
    fd['x'] = fileID.variables['x'][:].copy()
    fd['y'] = fileID.variables['y'][:].copy()
    xg,yg = np.meshgrid(fd['x'],fd['y'])
    # close the NetCDF files
    fileID.close()

    # calculate grid areas (read file or assume fully ice covered)
    fd['area'] = np.ma.zeros((ny,nx), fill_value=fv)
    fd['area'].mask = np.zeros((ny,nx), dtype=bool)
    if AREA is None:
        dx = np.abs(fd['x'][1] - fd['x'][0])
        dy = np.abs(fd['y'][1] - fd['y'][0])
        fd['area'].data[:,:] = dx*dy
    else:
        # read area file (km^2)
        logging.info(AREA)
        fileID = netCDF4.Dataset(AREA, 'r')
        fd['area'][:,:] = 1e6*fileID.variables['area'][:]
        fileID.close()

    # create mask object for reducing data
    if not MASKS:
        fd['mask'] = np.ones((ny,nx),dtype=bool)
    else:
        fd['mask'] = np.zeros((ny,nx),dtype=bool)
    # read masks for reducing regions before converting to harmonics
    for mask_file in MASKS:
        logging.info(mask_file)
        fileID = netCDF4.Dataset(mask_file,'r')
        fd['mask'] |= fileID.variables['mask'][:].astype(bool)
        fileID.close()
    # indices of valid GEMB data
    fd['mask'] &= (fd[VARIABLE].data[0,:,:] != fv)
    fd['mask'] &= np.isfinite(fd[VARIABLE].data[0,:,:])
    fd['mask'] &= np.logical_not(fd['area'].mask)

    # pyproj transformer for converting to input coordinates (EPSG)
    MODEL_EPSG = set_projection(region)
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
    # semimajor axis of the ellipsoid [m]
    a_axis = ellipsoid_params.a_axis
    # ellipsoidal flattening
    flat = ellipsoid_params.flat
    # Average Radius of the Earth with equal surface area [m]
    rad_e = ellipsoid_params.rad_e
    # calculate geocentric latitude and convert to degrees
    latitude_geocentric = mdlhmc.spatial.geocentric_latitude(gridlon, gridlat,
        a_axis=a_axis, flat=flat)

    # reduce latitude and longitude to valid and masked points
    indy,indx = np.nonzero(fd['mask'])
    lon,lat = (gridlon[indy,indx],latitude_geocentric[indy,indx])
    # scaled areas
    ps_scale = mdlhmc.spatial.scale_factors(gridlat[indy,indx], flat=flat,
        reference_latitude=reference_latitude)

    # unit parameters for each input variable type
    if (VARIABLE == 'accum_SMB'):
        # densities of meteoric ice [kg/m^3]
        rho_ice = 917.0
        # scaling factor to convert inputs from from meters ice eq to g
        scaling_factor = 1000.0*rho_ice*ps_scale*fd['area'][indy,indx]
        product_name = 'SMB'
        # use named point mass units code (grams)
        UNITS = 1
        # output spherical harmonic units
        harmonic_units = 'Geodesy_Normalization'
    elif (VARIABLE == 'dFAC'):
        # areas in terms of solid angle (steradians)
        scaling_factor = ps_scale*fd['area'][indy,indx]/(rad_e**2)
        product_name = 'FAC'
        # use custom UNITS to keep as inputs but use 4-pi norm
        UNITS = np.ones((LMAX+1))/(4.0*np.pi)
        # output spherical harmonic units
        harmonic_units = copy.copy(variable_units)

    # read load love numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE, FORMAT='class')
    # upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    # output string for both LMAX == MMAX and LMAX != MMAX cases
    order_str = 'M{MMAX:d}' if (MMAX != LMAX) else ''

    # allocate for output spherical harmonics
    Ylms = gravtk.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX+1,MMAX+1,nt-1))
    Ylms.slm = np.zeros((LMAX+1,MMAX+1,nt-1))
    Ylms.time = np.zeros((nt-1))
    Ylms.month = np.zeros((nt-1), dtype=np.int64)
    # for each time step
    for t in range(nt-1):
        # calculate date parameters for time step
        # dates are already set as mid-month values
        Ylms.time[t] = np.copy(fd['time'][t])
        Ylms.month[t] = gravtk.time.calendar_to_grace(Ylms.time[t])
        dpm = gravtk.time.calendar_days(np.floor(Ylms.time[t]))
        # calculate 2-month moving average
        # weighting by number of days in each month
        M1 = dpm[t % 12]*fd[VARIABLE][t,indy,indx]
        M2 = dpm[(t+1) % 12]*fd[VARIABLE][t+1,indy,indx]
        W = np.float64(dpm[(t+1) % 12] + dpm[t % 12])
        # reduce data for date and scale
        scaled = scaling_factor*(M1+M2)/W
        # convert to spherical harmonics
        YLMS = gravtk.gen_point_load(scaled, lon, lat,
            LMAX=LMAX, MMAX=MMAX, UNITS=UNITS, LOVE=LOVE)
        # copy harmonics for time step
        Ylms.clm[:,:,t] = YLMS.clm[:,:].copy()
        Ylms.slm[:,:,t] = YLMS.slm[:,:].copy()

    # output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')
    # attributes for output files
    attributes = {}
    attributes['institution'] = 'NASA Jet Propulsion Laboratory (JPL)'
    attributes['project'] = 'Glacier Energy and Mass Balance (GEMB)'
    attributes['product_region'] = region
    attributes['product_version'] = version.replace('_','.')
    attributes['product_name'] = product_name
    attributes['product_type'] = 'gravity_field'
    # add attributes for earth parameters if converting from mass
    if (VARIABLE == 'accum_SMB'):
        attributes['earth_model'] = LOVE.model
        attributes['earth_love_numbers'] = LOVE.citation
        attributes['reference_frame'] = LOVE.reference
    # add attributes for maximum degree and order
    attributes['max_degree'] = LMAX
    attributes['max_order'] = MMAX
    attributes['lineage'] = model_file.name
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'
    # add attributes to output harmonics
    Ylms.attributes['ROOT'] = attributes
    # output spherical harmonic data file
    args = (version,region,product_name,LMAX,order_str,suffix[DATAFORM])
    FILE = 'GEMB_{0}_{1}_{2}_CLM_L{3:d}{4}.{5}'.format(*args)
    output_file = model_file.with_name(FILE)
    Ylms.to_file(output_file, format=DATAFORM, date=True)
    # change the permissions mode of the output file to MODE
    output_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Read GEMB variables and convert to spherical harmonics
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('infile',
        type=pathlib.Path,
        help='GEMB SMB file to run')
    # GEMB product to convert to spherical harmonics
    parser.add_argument('--product','-P',
        type=str, default='accum_SMB', choices=('accum_SMB','dFAC'),
        help='GEMB product to calculate')
    # mask file for reducing to regions
    parser.add_argument('--mask',
        type=pathlib.Path,
        nargs='+', default=[],
        help='netCDF4 masks file for reducing to regions')
    # area file for reducing to regions
    parser.add_argument('--area',
        type=pathlib.Path,
        help='netCDF4 area file for calculating mass')
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
    gemb_smb_harmonics(args.infile,
        PRODUCT=args.product,
        MASKS=args.mask,
        AREA=args.area,
        LMAX=args.lmax,
        MMAX=args.mmax,
        LOVE_NUMBERS=args.love,
        REFERENCE=args.reference,
        DATAFORM=args.format,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
