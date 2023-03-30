#!/usr/bin/env python
u"""
gemb_smb_harmonics.py
Written by Tyler Sutterley (03/2023)
Read GEMB SMB variables and convert to spherical harmonics
Shifts dates of SMB point masses to mid-month values to correspond with GRACE

CALLING SEQUENCE:
    python gemb_smb_harmonics.py --lmax 60 --verbose <path_to_gemb_file>

COMMAND LINE OPTIONS:
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
    constants.py: calculate reference parameters for common ellipsoids
    gen_point_load.py: calculates spherical harmonics from point masses
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: spatial data class for reading, writing and processing data

UPDATE HISTORY:
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
import os
import re
import pyproj
import logging
import netCDF4
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
    MASKS=None,
    AREA=None,
    LMAX=None,
    MMAX=None,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    DATAFORM=None,
    MODE=0o775):

    # GEMB directory
    DIRECTORY = os.path.dirname(model_file)
    # regular expression pattern for extracting parameters
    pattern = r'GEMB_(Greenland|Antarctica)_(.*?)_(v.*?).nc$'
    region, _, version = re.findall(pattern, model_file).pop()

    # Open the GEMB NetCDF file for reading
    fileID = netCDF4.Dataset(os.path.expanduser(model_file), 'r')

    # Output NetCDF file information
    logging.info(model_file)
    logging.info(list(fileID.variables.keys()))

    # Get data from each netCDF variable and remove singleton dimensions
    fd = {}
    # read reduced time variables
    fd['time'] = fileID.variables['time'][:].copy()
    # read surface mass balance variables
    fd['accum_SMB'] = fileID.variables['accum_SMB'][:,:,:].copy()
    # invalid data value
    try:
        fv = np.float64(fileID.variables['accum_SMB']._FillValue)
    except (ValueError,AttributeError):
        fv = np.nan
    # input shape of GEMB firn data
    nt,ny,nx = np.shape(fd['accum_SMB'])
    # extract x and y coordinate arrays
    fd['x'] = fileID.variables['x'][:].copy()
    fd['y'] = fileID.variables['y'][:].copy()
    xg,yg = np.meshgrid(fd['x'],fd['y'])
    # calculate grid areas (read file or assume fully ice covered)
    fd['area'] = np.ma.zeros((ny,nx), fill_value=fv)
    fd['area'].mask = np.zeros((ny,nx), dtype=bool)
    if AREA is None:
        dx = np.abs(fd['x'][1] - fd['x'][0])
        dy = np.abs(fd['y'][1] - fd['y'][0])
        fd['area'].data[:,:] = dx*dy
    else:
        # read area file (km^2)
        fileID = netCDF4.Dataset(AREA, 'r')
        fd['area'][:,:] = 1e6*fileID.variables['area'][:]
        fileID.close()
    # close the NetCDF files
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
    fd['mask'] &= (fd['accum_SMB'].data[0,:,:] != fv)
    fd['mask'] &= np.isfinite(fd['accum_SMB'].data[0,:,:])
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
    ellipsoid_params = mdlhmc.constants(ellipsoid='WGS84')
    # semimajor axis of the ellipsoid [m]
    a_axis = ellipsoid_params.a_axis
    # ellipsoidal flattening
    flat = ellipsoid_params.flat
    # calculate geocentric latitude and convert to degrees
    latitude_geocentric = mdlhmc.spatial.geocentric_latitude(gridlon, gridlat,
        a_axis=a_axis, flat=flat)

    # reduce latitude and longitude to valid and masked points
    indy,indx = np.nonzero(fd['mask'])
    lon,lat = (gridlon[indy,indx],latitude_geocentric[indy,indx])
    # scaled areas
    ps_scale = mdlhmc.spatial.scale_areas(gridlat[indy,indx], flat=flat,
        ref=reference_latitude)
    scaled_area = ps_scale*fd['area'][indy,indx]
    # read load love numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE, FORMAT='class')
    # upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    # output string for both LMAX == MMAX and LMAX != MMAX cases
    order_str = 'M{MMAX:d}' if (MMAX != LMAX) else ''
    # densities of meteoric ice
    rho_ice = 917.0

    # allocate for output spherical harmonics
    Ylms = gravtk.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX+1,MMAX+1,nt-1))
    Ylms.slm = np.zeros((LMAX+1,MMAX+1,nt-1))
    Ylms.time = np.zeros((nt-1))
    Ylms.month = np.zeros((nt-1),dtype=np.int64)
    # for each time step
    for t in range(nt-1):
        # calculate date parameters for time step
        # dates are already set as mid-month values
        Ylms.time[t] = np.copy(fd['time'][t])
        Ylms.month[t] = gravtk.time.calendar_to_grace(Ylms.time[t])
        dpm = gravtk.time.calendar_days(np.floor(Ylms.time[t]))
        # calculate 2-month moving average
        # weighting by number of days in each month
        M1 = dpm[t % 12]*fd['accum_SMB'][t,indy,indx]
        M2 = dpm[(t+1) % 12]*fd['accum_SMB'][t+1,indy,indx]
        W = np.float64(dpm[(t+1) % 12] + dpm[t % 12])
        # reduce data for date and convert to mass (g)
        GEMB_mass = 1000.0*rho_ice*scaled_area*(M1+M2)/W
        # convert to spherical harmonics
        GEMB_Ylms = gravtk.gen_point_load(GEMB_mass, lon, lat,
            LMAX=LMAX, MMAX=MMAX, UNITS=1, LOVE=LOVE)
        # copy harmonics for time step
        Ylms.clm[:,:,t] = GEMB_Ylms.clm[:,:].copy()
        Ylms.slm[:,:,t] = GEMB_Ylms.slm[:,:].copy()

    # output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')
    # attributes for output files
    attributes = {}
    attributes['institution'] = 'NASA Jet Propulsion Laboratory (JPL)'
    attributes['project'] = 'Glacier Energy and Mass Balance (GEMB)'
    attributes['product_region'] = region
    attributes['product_version'] = version.replace('_','.')
    attributes['product_name'] = 'SMB'
    attributes['product_type'] = 'gravity_field'
    # add attributes for earth parameters
    attributes['earth_model'] = LOVE.model
    attributes['earth_love_numbers'] = LOVE.citation
    attributes['reference_frame'] = LOVE.reference
    # add attributes for maximum degree and order
    attributes['max_degree'] = LMAX
    attributes['max_order'] = MMAX
    attributes['lineage'] = os.path.basename(model_file)
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'
    # add attributes to output harmonics
    Ylms.attributes['ROOT'] = attributes
    # output spherical harmonic data file
    args = (version,region,'SMB',LMAX,order_str,suffix[DATAFORM])
    FILE = 'GEMB_{0}_{1}_{2}_CLM_L{3:d}{4}.{5}'.format(*args)
    Ylms.to_file(os.path.join(DIRECTORY,FILE), format=DATAFORM, date=True)
    # change the permissions mode of the output file to MODE
    os.chmod(os.path.join(DIRECTORY,FILE),MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Read GEMB SMB variables and convert to spherical harmonics
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='GEMB SMB file to run')
    # mask file for reducing to regions
    parser.add_argument('--mask',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        nargs='+', default=[],
        help='netCDF4 masks file for reducing to regions')
    # area file for reducing to regions
    parser.add_argument('--area',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
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