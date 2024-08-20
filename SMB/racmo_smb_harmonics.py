#!/usr/bin/env python
u"""
racmo_smb_harmonics.py
Written by Tyler Sutterley (03/2023)
Read RACMO surface mass balance products and converts to spherical harmonics
Shifts dates of SMB point masses to mid-month values to correspond with GRACE

CALLING SEQUENCE:
    python racmo_smb_harmonics.py --product smb --verbose <path_to_racmo_file>

COMMAND LINE OPTIONS:
    -P X, --product X: RACMO SMB product to calculate
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
    -F X, --format X: Input and output data format
        ascii
        netcdf
        HDF5
    -M X, --mode X: Permission mode of directories and files
    -V, --verbose: Output information for each output file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://h5py.org

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    utilities.py: download and management utilities for files
    datum.py: calculate reference parameters for common ellipsoids
    gen_point_load.py: calculates spherical harmonics from point masses
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors

UPDATE HISTORY:
    Updated 03/2023: add root attributes to output netCDF4 and HDF5 files
        use spatial function for calculating geocentric latitude
    Updated 02/2023: use love numbers class with additional attributes
    Updated 12/2022: single implicit import of spherical harmonic tools
        use constants class in place of geoid-toolkit ref_ellipsoid
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: use wrapper function for reading load Love numbers
        deprecation fixes for regular expressions
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 11/2021: complete rewrite of program
        dropped old RACMO ascii file read portions
    Updated 06/2018: using python3 compatible octal and input
    Updated 04/2015: minor code update (os)
    Updated 10/2014: code update to improve computational times
    Updated 06/2014 added main definition
        can now run from command line or parameter file
    Updated 02/2014: more general code updates
    Written 10/2011
"""
from __future__ import print_function

import sys
import re
import copy
import gzip
import uuid
import logging
import netCDF4
import pathlib
import argparse
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc

# PURPOSE: convert RACMO surface mass balance products to spherical harmonics
def racmo_smb_harmonics(model_file, VARIABLE,
    MASKS=None,
    LMAX=None,
    MMAX=None,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    DATAFORM=None,
    GZIP=False,
    MODE=0o775):

    # get reference parameters for ellipsoid
    ellipsoid_params = mdlhmc.datum(ellipsoid='WGS84')
    # semimajor axis of ellipsoid [cm]
    a_axis = ellipsoid_params.a_axis
    # ellipsoidal flattening
    flat = ellipsoid_params.flat

    # RACMO SMB model file
    model_file = pathlib.Path(model_file).expanduser().absolute()
    # try to extract region and version from filename
    R1 = re.compile(r'[XF]?(ANT27|GRN11|GRN055|PEN055|ASE055)',re.VERBOSE)
    R2 = re.compile(r'(RACMO\d+(\.\d+)?(p\d+)?)',re.VERBOSE)
    REGION = R1.search(model_file.name).group(0)
    VERSION = R2.search(model_file.name).group(0)
    # RACMO products
    racmo_products = {}
    racmo_products['precip'] = 'Precipitation'
    racmo_products['rainfall'] = 'Rainfall'
    racmo_products['refreeze'] = 'Meltwater Refreeze'
    racmo_products['runoff'] = 'Meltwater Runoff'
    racmo_products['smb'] = 'Surface Mass Balance'
    racmo_products['sndiv'] = 'Snow Drift Erosion'
    racmo_products['snowfall'] = 'Snowfall'
    racmo_products['snowmelt'] = 'Snowmelt'
    racmo_products['subl'] = 'Sublimation'

    # Open the RACMO SMB NetCDF file for reading
    if GZIP:
        # read as in-memory (diskless) netCDF4 dataset
        with gzip.open(str(model_file), mode='r') as f:
            fileID = netCDF4.Dataset(uuid.uuid4().hex, memory=f.read())
    else:
        # read netCDF4 dataset
        fileID = netCDF4.Dataset(model_file, mode='r')

    # Output NetCDF file information
    logging.info(str(model_file))
    logging.info(list(fileID.variables.keys()))

    # Get data from each netCDF variable and remove singleton dimensions
    fd = {}
    # extract data variable
    fd[VARIABLE] = np.squeeze(fileID.variables[VARIABLE][:].copy())
    # read latitude, longitude and rotated pole coordinates
    fd['lon'] = fileID.variables['lon'][:,:].copy()
    gridlat = fileID.variables['lat'][:,:].copy()
    rlon = fileID.variables['rlon'][:].copy()
    rlat = fileID.variables['rlat'][:].copy()
    # time within netCDF files is days since epoch
    TIME = fileID.variables['time'][:].copy()
    time_string = fileID.variables['time'].units
    epoch1,to_secs = gravtk.time.parse_date_string(time_string)
    # calculate Julian day by converting to MJD and adding offset
    JD = gravtk.time.convert_delta_time(TIME*to_secs,
        epoch1=epoch1, epoch2=(1858,11,17,0,0,0),
        scale=1.0/86400.0) + 2400000.5
    # convert from Julian days to calendar dates
    YY,MM,DD,hh,mm,ss = gravtk.time.convert_julian(JD,
        FORMAT='tuple')
    # convert from calendar dates to year-decimal
    fd['time'] = gravtk.time.convert_calendar_decimal(YY,MM,
        day=DD,hour=hh,minute=mm,second=ss)
    # invalid data value
    fv = np.float64(fileID.variables[VARIABLE]._FillValue)
    # input shape of RACMO SMB firn data
    nt,ny,nx = np.shape(fd[VARIABLE])
    # calculate grid areas (assume fully ice covered)
    dlon = np.pi*np.abs(rlon[1] - rlon[0])/180.0
    dlat = np.pi*np.abs(rlat[1] - rlat[0])/180.0
    _,gridrlat = np.meshgrid(rlon,rlat)
    fd['area'] = np.cos(gridrlat*np.pi/180.0)*(a_axis**2)*dlon*dlat
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
        fd['mask'] |= fileID.variables['maskgrounded2d'][:].astype(bool)
        fileID.close()
    # indices of valid RACMO data
    fd['mask'] &= (fd[VARIABLE].data[0,:,:] != fv)

    # calculate geocentric latitude and convert to degrees
    fd['lat'] = mdlhmc.spatial.geocentric_latitude(fd['lon'], gridlat,
        a_axis=a_axis, flat=flat)

    # reduce latitude and longitude to valid and masked points
    indy,indx = np.nonzero(fd['mask'])
    lon,lat = (fd['lon'][indy,indx],fd['lat'][indy,indx])
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
    Ylms.month = np.zeros((nt-1),dtype=np.int64)
    # for each time step
    for t in range(nt-1):
        # calculate date parameters for time step
        Ylms.time[t] = np.mean([fd['time'][t],fd['time'][t+1]])
        Ylms.month[t] = gravtk.time.calendar_to_grace(Ylms.time[t])
        dpm = gravtk.time.calendar_days(np.floor(Ylms.time[t]))
        # calculate 2-month moving average
        # weighting by number of days in each month
        M1 = dpm[t % 12]*fd[VARIABLE][t,indy,indx]
        M2 = dpm[(t+1) % 12]*fd[VARIABLE][t+1,indy,indx]
        W = np.float64(dpm[(t+1) % 12] + dpm[t % 12])
        # reduce data for date and convert to mass (g)
        ptms = 1000.0*fd['area'][indy,indx]*(M1+M2)/W
        racmo_Ylms = gravtk.gen_point_load(ptms, lon, lat, LMAX=LMAX,
            MMAX=MMAX, UNITS=1, LOVE=LOVE)
        # copy harmonics for time step
        Ylms.clm[:,:,t] = racmo_Ylms.clm[:,:].copy()
        Ylms.slm[:,:,t] = racmo_Ylms.slm[:,:].copy()
    # adjust months to be consistent
    Ylms.month = gravtk.time.adjust_months(Ylms.month)

    # output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')
    # attributes for output files
    attributes = {}
    attributes['project'] = 'Regional Atmospheric and Climate Model (RACMO)'
    attributes['title'] = copy.copy(racmo_products[VARIABLE])
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
    attributes['lineage'] = model_file.name
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'
    # add attributes to output harmonics
    Ylms.attributes['ROOT'] = attributes
    # output spherical harmonic data file
    args = (VERSION,REGION,VARIABLE.upper(),LMAX,order_str,suffix[DATAFORM])
    FILE = '{0}_{1}_{2}_CLM_L{3:d}{4}.{5}'.format(*args)
    output_file = model_file.with_name(FILE)
    Ylms.to_file(output_file, format=DATAFORM, date=True)
    # change the permissions mode of the output file to MODE
    output_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Read RACMO surface mass balance products and
            converts to spherical harmonics
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('infile',
        type=pathlib.Path,
        help='RACMO SMB file to run')
    # products from SMB model
    choices = ['precip','rainfall','refreeze','runoff','smb',
        'sndiv','snowfall','snowmelt','subl']
    parser.add_argument('--product','-P',
        type=str, metavar='PRODUCT', default='smb', choices=choices,
        help='RACMO SMB product to calculate')
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
    racmo_smb_harmonics(args.infile, args.product,
        MASKS=args.mask,
        LMAX=args.lmax,
        MMAX=args.mmax,
        LOVE_NUMBERS=args.love,
        REFERENCE=args.reference,
        DATAFORM=args.format,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
