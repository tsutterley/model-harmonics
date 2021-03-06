#!/usr/bin/env python
u"""
merra_smb_harmonics.py
Written by Tyler Sutterley (03/2021)
Reads monthly MERRA-2 surface mass balance anomalies and
    converts to spherical harmonic coefficients

INPUTS:
    SMB: Surface Mass Balance
    PRECIP: Total Precipitation
    RUNOFF: Meltwater Runoff

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -m X, --mean X: Year range for mean
    -Y X, --year X: Years to run
    -R X, --region X: region name for subdirectory
    --mask X: netCDF4 mask files for reducing to regions
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -n X, --love X: Load Love numbers dataset
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    -r X, --reference X: Reference frame for load love numbers
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
    plm_holmes.py: computes fully-normalized associated Legendre polynomials
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    ref_ellipsoid.py: calculate reference parameters for common ellipsoids
    gen_stokes.py: converts a spatial field into a series of spherical harmonics
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
        destripe_harmonics.py: calculates the decorrelation (destriping) filter
            and filters the GRACE/GRACE-FO coefficients for striping errors
        ncdf_read_stokes.py: reads spherical harmonic netcdf files
        ncdf_stokes.py: writes output spherical harmonic data to netcdf
        hdf5_read_stokes.py: reads spherical harmonic HDF5 files
        hdf5_stokes.py: writes output spherical harmonic data to HDF5
    spatial.py: spatial data class for reading, writing and processing data
        ncdf_read.py: reads input spatial data from netCDF4 files
        hdf5_read.py: reads input spatial data from HDF5 files
        ncdf_write.py: writes output spatial data to netCDF4
        hdf5_write.py: writes output spatial data to HDF5
    time.py: utilities for calculating time operations
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 03/2021: automatically update years to run based on current time
    Updated 02/2021: can use multiple mask files to create a combined solution
        replaced numpy bool to prevent deprecation warning
    Updated 01/2021: added more love number options
        set spatial variables for both 025 and 10 cases
        using utilities from time module. added maximum harmonic order option
        harmonics object output from gen_stokes.py
    Updated 10/2020: use argparse to set command line parameters
    Updated 10/2019: changing Y/N flags to True/False
    Updated 08/2017: convert from geodetic coordinates to geocentric
    Updated 01/2017: can output different data products (SMB, PRECIP, RUNOFF)
    Updated 11/2016: changes to shapefile read for Antarctica
    Written 11/2016
"""
from __future__ import print_function

import sys
import os
import re
import netCDF4
import argparse
import numpy as np
import gravity_toolkit.time
from gravity_toolkit.read_love_numbers import read_love_numbers
from gravity_toolkit.harmonics import harmonics
from gravity_toolkit.spatial import spatial
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.gen_stokes import gen_stokes
from gravity_toolkit.utilities import get_data_path
from geoid_toolkit.ref_ellipsoid import ref_ellipsoid

#-- PURPOSE: read Merra-2 cumulative data and convert to spherical harmonics
def merra_smb_harmonics(ddir, PRODUCT, YEARS, RANGE=None, REGION=None,
    MASKS=None, LMAX=0, MMAX=None, LOVE_NUMBERS=0, REFERENCE=None,
    DATAFORM=None, VERBOSE=False, MODE=0o775):

    #-- setup subdirectories
    cumul_sub = '{0}.5.12.4.CUMUL.{1:d}.{2:d}'.format(PRODUCT,*RANGE)
    #-- Creating output subdirectory if it doesn't exist
    prefix = '{0}_'.format(REGION) if REGION else ''
    output_sub = '{0}{1}_5.12.4_CUMUL_CLM_L{2:d}'.format(prefix,PRODUCT,LMAX)
    if (not os.access(os.path.join(ddir,output_sub), os.F_OK)):
        os.makedirs(os.path.join(ddir,output_sub),MODE)
    #-- titles for each data product
    merra_products = {}
    merra_products['SMB'] = 'MERRA-2 Surface Mass Balance'
    merra_products['PRECIP'] = 'MERRA-2 Precipitation'
    merra_products['RUNOFF'] = 'MERRA-2 Meltwater Runoff'
    #-- output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    #-- upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    #-- output string for both LMAX == MMAX and LMAX != MMAX cases
    order_str = 'M{0:d}'.format(MMAX) if (MMAX != LMAX) else ''

    #-- output dimensions and extents
    nlat,nlon = (361,576)
    extent = [-180.0,179.375,-90.0,90.0]
    #-- grid spacing
    dlon,dlat = (0.625,0.5)
    #-- latitude and longitude
    glon = np.arange(extent[0],extent[1]+dlon,dlon)
    glat = np.arange(extent[2],extent[3]+dlat,dlat)
    #-- create mesh grid of latitude and longitude
    gridlon,gridlat = np.meshgrid(glon,glat)

    #-- read masks for reducing regions before converting to harmonics
    input_mask = np.zeros((nlat,nlon),dtype=bool)
    for mask_file in MASKS:
        fileID = netCDF4.Dataset(mask_file,'r')
        input_mask |= fileID.variables['mask'][:].astype(bool)
        fileID.close()

    #-- Earth Parameters
    ellipsoid_params = ref_ellipsoid('WGS84')
    #-- semimajor axis of ellipsoid [m]
    a_axis = ellipsoid_params['a']
    #--  first numerical eccentricity
    ecc1 = ellipsoid_params['ecc1']
    #-- convert from geodetic latitude to geocentric latitude
    #-- geodetic latitude in radians
    latitude_geodetic_rad = np.pi*gridlat/180.0
    #-- prime vertical radius of curvature
    N = a_axis/np.sqrt(1.0 - ecc1**2.*np.sin(latitude_geodetic_rad)**2.)
    #-- calculate X, Y and Z from geodetic latitude and longitude
    X = N * np.cos(latitude_geodetic_rad) * np.cos(np.pi*gridlon/180.0)
    Y = N * np.cos(latitude_geodetic_rad) * np.sin(np.pi*gridlon/180.0)
    Z = (N * (1.0 - ecc1**2.0)) * np.sin(latitude_geodetic_rad)
    #-- calculate geocentric latitude and convert to degrees
    latitude_geocentric = 180.0*np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))/np.pi
    #-- colatitude in radians
    theta = (90.0 - latitude_geocentric[:,0])*np.pi/180.0

    #-- read load love numbers
    LOVE = load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)

    #-- calculate Legendre polynomials
    PLM,dPLM = plm_holmes(LMAX,np.cos(theta))

    #-- find input files from merra_smb_cumulative.py
    regex_years = r'\d+' if (YEARS is None) else '|'.join(map(str,YEARS))
    args = (PRODUCT, regex_years, suffix[DATAFORM])
    regex_pattern = r'MERRA2_(\d+).tavgM_2d_{0}_cumul_Nx.({1})(\d+).{2}$'
    rx = re.compile(regex_pattern.format(*args), re.VERBOSE)
    #-- will be out of order for 2020 due to September reprocessing
    FILES = sorted([fi for fi in os.listdir(os.path.join(ddir,cumul_sub))
        if rx.match(fi)])

    #-- for each input file
    for t,fi in enumerate(FILES):
        #-- extract parameters from input flux file
        MOD,Y1,M1 = rx.findall(fi).pop()
        #-- read data file for data format
        if (DATAFORM == 'ascii'):
            #-- ascii (.txt)
            merra_data = spatial(spacing=[dlon,dlat],nlat=nlat,nlon=nlon,
                extent=extent).from_ascii(os.path.join(ddir,cumul_sub,fi))
        elif (DATAFORM == 'netCDF4'):
            #-- netCDF4 (.nc)
            merra_data = spatial().from_netCDF4(os.path.join(ddir,cumul_sub,fi),
                varname=PRODUCT, verbose=VERBOSE)
        elif (DATAFORM == 'HDF5'):
            #-- HDF5 (.H5)
            merra_data = spatial().from_HDF5(os.path.join(ddir,cumul_sub,fi),
                varname=PRODUCT, verbose=VERBOSE)

        #-- if reducing to a region of interest before converting to harmonics
        if np.any(input_mask):
            #-- replace fill value points and masked points with 0
            merra_data.replace_invalid(0.0, mask=np.logical_not(input_mask))
        else:
            #-- replace fill value points points with 0
            merra_data.replace_invalid(0.0)

        #-- convert to spherical harmonics from mm w.e.
        merra_Ylms = gen_stokes(merra_data.data, glon, latitude_geocentric[:,0],
            LMAX=LMAX, MMAX=MMAX, UNITS=3, PLM=PLM, LOVE=LOVE)
        #-- copy date information
        merra_Ylms.time = np.copy(merra_data.time)
        #-- calculate GRACE/GRACE-FO month
        merra_Ylms.month = np.int(12.0*(np.float(Y1) - 2002.0) + np.float(M1))
        #-- output spherical harmonic data file
        args=(MOD,PRODUCT,LMAX,order_str,merra_Ylms.month,suffix[DATAFORM])
        FILE='MERRA2_{0}_tavgM_2d_{1}_CLM_L{2:d}{3}_{4:03d}.{5}'.format(*args)

        #-- output data for month
        print(os.path.join(ddir,output_sub,FILE)) if VERBOSE else None
        if (DATAFORM == 'ascii'):
            #-- ascii (.txt)
            merra_Ylms.to_ascii(os.path.join(ddir,output_sub,FILE))
        elif (DATAFORM == 'netCDF4'):
            #-- netcdf (.nc)
            merra_Ylms.to_netCDF4(os.path.join(ddir,output_sub,FILE))
        elif (DATAFORM == 'HDF5'):
            #-- HDF5 (.H5)
            merra_Ylms.to_HDF5(os.path.join(ddir,output_sub,FILE))
        #-- change the permissions mode of the output file to MODE
        os.chmod(os.path.join(ddir,output_sub,FILE),MODE)

    #-- Output date ascii file
    output_date_file = 'MERRA2_{0}_DATES.txt'.format(PRODUCT)
    fid1 = open(os.path.join(ddir,output_sub,output_date_file), 'w')
    #-- date file header information
    print('{0:8} {1:^6} {2:^5}'.format('Mid-date','GRACE','Month'), file=fid1)
    #-- index file listing all output spherical harmonic files
    output_index_file = 'index.txt'
    fid2 = open(os.path.join(ddir,output_sub,output_index_file),'w')
    #-- find all available output files
    args = (PRODUCT, LMAX, order_str, suffix[DATAFORM])
    output_pattern = r'MERRA2_(\d+)_tavgM_2d_{0}_CLM_L{1:d}{2}_([-]?\d+).{3}'
    output_regex = re.compile(output_pattern.format(*args))
    #-- find all output ECCO OBP harmonic files (not just ones created in run)
    output_files=[fi for fi in os.listdir(os.path.join(ddir,output_sub))
        if re.match(output_regex,fi)]
    for fi in sorted(output_files):
        #-- full path to output file
        full_output_file = os.path.join(ddir,output_sub,fi)
        #-- extract GRACE month
        MOD,grace_month=np.array(re.findall(output_regex,fi).pop(),dtype=np.int)
        YY = 2002.0 + np.floor((grace_month-1)/12.0)
        MM = ((grace_month-1) % 12) + 1
        tdec, = gravity_toolkit.time.convert_calendar_decimal(YY, MM)
        #-- full path to output file
        full_output_file = os.path.join(ddir,output_sub,fi)
        #-- print date, GRACE month and calendar month to date file
        fid1.write('{0:11.6f} {1:03d} {2:02.0f}\n'.format(tdec,grace_month,MM))
        #-- print output file to index
        print(full_output_file.replace(os.path.expanduser('~'),'~'), file=fid2)
    #-- close the date and index files
    fid1.close()
    fid2.close()
    #-- set the permissions level of the output date and index files to MODE
    os.chmod(os.path.join(ddir,output_sub,output_date_file), MODE)
    os.chmod(os.path.join(ddir,output_sub,output_index_file), MODE)

#-- PURPOSE: read load love numbers for the range of spherical harmonic degrees
def load_love_numbers(LMAX, LOVE_NUMBERS=0, REFERENCE='CF'):
    """
    Reads PREM load Love numbers for the range of spherical harmonic degrees
    and applies isomorphic parameters

    Arguments
    ---------
    LMAX: maximum spherical harmonic degree

    Keyword arguments
    -----------------
    LOVE_NUMBERS: Load Love numbers dataset
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    REFERENCE: Reference frame for calculating degree 1 love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth

    Returns
    -------
    kl: Love number of Gravitational Potential
    hl: Love number of Vertical Displacement
    ll: Love number of Horizontal Displacement
    """
    #-- load love numbers file
    if (LOVE_NUMBERS == 0):
        #-- PREM outputs from Han and Wahr (1995)
        #-- https://doi.org/10.1111/j.1365-246X.1995.tb01819.x
        love_numbers_file = get_data_path(['data','love_numbers'])
        header = 2
        columns = ['l','hl','kl','ll']
    elif (LOVE_NUMBERS == 1):
        #-- PREM outputs from Gegout (2005)
        #-- http://gemini.gsfc.nasa.gov/aplo/
        love_numbers_file = get_data_path(['data','Load_Love2_CE.dat'])
        header = 3
        columns = ['l','hl','ll','kl']
    elif (LOVE_NUMBERS == 2):
        #-- PREM outputs from Wang et al. (2012)
        #-- https://doi.org/10.1016/j.cageo.2012.06.022
        love_numbers_file = get_data_path(['data','PREM-LLNs-truncated.dat'])
        header = 1
        columns = ['l','hl','ll','kl','nl','nk']
    #-- LMAX of load love numbers from Han and Wahr (1995) is 696.
    #-- from Wahr (2007) linearly interpolating kl works
    #-- however, as we are linearly extrapolating out, do not make
    #-- LMAX too much larger than 696
    #-- read arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = read_love_numbers(love_numbers_file, LMAX=LMAX, HEADER=header,
        COLUMNS=columns, REFERENCE=REFERENCE, FORMAT='tuple')
    #-- return a tuple of load love numbers
    return (hl,kl,ll)

#-- Main program that calls merra_monthly_harmonics()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Reads monthly MERRA-2 surface mass balance
            anomalies and converts to spherical harmonic coefficients
            """
    )
    #-- command line parameters
    parser.add_argument('product',
        type=str, nargs='+', choices=['SMB','PRECIP','RUNOFF'],
        help='MERRA-2 derived product')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- start and end years to run for mean
    parser.add_argument('--mean',
        metavar=('START','END'), type=int, nargs=2,
        default=[1980,1995],
        help='Start and end year range for mean')
    #-- years to run
    now = gravity_toolkit.time.datetime.datetime.now()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.year+1),
        help='Years of model outputs to run')
    #-- region name for subdirectory
    parser.add_argument('--region','-R',
        type=str, default=None,
        help='Region name for subdirectory')
    #-- mask file for reducing to regions
    parser.add_argument('--mask',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        nargs='+', help='netCDF4 masks file for reducing to regions')
    #-- maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
    #-- different treatments of the load Love numbers
    #-- 0: Han and Wahr (1995) values from PREM
    #-- 1: Gegout (2005) values from PREM
    #-- 2: Wang et al. (2012) values from PREM
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2],
        help='Treatment of the Load Love numbers')
    #-- option for setting reference frame for gravitational load love number
    #-- reference frame options (CF, CM, CE)
    parser.add_argument('--reference','-r',
        type=str.upper, default='CF', choices=['CF','CM','CE'],
        help='Reference frame for load Love numbers')
    #-- input and output data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input and output data format')
    #-- print information about each input and output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files')
    args = parser.parse_args()

    #-- run program for each input product
    for PRODUCT in args.product:
        #-- run program
        merra_smb_harmonics(args.directory, PRODUCT, args.year, RANGE=args.mean,
            REGION=args.region, MASKS=args.mask, LMAX=args.lmax, MMAX=args.mmax,
            LOVE_NUMBERS=args.love, REFERENCE=args.reference,
            DATAFORM=args.format, VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
