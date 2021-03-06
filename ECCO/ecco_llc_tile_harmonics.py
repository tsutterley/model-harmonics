#!/usr/bin/env python
u"""
ecco_llc_tile_harmonics.py
Written by Tyler Sutterley (03/2021)
Reads monthly ECCO ocean bottom pressure anomalies from LLC tiles
    and converts to spherical harmonic coefficients

ECCO Version 4, Release 4
    https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/nctiles_monthly

ECCO Version 5, Alpha
    https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/nctiles_monthly

INPUTS:
    ECCO version 4 or 5 models
        V4r4: Version 4, Revision 4
        V5alpha: Version 5, Alpha release

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -Y X, --year X: Years to run
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
    -F X, --format X: Output data format
        ascii
        netcdf
        HDF5
    -V, --verbose: Output information for each output file
    -M X, --mode X: Permission mode of directories and files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        https://www.h5py.org/

PROGRAM DEPENDENCIES:
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    ref_ellipsoid.py: calculate reference parameters for common ellipsoids
    norm_gravity.py: calculates the normal gravity for locations on an ellipsoid
    gen_point_pressure.py: converts pressure point values to spherical harmonics
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
    Written 02/2021
"""
from __future__ import print_function

import os
import re
import netCDF4
import argparse
import numpy as np
import gravity_toolkit.time
import gravity_toolkit.spatial
import gravity_toolkit.harmonics
from gravity_toolkit.utilities import get_data_path
from gravity_toolkit.read_love_numbers import read_love_numbers
from geoid_toolkit.ref_ellipsoid import ref_ellipsoid
from geoid_toolkit.norm_gravity import norm_gravity
from model_harmonics.gen_point_pressure import gen_point_pressure

#-- PURPOSE: convert monthly ECCO OBP data to spherical harmonics
def ecco_llc_tile_harmonics(ddir, MODEL, YEARS, LMAX=0, MMAX=None,
    LOVE_NUMBERS=0, REFERENCE=None, DATAFORM=None, VERBOSE=False,
    MODE=0o775):
    #-- input and output subdirectories
    input_dir = os.path.join(ddir,'ECCO_{0}_AveRmvd_OBP'.format(MODEL),
        'nctiles_monthly')
    output_sub = 'ECCO_{0}_AveRmvd_OBP_CLM_L{1:d}'.format(MODEL,LMAX)
    #-- upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    #-- output string for both LMAX == MMAX and LMAX != MMAX cases
    order_str = 'M{0:d}'.format(MMAX) if (MMAX != LMAX) else ''
    #-- output file format
    output_file_format = 'ECCO_{0}_AveRmvd_OBP_CLM_L{1:d}{2}_{3:03d}.{4}'
    #-- Creating subdirectory if it doesn't exist
    if (not os.access(os.path.join(ddir,output_sub), os.F_OK)):
        os.makedirs(os.path.join(ddir,output_sub),MODE)
    #-- output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    #-- input variable names for each model
    if (MODEL == 'V4r4'):
        LONNAME = 'XC'
        LATNAME = 'YC'
        ZNAME = 'Depth'
        VARNAME = 'PHIBOT'
        AREANAME = 'rA'
        Nt,Nj,Ni = (13,90,90)
    elif (MODEL == 'V5alpha'):
        LONNAME = 'XC'
        LATNAME = 'YC'
        ZNAME = 'Depth'
        VARNAME = 'PHIBOT'
        AREANAME = 'rA'
        Nt,Nj,Ni = (13,270,270)

    #-- read ECCO tile grid file
    input_invariant_file = os.path.join(ddir,'ECCO-{0}'.format(MODEL),
        'nctiles_monthly','ECCO-GRID.nc')
    invariant = ncdf_invariant(input_invariant_file,
        lon=LONNAME, lat=LATNAME, depth=ZNAME, area=AREANAME)
    #-- read ECCO tile geoid height file from ecco_geoid_llc_tiles.py
    input_geoid_file = os.path.join(ddir,'ECCO-{0}'.format(MODEL),
        'nctiles_monthly','ECCO-EGM2008.nc')
    geoid_undulation = ncdf_geoid(input_geoid_file)
    #-- use geoid and depth to calculate bathymetry
    bathymetry = geoid_undulation - invariant['depth']

    #-- Earth Parameters
    ellipsoid_params = ref_ellipsoid('WGS84')
    #-- semimajor axis of ellipsoid [m]
    a_axis = ellipsoid_params['a']
    #--  first numerical eccentricity
    ecc1 = ellipsoid_params['ecc1']
    #-- convert from geodetic latitude to geocentric latitude
    #-- geodetic latitude and longitude in radians
    latitude_geodetic_rad = np.pi*invariant['lat']/180.0
    longitude_rad = np.pi*invariant['lon']/180.0
    #-- prime vertical radius of curvature
    N = a_axis/np.sqrt(1.0 - ecc1**2.*np.sin(latitude_geodetic_rad)**2.)
    #-- calculate X, Y and Z from geodetic latitude and longitude
    X = (N+bathymetry)*np.cos(latitude_geodetic_rad)*np.cos(longitude_rad)
    Y = (N+bathymetry)*np.cos(latitude_geodetic_rad)*np.sin(longitude_rad)
    Z = (N * (1.0 - ecc1**2.0) + bathymetry) * np.sin(latitude_geodetic_rad)
    R = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
    #-- calculate geocentric latitude and convert to degrees
    latitude_geocentric = 180.0*np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))/np.pi
    #-- colatitude in radians
    theta = (90.0 - latitude_geocentric)*np.pi/180.0

    #-- calculate normal gravity at latitudes and bathymetry
    gamma_h,dgamma_dh = norm_gravity(latitude_geocentric,bathymetry,'WGS84')

    #-- read load love numbers
    LOVE = load_love_numbers(LMAX,LOVE_NUMBERS=LOVE_NUMBERS,REFERENCE=REFERENCE)

    #-- regular expression pattern to find files and extract dates
    regex_years = r'\d+' if (YEARS is None) else '|'.join(map(str,YEARS))
    args = (MODEL, regex_years, suffix[DATAFORM])
    rx = re.compile(r'ECCO_{0}_AveRmvd_OBP_({1})_(\d+).{2}$'.format(*args))

    #-- find input ECCO OBP files
    FILES = [fi for fi in os.listdir(input_dir) if rx.match(fi)]
    #-- for each input file
    for f in sorted(FILES):
        #-- extract dates from file
        year,month = np.array(rx.findall(f).pop(), dtype=np.int)
        #-- read input data file
        obp_data = {}
        with netCDF4.Dataset(os.path.join(input_dir,f),'r') as fileID:
            for key,val in fileID.variables.items():
                obp_data[key] = val[:].copy()
        #-- allocate for output spherical harmonics
        obp_Ylms = gravity_toolkit.harmonics(lmax=LMAX, mmax=MMAX)
        obp_Ylms.clm = np.zeros((LMAX+1,MMAX+1))
        obp_Ylms.slm = np.zeros((LMAX+1,MMAX+1))
        #-- copy date information
        obp_Ylms.time = np.copy(obp_data['time'])
        obp_Ylms.month = 12*(year - 2002) + month
        #-- calculate harmonics for each spatial tile
        #-- add to the total set of harmonics
        for k in range(Nt):
            #-- valid indices for for tile
            indj,indi = np.nonzero(~obp_data[VARNAME].mask[k,:,:])
            #-- reduce to valid
            obp_tile = obp_data[VARNAME].data[k,indj,indi]
            gamma_tile = gamma_h[k,indj,indi]
            rad_tile = R[k,indi,indj]
            #-- grid point areas (m^2)
            area = invariant['area'][k,indj,indi]
            #-- geocentric latitude and longitude
            lat = latitude_geocentric[k,indj,indi]
            lon = invariant['lon'][k,indj,indi]
            #-- calculate harmonics from pressure/gravity ratio
            Ylms = gen_point_pressure(obp_tile, gamma_tile, rad_tile,
                lon, lat, AREA=area, LMAX=LMAX, MMAX=MMAX, LOVE=LOVE)
            #-- add tile Ylms to total for month
            obp_Ylms.add(Ylms)
        #-- output spherical harmonic data file
        args = (MODEL, LMAX, order_str, obp_Ylms.month, suffix[DATAFORM])
        FILE = output_file_format.format(*args)
        #-- output data for month
        print(os.path.join(ddir,output_sub,FILE)) if VERBOSE else None
        if (DATAFORM == 'ascii'):
            #-- ascii (.txt)
            obp_Ylms.to_ascii(os.path.join(ddir,output_sub,FILE))
        elif (DATAFORM == 'netCDF4'):
            #-- netcdf (.nc)
            obp_Ylms.to_netCDF4(os.path.join(ddir,output_sub,FILE))
        elif (DATAFORM == 'HDF5'):
            #-- HDF5 (.H5)
            obp_Ylms.to_HDF5(os.path.join(ddir,output_sub,FILE))
        #-- change the permissions mode of the output file to MODE
        os.chmod(os.path.join(ddir,output_sub,FILE),MODE)

    #-- Output date ascii file
    output_date_file = 'ECCO_{0}_OBP_DATES.txt'.format(MODEL)
    fid1 = open(os.path.join(ddir,output_sub,output_date_file), 'w')
    #-- date file header information
    print('{0:8} {1:^6} {2:^5}'.format('Mid-date','GRACE','Month'), file=fid1)
    #-- index file listing all output spherical harmonic files
    output_index_file = 'index.txt'
    fid2 = open(os.path.join(ddir,output_sub,output_index_file),'w')
    #-- find all available output files
    args = (MODEL, LMAX, suffix[DATAFORM])
    output_regex=r'ECCO_{0}_AveRmvd_OBP_CLM_L{1:d}_([-]?\d+).{2}'.format(*args)
    #-- find all output ECCO OBP harmonic files (not just ones created in run)
    output_files = [fi for fi in os.listdir(os.path.join(ddir,output_sub))
        if re.match(output_regex,fi)]
    for fi in sorted(output_files):
        #-- full path to output file
        full_output_file = os.path.join(ddir,output_sub,fi)
        #-- extract GRACE month
        grace_month, = np.array(re.findall(output_regex,fi),dtype=np.int)
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

#-- PURPOSE: read ECCO invariant grid file
def ncdf_invariant(invariant_file,**kwargs):
    #-- output dictionary with invariant parameters
    invariant = {}
    #-- open netCDF4 file for reading
    with netCDF4.Dataset(os.path.expanduser(invariant_file),'r') as fileID:
        #-- extract latitude, longitude, depth, area and valid mask
        for key,val in kwargs.items():
            invariant[key] = fileID.variables[val][:].copy()
    #-- return the invariant parameters
    return invariant

#-- PURPOSE: read geoid height netCDF4 files
def ncdf_geoid(FILENAME):
    with netCDF4.Dataset(FILENAME,'r') as fileID:
        geoid_undulation = fileID.variables['geoid'][:,:,:].copy()
    return geoid_undulation

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

#-- Main program that calls ecco_llc_tile_harmonics()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Reads monthly ECCO ocean bottom pressure
            anomalies from LLC tiles and converts to spherical harmonic
            coefficients
            """
    )
    #-- command line parameters
    parser.add_argument('model',
        metavar='MODEL', type=str, nargs='+',
        default=['V5alpha'],
        choices=['V4r4','V5alpha'],
        help='ECCO Version 4 or 5 Model')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- years to run
    now = gravity_toolkit.time.datetime.datetime.now()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.year+1),
        help='Years of model outputs to run')
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
    #-- Output data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Output data format')
    #-- print information about each input and output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files')
    args = parser.parse_args()

    #-- for each ECCO model
    for MODEL in args.model:
        #-- run program
        ecco_llc_tile_harmonics(args.directory, MODEL, args.year,
            LMAX=args.lmax, MMAX=args.mmax, LOVE_NUMBERS=args.love,
            REFERENCE=args.reference, DATAFORM=args.format,
            VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
