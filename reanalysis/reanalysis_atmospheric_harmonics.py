#!/usr/bin/env python
u"""
reanalysis_atmospheric_harmonics.py
Written by Tyler Sutterley (12/2022)
Reads atmospheric geopotential heights fields from reanalysis and calculates
    sets of spherical harmonics using a 3D geometry

INPUTS:
    Reanalysis model to run
    ERA-Interim: http://apps.ecmwf.int/datasets/data/interim-full-moda
    ERA5: http://apps.ecmwf.int/data-catalogues/era5/?class=ea
    MERRA-2: https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -Y X, --year X: Years of model outputs to run
    --mean X: start and end year for mean
    --redistribute: uniformly redistribute values over the ocean
    -l X, --lmax=X: maximum spherical harmonic degree
    -m X, --mmax=X: maximum spherical harmonic order
    -n X, --love X: Load Love numbers dataset
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -F X, --format X: output data format
        ascii
        netCDF4
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
    associated_legendre.py: computes fully-normalized associated Legendre polynomials
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    gen_atmosphere_stokes.py: converts atmospheric fields to spherical harmonics
    constants.py: calculate reference parameters for common ellipsoids
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    time.py: utilities for calculating time operations
    units.py: class for converting spherical harmonic data to specific units
    utilities.py: download and management utilities for files

REFERENCES:
    JP Boy and B Chao, "Precise evaluation of atmospheric loading effects on
        Earth's time-variable gravity field", Journal of Geophysical Research:
        Solid Earth, 110(B8), (2005).
        https://doi.org/10.1029/2002JB002333

    S Swenson and J Wahr, "Estimated effects of the vertical structure of
        atmospheric mass on the time-variable geoid", Journal of Geophysical
        Research: Solid Earth, 107(B9), (2002).
        https://doi.org/10.1029/2000JB000024

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
        use constants class in place of geoid-toolkit ref_ellipsoid
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: use wrapper function for reading load Love numbers
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
        use output harmonic file wrapper routine to write to file
    Updated 09/2021: use GRACE/GRACE-FO month to calendar month converters
    Updated 07/2021: can use input files to define command line arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 03/2021: automatically update years to run based on current time
    Updated 01/2021: read from netCDF4 file in slices to reduce memory load
        separated gen_atmosphere_stokes to a separate function
    Updated 12/2020: using argparse to set command line options
        using time module for operations and for extracting time units
    Updated 05/2020: use harmonics class for spherical harmonic operations
    Updated 04/2020: set path to load love numbers file
    Updated 01/2020: iterate over dates to calculate for incomplete files
    Updated 10/2019: changing Y/N flags to True/False
    Updated 08/2019: adjust time scale variable for MERRA-2
    Updated 07/2018: added parameters for ERA5.  added find_new_files function
    Updated 05/2018: added uniform redistribution of oceanic values
    Updated 03/2018: simplified love number extrapolation if LMAX > 696
    Written 03/2018
"""
from __future__ import print_function

import sys
import os
import re
import logging
import netCDF4
import argparse
import datetime
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc

# PURPOSE: read atmospheric 3D geopotential height fields
# converts model geopotential heights to spherical harmonics
def reanalysis_atmospheric_harmonics(base_dir, MODEL, YEARS, RANGE=None,
    REDISTRIBUTE=False, LMAX=0, MMAX=None, LOVE_NUMBERS=0, REFERENCE=None,
    DATAFORM=None, MODE=0o775):

    # directory setup
    ddir = os.path.join(base_dir,MODEL)

    # set model specific parameters
    if (MODEL == 'ERA-Interim'):
        # invariant parameters file
        input_invariant_file = 'ERA-Interim-Invariant-Parameters.nc'
        # geoid file from read_gfz_geoid_grids.py
        input_geoid_file = 'ERA-Interim-EGM2008-geoid.nc'
        # input land-sea mask for ocean redistribution
        input_mask_file = 'ERA-Interim-Invariant-Parameters.nc'
        # regular expression pattern for finding files
        # calculated from calculate_geopotential_heights.py
        regex_pattern = r'ERA\-Interim\-GPH\-Levels\-({0})\.nc$'
        ZNAME = 'z'
        DIFFNAME = 'dp'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'time'
        LEVELNAME = 'lvl'
        ELLIPSOID = 'WGS84'
        # land-sea mask variable name and value of oceanic points
        MASKNAME = 'lsm'
        OCEAN = 0
        GRAVITY = 9.80665
    elif (MODEL == 'ERA5'):
        # invariant parameters file
        input_invariant_file = 'ERA5-Invariant-Parameters.nc'
        # geoid file from read_gfz_geoid_grids.py
        input_geoid_file = 'ERA5-EGM2008-geoid.nc'
        # input land-sea mask for ocean redistribution
        input_mask_file = 'ERA5-Invariant-Parameters.nc'
        # regular expression pattern for finding files
        # calculated from calculate_geopotential_heights.py
        regex_pattern = r'ERA5\-GPH\-Levels\-({0})\.nc$'
        ZNAME = 'z'
        DIFFNAME = 'dp'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'time'
        LEVELNAME = 'lvl'
        ELLIPSOID = 'WGS84'
        # land-sea mask variable name and value of oceanic points
        MASKNAME = 'lsm'
        OCEAN = 0
        GRAVITY = 9.80665
    elif (MODEL == 'MERRA-2'):
        # invariant parameters file
        input_invariant_file = 'MERRA2_101.const_2d_asm_Nx.00000000.nc4'
        # geoid file form read_gfz_geoid_grids.py
        input_geoid_file = 'MERRA2_101.EGM2008_Nx.00000000.nc4'
        # input land-sea mask for ocean redistribution
        input_mask_file = 'MERRA2_101.const_2d_asm_Nx.00000000.nc4'
        # regular expression pattern for finding files
        # calculated from calculate_geopotential_heights.py
        regex_pattern = r'MERRA2_\d{{3}}.GPH_levels.({0})(\d{{2}}).SUB.nc$'
        ZNAME = 'PHIS'
        DIFFNAME = 'dP'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'
        LEVELNAME = 'lev'
        ELLIPSOID = 'WGS84'
        # land-sea mask variable name and value of oceanic points
        MASKNAME = 'FROCEAN'
        OCEAN = 1
        GRAVITY = 9.80665

    # upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    # output string for both LMAX == MMAX and LMAX != MMAX cases
    order_str = 'M{MMAX:d}' if (MMAX != LMAX) else ''
    # if redistributing oceanic values to a mean value
    ocean_str = '_OCN' if REDISTRIBUTE else ''
    # output suffix for data formats
    suffix = dict(ascii='txt',netCDF4='nc',HDF5='H5')
    # output subdirectory
    args = (MODEL.upper(),LMAX,order_str,ocean_str)
    output_sub = '{0}_ATMOSPHERE_CLM_L{1:d}{2}{3}'.format(*args)
    if not os.access(os.path.join(ddir,output_sub),os.F_OK):
        os.makedirs(os.path.join(ddir,output_sub))
    # attributes for output files
    attributes = {}
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'

    # read model latitude and longitude from invariant parameters file
    with netCDF4.Dataset(os.path.join(ddir,input_invariant_file),'r') as fileID:
        lon = fileID.variables[LONNAME][:].copy()
        lat = fileID.variables[LATNAME][:].copy()
    # calculate colatitude
    theta = (90.0 - lat)*np.pi/180.0
    # calculate meshgrid from latitude and longitude
    gridlon,gridlat = np.meshgrid(lon,lat)
    gridphi = gridlon*np.pi/180.0
    gridtheta = (90.0 - gridlat)*np.pi/180.0

    # read load love numbers and calculate Legendre polynomials
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)
    PLM, dPLM = gravtk.plm_holmes(LMAX, np.cos(theta))
    # read geoid heights and grid step size
    geoid,gridstep = ncdf_geoid(os.path.join(ddir,input_geoid_file))
    nlat,nlon = np.shape(geoid)

    # get reference parameters for ellipsoid
    ellipsoid_params = mdlhmc.constants(ellipsoid=ELLIPSOID)
    # semimajor and semiminor axes of the ellipsoid [m]
    a_axis = ellipsoid_params.a_axis
    b_axis = ellipsoid_params.b_axis

    # step size in radians
    if (np.ndim(gridstep) == 0):
        dphi = np.pi*gridstep/180.0
        dth = np.pi*gridstep/180.0
    else: # dlon ne dlat
        dphi = np.pi*gridstep[0]/180.0
        dth = np.pi*gridstep[1]/180.0
    # calculate grid areas globally
    AREA = dphi*dth*np.sin(gridtheta)*np.sqrt((a_axis**2)*(b_axis**2) *
        ((np.sin(gridtheta)**2)*(np.cos(gridphi)**2) +
        (np.sin(gridtheta)**2)*(np.sin(gridphi)**2)) +
        (a_axis**4)*(np.cos(gridtheta)**2))

    # get indices of land-sea mask if redistributing oceanic points
    if REDISTRIBUTE:
        ii,jj = ncdf_landmask(os.path.join(ddir,input_mask_file),MASKNAME,OCEAN)
        # calculate total area of oceanic points
        TOTAL_AREA = np.sum(AREA[ii,jj])

    # read each reanalysis pressure field and convert to spherical harmonics
    regex_years = r'\d{4}' if (YEARS is None) else '|'.join(map(str,YEARS))
    rx = re.compile(regex_pattern.format(regex_years),re.VERBOSE)
    input_files = sorted([fi for fi in os.listdir(ddir) if rx.match(fi)])
    # open output date and index files
    output_date_file = f'{MODEL.upper()}_DATES.txt'
    fid1 = open(os.path.join(ddir,output_sub,output_date_file),
        mode='w', encoding='utf8')
    output_index_file = f'{MODEL}_index.txt'
    fid2 = open(os.path.join(ddir,output_sub,output_index_file),
        mode='w', encoding='utf8')
    # date file header information
    print('{0:8} {1:10}'.format('Month','Date'), file=fid1)
    # output file format for spherical harmonic data
    output_file_format = '{0}_CLM_L{1:d}{2}_{3:03d}.{4}'

    # read mean spherical harmonics
    args = (MODEL.upper(),LMAX,order_str,RANGE[0],RANGE[1],suffix[DATAFORM])
    mean_file = '{0}_MEAN_CLM_L{1:d}{2}_{3:4d}-{4:4d}.{5}'.format(*args)
    mean_Ylms = gravtk.harmonics().from_file(
        os.path.join(ddir,output_sub,mean_file),
        format=DATAFORM)
    # truncating mean spherical harmonics to d/o LMAX/MMAX
    mean_Ylms = mean_Ylms.truncate(lmax=LMAX,mmax=MMAX)

    # read each reanalysis data file and convert to spherical harmonics
    for fi in input_files:
        # read model level geopotential height data
        fileID = netCDF4.Dataset(os.path.join(ddir,fi),'r')
        # extract shape from geopotential variable
        ntime,nlevels,nlat,nlon = fileID.variables[ZNAME].shape
        # convert time to Modified Julian Days
        delta_time = np.copy(fileID.variables[TIMENAME][:])
        date_string = fileID.variables[TIMENAME].units
        epoch,to_secs = gravtk.time.parse_date_string(date_string)
        MJD = gravtk.time.convert_delta_time(delta_time*to_secs,
            epoch1=epoch, epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)
        # iterate over Julian days
        for t,JD in enumerate(MJD+2400000.5):
            # convert geopotential to geopotential height
            GPH = np.squeeze(fileID.variables[ZNAME][t,:,:,:])/GRAVITY
            # extract pressure difference for month
            PD = np.squeeze(fileID.variables[DIFFNAME][t,:,:,:])
            # if redistributing oceanic pressure values
            if REDISTRIBUTE:
                for p in range(nlevels):
                    PD[p,ii,jj]=np.sum(PD[p,ii,jj]*AREA[ii,jj])/TOTAL_AREA
            # calculate spherical harmonics for month
            Ylms = mdlhmc.gen_atmosphere_stokes(GPH, PD, lon, lat,
                LMAX=LMAX, MMAX=MMAX, ELLIPSOID=ELLIPSOID, GEOID=geoid,
                PLM=PLM, LOVE=LOVE)
            # remove mean harmonics from month
            Ylms.subtract(mean_Ylms)
            # convert julian dates to calendar then to year-decimal
            YY,MM,DD,hh,mm,ss = gravtk.time.convert_julian(JD,
                FORMAT='tuple')
            Ylms.time, = gravtk.time.convert_calendar_decimal(YY,
                MM, day=DD, hour=hh, minute=mm, second=ss)
            # calculate GRACE month from calendar dates
            Ylms.month = gravtk.time.calendar_to_grace(YY, MM)
            # output data to file
            args = (MODEL.upper(),LMAX,order_str,Ylms.month,suffix[DATAFORM])
            FILE = output_file_format.format(*args)
            Ylms.to_file(os.path.join(ddir,output_sub,FILE),
                format=DATAFORM, **attributes)
            # set the permissions level of the output file to MODE
            os.chmod(os.path.join(ddir,output_sub,FILE), MODE)
        # close the input netCDF4 file
        fileID.close()

    # output file format for spherical harmonic data
    args = (MODEL.upper(),LMAX,order_str,suffix[DATAFORM])
    output_regex = re.compile(r'{0}_CLM_L{1:d}{2}_(\d+).{3}'.format(*args))
    # find all output harmonic files (not just ones created in run)
    output_files = [fi for fi in os.listdir(os.path.join(ddir,output_sub))
        if re.match(output_regex,fi)]
    for fi in sorted(output_files):
        # extract GRACE month
        grace_month, = np.array(re.findall(output_regex,fi),dtype=np.int64)
        YY,MM = gravtk.time.grace_to_calendar(grace_month)
        tdec, = gravtk.time.convert_calendar_decimal(YY, MM)
        # full path to output file
        full_output_file = os.path.join(ddir,output_sub,fi)
        # print date, GRACE month and calendar month to date file
        fid1.write('{0:11.6f} {1:03d} {2:02.0f}\n'.format(tdec,grace_month,MM))
        # print output file to index
        print(full_output_file.replace(os.path.expanduser('~'),'~'), file=fid2)
    # close the date and index files
    fid1.close()
    fid2.close()
    # set the permissions level of the output date and index files to MODE
    os.chmod(os.path.join(ddir,output_sub,output_date_file), MODE)
    os.chmod(os.path.join(ddir,output_sub,output_index_file), MODE)

# PURPOSE: read geoid height netCDF4 files from read_gfz_geoid_grids.py
def ncdf_geoid(FILENAME):
    with netCDF4.Dataset(FILENAME,'r') as fileID:
        geoid_undulation = fileID.variables['geoid'][:].copy()
        gridstep = np.array(fileID.gridstep.split(','),dtype=np.float64)
    return (geoid_undulation,np.squeeze(gridstep))

# PURPOSE: read land sea mask to get indices of oceanic values
def ncdf_landmask(FILENAME,MASKNAME,OCEAN):
    with netCDF4.Dataset(FILENAME,'r') as fileID:
        landsea = np.squeeze(fileID.variables[MASKNAME][:].copy()).astype('f2')
    return np.nonzero(landsea == OCEAN)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads atmospheric geopotential heights
            fields from reanalysis and calculates sets of
            spherical harmonics using a 3D geometry
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    choices = ['ERA-Interim','ERA5','MERRA-2']
    parser.add_argument('model',
        type=str, nargs='+',
        default=['ERA5','MERRA-2'], choices=choices,
        help='Reanalysis Model')
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # years to run
    now = datetime.datetime.now()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.year+1),
        help='Years of model outputs to run')
    # start and end years to run for mean
    parser.add_argument('--mean',
        metavar=('START','END'), type=int, nargs=2,
        default=[2001,2002],
        help='Start and end year range for mean')
    # uniformly redistribute pressure values over the ocean
    parser.add_argument('--redistribute',
        default=False, action='store_true',
        help='Redistribute pressure values over the ocean')
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
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2],
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
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # for each reanalysis model
    for MODEL in args.model:
        # run program
        reanalysis_atmospheric_harmonics(args.directory, MODEL, args.year,
            RANGE=args.mean, REDISTRIBUTE=args.redistribute, LMAX=args.lmax,
            MMAX=args.mmax, LOVE_NUMBERS=args.love, REFERENCE=args.reference,
            DATAFORM=args.format, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
