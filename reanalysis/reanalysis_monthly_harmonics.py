#!/usr/bin/env python
u"""
reanalysis_monthly_harmonics.py
Written by Tyler Sutterley (03/2021)
Reads atmospheric surface pressure fields from reanalysis and calculates sets of
    spherical harmonics using a thin-layer 2D spherical geometry

INPUTS:
    Reanalysis model to run
    ERA-Interim: http://apps.ecmwf.int/datasets/data/interim-full-moda
    ERA5: http://apps.ecmwf.int/data-catalogues/era5/?class=ea
    MERRA-2: https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/
    NCEP-DOE-2: https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis2.html
    NCEP-CFSR: https://rda.ucar.edu/datasets/ds093.1/
    JRA-55: http://jra.kishou.go.jp/JRA-55/index_en.html

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -Y X, --year X: Years of model outputs to run
    --mean X: start and end year for mean
    --redistribute: uniformly redistribute pressure over oceanic values
    -l X, --lmax=X: maximum spherical harmonic degree
    -m X, --mmax=X: maximum spherical harmonic order
    -n X, --love X: Load Love numbers dataset
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    -r X, --reference X: Reference frame for load love numbers
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
    time.py: utilities for calculating time operations
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
    Updated 03/2021: automatically update years to run based on current time
    Updated 01/2021: harmonics object output from gen_stokes.py
    Updated 12/2020: using argparse to set command line options
        using time module for operations and for extracting time units
    Updated 05/2020: use harmonics class for spherical harmonic operations
    Updated 04/2020: set path to load love numbers file
    Updated 01/2020: iterate over dates to calculate for incomplete files
    Updated 10/2019: changing Y/N flags to True/False
    Updated 09/2019: modified regular expression pattern for MERRA-2
    Updated 08/2019: added parameters for NCEP-CFSR, time scale for MERRA-2
    Updated 07/2018: added parameters for ERA5.  added find_new_files function
    Updated 05/2018: added uniform redistribution of oceanic values
    Updated 03/2018: added portions to run different reanalysis model outputs
        simplified love number extrapolation if LMAX is greater than 696
        create an index file for use in least squares mascon program
    Written 03/2018
"""
from __future__ import print_function

import sys
import os
import re
import netCDF4
import argparse
import numpy as np
import gravity_toolkit.time
import gravity_toolkit.harmonics
from gravity_toolkit.read_love_numbers import read_love_numbers
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.gen_stokes import gen_stokes
from gravity_toolkit.utilities import get_data_path
from geoid_toolkit.ref_ellipsoid import ref_ellipsoid

#-- PURPOSE: read atmospheric surface pressure fields and convert to harmonics
def reanalysis_monthly_harmonics(base_dir, MODEL, YEARS, RANGE=None,
    REDISTRIBUTE=False, LMAX=0, MMAX=None, LOVE_NUMBERS=0, REFERENCE=None,
    DATAFORM=None, VERBOSE=False, MODE=0o775):
    #-- directory setup
    ddir = os.path.join(base_dir,MODEL)

    #-- set model specific parameters
    if (MODEL == 'ERA-Interim'):
        #-- mean file from calculate_mean_pressure.py
        input_mean_file = 'ERA-Interim-Mean-SP-{0:4d}-{1:4d}.nc'
        #-- input land-sea mask for ocean redistribution
        input_mask_file = 'ERA-Interim-Invariant-Parameters.nc'
        #-- regular expression pattern for finding files
        regex_pattern = r'ERA\-Interim\-Monthly\-SP\-({0})\.nc$'
        VARNAME = 'sp'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'time'
        ELLIPSOID = 'WGS84'
        #-- land-sea mask variable name and value of oceanic points
        MASKNAME = 'lsm'
        OCEAN = 0
        GRAVITY = 9.80665
    elif (MODEL == 'ERA5'):
        #-- mean file from calculate_mean_pressure.py
        input_mean_file = 'ERA5-Mean-SP-{0:4d}-{1:4d}.nc'
        #-- input land-sea mask for ocean redistribution
        input_mask_file = 'ERA5-Invariant-Parameters.nc'
        #-- regular expression pattern for finding files
        regex_pattern = r'ERA5\-Monthly\-SP\-({0})\.nc$'
        VARNAME = 'sp'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'time'
        ELLIPSOID = 'WGS84'
        #-- land-sea mask variable name and value of oceanic points
        MASKNAME = 'lsm'
        OCEAN = 0
        GRAVITY = 9.80665
    elif (MODEL == 'MERRA-2'):
        #-- mean file from calculate_mean_pressure.py
        input_mean_file = 'MERRA2.Mean_PS.{0:4d}-{1:4d}.nc'
        #-- input land-sea mask for ocean redistribution
        input_mask_file = 'MERRA2_101.const_2d_asm_Nx.00000000.nc4'
        #-- regular expression pattern for finding files
        regex_pattern = r'MERRA2_\d{{3}}.tavgM_2d_slv_Nx.({0})(\d{{2}}).(.*?).nc$'
        VARNAME = 'PS'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'
        ELLIPSOID = 'WGS84'
        #-- land-sea mask variable name and value of oceanic points
        MASKNAME = 'FROCEAN'
        OCEAN = 1
        GRAVITY = 9.80665
    elif (MODEL == 'NCEP-DOE-2'):
        #-- mean file from calculate_mean_pressure.py
        input_mean_file = 'pres.sfc.mean.{0:4d}-{1:4d}.nc'
        #-- input land-sea mask for ocean redistribution
        input_mask_file = 'land.nc'
        #-- regular expression pattern for finding files
        regex_pattern = r'pres.sfc.mon.mean.({0}).nc$'
        VARNAME = 'pres'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'
        ELLIPSOID = 'WGS84'
        #-- land-sea mask variable name and value of oceanic points
        MASKNAME = 'land'
        OCEAN = 0
        #-- NCEP-DOE-2 reanalysis geopotential heights are already in meters
        GRAVITY = 1.0
    elif (MODEL == 'NCEP-CFSR'):
        #-- mean file from calculate_mean_pressure.py
        input_mean_file = 'pgbh.mean.gdas.{0:4d}-{1:4d}.nc'
        #-- input land-sea mask for ocean redistribution
        input_mask_file = 'land.gdas.nc'
        #-- regular expression pattern for finding files
        regex_pattern = r'pgbh.gdas.({0}).nc$'
        VARNAME = 'PRES_L1_Avg'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'
        ELLIPSOID = 'WGS84'
        #-- land-sea mask variable name and value of oceanic points
        MASKNAME = 'LAND_L1'
        OCEAN = 0
        #-- NCEP-CFSR reanalysis geopotential heights are already in meters
        GRAVITY = 1.0
    elif (MODEL == 'JRA-55'):
        #-- mean file from calculate_mean_pressure.py
        input_mean_file = 'anl_surf.001_pres.mean.{0:4d}-{1:4d}.nc'
        #-- input land-sea mask for ocean redistribution
        input_mask_file = 'll125.081_land.2000.nc'
        #-- regular expression pattern for finding files
        regex_pattern = r'anl_surf125\.001_pres\.({0}).nc$'
        VARNAME = 'Pressure_surface'
        LONNAME = 'g0_lon_1'
        LATNAME = 'g0_lat_0'
        TIMENAME = 'time'
        ELLIPSOID = 'WGS84'
        #-- land-sea mask variable name and value of oceanic points
        MASKNAME = 'LAND_GDS0_SFC'
        OCEAN = 0
        GRAVITY = 9.80665

    #-- upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    #-- output string for both LMAX == MMAX and LMAX != MMAX cases
    order_str = 'M{0:d}'.format(MMAX) if (MMAX != LMAX) else ''
    #-- if redistributing oceanic values to a mean value
    ocean_str = '_OCN' if REDISTRIBUTE else ''
    #-- output suffix for data formats
    suffix = dict(ascii='txt',netCDF4='nc',HDF5='H5')
    #-- output subdirectory
    args = (MODEL.upper(),LMAX,order_str,ocean_str)
    output_sub = '{0}_CLM_L{1:d}{2}{3}'.format(*args)
    if not os.access(os.path.join(ddir,output_sub),os.F_OK):
        os.makedirs(os.path.join(ddir,output_sub))

    #-- read mean pressure field from calculate_mean_pressure.py
    mean_file = os.path.join(ddir,input_mean_file.format(RANGE[0],RANGE[1]))
    mean_pressure,lon,lat=ncdf_mean_pressure(mean_file,VARNAME,LONNAME,LATNAME)
    nlat,nlon = np.shape(mean_pressure)
    #-- calculate colatitude
    theta = (90.0 - lat)*np.pi/180.0
    #-- calculate meshgrid from latitude and longitude
    gridlon,gridlat = np.meshgrid(lon,lat)
    gridphi = gridlon*np.pi/180.0
    gridtheta = (90.0 - gridlat)*np.pi/180.0

    #-- read load love numbers and calculate Legendre polynomials
    LOVE = load_love_numbers(LMAX,LOVE_NUMBERS=LOVE_NUMBERS,REFERENCE=REFERENCE)
    PLM,dPLM = plm_holmes(LMAX,np.cos(theta))

    #-- Earth Parameters
    ellipsoid_params = ref_ellipsoid(ELLIPSOID)
    #-- semimajor and semiminor axes of ellipsoid [m]
    a_axis = ellipsoid_params['a']
    b_axis = ellipsoid_params['b']

    #-- step size in radians
    gridstep = np.zeros((2))
    gridstep[0] = np.abs(lon[1] - lon[0])
    gridstep[1] = np.abs(lat[1] - lat[0])
    dphi = np.pi*gridstep[0]/180.0
    dth = np.pi*gridstep[1]/180.0
    #-- calculate grid areas globally
    AREA = dphi*dth*np.sin(gridtheta)*np.sqrt((a_axis**2)*(b_axis**2) *
        ((np.sin(gridtheta)**2)*(np.cos(gridphi)**2) +
        (np.sin(gridtheta)**2)*(np.sin(gridphi)**2)) +
        (a_axis**4)*(np.cos(gridtheta)**2))

    #-- get indices of land-sea mask if redistributing oceanic points
    if REDISTRIBUTE:
        ii,jj = ncdf_landmask(os.path.join(ddir,input_mask_file),MASKNAME,OCEAN)
        #-- calculate total area of oceanic points
        TOTAL_AREA = np.sum(AREA[ii,jj])

    #-- read each reanalysis pressure field and convert to spherical harmonics
    regex_years = r'\d{4}' if (YEARS is None) else '|'.join(map(str,YEARS))
    rx = re.compile(regex_pattern.format(regex_years),re.VERBOSE)
    input_files = sorted([fi for fi in os.listdir(ddir) if rx.match(fi)])
    #-- open output date and index files
    output_date_file = '{0}_DATES.txt'.format(MODEL.upper())
    fid1 = open(os.path.join(ddir,output_sub,output_date_file),'w')
    output_index_file = '{0}_index.txt'.format(MODEL)
    fid2 = open(os.path.join(ddir,output_sub,output_index_file),'w')
    #-- date file header information
    print('{0:8} {1:10}'.format('Month','Date'), file=fid1)
    #-- output file format for spherical harmonic data
    output_file_format = '{0}_CLM_L{1:d}{2}_{3:03d}.{4}'

    #-- for each reanalysis file
    for fi in input_files:
        #-- read input data
        with netCDF4.Dataset(os.path.join(ddir,fi),'r') as fileID:
            pressure=np.copy(fileID.variables[VARNAME][:])
            #-- convert time to Modified Julian Days
            delta_time=np.copy(fileID.variables[TIMENAME][:])
            date_string=fileID.variables[TIMENAME].units
            epoch,to_secs=gravity_toolkit.time.parse_date_string(date_string)
            MJD=gravity_toolkit.time.convert_delta_time(delta_time*to_secs,
                epoch1=epoch, epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)
        #-- iterate over Julian days
        for t,JD in enumerate(MJD+2400000.5):
            #-- calculate pressure/gravity ratio for month
            PG = (pressure[t,:,:] - mean_pressure[:,:])/GRAVITY
            #-- if redistributing oceanic pressure values
            if REDISTRIBUTE:
                PG[ii,jj] = np.sum(PG[ii,jj]*AREA[ii,jj])/TOTAL_AREA
            #-- calculate pressure harmonics from pressure/gravity ratio
            Ylms = gen_stokes(PG, lon, lat, LMAX=LMAX, MMAX=MMAX, UNITS=3,
                PLM=PLM, LOVE=LOVE)
            #-- convert julian dates to calendar then to year-decimal
            YY,MM,DD,hh,mm,ss = gravity_toolkit.time.convert_julian(JD,
                FORMAT='tuple')
            Ylms.time,=gravity_toolkit.time.convert_calendar_decimal(YY,
                MM, day=DD, hour=hh, minute=mm, second=ss)
            #-- calculate GRACE month from calendar dates
            Ylms.month, = np.array([(YY - 2002)*12 + MM], dtype=np.int)
            #-- output data to file
            args = (MODEL.upper(),LMAX,order_str,Ylms.month,suffix[DATAFORM])
            FILE = output_file_format.format(*args)
            #-- output data for month
            print(os.path.join(ddir,output_sub,FILE)) if VERBOSE else None
            if (DATAFORM == 'ascii'):
                #-- ascii (.txt)
                Ylms.to_ascii(os.path.join(ddir,output_sub,FILE))
            elif (DATAFORM == 'netCDF4'):
                #-- netcdf (.nc)
                Ylms.to_netCDF4(os.path.join(ddir,output_sub,FILE))
            elif (DATAFORM == 'HDF5'):
                #-- HDF5 (.H5)
                Ylms.to_HDF5(os.path.join(ddir,output_sub,FILE))
            #-- set the permissions level of the output file to MODE
            os.chmod(os.path.join(ddir,output_sub,FILE), MODE)

    #-- output file format for spherical harmonic data
    args = (MODEL.upper(),LMAX,order_str,suffix[DATAFORM])
    output_regex = re.compile(r'{0}_CLM_L{1:d}{2}_(\d+).{3}'.format(*args))
    #-- find all output harmonic files (not just ones created in run)
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

#-- PURPOSE: read reanalysis mean pressure from calculate_mean_pressure.py
def ncdf_mean_pressure(FILENAME,VARNAME,LONNAME,LATNAME):
    with netCDF4.Dataset(FILENAME,'r') as fileID:
        mean_pressure = np.array(fileID.variables[VARNAME][:].squeeze())
        longitude = fileID.variables[LONNAME][:].squeeze()
        latitude = fileID.variables[LATNAME][:].squeeze()
    return (mean_pressure,longitude,latitude)

#-- PURPOSE: read land sea mask to get indices of oceanic values
def ncdf_landmask(FILENAME,MASKNAME,OCEAN):
    with netCDF4.Dataset(FILENAME,'r') as fileID:
        landsea = np.squeeze(fileID.variables[MASKNAME][:].copy())
    return np.nonzero(landsea == OCEAN)

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

#-- Main program that calls reanalysis_monthly_harmonics()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Reads atmospheric surface pressure fields from
            reanalysis and calculates sets of spherical harmonics
            using a thin-layer 2D spherical geometry
            """
    )
    #-- command line parameters
    choices = ['ERA-Interim','ERA5','MERRA-2','NCEP-DOE-2','NCEP-CFSR','JRA-55']
    parser.add_argument('model',
        metavar='MODEL', type=str, nargs='+',
        default=['ERA5','MERRA-2'], choices=choices,
        help='Reanalysis Model')
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
    #-- start and end years to run for mean
    parser.add_argument('--mean',
        metavar=('START','END'), type=int, nargs=2,
        default=[2001,2002],
        help='Start and end year range for mean')
    #-- uniformly redistribute pressure values over the ocean
    parser.add_argument('--redistribute',
        default=False, action='store_true',
        help='Redistribute pressure values over the ocean')
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

    #-- for each reanalysis model
    for MODEL in args.model:
        #-- run program
        reanalysis_monthly_harmonics(args.directory, MODEL, args.year,
            RANGE=args.mean, REDISTRIBUTE=args.redistribute, LMAX=args.lmax,
            MMAX=args.mmax, LOVE_NUMBERS=args.love, REFERENCE=args.reference,
            DATAFORM=args.format, VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
