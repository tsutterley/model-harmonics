#!/usr/bin/env python
u"""
reanalysis_mean_harmonics.py
Written by Tyler Sutterley (12/2020)
Reads atmospheric geopotential heights fields from reanalysis and calculates
    a multi-annual mean set of spherical harmonics using a 3D geometry

INPUTS:
    Reanalysis model to run
    ERA-Interim: http://apps.ecmwf.int/datasets/data/interim-full-moda
    ERA5: http://apps.ecmwf.int/data-catalogues/era5/?class=ea
    MERRA-2: https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    --mean X: start and end year for mean
    --redistribute: uniformly redistribute values over the ocean
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
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
        destripe_harmonics.py: calculates the decorrelation (destriping) filter
            and filters the GRACE/GRACE-FO coefficients for striping errors
        ncdf_read_stokes.py: reads spherical harmonic netcdf files
        ncdf_stokes.py: writes output spherical harmonic data to netcdf
        hdf5_read_stokes.py: reads spherical harmonic HDF5 files
        hdf5_stokes.py: writes output spherical harmonic data to HDF5
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
    Updated 12/2020: using argparse to set command line options
        using time module for operations and for extracting time units
    Updated 05/2020: use harmonics class for spherical harmonic operations
    Updated 04/2020: set path to load love numbers file
    Updated 01/2020: use EGM2008 geoid with ERA-Interim for consistency
        iterate over Julian days to calculate harmonics for consistency
    Updated 10/2019: changing Y/N flags to True/False
    Updated 08/2019: adjust time scale variable for MERRA-2
    Updated 09/2018: added common land-sea mask from create_common_masks.py
    Updated 07/2018: added parameters for ERA5
    Updated 05/2018: added uniform redistribution of oceanic values
    Updated 03/2018: simplified love number extrapolation if LMAX > 696
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
import gravity_toolkit.units
import gravity_toolkit.harmonics
from gravity_toolkit.read_love_numbers import read_love_numbers
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.utilities import get_data_path
from geoid_toolkit.ref_ellipsoid import ref_ellipsoid

#-- PURPOSE: read atmospheric surface pressure fields and convert to harmonics
def reanalysis_mean_harmonics(base_dir, MODEL, RANGE=None, REDISTRIBUTE=False,
    LMAX=0, MMAX=None,  LOVE_NUMBERS=0, REFERENCE=None, DATAFORM=None,
    VERBOSE=False, MODE=0o775):
    #-- directory setup
    ddir = os.path.join(base_dir,MODEL)

    #-- set model specific parameters
    if (MODEL == 'ERA-Interim'):
        #-- invariant parameters file
        input_invariant_file = 'ERA-Interim-Invariant-Parameters.nc'
        #-- geoid file from read_gfz_geoid_grids.py
        input_geoid_file = 'ERA-Interim-EGM2008-geoid.nc'
        #-- input land-sea mask for ocean redistribution
        input_mask_file = 'ERA-Interim-Invariant-Parameters.nc'
        #-- regular expression pattern for finding files
        #-- calculated from calculate_geopotential_heights.py
        regex_pattern = r'ERA\-Interim\-GPH\-Levels\-({0})\.nc$'
        ZNAME = 'z'
        DIFFNAME = 'dp'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'time'
        LEVELNAME = 'lvl'
        ELLIPSOID = 'WGS84'
        #-- land-sea mask variable name and value of oceanic points
        MASKNAME = 'lsm'
        OCEAN = 0
        GRAVITY = 9.80665
    elif (MODEL == 'ERA5'):
        #-- invariant parameters file
        input_invariant_file = 'ERA5-Invariant-Parameters.nc'
        #-- geoid file from read_gfz_geoid_grids.py
        input_geoid_file = 'ERA5-EGM2008-geoid.nc'
        #-- input land-sea mask for ocean redistribution
        input_mask_file = 'ERA5-Invariant-Parameters.nc'
        #-- regular expression pattern for finding files
        #-- calculated from calculate_geopotential_heights.py
        regex_pattern = r'ERA5\-GPH\-Levels\-({0})\.nc$'
        ZNAME = 'z'
        DIFFNAME = 'dp'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'time'
        LEVELNAME = 'lvl'
        ELLIPSOID = 'WGS84'
        #-- land-sea mask variable name and value of oceanic points
        MASKNAME = 'lsm'
        OCEAN = 0
        GRAVITY = 9.80665
    elif (MODEL == 'MERRA-2'):
        #-- invariant parameters file
        input_invariant_file = 'MERRA2_101.const_2d_asm_Nx.00000000.nc4'
        #-- geoid file form read_gfz_geoid_grids.py
        input_geoid_file = 'MERRA2_101.EGM2008_Nx.00000000.nc4'
        #-- input land-sea mask for ocean redistribution
        input_mask_file = 'MERRA2_101.const_2d_asm_Nx.00000000.nc4'
        #-- regular expression pattern for finding files
        #-- calculated from calculate_geopotential_heights.py
        regex_pattern = r'MERRA2_\d{{3}}.GPH_levels.({0})(\d{{2}}).SUB.nc$'
        ZNAME = 'PHIS'
        DIFFNAME = 'dP'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'
        LEVELNAME = 'lev'
        ELLIPSOID = 'WGS84'
        #-- land-sea mask variable name and value of oceanic points
        MASKNAME = 'FROCEAN'
        OCEAN = 1
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
    output_sub = '{0}_ATMOSPHERE_CLM_L{1:d}{2}{3}'.format(*args)
    if not os.access(os.path.join(ddir,output_sub),os.F_OK):
        os.makedirs(os.path.join(ddir,output_sub))

    #-- read model latitude and longitude from invariant parameters file
    with netCDF4.Dataset(os.path.join(ddir,input_invariant_file),'r') as fileID:
        lon = fileID.variables[LONNAME][:].copy()
        lat = fileID.variables[LATNAME][:].copy()
    #-- calculate colatitude
    theta = (90.0 - lat)*np.pi/180.0
    #-- calculate meshgrid from latitude and longitude
    gridlon,gridlat = np.meshgrid(lon,lat)
    gridphi = gridlon*np.pi/180.0
    gridtheta = (90.0 - gridlat)*np.pi/180.0

    #-- read load love numbers and calculate Legendre polynomials
    LOVE = load_love_numbers(LMAX,LOVE_NUMBERS=LOVE_NUMBERS,REFERENCE=REFERENCE)
    PLM,dPLM = plm_holmes(LMAX,np.cos(theta))
    #-- read geoid heights and grid step size
    geoid,gridstep = ncdf_geoid(os.path.join(ddir,input_geoid_file))
    nlat,nlon = np.shape(geoid)

    #-- Earth Parameters
    ellipsoid_params = ref_ellipsoid(ELLIPSOID)
    #-- semimajor and semiminor axes of ellipsoid [m]
    a_axis = ellipsoid_params['a']
    b_axis = ellipsoid_params['b']

    #-- step size in radians
    if (np.ndim(gridstep) == 0):
        dphi = np.pi*gridstep/180.0
        dth = np.pi*gridstep/180.0
    else: #-- dlon ne dlat
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
        # ii,jj = ncdf_landmask(os.path.join(ddir,common_mask_file),'mask',1)
        #-- calculate total area of oceanic points
        TOTAL_AREA = np.sum(AREA[ii,jj])

    #-- read each reanalysis pressure field and convert to spherical harmonics
    #-- then calculate mean of the spherical harmonic fields
    regex_years = '|'.join(['{0:4d}'.format(Y) for Y in range(SY,EY+1)])
    rx = re.compile(regex_pattern.format(regex_years),re.VERBOSE)
    input_files = sorted([fi for fi in os.listdir(ddir) if rx.match(fi)])

    #-- list of monthly spherical harmonics
    harmonics_list = []
    #-- for each reanalysis file
    for fi in input_files:
        #-- read geopotential height data from calculate_geopotential_heights.py
        with netCDF4.Dataset(os.path.join(ddir,fi),'r') as fileID:
            geopotential = np.array(fileID.variables[ZNAME][:].copy())
            pressure_difference = np.array(fileID.variables[DIFFNAME][:].copy())
            nlevels, = fileID.variables[LEVELNAME][:].shape
            #-- convert time to Modified Julian Days
            delta_time=np.copy(fileID.variables[TIMENAME][:])
            date_string=fileID.variables[TIMENAME].units
            epoch,to_secs=gravity_toolkit.time.parse_date_string(date_string)
            MJD=gravity_toolkit.time.convert_delta_time(delta_time*to_secs,
                epoch1=epoch, epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)
        #-- iterate over Julian days
        for t,JD in enumerate(MJD+2400000.5):
            #-- convert geopotential to geopotential height
            GPH = np.squeeze(geopotential[t,:,:,:])/GRAVITY
            #-- extract pressure difference for month
            PD = np.squeeze(pressure_difference[t,:,:,:])
            #-- if redistributing oceanic pressure values
            if REDISTRIBUTE:
                for p in range(nlevels):
                    PD[p,ii,jj] = np.sum(PD[p,ii,jj]*AREA[ii,jj])/TOTAL_AREA
            #-- calculate spherical harmonics for month
            Ylms=gen_atmosphere_stokes(GPH, PD, lon, lat, LMAX=LMAX, MMAX=MMAX,
                ELLIPSOID=ELLIPSOID, GEOID=geoid, PLM=PLM, LOVE=LOVE)
            #-- convert julian dates to calendar then to year-decimal
            YY,MM,DD,hh,mm,ss = gravity_toolkit.time.convert_julian(JD,
                FORMAT='tuple')
            Ylms.time, = gravity_toolkit.time.convert_calendar_decimal(YY,
                MM, day=DD, hour=hh, minute=mm, second=ss)
            #-- calculate GRACE month from calendar dates
            Ylms.month, = np.array([(YY - 2002)*12 + MM], dtype=np.int)
            #-- append to list of harmonics
            harmonics_list.append(Ylms)

    #-- divide by the number of dates to calculate mean from totals
    mean_Ylms = gravity_toolkit.harmonics().from_list(harmonics_list).mean()

    #-- output mean spherical harmonics file
    args = (MODEL.upper(), LMAX, order_str, RANGE[0], RANGE[1], suffix[DATAFORM])
    output_mean_file = '{0}_MEAN_CLM_L{1:d}{2}_{3:4d}-{4:4d}.{5}'.format(*args)
    #-- output data for month
    if (DATAFORM == 'ascii'):
        #-- ascii (.txt)
        mean_Ylms.to_ascii(os.path.join(ddir,output_sub,output_mean_file))
    elif (DATAFORM == 'netCDF4'):
        #-- netcdf (.nc)
        mean_Ylms.to_netCDF4(os.path.join(ddir,output_sub,output_mean_file))
    elif (DATAFORM == 'HDF5'):
        #-- HDF5 (.H5)
        mean_Ylms.to_HDF5(os.path.join(ddir,output_sub,output_mean_file))
    #-- set the permissions level of the output file to MODE
    os.chmod(os.path.join(ddir,output_sub,output_mean_file), MODE)

#-- PURPOSE: read geoid height netCDF4 files from read_gfz_geoid_grids.py
def ncdf_geoid(FILENAME):
    with netCDF4.Dataset(FILENAME,'r') as fileID:
        geoid_undulation = fileID.variables['geoid'][:].copy()
        gridstep = np.array(fileID.gridstep.split(','),dtype=np.float)
    return (geoid_undulation,np.squeeze(gridstep))

#-- PURPOSE: read land sea mask to get indices of oceanic values
def ncdf_landmask(FILENAME,MASKNAME,OCEAN):
    with netCDF4.Dataset(FILENAME,'r') as fileID:
        landsea = np.squeeze(fileID.variables[MASKNAME][:].copy()).astype('f2')
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

#-- PURPOSE: calculates spherical harmonic fields from atmospheric pressure
def gen_atmosphere_stokes(GPH, pressure, lon, lat, LMAX=0, MMAX=None,
    ELLIPSOID=None, GEOID=None, PLM=None, LOVE=None, METHOD='BC05'):
    #-- converting LMAX to integer
    LMAX = np.int(LMAX)
    #-- upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX

    #-- number of pressure levels, longitudes and latitudes
    nlevels,nlat,nlon = np.shape(GPH)
    #-- grid step
    dlon = np.abs(lon[1]-lon[0])
    dlat = np.abs(lat[1]-lat[0])
    #-- longitude degree spacing in radians
    dphi = dlon*np.pi/180.0
    #-- colatitude degree spacing in radians
    dth = dlat*np.pi/180.0

    #-- calculate longitudes and colatitudes in radians
    phi = lon*np.pi/180.0
    phi = np.squeeze(phi)[np.newaxis,:]
    th = (90.0 - np.squeeze(lat))*np.pi/180.0
    #-- calculate meshgrid from latitude and longitude
    gridlon,gridlat = np.meshgrid(lon,lat)
    gridphi = gridlon*np.pi/180.0
    gridtheta = (90.0 - gridlat)*np.pi/180.0

    #-- Earth Parameters
    ellipsoid_params = ref_ellipsoid(ELLIPSOID)
    #-- semimajor axis of ellipsoid [m]
    a_axis = ellipsoid_params['a']
    #-- ellipsoidal flattening
    flat = ellipsoid_params['f']
    #-- Average Radius of the Earth having the same volume [m]
    rad_e = ellipsoid_params['rad_e']
    #--  first numerical eccentricity
    ecc1 = ellipsoid_params['ecc1']
    #-- convert from geodetic latitude to geocentric latitude
    #-- prime vertical radius of curvature
    N = a_axis/np.sqrt(1.0 - ecc1**2.0*np.cos(gridtheta)**2.0)

    #-- Coefficient for calculating Stokes coefficients from pressure field
    #-- SH Degree dependent factors with indirect loading components
    factors = gravity_toolkit.units(lmax=LMAX,a_axis=100.0*a_axis,flat=flat)
    dfactor = factors.spatial(*LOVE).mmwe
    #-- Multiplying sin(th) with differentials of theta and phi
    #-- to calculate the integration factor at each latitude
    int_fact = np.sin(th)*dphi*dth

    #-- Calculating cos/sin of phi arrays
    #-- output [m,phi]
    m = np.arange(MMAX+1)
    ccos = np.cos(np.dot(m[:, np.newaxis],phi))
    ssin = np.sin(np.dot(m[:, np.newaxis],phi))

    #-- Fully-normalized Legendre Polynomials were computed outside
    #-- Multiplying by the units conversion factor (conv) to
    #-- Multiplying by integration factors [sin(theta)*dtheta*dphi]
    plm = np.zeros((LMAX+1,MMAX+1,nlat))
    for j in range(0,nlat):
        plm[:,m,j] = PLM[:,m,j]*int_fact[j]

    #-- gravitational acceleration at the Earth's mean spherical surface
    g0 = 9.80665
    #-- gravitational acceleration at the equator and at mean sea level
    ge = 9.780356
    #-- gravitational acceleration at the mean sea level over gridtheta
    gs = ge*(1.0+5.2885e-3*np.cos(gridtheta)**2-5.9e-6*np.cos(2.0*gridtheta)**2)
    #-- calculate radii and gravity for each pressure level
    R = np.zeros((nlevels,nlat,nlon))
    gamma_h = np.zeros((nlevels,nlat,nlon))
    for p in range(nlevels):
        #-- orthometric height from List (1958)
        #-- as described in Boy and Chao (2005)
        orthometric = (1.0 - 0.002644*np.cos(2.0*gridtheta))*GPH[p,:,:] + \
            (1.0 - 0.0089*np.cos(2.0*gridtheta))*(GPH[p,:,:]**2)/6.245e6
        #-- calculate X, Y and Z from geodetic latitude and longitude
        X = (N + GEOID + orthometric) * np.sin(gridtheta) * np.cos(gridphi)
        Y = (N + GEOID + orthometric) * np.sin(gridtheta) * np.sin(gridphi)
        Z = (N * (1.0 - ecc1**2.0) + GEOID + orthometric) * np.cos(gridtheta)
        #-- calculate radius of levelZ
        R[p,:,:] = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
        #-- calculate normal gravity at each height above mean sea level
        gamma_h[p,:,:] = gs*(1.0-2.0*(1.006803-0.06706*np.cos(gridtheta)**2)*
            (orthometric/rad_e) + 3.0*(orthometric/rad_e)**2)

    #-- Initializing preliminary spherical harmonic matrices
    yclm = np.zeros((LMAX+1,MMAX+1))
    yslm = np.zeros((LMAX+1,MMAX+1))
    #-- Initializing output spherical harmonic matrices
    Ylms = gravity_toolkit.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX+1,MMAX+1))
    Ylms.slm = np.zeros((LMAX+1,MMAX+1))
    for l in range(0,LMAX+1):#-- equivalent to 0:LMAX
        mm = np.min([MMAX,l])#-- truncate to MMAX if specified (if l > MMAX)
        m = np.arange(0,mm+1)#-- mm+1 elements between 0 and mm
        #-- total pressure factor
        pfactor = np.zeros((nlat,nlon))
        #-- iterate over pressure levels
        for p in range(nlevels):
            #-- if using Swenson and Wahr (2002) or Boy and Chao (2005)
            if (METHOD == 'SW02'):
                #-- calculate pressure change/gravity ratio
                PG = pressure[p,:,:]/g0
                #-- add to pressure factor (pfactor) to integrate over levels
                pfactor += PG*(rad_e/(rad_e-GPH[p,:,:])+(GEOID/rad_e))**(l+4)
            elif (METHOD == 'BC05'):
                #-- calculate pressure change/gravity ratio
                PG = pressure[p,:,:]/gamma_h[p,:,:]
                #-- add to pressure factor (pfactor) to integrate over levels
                pfactor += PG*(R[p,:,:]/rad_e)**(l+2)
        #-- Multiplying gridded data with sin/cos of m#phis
        #-- This will sum through all phis in the dot product
        #-- need to reform pfactor to lonXlat as is originally latXlon
        #-- output [m,theta]
        dcos = np.dot(ccos,-np.transpose(pfactor))
        dsin = np.dot(ssin,-np.transpose(pfactor))
        #-- Summing product of plms and data over all latitudes
        #-- axis=1 signifies the direction of the summation (colatitude (th))
        #-- ycos and ysin are the SH coefficients before normalizing
        yclm[l,m] = np.sum(plm[l,m,:]*dcos[m,:], axis=1)
        yslm[l,m] = np.sum(plm[l,m,:]*dsin[m,:], axis=1)
        #-- Multiplying by coefficients to normalize
        Ylms.clm[l,m] = dfactor[l]*yclm[l,m]
        Ylms.slm[l,m] = dfactor[l]*yslm[l,m]

    #-- return the harmonics object
    return Ylms

#-- Main program that calls reanalysis_mean_harmonics()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Reads atmospheric geopotential heights
            fields from reanalysis and calculates a multi-annual
            mean set of spherical harmonics using a 3D geometry
            """
    )
    #-- command line parameters
    choices = ['ERA-Interim','ERA5','MERRA-2']
    parser.add_argument('model',
        metavar='MODEL', type=str, nargs='+',
        default=['ERA5','MERRA-2'], choices=choices,
        help='Reanalysis Model')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- mean pressure field to calculate
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
        reanalysis_mean_harmonics(args.directory, MODEL, RANGE=args.mean,
            REDISTRIBUTE=args.redistribute, LMAX=args.lmax, MMAX=args.mmax,
            LOVE_NUMBERS=args.love, REFERENCE=args.reference,
            DATAFORM=args.format, VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
