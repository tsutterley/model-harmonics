#!/usr/bin/env python
u"""
reanalysis_mean_harmonics.py
Written by Tyler Sutterley (05/2023)
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
    constants.py: calculate reference parameters for common ellipsoids
    gen_atmosphere_stokes.py: converts atmospheric fields to spherical harmonics
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
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: add root attributes to output netCDF4 and HDF5 files
    Updated 02/2023: use love numbers class with additional attributes
    Updated 12/2022: single implicit import of spherical harmonic tools
        use constants class in place of geoid-toolkit ref_ellipsoid
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: use wrapper function for reading load Love numbers
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 09/2021: use GRACE/GRACE-FO month to calendar month converters
    Updated 07/2021: can use input files to define command line arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 01/2021: read from netCDF4 file in slices to reduce memory load
        separated gen_atmosphere_stokes to a separate function
    Updated 12/2020: using argparse to set command line options
        using time module for operations and for extracting time units
    Updated 05/2020: use harmonics class for spherical harmonic operations
    Updated 04/2020: set path to load love numbers file
    Updated 01/2020: use EGM2008 geoid with ERA-Interim for consistency
        iterate over Julian days to calculate harmonics for consistency
    Updated 10/2019: changing Y/N flags to True/False
    Updated 08/2019: adjust time scale variable for MERRA-2
    Updated 07/2018: added parameters for ERA5
    Updated 05/2018: added uniform redistribution of oceanic values
    Updated 03/2018: simplified love number extrapolation if LMAX > 696
    Written 03/2018
"""
from __future__ import print_function

import sys
import re
import logging
import netCDF4
import pathlib
import argparse
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc

# PURPOSE: read atmospheric surface pressure fields and convert to harmonics
def reanalysis_mean_harmonics(base_dir, MODEL, RANGE=None, REDISTRIBUTE=False,
    LMAX=0, MMAX=None,  LOVE_NUMBERS=0, REFERENCE=None, DATAFORM=None,
    VERBOSE=False, MODE=0o775):

    # create logger for verbosity level
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[VERBOSE])

    # directory setup
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    ddir = base_dir.joinpath(MODEL)

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
    output_dir = ddir.joinpath('{0}_ATMOSPHERE_CLM_L{1:d}{2}{3}'.format(*args))
    output_dir.mkdir(mode=MODE, parents=True, exist_ok=True)
    # attributes for output files
    attributes = {}
    attributes['project'] = MODEL
    attributes['product_type'] = 'gravity_field'
    # attributes for input files
    attributes['lineage'] = []
    attributes['lineage'].append(input_invariant_file)
    attributes['lineage'].append(input_geoid_file)
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'

    # read model latitude and longitude from invariant parameters file
    with netCDF4.Dataset(ddir.joinpath(input_invariant_file),'r') as fileID:
        lon = fileID.variables[LONNAME][:].copy()
        lat = fileID.variables[LATNAME][:].copy()
    # calculate colatitude
    theta = (90.0 - lat)*np.pi/180.0
    # calculate meshgrid from latitude and longitude
    gridlon,gridlat = np.meshgrid(lon,lat)
    gridphi = gridlon*np.pi/180.0
    gridtheta = (90.0 - gridlat)*np.pi/180.0

    # read load love numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE, FORMAT='class')
    # add attributes for earth parameters
    attributes['earth_model'] = LOVE.model
    attributes['earth_love_numbers'] = LOVE.citation
    attributes['reference_frame'] = LOVE.reference
    # add attributes for maximum degree and order
    attributes['max_degree'] = LMAX
    attributes['max_order'] = MMAX

    # calculate Legendre polynomials
    PLM, dPLM = gravtk.plm_holmes(LMAX, np.cos(theta))
    # read geoid heights and grid step size
    geoid,gridstep = ncdf_geoid(ddir.joinpath(input_geoid_file))

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
        ii,jj = ncdf_landmask(ddir.joinpath(input_mask_file),
            MASKNAME, OCEAN)
        # calculate total area of oceanic points
        TOTAL_AREA = np.sum(AREA[ii,jj])

    # read each reanalysis pressure field and convert to spherical harmonics
    # then calculate mean of the spherical harmonic fields
    regex_years = r'|'.join([rf'{Y:4d}' for Y in range(RANGE[0],RANGE[1]+1)])
    rx = re.compile(regex_pattern.format(regex_years),re.VERBOSE)
    input_files = sorted([f for f in ddir.iterdir() if rx.match(f.name)])

    # list of monthly spherical harmonics
    harmonics_list = []
    # for each reanalysis file
    for i,input_file in enumerate(input_files):
        # read model level geopotential height data
        logging.debug(str(input_file))
        fileID = netCDF4.Dataset(input_file, mode='r')
        # append file to lineage
        attributes['lineage'].append(input_file.name)
        # extract shape from geopotential variable
        ntime, nlevels, nlat, nlon = fileID.variables[ZNAME].shape
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
                    PD[p,ii,jj] = np.sum(PD[p,ii,jj]*AREA[ii,jj])/TOTAL_AREA
            # calculate spherical harmonics for month
            Ylms = mdlhmc.gen_atmosphere_stokes(GPH, PD, lon, lat,
                LMAX=LMAX, MMAX=MMAX, ELLIPSOID=ELLIPSOID, GEOID=geoid,
                PLM=PLM, LOVE=LOVE)
            # convert julian dates to calendar then to year-decimal
            YY,MM,DD,hh,mm,ss = gravtk.time.convert_julian(JD,
                FORMAT='tuple')
            Ylms.time, = gravtk.time.convert_calendar_decimal(YY,
                MM, day=DD, hour=hh, minute=mm, second=ss)
            # calculate GRACE month from calendar dates
            Ylms.month = gravtk.time.calendar_to_grace(YY, MM)
            # append to list of harmonics
            harmonics_list.append(Ylms)
        # close the input netCDF4 file
        fileID.close()

    # divide by the number of dates to calculate mean from totals
    mean_Ylms = gravtk.harmonics().from_list(harmonics_list).mean()
    # add attributes to output harmonics
    mean_Ylms.attributes['ROOT'] = attributes

    # output mean spherical harmonics file
    args = (MODEL.upper(), LMAX, order_str, RANGE[0], RANGE[1], suffix[DATAFORM])
    FILE = '{0}_MEAN_CLM_L{1:d}{2}_{3:4d}-{4:4d}.{5}'.format(*args)
    output_mean_file = output_dir.joinpath(FILE)
    mean_Ylms.to_file(output_mean_file, format=DATAFORM, verbose=VERBOSE)
    # set the permissions level of the output file to MODE
    output_mean_file.chmod(mode=MODE)

# PURPOSE: read geoid height netCDF4 files from read_gfz_geoid_grids.py
def ncdf_geoid(FILENAME):
    logging.debug(str(FILENAME))
    with netCDF4.Dataset(FILENAME, mode='r') as fileID:
        geoid_undulation = fileID.variables['geoid'][:].copy()
        gridstep = np.array(fileID.gridstep.split(','),dtype=np.float64)
    return (geoid_undulation,np.squeeze(gridstep))

# PURPOSE: read land sea mask to get indices of oceanic values
def ncdf_landmask(FILENAME, MASKNAME, OCEAN):
    logging.debug(str(FILENAME))
    with netCDF4.Dataset(FILENAME, mode='r') as fileID:
        landsea = np.squeeze(fileID.variables[MASKNAME][:].copy()).astype('f2')
    return np.nonzero(landsea == OCEAN)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads atmospheric geopotential heights
            fields from reanalysis and calculates a multi-annual
            mean set of spherical harmonics using a 3D geometry
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
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # mean pressure field to calculate
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

    # for each reanalysis model
    for MODEL in args.model:
        # run program
        reanalysis_mean_harmonics(args.directory, MODEL, RANGE=args.mean,
            REDISTRIBUTE=args.redistribute, LMAX=args.lmax, MMAX=args.mmax,
            LOVE_NUMBERS=args.love, REFERENCE=args.reference,
            DATAFORM=args.format, VERBOSE=args.verbose, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
