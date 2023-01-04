#!/usr/bin/env python
u"""
ecco_monthly_harmonics.py
Written by Tyler Sutterley (12/2022)
Reads monthly ECCO ocean bottom pressure anomalies and converts to
    spherical harmonic coefficients

INPUTS:
    ECCO Near Real-Time models
        kf080i: Kalman filter analysis
            https://ecco.jpl.nasa.gov/drive/files/NearRealTime/KalmanFilter/
        dr080i: RTS smoother analysis
            https://ecco.jpl.nasa.gov/drive/files/NearRealTime/Smoother/
    ECCO2 Cube92 models
        Cube92
    ECCO version 4 models
        V4r3: Version 4, Revision 3
        V4r4: Version 4, Revision 4

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -Y X, --year X: Years to run
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -n X, --love X: Load Love numbers dataset
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -F X, --format X: Input and output data format
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
    associated_legendre.py: computes fully-normalized associated Legendre polynomials
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    constants.py: calculate reference parameters for common ellipsoids
    ref_ellipsoid.py: calculate reference parameters for common ellipsoids
    norm_gravity.py: calculates the normal gravity for locations on an ellipsoid
    gen_pressure_stokes.py: converts a pressure field into spherical harmonics
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: spatial data class for reading, writing and processing data
    time.py: utilities for calculating time operations
    utilities.py: download and management utilities for files

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
    Updated 02/2021: separate inputs to gen_pressure_stokes
    Updated 01/2021: added Cube92 choice to input model types
        outputs from gen_pressure_stokes are now harmonics objects
    Updated 12/2020: use argparse to set command line parameters
        using spatial and harmonics modules for read/write operations
        added more love number options. using utilities from time module
    Updated 10/2019: changing Y/N flags to True/False
    Updated 06/2019: recommending kf080i for the Kalman filtered solution
    Updated 10/2018: separated gen_pressure_stokes into separate function
    Updated 07/2018: output index and date files in separate loop for all files
    Updated 03/2018: use realistic geometry from bathymetry and local gravity
        simplified love number extrapolation if LMAX is greater than 696
    Updated 01/2018: using getopt to set parameters
    Updated 08/2017: convert from geodetic coordinates to geocentric
    Updated 08/2016: fixed find_new_files function with previous updates
    Updated 06/2016: can use dr080g model, using __future__ print option
    Updated 05/2016: complete rewrite of program
    Written 05/2013
"""
from __future__ import print_function

import sys
import os
import re
import logging
import netCDF4
import argparse
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc
import geoid_toolkit as geoidtk

# PURPOSE: convert monthly ECCO OBP data to spherical harmonics
def ecco_monthly_harmonics(ddir, MODEL, YEARS, LMAX=0, MMAX=None,
    LOVE_NUMBERS=0, REFERENCE=None, DATAFORM=None, MODE=0o775):

    # input and output subdirectory
    input_sub = f'ECCO_{MODEL}_AveRmvd_OBP'
    output_sub = f'ECCO_{MODEL}_AveRmvd_OBP_CLM_L{LMAX:d}'

    # upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    # output string for both LMAX == MMAX and LMAX != MMAX cases
    order_str = 'M{MMAX:d}' if (MMAX != LMAX) else ''
    # output file format
    output_file_format = 'ECCO_{0}_AveRmvd_OBP_CLM_L{1:d}{2}_{3:03d}.{4}'
    # Creating subdirectory if it doesn't exist
    if (not os.access(os.path.join(ddir,output_sub), os.F_OK)):
        os.makedirs(os.path.join(ddir,output_sub),MODE)
    # input/output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # parameters for each model
    if MODEL in ('kf080i','dr080i'):
        # grid step size
        dlon,dlat = (1.0,1.0)
        # grid extent
        LAT_MAX = 78.5
        extent = [0.5,359.5,-LAT_MAX,LAT_MAX]
        input_depth_file = os.path.join(ddir,'depth.nc')
        input_geoid_file = os.path.join(ddir,'egm_2008.nc')
        # indices to read
        indices = np.arange(1,2*LAT_MAX+2).astype(np.int64)
    elif MODEL in ('Cube92',):
        # grid step size
        dlon,dlat = (0.25,0.25)
        # grid extent
        extent = [0.125,359.875,-89.875,89.875]
        input_depth_file = os.path.join(ddir,'DEPTH.2020.1440x720.nc')
        input_geoid_file = os.path.join(ddir,'EGM_2008.1440x720.nc')
        # indices to read (all)
        indices = Ellipsis
    elif MODEL in ('V4r3','V4r4'):
        # grid step size
        dlon,dlat = (0.5,0.5)
        # grid extent
        extent = [-179.75,179.75,-89.75,89.75]
        input_depth_file = os.path.join(ddir,'DEPTH.2020.720x360.nc')
        input_geoid_file = os.path.join(ddir,'EGM_2008.720x360.nc')
        # indices to read (all)
        indices = Ellipsis
    # attributes for output files
    attributes = {}
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'

    # input grid dimensions
    glon = np.arange(extent[0],extent[1]+dlon,dlon)
    glat = np.arange(extent[2],extent[3]+dlat,dlat)
    # create mesh grids of datasets
    gridlon,gridlat = np.meshgrid(glon,glat)

    # read geoid and depth to calculate bathymetry
    depth = ncdf_depth(input_depth_file, indices=indices)
    geoid_undulation,gridstep = ncdf_geoid(input_geoid_file, indices=indices)
    bathymetry = geoid_undulation - depth

    # Earth Parameters
    ellipsoid_params = mdlhmc.constants(ellipsoid='WGS84')
    # semimajor axis of ellipsoid [m]
    a_axis = ellipsoid_params.a_axis
    # first numerical eccentricity
    ecc1 = ellipsoid_params.ecc1
    # convert from geodetic latitude to geocentric latitude
    # geodetic latitude in radians
    latitude_geodetic_rad = np.pi*gridlat/180.0
    # prime vertical radius of curvature
    N = a_axis/np.sqrt(1.0 - ecc1**2.*np.sin(latitude_geodetic_rad)**2.)
    # calculate X, Y and Z from geodetic latitude and longitude
    X = (N+bathymetry)*np.cos(latitude_geodetic_rad)*np.cos(np.pi*gridlon/180.0)
    Y = (N+bathymetry)*np.cos(latitude_geodetic_rad)*np.sin(np.pi*gridlon/180.0)
    Z = (N * (1.0 - ecc1**2.0) + bathymetry) * np.sin(latitude_geodetic_rad)
    R = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
    # calculate geocentric latitude and convert to degrees
    latitude_geocentric = 180.0*np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))/np.pi
    # colatitude in radians
    theta = (90.0 - latitude_geocentric)*np.pi/180.0

    # calculate normal gravity at latitudes and bathymetry
    gamma_h,dgamma_dh = geoidtk.norm_gravity(latitude_geocentric,
        bathymetry, 'WGS84')

    # read load love numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)
    # calculate Legendre polynomials
    PLM, dPLM = gravtk.plm_holmes(LMAX,np.cos(theta[:,0]))

    # regular expression pattern to find files and extract dates
    regex_years = r'\d+' if (YEARS is None) else '|'.join(map(str,YEARS))
    args = (MODEL, regex_years, suffix[DATAFORM])
    rx = re.compile(r'ECCO_{0}_AveRmvd_OBP_({1})_(\d+).{2}$'.format(*args))

    # find input ECCO OBP files
    FILES = [fi for fi in os.listdir(os.path.join(ddir,input_sub)) if rx.match(fi)]

    # for each input file
    for f in sorted(FILES):
        # extract dates from file
        year,month = np.array(rx.findall(f).pop(), dtype=np.int64)
        # read input data file
        if (DATAFORM == 'ascii'):
            obp_data = gravtk.spatial(spacing=[dlon,dlat],
                nlat=150,nlon=360,extent=extent).from_ascii(
                os.path.join(ddir,input_sub,f))
        elif (DATAFORM == 'netCDF4'):
            obp_data = gravtk.spatial().from_netCDF4(
                os.path.join(ddir,input_sub,f))
        elif (DATAFORM == 'HDF5'):
            obp_data = gravtk.spatial().from_HDF5(
                os.path.join(ddir,input_sub,f))
        # replace fill value points with 0
        obp_data.replace_invalid(0.0)
        # calculate spherical harmonics from pressure/gravity ratio
        obp_Ylms = mdlhmc.gen_pressure_stokes(obp_data.data, gamma_h, R,
            glon, latitude_geocentric[:,0], LMAX=LMAX, MMAX=MMAX,
            PLM=PLM, LOVE=LOVE)
        obp_Ylms.time = np.copy(obp_data.time)
        obp_Ylms.month = gravtk.time.calendar_to_grace(year,month)
        # output spherical harmonic data file
        args = (MODEL, LMAX, order_str, obp_Ylms.month, suffix[DATAFORM])
        FILE = output_file_format.format(*args)
        obp_Ylms.to_file(os.path.join(ddir,output_sub,FILE),
            format=DATAFORM, **attributes)
        # change the permissions mode of the output file to MODE
        os.chmod(os.path.join(ddir,output_sub,FILE),MODE)

    # Output date ascii file
    output_date_file = f'ECCO_{MODEL}_OBP_DATES.txt'
    fid1 = open(os.path.join(ddir,output_sub,output_date_file),
        mode='w', encoding='utf8')
    # date file header information
    print('{0:8} {1:^6} {2:^5}'.format('Mid-date','GRACE','Month'), file=fid1)
    # index file listing all output spherical harmonic files
    output_index_file = 'index.txt'
    fid2 = open(os.path.join(ddir,output_sub,output_index_file),
        mode='w', encoding='utf8')
    # find all available output files
    args = (MODEL, LMAX, suffix[DATAFORM])
    output_regex=r'ECCO_{0}_AveRmvd_OBP_CLM_L{1:d}_([-]?\d+).{2}'.format(*args)
    # find all output ECCO OBP harmonic files (not just ones created in run)
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

# PURPOSE: read ECCO2 depth file
# ftp://mit.ecco-group.org/ecco_for_las/grid_fields/
def ncdf_depth(FILENAME, indices=Ellipsis):
    with netCDF4.Dataset(FILENAME,'r') as fileID:
        depth = np.array(fileID.variables['depth'][indices,:])
        fill_value = fileID.variables['depth']._FillValue
        depth[depth == fill_value] = 0.0
    return depth

# PURPOSE: read geoid height netCDF4 files from read_gfz_geoid_grids.py
def ncdf_geoid(FILENAME, indices=Ellipsis):
    with netCDF4.Dataset(FILENAME,'r') as fileID:
        geoid_undulation = np.array(fileID.variables['geoid'][indices,:])
        gridstep = [float(s) for s in fileID.gridstep.split(',')]
    return (geoid_undulation,np.squeeze(gridstep))

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads monthly ECCO ocean bottom pressure
            anomalies and converts to spherical harmonic coefficients
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('model',
        type=str, nargs='+',
        default=['kf080i','dr080i'],
        choices=['kf080i','dr080i','Cube92','V4r3','V4r4'],
        help='ECCO Model')
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # years to run
    now = gravtk.time.datetime.datetime.now()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.year+1),
        help='Years of model outputs to run')
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

    # for each ECCO model
    for MODEL in args.model:
        # run program
        ecco_monthly_harmonics(args.directory, MODEL, args.year,
            LMAX=args.lmax, MMAX=args.mmax, LOVE_NUMBERS=args.love,
            REFERENCE=args.reference, DATAFORM=args.format,
            MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
