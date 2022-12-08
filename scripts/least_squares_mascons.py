#!/usr/bin/env python
u"""
least_squares_mascons.py
Written by Tyler Sutterley (12/2022)

Calculates regional mass anomalies through a least-squares mascon procedure
    from an index of spherical harmonic coefficient files

INPUTS:
    Input index file for model harmonics

COMMAND LINE OPTIONS:
    --help: list the command line options
    -O X, --output-directory X: output directory for mascon files
    -P X, --file-prefix X: prefix string for mascon files
    -D, --date: Model harmonics are a time series
    --lmin X: minimum spherical harmonic degree
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -R X, --radius X: Gaussian smoothing radius (km)
    -d, --destripe: use decorrelation filter (destriping filter)
    -n X, --love X: Load Love numbers dataset
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -F X, --format X: input data format
        ascii
        netCDF4
        HDF5
    --mask X: Land-sea mask for redistributing mascon mass and land water flux
    --mascon-file X: index file of mascons spherical harmonics
    --mascon-format X: input format for mascon files
    --redistribute-mascons: redistribute mascon mass over the ocean
    --redistribute-mass: redistribute input mass fields over the ocean
    --harmonic-errors: input spherical harmonic fields are data errors
    --fit-method X: method for fitting sensitivity kernel to harmonics
        1: mass coefficients
        2: geoid coefficients
    --log: Output log of files created for each job
    -V, --verbose: Verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        http://www.h5py.org/
    future: Compatibility layer between Python 2 and Python 3
        http://python-future.org/

PROGRAM DEPENDENCIES:
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    gen_stokes.py: converts a spatial field into spherical harmonic coefficients
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    ocean_stokes.py: reads a land-sea mask and converts to spherical harmonics
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    utilities.py: download and management utilities for files

REFERENCES:
    I Velicogna, T C Sutterley and M R van den Broeke. "Regional acceleration
        in ice mass loss from Greenland and Antarctica using GRACE
        time-variable gravity data". Geophysical Research Letters,
        41(22):8130-8137, 2014. https://doi.org/10.1002/2014GL061052

    T Jacob, J Wahr, W Pfeffer, and S C Swenson "Recent contributions of
        glaciers and ice caps to sea level rise". Nature, 482, 514-518 (2012).
        https://doi.org/10.1038/nature10847

    V M Tiwari, J Wahr, S and Swenson, "Dwindling groundwater resources in
        northern India, from satellite gravity observations",
        Geophysical Research Letters, 36(18), L18401, (2009).
        https://doi.org/10.1029/2009GL039401

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: use wrapper function for reading load Love numbers
        include utf-8 encoding in reads to be windows compliant
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 08/2021: add option for setting input format of the mascon files
    Updated 07/2021: switch from parameter files to argparse arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: add parser object for removing commented or empty lines
    Updated 01/2021: harmonics object output from gen_stokes.py/ocean_stokes.py
    Updated 12/2020: added more love number options
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 04/2020: using the harmonics class for spherical harmonic operations
        updated load love numbers read function
    Updated 03/2020: switched to destripe_harmonics for filtering harmonics
    Updated 10/2019: changing Y/N flags to True/False
    Updated 10/2018: verify integers for python3 compatibility
    Updated 06/2018: using python3 compatible octal and input
    Updated 03/2018: include order_str denoting if MMAX != LMAX
        added extrapolation of load love numbers if LMAX > 696
    Updated 01/2018: recursively make output directories
    Updated 09/2017: use a different land-sea mask for calculating ocean_Ylms
        use rcond=-1 in numpy least-squares algorithm
    Updated 02/2017: input spherical harmonic mapping indices as integers
        clean up ocean redistribution read function
    Updated 01/2017: ocean_stokes for reading/converting ocean function
    Updated 05-07/2016: output full mass of input harmonics for test cases
        using __future__ print function
    Updated 02/2016: direct calculation of number of harmonics n_harm
        use getopt parameters to set number of PROCESSES to run in parallel,
            whether or not to output a log file, added new help module
    Updated 11/2015: create unique log filenames
    Updated 06/2015: added output_files for log files
    Updated 05/2015: added parameter MMAX for LMAX != MMAX
    Updated 03/2015: updated code with updates from calc_mascon program
         error handling with traceback, ocean redistribution, multiple data types
    Updated 11/2014: general program for non-grace fields
    Updated 10/2014: forked from R and converted to python
        distributed computing of tasks with the multiprocessing module
    Updated 05/2014: added import functions
    Updated 02/2014: updated comments and added os.path.joins for connecting
        directories and files (generalizing code)
        some general updates to the program code
    Updated 09/2013: saving GRACE DELTA file (won't calculate each time)
        added option to remove RACMO data
    Updated 08/2013: general updates to inputting data
        wrote grace_find_months, grace_input_months, gia_input
        to input spherical harmonics similar to python programs
     Updated 03/2012: edited to use new gen_stokes time-series option
     Updated 02/2012: Added sensitivity kernels
     Written 02/2012
"""
from __future__ import print_function, division

import sys
import os
import re
import time
import logging
import argparse
import numpy as np
import traceback
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(os.path.basename(sys.argv[0]))
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: calculate a regional time-series through a least
# squares mascon process
def least_squares_mascons(input_file, LMAX, RAD,
    LMIN=None,
    MMAX=None,
    DESTRIPE=False,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    DATAFORM=None,
    MASCON_FILE=None,
    MASCON_FORMAT=None,
    REDISTRIBUTE_MASCONS=False,
    FIT_METHOD=0,
    REDISTRIBUTE=False,
    DATA_ERROR=False,
    LANDMASK=None,
    OUTPUT_DIRECTORY=None,
    FILE_PREFIX=None,
    DATE=False,
    MODE=0o775):

    # Recursively create output directory if not currently existing
    if (not os.access(OUTPUT_DIRECTORY,os.F_OK)):
        os.makedirs(OUTPUT_DIRECTORY, mode=MODE, exist_ok=True)

    # list object of output files for file logs (full path)
    output_files = []
    # file parser for reading index files
    # removes commented lines (can comment out files in the index)
    # removes empty lines (if there are extra empty lines)
    parser = re.compile(r'^(?!\#|\%|$)', re.VERBOSE)

    # read arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)

    # Earth Parameters
    factors = gravtk.units(lmax=LMAX).harmonic(hl,kl,ll)
    # Average Density of the Earth [g/cm^3]
    rho_e = factors.rho_e
    # Average Radius of the Earth [cm]
    rad_e = factors.rad_e

    # Calculating the Gaussian smoothing for radius RAD
    if (RAD != 0):
        wt = 2.0*np.pi*gravtk.gauss_weights(RAD,LMAX)
        gw_str = f'_r{RAD:0.0f}km'
    else:
        # else = 1
        wt = np.ones((LMAX+1))
        gw_str = ''

    # output string for both LMAX==MMAX and LMAX != MMAX cases
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    order_str = 'M{MMAX:d}' if (MMAX != LMAX) else ''
    # output string for destriped harmonics
    ds_str = '_FL' if DESTRIPE else ''

    # Read Ocean function and convert to Ylms for redistribution
    if (REDISTRIBUTE_MASCONS | REDISTRIBUTE):
        # read Land-Sea Mask and convert to spherical harmonics
        ocean_Ylms = gravtk.ocean_stokes(LANDMASK, LMAX, MMAX=MMAX,
            LOVE=(hl,kl,ll))
        ocean_str = '_OCN'
    else:
        # not distributing uniformly over ocean
        ocean_str = ''

    # read spherical harmonics file in data format
    if DATAFORM in ('ascii','netCDF4','HDF5'):
        # ascii (.txt)
        # netCDF4 (.nc)
        # HDF5 (.H5)
        data_Ylms = gravtk.harmonics().from_file(input_file,
            format=DATAFORM, date=DATE)
    elif DATAFORM in ('index-ascii','index-netCDF4','index-HDF5'):
        # read from index file
        _,dataform = DATAFORM.split('-')
        data_Ylms = gravtk.harmonics().from_index(input_file,
            format=dataform, date=DATE)
    # number of files within the index
    n_files = data_Ylms.shape[-1]
    # truncate to degree and order
    data_Ylms = data_Ylms.truncate(lmax=LMAX,mmax=MMAX)
    # Total mass in the simulated input data in gigatonnes
    rmass = 4.0*np.pi*(rad_e**3.0)*rho_e*data_Ylms.clm[0,0,:]/3.0/1e15
    # distribute Ylms uniformly over the ocean
    if REDISTRIBUTE:
        # calculate ratio between total removed mass and
        # a uniformly distributed cm of water over the ocean
        ratio = data_Ylms.clm[0,0,:]/ocean_Ylms.clm[0,0]
        # for each spherical harmonic
        for m in range(0,MMAX+1):# MMAX+1 to include MMAX
            for l in range(m,LMAX+1):# LMAX+1 to include LMAX
                # remove the ratio*ocean Ylms from Ylms
                # note: x -= y is equivalent to x = x - y
                data_Ylms.clm[l,m,:] -= ratio*ocean_Ylms.clm[l,m]
                data_Ylms.slm[l,m,:] -= ratio*ocean_Ylms.slm[l,m]
    # filter data coefficients
    if DESTRIPE:
        data_Ylms = data_Ylms.destripe()

    # input mascon spherical harmonic datafiles
    with open(MASCON_FILE, mode='r', encoding='utf8') as f:
        mascon_files = [l for l in f.read().splitlines() if parser.match(l)]
    # number of mascons
    n_mas = len(mascon_files)
    # spatial area of the mascon
    area_tot = np.zeros((n_mas))
    # name of each mascon
    mascon_name = []
    # for each valid file in the index (iterate over mascons)
    mascon_list = []
    for k in range(n_mas):
        # read mascon spherical harmonics
        Ylms = gravtk.harmonics().from_file(mascon_files[k],
            format=MASCON_FORMAT, date=False)
        # Calculating the total mass of each mascon (1 cmH2O uniform)
        area_tot[k] = 4.0*np.pi*(rad_e**3)*rho_e*Ylms.clm[0,0]/3.0
        # distribute MASCON mass uniformly over the ocean
        if REDISTRIBUTE_MASCONS:
            # calculate ratio between total mascon mass and
            # a uniformly distributed cm of water over the ocean
            ratio = Ylms.clm[0,0]/ocean_Ylms.clm[0,0]
            # for each spherical harmonic
            for m in range(0,MMAX+1):# MMAX+1 to include MMAX
                for l in range(m,LMAX+1):# LMAX+1 to include LMAX
                    # remove ratio*ocean Ylms from mascon Ylms
                    # note: x -= y is equivalent to x = x - y
                    Ylms.clm[l,m] -= ratio*ocean_Ylms.clm[l,m]
                    Ylms.slm[l,m] -= ratio*ocean_Ylms.slm[l,m]
        # truncate mascon spherical harmonics to d/o LMAX/MMAX and add to list
        mascon_list.append(Ylms.truncate(lmax=LMAX, mmax=MMAX))
        # mascon base is the file without directory or suffix
        mascon_base = os.path.basename(mascon_files[k])
        mascon_base = os.path.splitext(mascon_base)[0]
        # if lower case, will capitalize
        mascon_base = mascon_base.upper()
        # if mascon name contains degree and order info, remove
        mascon_name.append(mascon_base.replace(f'_L{LMAX:d}', ''))
    # create single harmonics object from list
    mascon_Ylms = gravtk.harmonics().from_list(mascon_list, date=False)
    # clear mascon list variable
    del mascon_list

    # Calculating the number of cos and sin harmonics between LMIN and LMAX
    # taking into account MMAX (if MMAX == LMAX then LMAX-MMAX=0)
    n_harm=np.int64(LMAX**2 - LMIN**2 + 2*LMAX + 1 - (LMAX-MMAX)**2 - (LMAX-MMAX))

    # Initialing harmonics for least squares fitting
    # mascon kernel
    M_lm = np.zeros((n_harm,n_mas))
    # mascon kernel converted to output unit
    MA_lm = np.zeros((n_harm,n_mas))
    # corrected clm and slm
    Y_lm = np.zeros((n_harm,n_files))
    # sensitivity kernel
    A_lm = np.zeros((n_harm,n_mas))
    # Initializing output Mascon time-series
    mascon = np.zeros((n_mas,n_files))
    # Initializing conversion factors
    # factor for converting to coefficients of mass
    fact = np.zeros((n_harm))
    # smoothing factor
    wt_lm = np.zeros((n_harm))

    # ii is a counter variable for building the mascon column array
    ii = 0
    # Creating column array of clm/slm coefficients
    # Order is [C00...C6060,S11...S6060]
    coeff = rho_e*rad_e/3.0
    # Switching between Cosine and Sine Stokes
    for cs,csharm in enumerate(['clm','slm']):
        # copy cosine and sin harmonics
        mascon_harm = getattr(mascon_Ylms, csharm)
        data_harm = getattr(data_Ylms, csharm)
        # for each spherical harmonic degree
        # +1 to include LMAX
        for l in range(LMIN,LMAX+1):
            # for each spherical harmonic order
            # Sine Stokes for (m=0) = 0
            mm = np.min([MMAX,l])
            # +1 to include l or MMAX (whichever is smaller)
            for m in range(cs,mm+1):
                # Mascon Spherical Harmonics
                M_lm[ii,:] = np.copy(mascon_harm[l,m,:])
                # Data Spherical Harmonics
                Y_lm[ii,:] = np.copy(data_harm[l,m,:])
                # degree dependent factor to convert to mass
                fact[ii] = (2.0*l+1.0)/(1.0 + kl[l])
                # degree dependent smoothing
                wt_lm[ii] = np.copy(wt[l])
                # add 1 to counter
                ii += 1
    # free up memory from data harmonics
    data_Ylms.clm = None
    data_Ylms.slm = None

    # Converting mascon coefficients to fit method
    if (FIT_METHOD == 1):
        # Fitting Sensitivity Kernel as mass coefficients
        # converting M_lm to mass coefficients of the kernel
        for i in range(n_harm):
            MA_lm[i,:] = wt_lm[i]*M_lm[i,:]*fact[i]
        fit_factor = wt_lm*fact
    else:
        # Fitting Sensitivity Kernel as geoid coefficients
        MA_lm[:,:] = np.copy(M_lm)
        fit_factor = wt_lm*np.ones((n_harm))

    # Fitting the sensitivity kernel from the input kernel
    for i in range(n_harm):
        # setting kern_i equal to 1 for d/o
        kern_i = np.zeros((n_harm))
        # converting to mass coefficients if specified
        kern_i[i] = 1.0*fit_factor[i]
        # spherical harmonics solution for the
        # mascon sensitivity kernels
        # Least Squares Solutions: Inv(X'.X).(X'.Y)
        kern_lm = np.linalg.lstsq(MA_lm,kern_i,rcond=-1)[0]
        for k in range(n_mas):
            A_lm[i,k] = kern_lm[k]*area_tot[k]

    # output formatting string if containing date variables
    if DATE:
        formatting_string = '{0:03d} {1:12.4f} {2:16.10f} {3:16.5f} {4:16.10f}'
    else:
        formatting_string = '{0:16.10f} {1:16.5f} {2:16.10f}'

    # for each mascon
    for k in range(n_mas):
        # output filename format:
        # mascon name, LMAX, Gaussian smoothing radii, filter flag
        fargs = (FILE_PREFIX, mascon_name[k], LMAX, order_str,
            gw_str, ds_str, ocean_str)
        file_format = '{0}{1}_L{2:d}{3}{4}{5}{6}.txt'
        output_file = os.path.join(OUTPUT_DIRECTORY,file_format.format(*fargs))

        # Output mascon datafiles
        # Will output each mascon mass series
        # if with dates:
        # month, date, mascon mass, mascon area
        # else:
        # mascon mass, mascon area
        # open output mascon time-series file
        fid = open(output_file, mode='w', encoding='utf8')
        # for each date
        for f in range(n_files):
            # Summing over all spherical harmonics for mascon k, and time t
            # multiplies by the degree dependent factor to convert
            # the harmonics into mass coefficients
            # Converting mascon mass time-series from g to gigatonnes
            if DATA_ERROR:
                # if input data are errors (i.e. estimated dealiasing errors)
                mascon[k,f] = np.sqrt(np.sum((A_lm[:,k]*Y_lm[:,f])**2))/1e15
            else:
                mascon[k,f] = np.sum(A_lm[:,k]*Y_lm[:,f])/1e15
            # output to file
            if DATE:
                # if files contain date information
                args = (data_Ylms.month[f], data_Ylms.time[f],
                    mascon[k,f], area_tot[k]/1e10, rmass[f])
                print(formatting_string.format(*args), file=fid)
            else:
                # just print the time-series and mascon areas
                args = (mascon[k,f], area_tot[k]/1e10, rmass[f])
                print(formatting_string.format(*args), file=fid)
        # close the output file
        fid.close()
        # change the permissions mode of the output file
        os.chmod(output_file, MODE)
        # add output files to list object
        output_files.append(output_file)

    # return the list of output files
    return output_files

# PURPOSE: print a file log for the mascon analysis
def output_log_file(input_arguments, output_files):
    # format: calc_mascon_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'calc_mascon_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = os.path.expanduser(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(os.path.join(DIRECTORY,LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    # print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(input_arguments).items()):
        logging.info(f'{arg}: {value}')
    # print output files
    logging.info('\n\nOUTPUT FILES:')
    for f in output_files:
        logging.info(f)
    # close the log file
    fid.close()

# PURPOSE: print a error file log for the mascon analysis
def output_error_log_file(input_arguments):
    # format: calc_mascon_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'calc_mascon_failed_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = os.path.expanduser(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(os.path.join(DIRECTORY,LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    # print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(input_arguments).items()):
        logging.info(f'{arg}: {value}')
    # print traceback error
    logging.info('\n\nTRACEBACK ERROR:')
    traceback.print_exc(file=fid)
    # close the log file
    fid.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates a time-series of regional mass anomalies
            through a least-squares mascon procedure procedure from an index
            of spherical harmonic coefficient files
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Input index file with spherical harmonic data files')
    # working data directory
    parser.add_argument('--output-directory','-O',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Output directory for mascon files')
    parser.add_argument('--file-prefix','-P',
        type=str,
        help='Prefix string for mascon files')
    parser.add_argument('--date','-D',
        default=False, action='store_true',
        help='Model harmonics are a time series')
    # minimum spherical harmonic degree
    parser.add_argument('--lmin',
        type=int, default=1,
        help='Minimum spherical harmonic degree')
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
    # Gaussian smoothing radius (km)
    parser.add_argument('--radius','-R',
        type=float, default=0,
        help='Gaussian smoothing radius (km)')
    # Use a decorrelation (destriping) filter
    parser.add_argument('--destripe','-d',
        default=False, action='store_true',
        help='Use decorrelation (destriping) filter')
    # input data format (ascii, netCDF4, HDF5)
    choices = []
    choices.extend(['ascii','netCDF4','HDF5'])
    choices.extend(['index-ascii','index-netCDF4','index-HDF5'])
    parser.add_argument('--format','-F',
        metavar='FORMAT', type=str,
        default='netCDF4', choices=choices,
        help='Input data format')
    # mascon index file and parameters
    parser.add_argument('--mascon-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Index file of mascons spherical harmonics')
    parser.add_argument('--mascon-format',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input data format for mascon files')
    parser.add_argument('--redistribute-mascons',
        default=False, action='store_true',
        help='Redistribute mascon mass over the ocean')
    # 1: mass coefficients
    # 2: geoid coefficients
    parser.add_argument('--fit-method',
        type=int, default=1, choices=(1,2),
        help='Method for fitting sensitivity kernel to harmonics')
    # redistribute total mass over the ocean
    parser.add_argument('--redistribute-mass',
        default=False, action='store_true',
        help='Redistribute total mass over the ocean')
    # calculating with data errors
    parser.add_argument('--harmonic-errors',
        default=False, action='store_true',
        help='Input spherical harmonic fields are data errors')
    # land-sea mask for redistributing mascon mass and land water flux
    lsmask = gravtk.utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), default=lsmask,
        help='Land-sea mask for redistributing mascon mass and land water flux')
    # Output log file for each job in forms
    # calc_mascon_run_2002-04-01_PID-00000.log
    # calc_mascon_failed_run_2002-04-01_PID-00000.log
    parser.add_argument('--log',
        default=False, action='store_true',
        help='Output log file for each job')
    # print information about processing run
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output files')
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

    # try to run the analysis with listed parameters
    try:
        info(args)
        # run least_squares_mascons algorithm with parameters
        output_files = least_squares_mascons(
            args.infile,
            args.lmax,
            args.radius,
            LMIN=args.lmin,
            MMAX=args.mmax,
            DESTRIPE=args.destripe,
            LOVE_NUMBERS=args.love,
            REFERENCE=args.reference,
            DATAFORM=args.format,
            MASCON_FILE=args.mascon_file,
            MASCON_FORMAT=args.mascon_format,
            REDISTRIBUTE_MASCONS=args.redistribute_mascons,
            FIT_METHOD=args.fit_method,
            REDISTRIBUTE=args.redistribute_mass,
            DATA_ERROR=args.harmonic_errors,
            LANDMASK=args.mask,
            OUTPUT_DIRECTORY=args.output_directory,
            FILE_PREFIX=args.file_prefix,
            DATE=args.date,
            MODE=args.mode)
    except Exception as e:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())
        if args.log:# write failed job completion log file
            output_error_log_file(args)
    else:
        if args.log:# write successful job completion log file
            output_log_file(args,output_files)

# run main program
if __name__ == '__main__':
    main()
