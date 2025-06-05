#!/usr/bin/env python
u"""
seismic_monthly_harmonics.py
Written by Tyler Sutterley (05/2023)

Reads monthly gravity changes due to earthquakes and
    converts to spherical harmonic coefficients

INPUTS:
    Seismic model:
        SHdegree40_M3_microGal
        SHdegree40_M5_microGal
        SHdegree60_M3_microGal
        SHdegree60_M5_microGal

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
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
    associated_legendre.py: computes fully-normalized associated Legendre polynomials
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    gen_stokes.py: converts a spatial field into a series of spherical harmonics
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: spatial data class for reading, writing and processing data
    time.py: utilities for calculating time operations

REFERENCES:
    T Jeon, K-W Seo, S-C Han, "Impact of the solid Earth mass adjustment
        by the 2011 Tohoku--Oki earthquake on the regional sea level
        and hydrological mass change recovery from GRACE",
        Geophysical Journal International, 235(2), 1373--1383, 2023.

UPDATE HISTORY:
    Written 10/2023
"""
from __future__ import print_function

import sys
import io
import re
import pathlib
import logging
import argparse
import datetime
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: convert seismic data to spherical harmonics
def seismic_monthly_harmonics(base_dir, MODEL, LMAX=60, MMAX=None,
    LOVE_NUMBERS=0, REFERENCE=None, DATAFORM=None, MODE=0o775):

    # directory setup
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    ddir = base_dir.joinpath(MODEL)
    mv = re.findall(r'SHdegree\d{2}_M(\d)_microGal', MODEL).pop()

    # regular expression pattern for finding files
    rx = re.compile(r'^T\d{2}_(\w{3}\d{4})\.txt$', re.IGNORECASE)
    # find files within directory
    input_files = [f for f in ddir.iterdir() if rx.match(f.name)]

    # output file attributes
    attributes = {}
    attributes['product_type'] = 'gravity_field'
    attributes['elastic_layer'] = '0\u201365 km in depth'
    attributes['asthenosphere'] = '65\u2013220 km in depth'
    attributes['Maxwell_viscosity'] = f'{mv}\u00d7 18 Pa s'
    attributes['Kelvin_viscosity'] = '0 (Maxwell body)'
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'

    # read load love numbers and calculate Legendre polynomials
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE, FORMAT='class')
    # add attributes for earth parameters
    attributes['earth_model'] = LOVE.model
    attributes['earth_love_numbers'] = LOVE.citation
    attributes['reference_frame'] = LOVE.reference
    # upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    # output string for both LMAX == MMAX and LMAX != MMAX cases
    order_str = 'M{MMAX:d}' if (MMAX != LMAX) else ''
    # add attributes for maximum degree and order
    attributes['max_degree'] = LMAX
    attributes['max_order'] = MMAX
    # add attributes for earth parameters
    factors = gravtk.units(lmax=LMAX).spatial(*LOVE)
    attributes['earth_radius'] = f'{factors.rad_e:0.3f} cm'
    attributes['earth_density'] = f'{factors.rho_e:0.3f} g/cm^3'
    attributes['earth_gravity_constant'] = f'{factors.GM:0.3f} cm^3/s^2'
    # output suffix for data formats
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')[DATAFORM]

    # for each input file
    for seismic_file in sorted(input_files):
        logging.info(str(seismic_file))
        # add attributes for input file
        attributes['lineage'] = seismic_file.name
        # extract date from file name
        calendar_date = rx.findall(seismic_file.name).pop()
        d = datetime.datetime.strptime(calendar_date, '%b%Y')
        # seismic files have a decoding issue
        with seismic_file.open('rb') as f_in:
            contents = f_in.read().splitlines()
            fid = io.StringIO()
            fid.write(b'\n'.join(contents[12:]).decode('utf-8'))
            fid.seek(0)
        # read data as spatial object
        dinput = gravtk.spatial().from_ascii(fid,
            compression='bytes', date=False,
            spacing=[1.0, 1.0], nlon=361, nlat=181,
            columns=['lon','lat','data'])
        # convert to spherical harmonics from microGal
        Ylms = gravtk.gen_stokes(dinput.data, dinput.lon, dinput.lat,
            LMAX=LMAX, MMAX=MMAX, UNITS=factors.microGal)
        # add time and grace month variables
        Ylms.time = gravtk.time.convert_calendar_decimal(d.year, d.month)
        Ylms.month = gravtk.time.calendar_to_grace(d.year, d.month)
        # add attributes to output harmonics
        Ylms.attributes['ROOT'] = attributes
        # output file name
        f = f'seismic_mv{mv}_CLM_L{LMAX:d}{order_str}_{Ylms.month:03d}.{suffix}'
        output_file = ddir.joinpath(f)
        # write spherical harmonics to file
        Ylms.to_file(output_file, format=DATAFORM)
        # set the permissions level of the output file to MODE
        output_file.chmod(mode=MODE)

    # open output date and index files
    output_date_file = ddir.joinpath(f'SEISMIC_mv{mv}_DATES.txt')
    fid1 = output_date_file.open(mode='w', encoding='utf8')
    output_index_file = ddir.joinpath(f'index.txt')
    fid2 = output_index_file.open(mode='w', encoding='utf8')
    # date file header information
    print('{0:8} {1:10}'.format('Month','Date'), file=fid1)

    # output file format for spherical harmonic data
    args = (mv, LMAX, order_str, suffix)
    regex_pattern = r'seismic_mv{0}_CLM_L{1:d}{2}_(\d+).{3}'
    output_regex = re.compile(regex_pattern.format(*args))
    # find all output harmonic files (not just ones created in run)
    output_files = [f for f in ddir.iterdir()
        if re.match(output_regex, f.name)]
    for fi in sorted(output_files):
        # extract GRACE month
        grace_month, = np.array(re.findall(output_regex,fi.name), dtype=int)
        YY,MM = gravtk.time.grace_to_calendar(grace_month)
        tdec, = gravtk.time.convert_calendar_decimal(YY, MM)
        # print date, GRACE month and calendar month to date file
        fid1.write(f'{tdec:11.6f} {grace_month:03d} {MM:02.0f}\n')
        # print output file to index
        full_output_file = gravtk.spatial().compressuser(fi)
        print(full_output_file, file=fid2)
    # close the date and index files
    fid1.close()
    fid2.close()
    # set the permissions level of the output date and index files to MODE
    output_date_file.chmod(mode=MODE)
    output_index_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads monthly gravity changes due to earthquakes and
            converts to spherical harmonic coefficients
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    choices = ['SHdegree40_M3_microGal','SHdegree40_M5_microGal',
        'SHdegree60_M3_microGal','SHdegree60_M5_microGal']
    parser.add_argument('model',
        type=str, metavar='MODEL', nargs='+', choices=choices,
        help='Seismic Model')
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
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

    # for each reanalysis model
    for MODEL in args.model:
        # run program
        seismic_monthly_harmonics(args.directory, MODEL,
            LMAX=args.lmax, MMAX=args.mmax,
            LOVE_NUMBERS=args.love, REFERENCE=args.reference,
            DATAFORM=args.format, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
