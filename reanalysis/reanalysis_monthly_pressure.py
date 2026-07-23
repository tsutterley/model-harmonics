#!/usr/bin/env python
"""
reanalysis_monthly_pressure.py
Written by Tyler Sutterley (05/2023)
Reads daily atmospheric pressure fields from reanalysis and outputs monthly averages

INPUTS:
    Reanalysis model to run
    NCEP-DOE-2: https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis2.html

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -Y X, --year X: years to run
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
    spatial.py: spatial data class for reading, writing and processing data
    time.py: utilities for calculating time operations
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 07/2026: output using gravity_toolkit spatial class
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: use full path to output file in verbose logging
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: can use input files to define command line arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 03/2021: automatically update years to run based on current time
    Updated 12/2020: using argparse to set command line options
        using spatial module for operations
    Updated 01/2020: changed year option to be specific years to run
    Written 03/2018
"""

from __future__ import print_function

import re
import logging
import pathlib
import argparse
import datetime
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc


# PURPOSE: read atmospheric surface pressure fields and calculate monthly mean
def reanalysis_monthly_pressure(base_dir, MODEL, YEARS, MODE=0o775):
    # directory setup
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    ddir = base_dir.joinpath(MODEL)

    # field mapping for different products
    field_mapping = {}
    # set model specific parameters
    if MODEL == 'NCEP-DOE-2':
        # regular expression pattern for finding files
        regex_pattern = r'pres.sfc.({0:4d}).nc$'
        FILL_VALUE = 'missing_value'
        # output file format
        output_file_format = 'pres.sfc.mon.mean.{0:4d}.nc'
        field_mapping['data'] = 'pres'
        field_mapping['lon'] = 'lon'
        field_mapping['lat'] = 'lat'
        field_mapping['time'] = 'time'

    # for each year
    for YEAR in YEARS:
        # days per month (check if year is a leap year)
        dpm = [0.0] + list(gravtk.time.calendar_days(YEAR))
        cumulative_days = np.cumsum(dpm)
        # list of spatial data
        p_mean = []
        # read each reanalysis pressure field and calculate mean
        rx = re.compile(regex_pattern.format(YEAR), re.VERBOSE)
        input_files = sorted([f for f in ddir.iterdir() if rx.match(f.name)])
        # for each input file
        for i, input_file in enumerate(input_files):
            # read input data
            p = gravtk.spatial().from_netCDF4(
                input_file,
                field_mapping=field_mapping,
            )
            # reorder dimensions to match the required order
            p = p.transpose(axes=(1, 2, 0))
            p.fill_value = p.attributes['data'][FILL_VALUE]
            # iterate over months
            for m in range(0, 12):
                # for each day in the month
                indices = np.arange(cumulative_days[m], cumulative_days[m + 1])
                try:
                    p_mean.append(p.mean(indices=indices.astype(np.int64)))
                except:
                    break
                else:
                    p_mean[m].month = m + 1

        # mean pressure for each month
        p_month = gravtk.spatial().from_list(p_mean)
        # copy global attributes from input file
        p_month.attributes = p.attributes
        # copy attributes for each output field
        for field, key in field_mapping.items():
            p_month.attributes[key] = p.attributes.get(field)
        # save to file
        filename = output_file_format.format(YEAR)
        output = ddir.joinpath(filename)
        logging.info(f'\t{output}')
        # write monthly data to netCDF4 file
        p_month.to_netCDF4(
            output,
            field_mapping=field_mapping,
            attributes=p_month.attributes,
            clobber=True,
        )
        # set the permissions level of the output file to MODE
        output.chmod(mode=MODE)


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads daily atmospheric pressure fields
            from reanalysis and outputs monthly averages
            """,
        fromfile_prefix_chars='@',
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    choices = ['NCEP-DOE-2']
    parser.add_argument(
        'model',
        type=str,
        nargs='+',
        default=['NCEP-DOE-2'],
        choices=choices,
        help='Reanalysis Model',
    )
    # working data directory
    parser.add_argument(
        '--directory',
        '-D',
        type=pathlib.Path,
        default=pathlib.Path.cwd(),
        help='Working data directory',
    )
    # years to run
    now = datetime.datetime.now()
    parser.add_argument(
        '--year',
        '-Y',
        type=int,
        nargs='+',
        default=range(2000, now.year + 1),
        help='Years of model outputs to run',
    )
    # print information about each input and output file
    parser.add_argument(
        '--verbose',
        '-V',
        action='count',
        default=0,
        help='Verbose output of processing run',
    )
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument(
        '--mode',
        '-M',
        type=lambda x: int(x, base=8),
        default=0o775,
        help='Permission mode of directories and files',
    )
    # return the parser
    return parser


# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args, _ = parser.parse_known_args()

    # create logger
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # for each reanalysis model
    for MODEL in args.model:
        # run program
        reanalysis_monthly_pressure(
            args.directory, MODEL, args.year, MODE=args.mode
        )


# run main program
if __name__ == '__main__':
    main()
