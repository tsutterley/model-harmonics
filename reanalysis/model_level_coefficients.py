#!/usr/bin/env python
"""
model_level_coefficients.py
Written by Tyler Sutterley (07/2026)
Creates a netCDF4 file of reanalysis A and B coefficients for model levels
Model level coefficients are obtained using equation 3.17 of
    Simmons and Burridge (1981) and the methodology of Trenberth et al (1993)

ERA-Interim coefficients:
https://rda.ucar.edu/datasets/ds627.1/docs/Eta_coordinate/
https://rda.ucar.edu/datasets/ds627.1/docs/Eta_coordinate/ERA-Interim_coordvars.nc

ERA5 coefficients:
https://doi.org/10.1175/1520-0493(1981)109<0758:AEAAMC>2.0.CO;2
https://www.ecmwf.int/en/forecasts/documentation-and-support/137-model-levels
https://confluence.ecmwf.int/spaces/UDOC/pages/108117123/L137+model+level+definitions

MERRA-2 coefficients:
https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf
http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_vertical_grids
https://geos-chem.readthedocs.io/en/latest/supplemental-guides/vertical-grids.html

JRA-3Q coefficients:
https://www.data.jma.go.jp/jra/html/JRA-3Q/document/JRA-3Q_TL479_format_v1_en.pdf
https://osdata.gdex.ucar.edu/web/datasets/d640000/docs/JRA-3Q_TL479_format_en.pdf

Reanalysis models:
    ERA-Interim:
        http://apps.ecmwf.int/datasets/data/interim-full-moda
    ERA5:
        http://apps.ecmwf.int/data-catalogues/era5/?class=ea
    MERRA-2:
        https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/
    JRA-3Q:
        https://gdex.ucar.edu/datasets/d640000/
        https://www.data.jma.go.jp/jra/html/JRA-3Q/index_en.html

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -M X, --mode X: Permission mode of directories and files

REFERENCES:
    A. J. Simmons, and D. M. Burridge, "An energy and angular-momentum
        conserving finite-difference scheme and hybrid vertical coordinates."
        Monthly Weather Review., 109, 758-766, (1981).
        https://doi.org/10.1175/1520-0493(1981)109<0758:AEAAMC>2.0.CO;2

    K. E. Trenberth, J. C. Berry, and L. E. Buja, "Vertical interpolation
        and truncation of model-coordinate data."
        NCAR Technical Note NCAR/TN-396+STR, 54 pp., (1993).
        https://doi.org/10.5065/D6HX19NH

UPDATE HISTORY:
    Updated 07/2026: use struct dictionary to define netCDF4 parameters
        add JRA-3Q reanalysis to list of optional models to run
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 12/2020: using argparse to set command line options
    Written 03/2018
"""

from __future__ import print_function

import sys
import pathlib
import argparse
import numpy as np
from model_harmonics.utilities import get_data_path
import model_harmonics as mdlhmc


# PURPOSE: create netCDF4 file with the model level A and B coefficients
def model_level_coefficients(base_dir, MODEL, MODE=0o775):
    # directory setup
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    ddir = base_dir.joinpath(MODEL)
    ddir.mkdir(mode=MODE, parents=True, exist_ok=True)

    # create output dictionary with variables
    output = {}
    # model parameters
    # input and output files
    # output name for half-levels and interfaces
    # extract A and B coefficients
    if MODEL in ('ERA-Interim', 'ERA5'):
        # input and output coordinate files
        input_file = get_data_path(['data', f'{MODEL}_coordvars.csv'])
        filename = f'{MODEL}_coordvars.nc'
        # read input file
        dinput = np.loadtxt(input_file, delimiter=',', skiprows=1)
        # output level names
        INTERFACE = 'intf'
        LEVELNAME = 'lvl'
        ANAME, BNAME = ('a_half', 'b_half')
        AINTERFACE, BINTERFACE = ('a_interface', 'b_interface')
        # extract levels
        output[INTERFACE] = dinput[:, 0]
        output[LEVELNAME] = 0.5 + dinput[0:-1, 0]
        # extract A and B coefficients
        output[AINTERFACE] = dinput[:, 1]
        output[BINTERFACE] = dinput[:, 2]
        # calculate half-level A and B coefficients from interfaces
        output[ANAME] = (dinput[1:, 1] + dinput[:-1, 1]) / 2.0
        output[BNAME] = (dinput[1:, 2] + dinput[:-1, 2]) / 2.0
    elif MODEL == 'MERRA-2':
        # input and output coordinate files
        input_file = get_data_path(['data', 'GMAO_72-level_vertical_grid.csv'])
        filename = 'MERRA2_101.const_3d_coords_Nx.00000000.nc4'
        # read input file
        dinput = np.loadtxt(input_file, delimiter=',', skiprows=1)
        # invert so top-of-atmosphere == layer 1
        dinput = np.flipud(dinput)
        # output level names
        INTERFACE = 'intf'
        LEVELNAME = 'lev'
        ANAME, BNAME = ('a_half', 'b_half')
        AINTERFACE, BINTERFACE = ('a_interface', 'b_interface')
        # extract levels
        output[INTERFACE] = dinput[:, 0]
        output[LEVELNAME] = 0.5 + dinput[0:-1, 0]
        # extract A and B coefficients
        # convert units from millibars to pascals
        # Ap [pascals] for 72 levels (73 edges)
        # Bp [unitless] for 72 levels (73 edges)
        output[AINTERFACE] = 100.0 * dinput[:, 1]
        output[BINTERFACE] = dinput[:, 2]
        # calculate half-level A and B coefficients from interfaces
        output[ANAME] = 100.0 * (dinput[1:, 1] + dinput[:-1, 1]) / 2.0
        output[BNAME] = (dinput[1:, 2] + dinput[:-1, 2]) / 2.0
    elif MODEL == 'JRA-3Q':
        # input and output coordinate files
        input_file = get_data_path(['data', 'jra3q.mdl-hyb.coefficients.csv'])
        filename = 'jra3q.mdl-hyb.coefficients.nc'
        # read input file
        dinput = np.loadtxt(input_file, delimiter=',', skiprows=1)
        # invert so top-of-atmosphere == layer 1
        dinput = np.flipud(dinput)
        # output level names
        INTERFACE = 'hybrid_level'
        LEVELNAME = 'hybrid_half_level'
        ANAME, BNAME = ('a_hybrid_half_level', 'b_hybrid_half_level')
        AINTERFACE, BINTERFACE = ('a_hybrid_level', 'b_hybrid_level')
        # extract levels
        output[INTERFACE] = dinput[:-1, 5]
        output[LEVELNAME] = dinput[:, 3]
        # extract  half-level A and B coefficients
        output[ANAME] = dinput[:, 1]
        output[BNAME] = dinput[:, 2]
        # calculate A and B coefficients from half-levels
        output[AINTERFACE] = (dinput[1:, 1] + dinput[:-1, 1]) / 2.0
        output[BINTERFACE] = (dinput[1:, 2] + dinput[:-1, 2]) / 2.0

    # dictionary defining output structure
    struct = dict(
        dimensions=(LEVELNAME, INTERFACE),
        variables={
            ANAME: (LEVELNAME,),
            BNAME: (LEVELNAME,),
            AINTERFACE: (INTERFACE,),
            BINTERFACE: (INTERFACE,),
        },
    )

    # dictionary defining file-level and variable attributes
    attributes = dict(ROOT={})
    # reference attribute
    REFERENCE = f'Output from {pathlib.Path(sys.argv[0]).name}'
    attributes['ROOT']['reference'] = REFERENCE
    # Defining attributes for model levels
    attributes[LEVELNAME] = dict(
        long_name='Model Level Number',
        units='1',
    )
    attributes[INTERFACE] = dict(
        long_name='Model Level Interfaces',
        units='1',
    )
    attributes[ANAME] = dict(
        long_name='A coefficients for model levels',
        units='Pa',
    )
    attributes[BNAME] = dict(
        long_name='B coefficients for model levels',
        units='1',
    )
    attributes[AINTERFACE] = dict(
        long_name='A coefficients at model level interfaces',
        units='Pa',
    )
    attributes[BINTERFACE] = dict(
        long_name='B coefficients at model level interfaces',
        units='1',
    )

    # write structured data to netCDF4 file
    output_file = ddir.joinpath(filename)
    mdlhmc.spatial.to_netCDF4(output_file, output, attributes, struct)
    # change the permissions level to MODE
    output_file.chmod(mode=MODE)


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Creates a netCDF4 file of reanalysis
            A and B coefficients for model levels
            """
    )
    # command line parameters
    choices = ['ERA-Interim', 'ERA5', 'JRA-3Q', 'MERRA-2']
    parser.add_argument(
        'model',
        type=str,
        nargs='+',
        metavar='MODEL',
        default=['ERA5', 'MERRA-2'],
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
    # permissions mode of the directories and files retrieved
    parser.add_argument(
        '--mode',
        '-M',
        type=lambda x: int(x, base=8),
        default=0o775,
        help='Permission mode of directories and files retrieved',
    )
    # return the parser
    return parser


# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args, _ = parser.parse_known_args()

    # run program
    for MODEL in args.model:
        model_level_coefficients(args.directory, MODEL, MODE=args.mode)


# run main program
if __name__ == '__main__':
    main()
