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

Reanalysis models:
    ERA-Interim:
        http://apps.ecmwf.int/datasets/data/interim-full-moda
    ERA5:
        http://apps.ecmwf.int/data-catalogues/era5/?class=ea
    MERRA-2:
        https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/

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
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 12/2020: using argparse to set command line options
    Written 03/2018
"""

from __future__ import print_function

import time
import netCDF4
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
        # extract A and B coefficients
        Ap = dinput[:, 1]
        Bp = dinput[:, 2]
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
        # extract A and B coefficients
        # convert units from millibars to pascals
        # Ap [pascals] for 72 levels (73 edges)
        # Bp [unitless] for 72 levels (73 edges)
        Ap = 100.0 * dinput[:, 1]
        Bp = dinput[:, 2]

    # create output dictionary with variables
    output = {}
    # interfaces and half-levels
    output[INTERFACE] = dinput[:, 0]
    output[LEVELNAME] = 0.5 + dinput[0:-1, 0]
    # add A and B coefficients to output dictionary
    output['a_interface'] = Ap
    output['b_interface'] = Bp
    output['a_half'] = (Ap[1:] + Ap[:-1]) / 2.0
    output['b_half'] = (Bp[1:] + Bp[:-1]) / 2.0

    # dictionary defining output structure
    struct = dict(
        dimensions=(LEVELNAME, INTERFACE),
        variables={
            'a_half': (LEVELNAME,),
            'b_half': (LEVELNAME,),
            'a_interface': (INTERFACE,),
            'b_interface': (INTERFACE,),
        },
    )

    # dictionary defining file-level and variable attributes
    attributes = dict(ROOT={})
    # Defining attributes for model levels
    attributes[LEVELNAME] = dict(
        long_name='Model Level Number',
        units='1',
    )
    attributes[INTERFACE] = dict(
        long_name='Model Level Interfaces',
        units='1',
    )
    attributes['a_half'] = dict(
        long_name='A coefficients for model levels',
        units='Pa',
    )
    attributes['b_half'] = dict(
        long_name='B coefficients for model levels',
        units='1',
    )
    attributes['a_interface'] = dict(
        long_name='A coefficients at model level interfaces',
        units='Pa',
    )
    attributes['b_interface'] = dict(
        long_name='B coefficients at model level interfaces',
        units='1',
    )

    # output coefficients to netCDF4 file
    output_file = ddir.joinpath(filename)
    ncdf_model_levels(output, attributes, struct, FILENAME=output_file)
    # change the permissions level to MODE
    output_file.chmod(mode=MODE)


# PURPOSE: write output model levels to file
def ncdf_model_levels(output, attributes, struct, FILENAME=None):
    # opening NetCDF file for writing
    FILENAME = pathlib.Path(FILENAME).expanduser().absolute()
    fileID = netCDF4.Dataset(FILENAME, 'w', format='NETCDF4')
    # dictionary with NetCDF4 variable objects
    nc = {}
    # defining the NetCDF4 dimensions
    for dim in struct['dimensions']:
        fileID.createDimension(dim, len(output[dim]))
        nc[dim] = fileID.createVariable(dim, output[dim].dtype, (dim,))
        # add data to NetCDF4 dimension variable
        nc[dim][:] = output[dim].copy()
        # set netCDF4 attributes for dimensions
        for att_name, att_val in attributes[dim].items():
            nc[dim].setncattr(att_name, att_val)

    # defining the NetCDF4 variables
    for var, dimensions in struct['variables'].items():
        nc[var] = fileID.createVariable(
            var,
            output[var].dtype,
            dimensions,
        )
        # add data to NetCDF4 variable
        nc[var][:] = output[var].copy()
        # set netCDF4 attributes for variables
        for att_name, att_val in attributes[var].items():
            nc[var].setncattr(att_name, att_val)

    # Defining file-level attributes
    for att_name, att_val in attributes['ROOT'].items():
        fileID.setncattr(att_name, att_val)

    # add software information
    fileID.software_reference = mdlhmc.version.project_name
    fileID.software_version = mdlhmc.version.full_version
    # date created
    fileID.date_created = time.strftime('%Y-%m-%d', time.localtime())
    # close the netCDF4 file
    fileID.close()
    # clear nc dictionary variable
    nc = None


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Creates a netCDF4 file of reanalysis
            A and B coefficients for model levels
            """
    )
    # command line parameters
    choices = ['ERA5', 'MERRA-2']
    parser.add_argument(
        'model',
        type=str,
        nargs='+',
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
