#!/usr/bin/env python
u"""
model_level_coefficients.py
Written by Tyler Sutterley (12/2022)
Creates a netCDF4 file of reanalysis A and B coefficients for model levels
Model level coefficients are obtained using equation 3.17 of
    Simmons and Burridge (1981) and the methodology of Trenberth et al (1993)

ERA-Interim coefficients:
https://rda.ucar.edu/datasets/ds627.1/docs/Eta_coordinate/
https://rda.ucar.edu/datasets/ds627.1/docs/Eta_coordinate/ERA-Interim_coordvars.nc

ERA5 coefficients:
https://doi.org/10.1175/1520-0493(1981)109<0758:AEAAMC>2.0.CO;2
https://www.ecmwf.int/en/forecasts/documentation-and-support/137-model-levels

MERRA-2 coefficients:
https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf
http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_vertical_grids

INPUTS:
    Reanalysis model to run
    ERA-Interim: http://apps.ecmwf.int/datasets/data/interim-full-moda
    ERA5: http://apps.ecmwf.int/data-catalogues/era5/?class=ea
    MERRA-2: https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/

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
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 12/2020: using argparse to set command line options
    Written 03/2018
"""
from __future__ import print_function

import sys
import os
import time
import netCDF4
import argparse
import numpy as np
import model_harmonics as mdlhmc

# PURPOSE: create netCDF4 file with the model level A and B coefficients
def model_level_coefficients(base_dir, MODEL, MODE=0o775):
    # directory setup
    ddir = os.path.join(base_dir,MODEL)
    if (MODEL == 'ERA5'):
        # output netCDF4 file
        output_coordinate_file = 'ERA5_coordvars.nc'
        # read input file
        dinput = np.loadtxt(os.path.join(ddir,'ERA5_coordvars.txt'))
        # create output dictionary with variables
        output = {}
        # interfaces
        output['intf'] = dinput[:,0]
        # half levels
        output['lvl'] = 0.5 + dinput[0:-1,0]
        # extract A and B coefficients
        output['a_interface'] = dinput[:,1]
        output['b_interface'] = dinput[:,2]
        output['a_half']=(output['a_interface'][1:]+output['a_interface'][0:-1])/2.0
        output['b_half']=(output['b_interface'][1:]+output['b_interface'][0:-1])/2.0
    elif (MODEL == 'MERRA-2'):
        # output netCDF4 file
        output_coordinate_file = 'MERRA2_101.Coords_Nx.00000000.nc'
        # python dictionary with output variables
        output = {}
        # Ap [millibars] for 72 levels (73 edges)
        Ap = np.array([0.000000e00, 4.804826e-02, 6.593752e00, 1.313480e01,
            1.961311e01, 2.609201e01, 3.257081e01, 3.898201e01,
            4.533901e01, 5.169611e01, 5.805321e01, 6.436264e01,
            7.062198e01, 7.883422e01, 8.909992e01, 9.936521e01,
            1.091817e02, 1.189586e02, 1.286959e02, 1.429100e02,
            1.562600e02, 1.696090e02, 1.816190e02, 1.930970e02,
            2.032590e02, 2.121500e02, 2.187760e02, 2.238980e02,
            2.243630e02, 2.168650e02, 2.011920e02, 1.769300e02,
            1.503930e02, 1.278370e02, 1.086630e02, 9.236572e01,
            7.851231e01, 6.660341e01, 5.638791e01, 4.764391e01,
            4.017541e01, 3.381001e01, 2.836781e01, 2.373041e01,
            1.979160e01, 1.645710e01, 1.364340e01, 1.127690e01,
            9.292942e00, 7.619842e00, 6.216801e00, 5.046801e00,
            4.076571e00, 3.276431e00, 2.620211e00, 2.084970e00,
            1.650790e00, 1.300510e00, 1.019440e00, 7.951341e-01,
            6.167791e-01, 4.758061e-01, 3.650411e-01, 2.785261e-01,
            2.113490e-01, 1.594950e-01, 1.197030e-01, 8.934502e-02,
            6.600001e-02, 4.758501e-02, 3.270000e-02, 2.000000e-02,
            1.000000e-02])
        # invert so top-of-atmosphere == layer 1
        # convert units from millibars to pascals
        output['a_interface'] = 100.0*Ap[::-1]
        # Ap at half levels
        output['a_half']=(output['a_interface'][1:]+output['a_interface'][0:-1])/2.0
        nlevels = len(Ap)
        # half levels
        output['lev'] = 0.5 + np.arange(nlevels-1)
        # interfaces
        output['intf'] = np.arange(nlevels) + 1
        # Bp [unitless] for 72 levels (73 edges)
        Bp = np.array([1.000000e00, 9.849520e-01, 9.634060e-01, 9.418650e-01,
            9.203870e-01, 8.989080e-01, 8.774290e-01, 8.560180e-01,
            8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,
            7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01,
            6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01,
            4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,
            2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01,
            6.372006e-02, 2.801004e-02, 6.960025e-03, 8.175413e-09,
            0.000000e00, 0.000000e00, 0.000000e00, 0.000000e00,
            0.000000e00, 0.000000e00, 0.000000e00, 0.000000e00,
            0.000000e00, 0.000000e00, 0.000000e00, 0.000000e00,
            0.000000e00, 0.000000e00, 0.000000e00, 0.000000e00,
            0.000000e00, 0.000000e00, 0.000000e00, 0.000000e00,
            0.000000e00, 0.000000e00, 0.000000e00, 0.000000e00,
            0.000000e00, 0.000000e00, 0.000000e00, 0.000000e00,
            0.000000e00, 0.000000e00, 0.000000e00, 0.000000e00,
            0.000000e00, 0.000000e00, 0.000000e00, 0.000000e00,
            0.000000e00, 0.000000e00, 0.000000e00, 0.000000e00,
            0.000000e00])
        # invert so top-of-atmosphere == layer 1
        output['b_interface'] = Bp[::-1]
        # Bp at half levels
        output['b_half']=(output['b_interface'][1:]+output['b_interface'][0:-1])/2.0

    # output coefficients to netCDF4 file
    fileID = netCDF4.Dataset(os.path.join(ddir,output_coordinate_file),'w')
    # Defining the NetCDF4 dimensions and creating dimension variables
    nc = {}
    for key in ['lvl','intf']:
        fileID.createDimension(key, len(output[key]))
        nc[key] = fileID.createVariable(key, output[key].dtype, (key,))
    # creating the half-layer NetCDF4 variables
    for key in ['a_half','b_half']:
        nc[key] = fileID.createVariable(key, output[key].dtype, ('lvl',))
    # creating the interface NetCDF4 variables
    for key in ['a_interface','b_interface']:
        nc[key] = fileID.createVariable(key, output[key].dtype, ('intf',))
    # filling NetCDF4 variables
    for key,val in output.items():
        nc[key][:] = np.copy(val)
    # add software information
    fileID.software_reference = mdlhmc.version.project_name
    fileID.software_version = mdlhmc.version.full_version
    # date created
    fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())
    # close the netCDF4 file
    fileID.close()
    # change the permissions level to MODE
    os.chmod(os.path.join(ddir,output_coordinate_file), MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Creates a netCDF4 file of reanalysis
            A and B coefficients for model levels
            """
    )
    # command line parameters
    choices = ['ERA5','MERRA-2']
    parser.add_argument('model',
        type=str, nargs='+',
        default=['ERA5','MERRA-2'], choices=choices,
        help='Reanalysis Model')
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # permissions mode of the directories and files retrieved
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files retrieved')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # run program
    for MODEL in args.model:
        model_level_coefficients(args.directory, MODEL, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
