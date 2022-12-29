#!/usr/bin/env python
u"""
merra_hybrid_regrid.py
Written by Tyler Sutterley (12/2022)
Read and regrid MERRA-2 hybrid variables
MERRA-2 Hybrid firn model outputs provided by Brooke Medley at GSFC

CALLING SEQUENCE:
    python merra_hybrid_regrid.py --directory <path> --region gris \
        --version v1.1 --product SMB_a --verbose

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -R X, --region X: Region to calculate (gris, ais)
    -v X, --version X: Version of firn model to calculate
        v0
        v1
        v1.0
        v1.1
        v1.2
        v1.2.1
    -P X, --product X: MERRA-2 hybrid product to calculate
    -Y X, --year X: Years to run
    --mask X: netCDF4 mask files for reducing to regions
    -S X, --spacing X: spatial resolution of output data (dlon,dlat)
    -I X, --interval X: output grid interval
        1: (0:360, 90:-90)
        2: (degree spacing/2)
        3: non-global grid (set with defined bounds)
    -B X, --bounds X: non-global grid bounding box (minlon,maxlon,minlat,maxlat)
    -G, --gzip: input netCDF4 file is gzip compressed
    -V, --verbose: Output information for each output file
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    utilities.py: download and management utilities for files
    constants.py: calculate reference parameters for common ellipsoids
    spatial.py: spatial data class for reading, writing and processing data

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 10/2022: move polar stereographic scaling function to spatial
        add Greenland and Antarctic versions v1.2.1
    Updated 06/2022: change default variables to include firn height anomaly
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 12/2021: open MERRA-2 hybrid product command line options
        added GSFC MERRA-2 Hybrid Greenland v1.2
        can use variable loglevels for verbose output
    Updated 10/2021: add pole case in stereographic area scale calculation
        using python logging for handling verbose output
    Updated 09/2021: use original FDM file for ais products
        use original FDM file for non-mass variables (height and FAC)
    Written 09/2021
"""
from __future__ import print_function

import sys
import os
import copy
import gzip
import uuid
import time
import pyproj
import logging
import netCDF4
import argparse
import warnings
import numpy as np
import sklearn.neighbors
import scipy.interpolate
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc
# ignore pyproj and divide by zero warnings
warnings.filterwarnings("ignore")

# PURPOSE: set the projection parameters based on the region name
def set_projection(REGION):
    if (REGION == 'ais'):
        projection_flag = 'EPSG:3031'
    elif (REGION == 'gris'):
        projection_flag = 'EPSG:3413'
    return projection_flag

# PURPOSE: read and regrid MERRA-2 hybrid SMB estimates
def merra_hybrid_regrid(base_dir, REGION, VARIABLE, YEARS,
    VERSION='v1',
    MASKS=None,
    DDEG=None,
    INTERVAL=None,
    BOUNDS=None,
    GZIP=False,
    MODE=0o775):

    # MERRA-2 hybrid directory
    DIRECTORY = os.path.join(base_dir,VERSION)
    # suffix if compressed
    suffix = '.gz' if GZIP else ''
    # set the input netCDF4 file for the variable of interest
    if VARIABLE in ('cum_smb_anomaly',):
        FILE_VERSION = copy.copy(VERSION)
        args = (VERSION,REGION.lower(),suffix)
        hybrid_file = 'gsfc_fdm_{0}_{1}.nc{2}'.format(*args)
    elif VARIABLE in ('FAC','height','h_a') or (REGION.lower() == 'ais'):
        FILE_VERSION = VERSION.replace('.','_')
        args = (FILE_VERSION,REGION.lower(),suffix)
        hybrid_file = 'gsfc_fdm_{0}_{1}.nc{2}'.format(*args)
    elif VARIABLE in ('Me_a','Ra_a','Ru_a','Sn-Ev_a','SMB_a'):
        FILE_VERSION = VERSION.replace('.','_')
        args = (FILE_VERSION,REGION.lower(),suffix)
        hybrid_file = 'gsfc_fdm_smb_cumul_{0}_{1}.nc{2}'.format(*args)
    else:
        raise ValueError(f'Unknown variable {VARIABLE}')

    # reference for GSFC-FDM model outputs
    reference = ("Medley, B., Neumann, T. A., Zwally, H. J., "
        "Smith, B. E., and Stevens, C. M.: Simulations of Firn Processes "
        "over the Greenland and Antarctic Ice Sheets: 1980--2021, "
        "The Cryosphere, https://doi.org/10.5194/tc-2020-266, 2022.")
    # Open the MERRA-2 Hybrid NetCDF file for reading
    if GZIP:
        # read as in-memory (diskless) netCDF4 dataset
        with gzip.open(os.path.join(DIRECTORY,hybrid_file),'r') as f:
            fileID = netCDF4.Dataset(uuid.uuid4().hex, memory=f.read())
    else:
        # read netCDF4 dataset
        fileID = netCDF4.Dataset(os.path.join(DIRECTORY,hybrid_file), 'r')

    # Output NetCDF file information
    logging.info(os.path.join(DIRECTORY,hybrid_file))
    logging.info(list(fileID.variables.keys()))

    # Get data from each netCDF variable and remove singleton dimensions
    fd = {}
    # time is year decimal at time step 5 days
    time_step = 5.0/365.25
    # reduce grids to time period of input buffered by time steps
    tmin = np.min(YEARS) - 2.0*time_step
    tmax = np.max(YEARS) + 2.0*time_step
    # find indices to times
    nt, = fileID.variables['time'].shape
    f = scipy.interpolate.interp1d(fileID.variables['time'][:],
        np.arange(nt), kind='nearest', bounds_error=False,
        fill_value=(0,nt))
    imin,imax = f((tmin,tmax)).astype(np.int64)
    # read reduced time variables
    fd['time'] = fileID.variables['time'][imin:imax+1].copy()
    # read reduced dataset and remove singleton dimensions
    fd[VARIABLE] = np.squeeze(fileID.variables[VARIABLE][imin:imax+1,:,:])
    # invalid data value
    fv = np.float64(fileID.variables[VARIABLE]._FillValue)
    # input shape of MERRA-2 Hybrid firn data
    nt,nx,ny = np.shape(fd[VARIABLE])
    # extract x and y coordinate arrays from grids if applicable
    # else create meshgrids of coordinate arrays
    if (np.ndim(fileID.variables['x'][:]) == 2):
        xg = fileID.variables['x'][:].copy()
        yg = fileID.variables['y'][:].copy()
        fd['x'],fd['y'] = (xg[:,0],yg[0,:])
    else:
        fd['x'] = fileID.variables['x'][:].copy()
        fd['y'] = fileID.variables['y'][:].copy()
        xg,yg = np.meshgrid(fd['x'],fd['y'],indexing='ij')
    # extract area of each grid cell if applicable
    # calculate using dimensions if not possible
    fd['area'] = np.zeros((nx,ny))
    try:
        # ice covered area
        fd['area'][:,:] = fileID.variables['iArea'][:,:].copy()
    except:
        # calculate grid areas (assume fully ice covered)
        dx = np.abs(fd['x'][1] - fd['x'][0])
        dy = np.abs(fd['y'][1] - fd['y'][0])
        fd['area'][:,:] = dx*dy
    # close the NetCDF files
    fileID.close()

    # create mask object for reducing data
    if not MASKS:
        fd['mask'] = np.ones((nx,ny),dtype=bool)
    else:
        fd['mask'] = np.zeros((nx,ny),dtype=bool)
    # read masks for reducing regions before regridding
    for mask_file in MASKS:
        fileID = netCDF4.Dataset(mask_file,'r')
        fd['mask'] |= fileID.variables['mask'][:].astype(bool)
        fileID.close()
    # indices of valid MERRA hybrid data
    fd['mask'] &= (fd[VARIABLE].data[0,:,:] != fv)

    # pyproj transformer for converting to input coordinates (EPSG)
    MODEL_EPSG = set_projection(REGION)
    crs1 = pyproj.CRS.from_string('EPSG:4326')
    crs2 = pyproj.CRS.from_string(MODEL_EPSG)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    direction = pyproj.enums.TransformDirection.INVERSE
    # convert projection from model coordinates
    gridlon,gridlat = transformer.transform(xg, yg, direction=direction)
    # polar stereographic standard parallel (latitude of true scale)
    reference_latitude = crs2.to_dict().pop('lat_ts')

    # get reference parameters for ellipsoid
    ellipsoid_params = mdlhmc.constants(ellipsoid='WGS84', units='CGS')
    # semimajor axis of ellipsoid [cm]
    a_axis = ellipsoid_params.a_axis
    # first numerical eccentricity
    ecc1 = ellipsoid_params.ecc1
    # flattening of the ellipsoid
    flat = ellipsoid_params.flat
    # Average Radius of the Earth with equal surface area [cm]
    rad_e = ellipsoid_params.rad_e
    # convert from geodetic latitude to geocentric latitude
    # geodetic latitude in radians
    latitude_geodetic_rad = np.pi*gridlat/180.0
    # prime vertical radius of curvature
    N = a_axis/np.sqrt(1.0 - ecc1**2.*np.sin(latitude_geodetic_rad)**2.)
    # calculate X, Y and Z from geodetic latitude and longitude
    X = N * np.cos(latitude_geodetic_rad) * np.cos(np.pi*gridlon/180.0)
    Y = N * np.cos(latitude_geodetic_rad) * np.sin(np.pi*gridlon/180.0)
    Z = (N * (1.0 - ecc1**2.0)) * np.sin(latitude_geodetic_rad)
    # calculate geocentric latitude and convert to degrees
    latitude_geocentric = 180.0*np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))/np.pi

    # reduce latitude and longitude to valid and masked points
    indx,indy = np.nonzero(fd['mask'])
    lon,lat = (gridlon[indx,indy],latitude_geocentric[indx,indy])
    # scaled areas
    ps_scale = mdlhmc.spatial.scale_areas(gridlat[indx,indy], flat=flat,
        ref=reference_latitude)
    scaled_area = ps_scale*fd['area'][indx,indy]
    npts = len(scaled_area)
    # densities of meteoric ice [kg/m^3]
    rho_ice = 917.0

    # Output spatial data
    grid = gravtk.spatial(fill_value=fv)
    grid.time = np.zeros((nt))
    grid.month = np.zeros((nt),dtype=np.int64)

    # Output Degree Spacing
    dlon,dlat = (DDEG[0],DDEG[0]) if (len(DDEG) == 1) else (DDEG[0],DDEG[1])
    # Output Degree Interval
    if (INTERVAL == 1):
        # (0:360,90:-90)
        nlon = np.int64((360.0/dlon)+1.0)
        nlat = np.int64((180.0/dlat)+1.0)
        grid.lon = dlon*np.arange(0,nlon)
        grid.lat = 90.0 - dlat*np.arange(0,nlat)
    elif (INTERVAL == 2):
        # (Degree spacing)/2
        grid.lon = np.arange(dlon/2.0,360+dlon/2.0,dlon)
        grid.lat = np.arange(90.0-dlat/2.0,-90.0-dlat/2.0,-dlat)
        nlon = len(grid.lon)
        nlat = len(grid.lat)
    elif (INTERVAL == 3):
        # non-global grid set with BOUNDS parameter
        minlon,maxlon,minlat,maxlat = BOUNDS.copy()
        grid.lon = np.arange(minlon+dlon/2.0,maxlon+dlon/2.0,dlon)
        grid.lat = np.arange(maxlat-dlat/2.0,minlat-dlat/2.0,-dlat)
        nlon = len(grid.lon)
        nlat = len(grid.lat)
    # degree spacing in radians
    dphi = np.pi*dlon/180.0
    dth = np.pi*dlat/180.0
    # meshgrid of output data coordinates
    gridlon,gridlat = np.meshgrid(grid.lon, grid.lat)
    # area of each new grid cells
    grid.area = (rad_e**2)*dphi*dth*np.cos(gridlat*np.pi/180.0)
    # allocate for output data and mask
    grid.data = np.zeros((nlat,nlon,nt))
    grid.mask = np.zeros((nlat,nlon,nt),dtype=bool)
    # update attributes
    grid.update_spacing()
    grid.update_extents()
    grid.update_dimensions()

    # Finding bin indices by finding the closest lat/lon
    # calculates the great circle distance between the
    # model lat/lon and the global lat/lon grids
    INDICE_FILE = 'gsfc_fdm_global_indices.nc'
    if (not os.access(os.path.join(DIRECTORY, INDICE_FILE),os.F_OK)):
        # convert latitude and longitude into radians
        phi,theta = (np.pi*lon/180.0,np.pi*lat/180.0)
        # allocate for output variables (indice maps and grid mask)
        ilon = np.zeros((npts), dtype=np.int32)
        ilat = np.zeros((npts), dtype=np.int32)
        point_mask = np.zeros((nlat,nlon),dtype=bool)
        # creating column arrays of the global grid to quickly
        # compute the great circle distance (all at once for each model pt)
        xphi,xth = (np.pi*gridlon.flatten()/180.0,np.pi*gridlat.flatten()/180.0)
        # create Ball search tree for finding nearest neighbors (great-circle)
        xpoints = np.concatenate((xphi[:,None],xth[:,None]),axis=1)
        tree = sklearn.neighbors.BallTree(xpoints, metric='haversine')
        # find indices of closest grid point for each lat/lon value
        points = np.concatenate((phi[:,None],theta[:,None]),axis=1)
        indices = np.squeeze(tree.query(points,k=1,return_distance=False))
        # longitude indice will be the remainder of the indices/number of lons
        ilon[:] = indices % nlon
        # ilat indice will be number of longitude loops
        ilat[:] = (indices - ilon)/nlon
        # set point mask for lon and lat to True
        point_mask[ilat,ilon] = True

        # Create the NetCDF file for saving indices and points
        fileID=netCDF4.Dataset(os.path.join(DIRECTORY,INDICE_FILE),'w')
        # Defining the netCDF dimensions
        fileID.createDimension('points', npts)
        fileID.createDimension('nlon', nlon)
        fileID.createDimension('nlat', nlat)
        # lat, lon and mask
        nc = {}
        nc['ilon'] = fileID.createVariable('ilon', 'i', ('points',))
        nc['ilat'] = fileID.createVariable('ilat', 'i', ('points',))
        nc['lon'] = fileID.createVariable('lon', 'f', ('nlon',))
        nc['lat'] = fileID.createVariable('lat', 'f', ('nlat',))
        nc['mask'] = fileID.createVariable('mask', 'b', ('nlat','nlon',))
        nc['area'] = fileID.createVariable('area', 'f', ('nlat','nlon',))
        # filling netCDF variables
        nc['ilon'][:] = ilon.copy()
        nc['ilat'][:] = ilat.copy()
        nc['lat'][:] = grid.lat.copy()
        nc['lon'][:] = grid.lon.copy()
        nc['mask'][:,:] = point_mask.copy()
        nc['area'][:,:] = grid.area/1e4
        # Defining attributes for each variable
        nc['ilon'].long_name = 'longitude_indices'
        nc['ilat'].long_name = 'latitude_indices'
        nc['lat'].long_name = 'latitude'
        nc['lat'].units = 'Degrees North'
        nc['lon'].long_name = 'longitude'
        nc['lon'].units = 'Degrees East'
        nc['mask'].long_name = 'point_mask'
        nc['area'].long_name = 'Area'
        nc['area'].units = 'meters squared'
        # global variable of netCDF file
        fileID.title = 'GSFC-FDM remapping'
        fileID.source = f'version {VERSION}'
        fileID.reference = copy.copy(reference)
        # add software information
        fileID.software_reference = mdlhmc.version.project_name
        fileID.software_version = mdlhmc.version.full_version
        # date created
        fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())
        # Closing the netCDF file
        fileID.close()
    else:
        # Open the NetCDF file for reading
        fileID = netCDF4.Dataset(os.path.join(DIRECTORY,INDICE_FILE), 'r')
        # Getting the data from each netCDF variable
        # filling numpy arrays with NetCDF objects
        ilon = np.array(fileID.variables['ilon'][:].copy())
        ilat = np.array(fileID.variables['ilat'][:].copy())
        point_mask=np.array(fileID.variables['mask'][:,:].copy(),dtype=bool)
        # Closing the netCDF file
        fileID.close()

    # for each time step
    for t in range(nt):
        # reduce data for date and convert to mass (g)
        ptms = 1000.0*rho_ice*scaled_area*fd[VARIABLE][t,indx,indy]
        # Sums model point masses nearest to the regular grid lat/lon
        # nearest lat/lon point was calculated using the minimum great circle
        # distance between the model lat/lon and the global lat/lon
        mass_grid = np.zeros((nlat,nlon))
        for i in range(npts):
            mass_grid[ilat[i],ilon[i]] += ptms[i]
        # new surface density in cm w.e.
        grid.data[:,:,t] = mass_grid/grid.area
        # set grid data mask
        grid.mask[:,:,t] = np.logical_not(point_mask)
        # copy date parameters for time step
        grid.time[t] = fd['time'][t].copy()
        grid.month[t] = gravtk.time.calendar_to_grace(grid.time[t])

    # update mask
    grid.update_mask()
    # output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')
    # output regridded data file
    FILE='gsfc_fdm_{0}_{1}_{2}.nc'.format(FILE_VERSION,REGION.lower(),VARIABLE)
    # output grid to netCDF4 (.nc)
    grid.to_netCDF4(os.path.join(DIRECTORY,FILE), varname=VARIABLE, units='cmwe',
        longname='Equivalent Water Thickness', reference=reference)
    # change the permissions mode
    os.chmod(os.path.join(DIRECTORY,FILE), MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Read and regrid MERRA-2 hybrid variables
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # region of firn model
    parser.add_argument('--region','-R',
        type=str, default='gris', choices=['gris','ais'],
        help='Region of firn model to calculate')
    # version of firn model
    versions = ['v0','v1','v1.0','v1.1','v1.2','v1.2.1']
    parser.add_argument('--version','-v',
        type=str, default='v1.2.1', choices=versions,
        help='Version of firn model to calculate')
    # products from firn model
    parser.add_argument('--product','-P',
        type=str, default=('SMB_a','h_a'), nargs='+',
        help='MERRA-2 hybrid product to calculate')
    # years to run
    now = time.gmtime()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.tm_year+1),
        help='Years of model outputs to run')
    # mask file for reducing to regions
    parser.add_argument('--mask',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        nargs='+', default=[],
        help='netCDF4 masks file for reducing to regions')
    # output grid parameters
    parser.add_argument('--spacing','-S',
        type=float, nargs='+', default=[0.5,0.5], metavar=('dlon','dlat'),
        help='Spatial resolution of output data')
    parser.add_argument('--interval','-I',
        type=int, default=2, choices=[1,2,3],
        help=('Output grid interval '
            '(1: global, 2: centered global, 3: non-global)'))
    parser.add_argument('--bounds','-B',
        type=float, nargs=4, metavar=('lon_min','lon_max','lat_min','lat_max'),
        help='Bounding box for non-global grid')
    # print information about each input and output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
    # netCDF4 files are gzip compressed
    parser.add_argument('--gzip','-G',
        default=False, action='store_true',
        help='netCDF4 file is locally gzip compressed')
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

    # run program
    for VARIABLE in args.product:
        merra_hybrid_regrid(args.directory, args.region, VARIABLE, args.year,
            VERSION=args.version,
            MASKS=args.mask,
            DDEG=args.spacing,
            INTERVAL=args.interval,
            BOUNDS=args.bounds,
            GZIP=args.gzip,
            MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()