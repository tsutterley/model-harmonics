#!/usr/bin/env python
u"""
spatial.py
Written by Tyler Sutterley (06/2024)
Functions for reading, writing and processing spatial data
Extends gravity_toolkit spatial module adding raster support

PYTHON DEPENDENCIES:
    spatial.py: spatial data class for reading, writing and processing data

UPDATE HISTORY:
    Updated 06/2024: added function for calculating latitude and longitude
    Updated 04/2024: changed polar stereographic area function to scale_factors
    Updated 11/2023: add class for creating spatial mosaics
    Updated 08/2023: add function for flipping raster object
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: convert spacing and extent to raster class properties
        improve typing for variables in docstrings
        add function for calculating geocentric latitude from geodetic
    Updated 02/2023: geotiff read and write to inheritance of spatial class
    Written 10/2022
"""
import copy
import uuid
import logging
import pathlib
import warnings
import numpy as np
import gravity_toolkit.spatial
from model_harmonics.datum import datum

# attempt imports
try:
    import osgeo.gdal, osgeo.osr, osgeo.gdalconst
except (AttributeError, ImportError, ModuleNotFoundError) as exc:
    warnings.filterwarnings("module")
    warnings.warn("GDAL not available", ImportWarning)
try:
    import pyproj
except (AttributeError, ImportError, ModuleNotFoundError) as exc:
    warnings.filterwarnings("module")
    warnings.warn("pyproj not available", ImportWarning)
# ignore warnings
warnings.filterwarnings("ignore")

# PURPOSE: additional routines for the spatial module
# for reading and writing raster data
class raster(gravity_toolkit.spatial):
    """
    Inheritance of ``spatial`` class for reading and writing
    raster data

    Attributes
    ----------
    data: np.ndarray
        spatial grid data
    mask: np.ndarray
        spatial grid mask
    x: np.ndarray
        x-coordinate array
    y: np.ndarray
        y-coordinate array
    lon: np.ndarray
        grid longitudes
    lat: np.ndarray
        grid latitudes
    fill_value: float or NoneType, default None
        invalid value for spatial grid data
    attributes: dict
        attributes of spatial variables
    extent: list, default [None,None,None,None]
        spatial grid bounds
        ``[minimum x, maximum x, minimum y, maximum y]``
    spacing: list, default [None,None]
        grid step size ``[x, y]``
    shape: tuple
        dimensions of spatial object
    ndim: int
        number of dimensions of spatial object
    filename: str
        input or output filename

    """
    np.seterr(invalid='ignore')
    # inherit spatial class to read more data types
    def __init__(self, projection=None, **kwargs):
        super().__init__(**kwargs)
        self.x = None
        self.y = None
        self.projection = projection

    def from_geotiff(self, filename, **kwargs):
        """
        Read data from a geotiff file

        Parameters
        ----------
        filename: str
            full path of input geotiff file
        compression: str or NoneType, default None
            file compression type
        bounds: list or NoneType, default bounds
            extent of the file to read:
            ``[minimum x, maximum x, minimum y, maximum y]``
        """
        # set filename
        self.case_insensitive_filename(filename)
        # set default keyword arguments
        kwargs.setdefault('compression',None)
        kwargs.setdefault('bounds',None)
        # Open the geotiff file for reading
        logging.info(self.filename)
        if (kwargs['compression'] == 'gzip'):
            # read as GDAL gzip virtual geotiff dataset
            mmap_name = f"/vsigzip/{str(self.filename)}"
            ds = osgeo.gdal.Open(mmap_name)
        elif (kwargs['compression'] == 'bytes'):
            # read as GDAL memory-mapped (diskless) geotiff dataset
            mmap_name = f"/vsimem/{uuid.uuid4().hex}"
            osgeo.gdal.FileFromMemBuffer(mmap_name, self.filename.read())
            ds = osgeo.gdal.Open(mmap_name)
        else:
            # read geotiff dataset
            ds = osgeo.gdal.Open(str(self.filename),
                osgeo.gdalconst.GA_ReadOnly)
        # get the spatial projection reference information
        srs = ds.GetSpatialRef()
        self.attributes['projection'] = srs.ExportToProj4()
        self.attributes['wkt'] = srs.ExportToWkt()
        # get dimensions
        xsize = ds.RasterXSize
        ysize = ds.RasterYSize
        bsize = ds.RasterCount
        # get geotiff info
        info_geotiff = ds.GetGeoTransform()
        # calculate image extents
        xmin = info_geotiff[0]
        ymax = info_geotiff[3]
        xmax = xmin + (xsize-1)*info_geotiff[1]
        ymin = ymax + (ysize-1)*info_geotiff[5]
        # x and y pixel center coordinates (converted from upper left)
        x = xmin + info_geotiff[1]/2.0 + np.arange(xsize)*info_geotiff[1]
        y = ymax + info_geotiff[5]/2.0 + np.arange(ysize)*info_geotiff[5]
        # if reducing to specified bounds
        if kwargs['bounds'] is not None:
            # reduced x and y limits
            xlimits = (kwargs['bounds'][0],kwargs['bounds'][1])
            ylimits = (kwargs['bounds'][2],kwargs['bounds'][3])
            # Specify offset and rows and columns to read
            xoffset = int((xlimits[0] - xmin)/info_geotiff[1])
            yoffset = int((ymax - ylimits[1])/np.abs(info_geotiff[5]))
            xcount = int((xlimits[1] - xlimits[0])/info_geotiff[1]) + 1
            ycount = int((ylimits[1] - ylimits[0])/np.abs(info_geotiff[5])) + 1
            # reduced x and y pixel center coordinates
            self.x = x[slice(xoffset, xoffset + xcount, None)]
            self.y = y[slice(yoffset, yoffset + ycount, None)]
            # read reduced image with GDAL
            self.data = ds.ReadAsArray(xoff=xoffset, yoff=yoffset,
                xsize=xcount, ysize=ycount)
            # reduced image extent (converted back to upper left)
            xmin = np.min(self.x) - info_geotiff[1]/2.0
            xmax = np.max(self.x) - info_geotiff[1]/2.0
            ymin = np.min(self.y) - info_geotiff[5]/2.0
            ymax = np.max(self.y) - info_geotiff[5]/2.0
        else:
            # x and y pixel center coordinates
            self.x = np.copy(x)
            self.y = np.copy(y)
            # read full image with GDAL
            self.data = ds.ReadAsArray()
        # image extent
        self.attributes['extent'] = (xmin, xmax, ymin, ymax)
        # check if image has fill values
        self.mask = np.zeros_like(self.data, dtype=bool)
        # get invalid value (0 is falsy)
        self.fill_value = ds.GetRasterBand(1).GetNoDataValue()
        if self.fill_value or (self.fill_value == 0):
            # mask invalid values
            self.mask[:] = (self.data == self.fill_value)
        # close the dataset
        ds = None
        self.update_mask()
        return self

    def to_geotiff(self, filename, **kwargs):
        """
        Write a spatial object to a geotiff file

        Parameters
        ----------
        filename: str
            full path of output geotiff file
        driver: str, default 'cog'
            GDAL driver

                - ``'GTiff'``: GeoTIFF
                - ``'cog'``: Cloud Optimized GeoTIFF
        dtype: obj, default osgeo.gdal.GDT_Float64
            GDAL data type
        options: list, default ['COMPRESS=LZW']
            GDAL driver creation options
        """
        # set filename
        self.filename = pathlib.Path(filename).expanduser().absolute()
        # set default keyword arguments
        kwargs.setdefault('driver', 'GTiff')
        kwargs.setdefault('dtype', osgeo.gdal.GDT_Float64)
        kwargs.setdefault('options', ['COMPRESS=LZW'])
        # verify grid dimensions to be iterable
        self.expand_dims()
        # grid shape
        ny,nx,nband = np.shape(self.data)
        # output as geotiff or specified driver
        driver = osgeo.gdal.GetDriverByName(kwargs['driver'])
        # set up the dataset with creation options
        ds = driver.Create(str(self.filename), nx, ny, nband,
            kwargs['dtype'], kwargs['options'])
        # top left x, w-e pixel resolution, rotation
        # top left y, rotation, n-s pixel resolution
        xmin, xmax, ymin, ymax = self.attributes['extent']
        dx, dy = self.spacing
        ds.SetGeoTransform([xmin, dx, 0, ymax, 0, dy])
        # set the spatial projection reference information
        srs = osgeo.osr.SpatialReference()
        srs.ImportFromWkt(self.attributes['wkt'])
        # export
        ds.SetProjection( srs.ExportToWkt() )
        # for each band
        for band in range(nband):
            # set fill value for band (0 is falsy)
            if self.fill_value or (self.fill_value == 0):
                ds.GetRasterBand(band+1).SetNoDataValue(self.fill_value)
            # write band to geotiff array
            ds.GetRasterBand(band+1).WriteArray(self.data[:,:,band])
        # print filename if verbose
        logging.info(self.filename)
        # close dataset
        ds.FlushCache()

    def get_latlon(self, srs_proj4=None, srs_wkt=None, srs_epsg=None):
        """
        Get the latitude and longitude of grid cells

        Parameters
        ----------
        srs_proj4: str or NoneType, default None
            PROJ4 projection string
        srs_wkt: str or NoneType, default None
            Well-Known Text (WKT) projection string
        srs_epsg: int or NoneType, default None
            EPSG projection code

        Returns
        -------
        longitude: np.ndarray
            longitude coordinates of grid cells
        latitude: np.ndarray
            latitude coordinates of grid cells
        """
        # set the spatial projection reference information
        if srs_proj4 is not None:
            source = pyproj.CRS.from_proj4(srs_proj4)
        elif srs_wkt is not None:
            source = pyproj.CRS.from_wkt(srs_wkt)
        elif srs_epsg is not None:
            source = pyproj.CRS.from_epsg(srs_epsg)
        else:
            source = pyproj.CRS.from_string(self.projection)
        # target spatial reference (WGS84 latitude and longitude)
        target = pyproj.CRS.from_epsg(4326)
        # create transformation
        transformer = pyproj.Transformer.from_crs(source, target,
            always_xy=True)
        # create meshgrid of points in original projection
        x, y = np.meshgrid(self.x, self.y)
        # convert coordinates to latitude and longitude
        self.lon, self.lat = transformer.transform(x, y)
        return self

    @property
    def spacing(self):
        """Step size of ``raster`` object ``[x, y]``
        """
        return (self.x[1] - self.x[0], self.y[1] - self.y[0])

    @property
    def extent(self):
        """Bounds of ``raster`` object
        ``[minimum x, maximum x, minimum y, maximum y]``
        """
        xmin = np.min(self.x)
        xmax = np.max(self.x)
        ymin = np.min(self.y)
        ymax = np.max(self.y)
        return [xmin, xmax, ymin, ymax]

    def copy(self):
        """
        Copy a ``raster`` object to a new ``raster`` object
        """
        temp = raster(fill_value=self.fill_value)
        # copy attributes or update attributes dictionary
        if isinstance(self.attributes, list):
            setattr(temp,'attributes',self.attributes)
        elif isinstance(self.attributes, dict):
            temp.attributes.update(self.attributes)
        # assign variables to self
        var = ['x','y','data','mask','error','time','month','filename']
        for key in var:
            try:
                val = getattr(self, key)
                setattr(temp, key, copy.copy(val))
            except AttributeError:
                pass
        # update mask
        temp.replace_masked()
        return temp

    def expand_dims(self):
        """
        Add a singleton dimension to a spatial object if non-existent
        """
        # output spatial with a third dimension
        if (np.ndim(self.data) == 2):
            self.data = self.data[:,:,None]
            # try expanding mask variable
            try:
                self.mask = self.mask[:,:,None]
            except Exception as exc:
                pass
        # get spacing and dimensions
        self.update_mask()
        return self

    def flip(self, axis=0):
        """
        Reverse the order of data and dimensions along an axis

        Parameters
        ----------
        axis: int, default 0
            axis to reorder
        """
        # output spatial object
        temp = self.copy()
        # copy dimensions and reverse order
        if (axis == 0):
            temp.y = temp.y[::-1].copy()
        elif (axis == 1):
            temp.x = temp.x[::-1].copy()
        # attempt to reverse possible data variables
        for key in ['data','mask','error']:
            try:
                setattr(temp, key, np.flip(getattr(self, key), axis=axis))
            except Exception as exc:
                pass
        # update mask
        temp.update_mask()
        return temp

    def __str__(self):
        """String representation of the ``raster`` object
        """
        properties = ['model_harmonics.raster']
        extent = ', '.join(map(str, self.extent))
        properties.append(f"    extent: {extent}")
        spacing = ', '.join(map(str, self.spacing))
        properties.append(f"    spacing: {spacing}")
        shape = ', '.join(map(str, self.shape))
        properties.append(f"    shape: {shape}")
        if self.month:
            properties.append(f"    start_month: {min(self.month)}")
            properties.append(f"    end_month: {max(self.month)}")
        return '\n'.join(properties)

class mosaic:
    """Utility for creating spatial mosaics
    """
    def __init__(self, **kwargs):
        self.extent = [np.inf,-np.inf,np.inf,-np.inf]
        self.spacing = [None,None]
        self.fill_value = np.nan

    def update_spacing(self, x, y):
        """
        update the step size of mosaic
        """
        try:
            self.spacing = (x[1] - x[0], y[1] - y[0])
        except:
            pass
        return self

    def update_bounds(self, x, y):
        """
        update the bounds of mosaic
        """
        # check that there is data
        if not np.any(x) or not np.any(y):
            return self
        # get extent of new data
        extent = [x.min(), x.max(), y.min(), y.max()]
        if (extent[0] < self.extent[0]):
            self.extent[0] = np.copy(extent[0])
        if (extent[1] > self.extent[1]):
            self.extent[1] = np.copy(extent[1])
        if (extent[2] < self.extent[2]):
            self.extent[2] = np.copy(extent[2])
        if (extent[3] > self.extent[3]):
            self.extent[3] = np.copy(extent[3])
        return self

    def image_coordinates(self, x, y):
        """
        get the image coordinates
        """
        # check that there is data
        if not np.any(x) or not np.any(y):
            return (None, None)
        # get the image coordinates
        iy = np.array((y[:,None] - self.extent[2])/self.spacing[1], dtype=np.int64)
        ix = np.array((x[None,:] - self.extent[0])/self.spacing[0], dtype=np.int64)
        return (iy, ix)

    @property
    def dimensions(self):
        """Dimensions of the mosaic"""
        dims = [None, None]
        # calculate y dimensions with new extents
        dims[0] = np.int64((self.extent[3] - self.extent[2])/self.spacing[1]) + 1
        # calculate x dimensions with new extents
        dims[1] = np.int64((self.extent[1] - self.extent[0])/self.spacing[0]) + 1
        return dims

    @property
    def shape(self):
        """Shape of the mosaic"""
        return (self.dimensions[0], self.dimensions[1], )

    @property
    def x(self):
        """X-coordinates of the mosaic"""
        return self.extent[0] + self.spacing[0]*np.arange(self.dimensions[1])

    @property
    def y(self):
        """Y-coordinates of the mosaic"""
        return self.extent[2] + self.spacing[1]*np.arange(self.dimensions[0])

# get WGS84 parameters in CGS (centimeters, grams, seconds)
_wgs84 = datum(ellipsoid='WGS84', units='CGS')

# PURPOSE: calculate the geocentric latitudes
def geocentric_latitude(
        lon: np.ndarray,
        lat: np.ndarray,
        a_axis: float = _wgs84.a_axis,
        flat: float = _wgs84.flat,
    ):
    """
    Converts from geodetic latitude to geocentric latitude for an ellipsoid
    :cite:p:`Snyder:1982gf`

    Parameters
    ----------
    lon: np.ndarray,
        longitude (degrees east)
    lat: np.ndarray,
        geodetic latitude (degrees north)
    a_axis: float, default 6378137.0
        semimajor axis of the ellipsoid
    flat: float, default 1.0/298.257223563
        ellipsoidal flattening

    Returns
    -------
    geocentric_latitude: np.ndarray
        latitude intersecting the center of the Earth (degrees north)
    """
    # first numerical eccentricity
    ecc1 = np.sqrt((2.0*flat - flat**2)*a_axis**2)/a_axis
    # geodetic latitude in radians
    latitude_geodetic_rad = np.pi*lat/180.0
    # prime vertical radius of curvature
    N = a_axis/np.sqrt(1.0 - ecc1**2.*np.sin(latitude_geodetic_rad)**2.)
    # calculate X, Y and Z from geodetic latitude and longitude
    X = N * np.cos(latitude_geodetic_rad) * np.cos(np.pi*lon/180.0)
    Y = N * np.cos(latitude_geodetic_rad) * np.sin(np.pi*lon/180.0)
    Z = (N * (1.0 - ecc1**2.0)) * np.sin(latitude_geodetic_rad)
    # calculate geocentric latitude and convert to degrees
    return 180.0*np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))/np.pi

def scale_factors(
        lat: np.ndarray,
        flat: float = _wgs84.flat,
        reference_latitude: float = 70.0,
        metric: str = 'area'
    ):
    """
    Calculates scaling factors to account for polar stereographic
    distortion including special case of at the exact pole
    :cite:p:`Snyder:1982gf`

    Parameters
    ----------
    lat: np.ndarray
        latitude (degrees north)
    flat: float, default 1.0/298.257223563
        ellipsoidal flattening
    reference_latitude: float, default 70.0
        reference latitude (true scale latitude)
    metric: str, default 'area'
        metric to calculate scaling factors

            - ``'distance'``: scale factors for distance
            - ``'area'``: scale factors for area

    Returns
    -------
    scale: np.ndarray
        scaling factors at input latitudes
    """
    assert metric.lower() in ['distance', 'area'], 'Unknown metric'
    # convert latitude from degrees to positive radians
    theta = np.abs(lat)*np.pi/180.0
    # convert reference latitude from degrees to positive radians
    theta_ref = np.abs(reference_latitude)*np.pi/180.0
    # square of the eccentricity of the ellipsoid
    # ecc2 = (1-b**2/a**2) = 2.0*flat - flat^2
    ecc2 = 2.0*flat - flat**2
    # eccentricity of the ellipsoid
    ecc = np.sqrt(ecc2)
    # calculate ratio at input latitudes
    m = np.cos(theta)/np.sqrt(1.0 - ecc2*np.sin(theta)**2)
    t = np.tan(np.pi/4.0 - theta/2.0)/((1.0 - ecc*np.sin(theta)) / \
        (1.0 + ecc*np.sin(theta)))**(ecc/2.0)
    # calculate ratio at reference latitude
    mref = np.cos(theta_ref)/np.sqrt(1.0 - ecc2*np.sin(theta_ref)**2)
    tref = np.tan(np.pi/4.0 - theta_ref/2.0)/((1.0 - ecc*np.sin(theta_ref)) / \
        (1.0 + ecc*np.sin(theta_ref)))**(ecc/2.0)
    # distance scaling
    k = (mref/m)*(t/tref)
    kp = 0.5*mref*np.sqrt(((1.0+ecc)**(1.0+ecc))*((1.0-ecc)**(1.0-ecc)))/tref
    if (metric.lower() == 'distance'):
        # distance scaling
        scale = np.where(np.isclose(theta, np.pi/2.0), 1.0/kp, 1.0/k)
    elif (metric.lower() == 'area'):
        # area scaling
        scale = np.where(np.isclose(theta, np.pi/2.0), 1.0/(kp**2), 1.0/(k**2))
    return scale

def get_latlon(x, y, srs_proj4=None, srs_wkt=None, srs_epsg=None):
    """
    Get the latitude and longitude of grid cells

    Parameters
    ----------
    srs_proj4: str or NoneType, default None
        PROJ4 projection string
    srs_wkt: str or NoneType, default None
        Well-Known Text (WKT) projection string
    srs_epsg: int or NoneType, default None
        EPSG projection code

    Returns
    -------
    longitude: np.ndarray
        longitude coordinates of grid cells
    latitude: np.ndarray
        latitude coordinates of grid cells
    """
    # set the spatial projection reference information
    if srs_proj4 is not None:
        source = pyproj.CRS.from_proj4(srs_proj4)
    elif srs_wkt is not None:
        source = pyproj.CRS.from_wkt(srs_wkt)
    elif srs_epsg is not None:
        source = pyproj.CRS.from_epsg(srs_epsg)
    else:
        raise ValueError('No projection information provided')
    # target spatial reference (WGS84 latitude and longitude)
    target = pyproj.CRS.from_epsg(4326)
    # create transformation
    transformer = pyproj.Transformer.from_crs(source, target,
        always_xy=True)
    # create meshgrid of points in original projection
    gridx, gridy = np.meshgrid(x, y)
    # convert coordinates to latitude and longitude
    longitude, latitude = transformer.transform(gridx, gridy)
    return (longitude, latitude)
