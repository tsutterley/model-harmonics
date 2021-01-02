gldas_monthly_harmonics.py
==========================

- Reads monthly GLDAS total water storage anomalies and converts to spherical harmonic coefficients

#### Calling Sequence
```bash
python gldas_monthly_harmonics.py --directory <path_to_directory> CLSM NOAH VIC
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/gldas_monthly_harmonics.py)

#### Inputs
- `'CLM'`: GLDAS Common Land Model
- `'CLSM'`: GLDAS Catchment Land Surface Model
- `'MOS'`: GLDAS Mosaic model
- `'NOAH'`: GLDAS Noah model
- `'VIC'`: GLDAS Variable Infiltration Capacity model

#### Command Line Options
- `-D X`, `--directory X`: working data directory
- `-Y X`, `--year X`: Years to run
- `-S X`, `--spacing X`: Spatial resolution of models to run
    * `'10'`: 1.0 degrees latitude/longitude
    * `'025'`: 0.25 degrees latitude/longitude
- `-v X`, `--version X`: GLDAS model version to run
- `-l X`, `--lmax X`: maximum spherical harmonic degree
- `-m X`, `--mmax X`: maximum spherical harmonic order
- `-n X`, `--love X`: Load Love numbers dataset
    * `0`: Han and Wahr (1995) values from PREM
    * `1`: Gegout (2005) values from PREM
    * `2`: Wang et al. (2012) values from PREM
- `-r X`, `--reference X`: Reference frame for load love numbers
    * `'CF'`: Center of Surface Figure (default)
    * `'CM'`: Center of Mass of Earth System
    * `'CE'`: Center of Mass of Solid Earth
- `-F X`, `--format X`: input and output data format
    * `'ascii'`
    * `'netCDF4'`
    * `'HDF5'`
- `-V`, `--verbose`: verbose output of processing run
- `-M X`, `--mode X`: Permissions mode of the files created
