merra_smb_harmonics.py
======================

- Reads monthly MERRA-2 surface mass balance anomalies and converts to spherical harmonic coefficients

#### Calling Sequence
```bash
python merra_smb_harmonics.py --directory <path_to_directory> SMB PRECIP RUNOFF
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/SMB/merra_smb_harmonics.py)

#### Inputs
- `'SMB'`: Surface Mass Balance
- `'PRECIP'`: Precipitation
- `'RUNOFF'`: Meltwater Runoff

#### Command Line Options
- `-D X`, `--directory X`: working data directory
- `-m X`, `--mean X`: Year range for mean
- `-Y X`, `--year X`: Years to run
- `-R X`, `--region X`: region name for subdirectory
- `--mask X`: netCDF4 mask file for reducing to regions
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
