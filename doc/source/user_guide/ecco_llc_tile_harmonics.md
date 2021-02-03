ecco_llc_tile_harmonics.py
==========================

- Reads monthly ECCO ocean bottom pressure anomalies from LLC tiles and converts to spherical harmonic coefficients

#### Calling Sequence
```bash
python ecco_llc_tile_harmonics.py --directory <path_to_directory> V5alpha
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/ECCO/ecco_llc_tile_harmonics.py)

#### Inputs
- ECCO Version 4 or 5 models
    * `'V4r4'`: Version 4, Revision 4
    * `'V5alpha'`: ECCO Version 5, Alpha release

#### Command Line Options
- `-D X`, `--directory X`: working data directory
- `-Y X`, `--year X`: Years to run
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
- `-F X`, `--format X`: output data format
    * `'ascii'`
    * `'netCDF4'`
    * `'HDF5'`
- `-V`, `--verbose`: verbose output of processing run
- `-M X`, `--mode X`: Permissions mode of the files created
