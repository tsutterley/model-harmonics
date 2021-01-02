reanalysis_monthly_harmonics.py
===============================

- Reads atmospheric surface pressure fields from reanalysis and calculates sets of spherical harmonics using a thin-layer 2D spherical geometry

#### Calling Sequence
```bash
python reanalysis_monthly_harmonics.py --directory <path_to_directory> ERA5 MERRA-2
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/reanalysis_monthly_harmonics.py)

#### Inputs
- [ERA-Interim](http://apps.ecmwf.int/datasets/data/interim-full-moda)
- [ERA5](http://apps.ecmwf.int/data-catalogues/era5/?class=ea)
- [MERRA-2](https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/)
- [NCEP-DOE-2](https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis2.html)
- [NCEP-CFSR](https://rda.ucar.edu/datasets/ds093.1/)
- [JRA-55](http://jra.kishou.go.jp/JRA-55/index_en.html)

#### Command Line Options
- `-D X`, `--directory X`: working data directory
- `-Y X`, `--year X`: Years of model outputs to run
- `--mean X`: Start and end year for mean
- `--redistribute`: Uniformly redistribute values over the ocean
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
- `-V`, `--verbose`:  Output information for each output file
- `-M X`, `--mode X`: Permissions mode of the files created
