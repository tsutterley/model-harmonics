ecco_monthly_harmonics.py
=========================

- Reads monthly ECCO ocean bottom pressure anomalies and converts to spherical harmonic coefficients

#### Calling Sequence
```bash
python ecco_monthly_harmonics.py --directory <path_to_directory> kf080i dr080i
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/ECCO/ecco_monthly_harmonics.py)

#### Inputs
- `'kf080i'`: ECCO Near Real-Time Kalman filter analysis
- `'dr080i'`: ECCO Near Real-Time RTS smoother analysis
- `'Cube92'`: ECCO2 Cube92 models
- `'V4r3'`: ECCO Version 4, Revision 3 models
- `'V4r4'`: ECCO Version 4, Revision 4 models

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
- `-F X`, `--format X`: input and output data format
    * `'ascii'`
    * `'netCDF4'`
    * `'HDF5'`
- `-V`, `--verbose`: verbose output of processing run
- `-M X`, `--mode X`: Permissions mode of the files created
