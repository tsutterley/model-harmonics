gldas_read_monthly.py
=====================

- Reads GLDAS monthly datafiles to calculate anomalies in total water storage from soil moisture, snow water equivalent and total canopy storage

#### Calling Sequence
```bash
python gldas_read_monthly.py --directory <path_to_directory> --version 2.1 CLSM NOAH VIC
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/gldas_read_monthly.py)

#### Inputs
- `'CLM'`: GLDAS Common Land Model
- `'CLSM'`: GLDAS Catchment Land Surface Model
- `'MOS'`: GLDAS Mosaic model
- `'NOAH'`: GLDAS Noah model
- `'VIC'`: GLDAS Variable Infiltration Capacity model

#### Command Line Options
- `-D X`, `--directory X`: Working data directory
- `-Y X,` `--year X`: Years to run
- `-M X`, `--mean X`: Year range for mean
- `-S X`, `--spacing X`: Spatial resolution of models to run
    * `'10'`: 1.0 degrees latitude/longitude
    * `'025'`: 0.25 degrees latitude/longitude
- `-v X`, `--version X`: GLDAS model version to run
- `-F X`, `--format=X`: input and output data format
    * `'ascii'`
    * `'netcdf'`
    * `'HDF5'`
- `-M X`, `--mode X`: Permission mode of directories and files
- `-V`, `--verbose`: Output information for each output file
