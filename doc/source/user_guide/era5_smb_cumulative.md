era5_smb_cumulative.py
======================

- Reads ERA5 datafiles to calculate monthly cumulative anomalies in derived surface mass balance products

#### Calling Sequence
```bash
python era5_smb_cumulative.py --directory <path_to_directory>
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/SMB/era5_smb_cumulative.py)

#### Command Line Options
- `-D X`, `--directory X`: Working data directory
- `-M X`, `--mean X`: Year range for mean
- `-F X`, `--format X`: input and output data format
    * `'ascii'`
    * `'netcdf'`
    * `'HDF5'`
- `-M X`, `--mode X`: Permission mode of directories and files
- `-V`, `--verbose`: Output information for each output file
