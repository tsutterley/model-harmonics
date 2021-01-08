merra_smb_mean.py
=================

- Reads MERRA-2 datafiles to calculate multi-annual means of derived surface mass balance products

#### Calling Sequence
```bash
python merra_smb_mean.py --directory <path_to_directory> SMB PRECIP RUNOFF
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/SMB/merra_smb_mean.py)

#### Inputs
- `'SMB'`: Surface Mass Balance
- `'PRECIP'`: Precipitation
- `'RUNOFF'`: Meltwater Runoff

#### Command Line Options
- `-D X`, `--directory X`: Working data directory
- `-M X`, `--mean X`: Year range for mean
- `-F X`, `--format X`: input and output data format
    * `'ascii'`
    * `'netcdf'`
    * `'HDF5'`
- `-M X`, `--mode X`: Permission mode of directories and files
- `-V`, `--verbose`: Output information for each output file
