ecco_read_version4.py
=====================

- Reads monthly ECCO ocean bottom pressure data from [Version 4 models](https://ecco-group.org/products-ECCO-V4r4.htm) and calculates monthly anomalies
- Global area average of each ocean bottom pressure map is removed [(Greatbatch, 1994)](https://doi.org/10.1029/94JC00847)

#### Calling Sequence
```bash
python ecco_read_version4.py --directory <path_to_directory> V4r3 V4r4
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/ECCO/ecco_read_version4.py)

#### Inputs
- ECCO near-real time models
    * `'V4r3'`: Version 4, Revision 3
    * `'V4r4'`: Version 4, Revision 4

#### Command Line Options
- `-D X`, `--directory X`: working data directory
- `-Y X`, `--year X`: Years to run
- `-m X`, `--mean X`: Year range for mean
- `-F X`, `--format X`: input and output data format
    * `'ascii'`
    * `'netcdf'`
    * `'HDF5'`
- `-M X`, `--mode X`: Permission mode of directories and files
- `-V`, `--verbose`: Output information for each output file
