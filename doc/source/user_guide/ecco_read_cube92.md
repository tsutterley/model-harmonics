ecco_read_cube92.py
===================

- Reads monthly ECCO2 ocean bottom pressure data from [Cube92 models](https://ecco-group.org/products-ECCO-V4r4.htm) and calculates monthly anomalies
- Global area average of each ocean bottom pressure map is removed [(Greatbatch, 1994)](https://doi.org/10.1029/94JC00847)

#### Calling Sequence
```bash
python ecco_read_cube92.py --directory <path_to_directory>
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/ECCO/ecco_read_cube92.py)

#### Command Line Options
- `-D X`, `--directory X`: working data directory
- `-Y X`, `--year X`: Years to run
- `-m X`, `--mean X`: Year range for mean
- `-F X`, `--format=X`: input and output data format
    * `'ascii'`
    * `'netcdf'`
    * `'HDF5'`
- `-M X`, `--mode X`: Permission mode of directories and files
- `-V`, `--verbose`: Output information for each output file
