ecco_mean_realtime.py
=====================

- Reads 12-hour [ECCO ocean bottom pressure data from JPL](https://ecco.jpl.nasa.gov/drive/files/NearRealTime/Readme)
- Global area average of each ocean bottom pressure map is removed [(Greatbatch, 1994)](https://doi.org/10.1029/94JC00847)
- Calculates multi-annual means on an equirectangular grid

#### Calling Sequence
```bash
python ecco_mean_realtime.py --directory <path_to_directory> kf080i dr080i
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/ECCO/ecco_mean_realtime.py)

#### Inputs
- ECCO Near Real-Time models
    * `'kf080i'`: Kalman filter analysis
    * `'dr080i'`: RTS smoother analysis

#### Command Line Options
- `-D X`, `--directory X`: working data directory
- `-m X`, `--mean X`: Year range for mean
- `-F X`, `--format X`: input and output data format
    * `'ascii'`
    * `'netcdf'`
    * `'HDF5'`
- `-M X`, `--mode X`: Permission mode of directories and files
- `-V`, `--verbose`: Output information for each output file