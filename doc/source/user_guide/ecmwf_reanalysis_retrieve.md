ecmwf_reanalysis_retrieve.py
============================

- Retrieves ERA-Interim reanalysis netCDF4 datasets from the ECMWF Web API

#### Calling Sequence
```bash
python ecmwf_reanalysis_retrieve.py --directory <path_to_directory>
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/ecmwf_reanalysis_retrieve.py)

#### Command Line Options

- `-D X`, `--directory X`: Working data directory
- `-Y X`, `--year X`: Year to retrieve
- `-I`, `--invariant`: Retrieve the model invariant parameters
- `-M X`, `--mode X`: Permissions mode of the directories and files
