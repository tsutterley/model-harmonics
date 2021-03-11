ecmwf_reanalysis_retrieve.py
============================

- Retrieves ERA-Interim reanalysis netCDF4 datasets from the ECMWF Web API
    * 2-metre Temperature (t2m)
    * Surface Pressure (ps)
    * Mean Sea Level Pressure (msl)
    * Temperature (t) and Specific Humidity (q) on Model Levels
    * Invariant Parameters

#### Calling Sequence
```bash
python ecmwf_reanalysis_retrieve.py --directory <path_to_directory>
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/ecmwf_reanalysis_retrieve.py)

#### Command Line Options
- `-U X`, `--api-url`: ECMWF api url
- `-K X`, `--api-key`: ECMWF api key
- `-E X,` `--api-email`: ECMWF api email
- `-D X`, `--directory X`: Working data directory
- `-Y X`, `--year X`: Year to retrieve
- `-I`, `--invariant`: Retrieve the model invariant parameters
- `-M X`, `--mode X`: Permissions mode of the directories and files
