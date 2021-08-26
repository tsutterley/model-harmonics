merra_hybrid_cumulative.py
==========================

- Reads MERRA-2 hybrid datafiles to calculate cumulative anomalies in derived surface mass balance products

#### Calling Sequence
```bash
python merra_hybrid_cumulative.py --directory <path_to_directory> --region gris --mean 1980 1995
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/SMB/merra_hybrid_cumulative.py)

#### Command Line Options
- `-D X`, `--directory X`: Working data directory
- `-R X`, `--region X`: Region to calculate
    * `gris`
    * `ais`
- `-v X`, `--version X`: Version of firn model to calculate
    * `v0`
    * `v1`
    * `v1.0`
    * `v1.1`
- `--mean`: Start and end year of mean
- `-G`, `--gzip`: netCDF4 file is locally gzip compressed
- `-M X`, `--mode X`: Local permissions mode of the directories and files

