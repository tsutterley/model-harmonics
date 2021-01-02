gesdisc_gldas_sync.py
=====================

- Syncs GLDAS monthly datafiles from the Goddard Earth Sciences Data and Information Server Center (GES DISC)

#### Calling Sequence
```bash
python gesdisc_gldas_sync.py --directory <path_to_directory>
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/gesdisc_gldas_sync.py)

#### Inputs
- `'CLM'`: GLDAS Common Land Model
- `'CLSM'`: GLDAS Catchment Land Surface Model
- `'MOS'`: GLDAS Mosaic model
- `'NOAH'`: GLDAS Noah model
- `'VIC'`: GLDAS Variable Infiltration Capacity model

#### Command Line Options
- `-U X`, `--user X`: username for NASA Earthdata Login
- `-N X`, `--netrc X`: path to .netrc file for authentication
- `-D X`, `--directory X`: working data directory
- `-Y X`, `--year X`: years to sync
- `-S X`, `--spacing X`: spatial resolution of models to sync
    * `'10'`: 1.0 degrees latitude/longitude
    * `'025'`: 0.25 degrees latitude/longitude
- `-T X`, `--temporal X`: temporal resolution of models to sync
    * `'M'`: Monthly
    * `'3H'`: 3-hourly
- `-v X`, `--version X`: GLDAS model version to sync
- `-e`, `--early`: Sync GLDAS early products
- `--log`: output log of files downloaded
- `--list`: print files to be transferred, but do not execute transfer
- `--clobber`: Overwrite existing data in transfer
- `-M X`, `--mode X`: permissions mode of the directories and files synced
