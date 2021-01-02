noaa_cdc_ncep_ftp.py
====================

- Syncs NOAA-DOE-2 surface reanalysis outputs with the [NOAA CDC ftp server](ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis2.dailyavgs/surface/)

#### Calling Sequence
```bash
python noaa_cdc_ncep_ftp.py --directory <path_to_directory>
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/noaa_cdc_ncep_ftp.py)

#### Command Line Options
- `-D X`, `--directory X`: Working data directory
- `-Y X`, `--year X`: Years to download
- `--mask`: Download land-sea mask file (`land.nc`)
- `-I`, `--invariant`: Download invariant parameters file (`hgt.sfc.nc`)
- `-l`, `--log`: output log of files downloaded
- `-M X`, `--mode X`: Permissions mode of the directories and files downloaded
