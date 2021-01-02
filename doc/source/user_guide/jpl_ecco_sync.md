jpl_ecco_sync.py
================

- Syncs ECCO Ocean Bottom Pressure outputs from the [NASA JPL ECCO Drive server](https://ecco.jpl.nasa.gov/drive/files/NearRealTime/Readme)

#### Calling Sequence
```bash
python jpl_ecco_sync.py --directory <path_to_directory> kf080i dr080i
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/ECCO/jpl_ecco_sync.py)

#### Inputs
- ECCO near-real time models
    * `'kf080i'`: Kalman filter analysis
    * `'dr080i'`: RTS smoother analysis

#### Command Line Options
- `-U X`, `--user X`: username for NASA Earthdata Login
- `-N X`, `--netrc X`: path to .netrc file for authentication
- `-D X`, `--directory X`: working data directory
- `-Y X`, `--year X`: Years to sync
- `-L`, `--list`: print files to be transferred, but do not execute transfer
- `-l`, `--log`: output log of files downloaded
- `-C`, `--clobber:` Overwrite existing data in transfer
- `--checksum`: compare hashes to check if overwriting existing data
- `-M X`, `--mode X`: Local permissions mode of the directories and files synced
