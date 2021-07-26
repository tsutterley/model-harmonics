jpl_ecco_v4_sync.py
===================

- Syncs ECCO Version 4 model outputs from the [NASA JPL ECCO Drive server](https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/interp_monthly/README)
- Errata document for [Version 4, Revision 4 Atmospheric Pressure Forcing](https://ecco-group.org/docs/ECCO_V4r4_errata.pdf)

#### Calling Sequence
```bash
python jpl_ecco_v4_sync.py --directory <path_to_directory> V4r3 V4r4
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/ECCO/jpl_ecco_v4_sync.py)

#### Inputs
- ECCO Version 4 models
    * `'V4r3'`: Version 4, Revision 3
    * `'V4r4'`: Version 4, Revision 4

#### Command Line Options
- `-U X`, `--user X`: username for NASA Earthdata Login
- `-W X`, `--webdav X`: WebDAV password for JPL ECCO Drive Login
- `-N X`, `--netrc X`: path to .netrc file for authentication
- `-D X`, `--directory X`: working data directory
- `-Y X`, `--year X`: Years to sync
- `-P X`, `--product X`: Product to sync
- `-t X`, `--timeout X`: Timeout in seconds for blocking operations
- `-l`, `--log`: output log of files downloaded
- `-L`, `--list`: print files to be transferred, but do not execute transfer
- `-C`, `--clobber:` Overwrite existing data in transfer
- `--checksum`: compare hashes to check if overwriting existing data
- `-M X`, `--mode X`: Local permissions mode of the directories and files synced

