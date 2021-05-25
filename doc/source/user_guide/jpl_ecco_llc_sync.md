jpl_ecco_llc_sync.py
===================

- Syncs ECCO Version 4 and 5 LLC tile model outputs from the [NASA JPL ECCO Drive server](https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/nctiles_monthly)

#### Calling Sequence
```bash
python jpl_ecco_llc_sync.py --directory <path_to_directory> V4r4 V5alpha
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/ECCO/jpl_ecco_llc_sync.py)

#### Inputs
- ECCO Version 4 or 5 models
    * `'V4r4'`: Version 4, Revision 4
    * `'V5alpha'`: Version 5, Alpha release

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
