jpl_ecco_cube92_sync.py
=======================

- Converts ECCO2 Cube92 daily model outputs from the [NASA JPL ECCO2 server](https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/readme.txt) into monthly averages

#### Calling Sequence
```bash
python jpl_ecco_cube92_sync.py --directory <path_to_directory>
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/ECCO/jpl_ecco_cube92_sync.py)

#### Command Line Options
- `-U X`, `--user X`: username for NASA Earthdata Login
- `-W X`, `--webdav X`: WebDAV password for JPL ECCO Drive Login
- `-N X`, `--netrc X`: path to .netrc file for authentication
- `-D X`, `--directory X`: working data directory
- `-Y X`, `--year X`: Years to sync
- `-P X`, `--product X`: Product to sync
- `-l`, `--log`: output log of files downloaded
- `-V`, `--verbose:`:  Output information for each output file
- `-M X`, `--mode X`: Local permissions mode of the directories and files synced
