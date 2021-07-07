ucar_rda_cfsr_surface.py
========================

- Downloads NCEP-CFSR products using a links list provided by the [NCAR/UCAR Research Data Archive (RDA)](https://rda.ucar.edu/)
- For CFSR version 2 files: combines monthly files into yearly to standardize between versions

#### Calling Sequence
```bash
python ucar_rda_cfsr_surface.py --directory <path_to_directory> links_list_file
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/ucar_rda_cfsr_surface.py)

#### Inputs
- `links_list_file`: UCAR links list file (`csh`)

#### Command Line Options
- `-D X`, `--directory X:` Working data directory
- `-U X`, `--user X`: Username for UCAR/NCAR RDA login
- `-P X`, `--password X`: Password for UCAR/NCAR RDA login
- `-N X`, `--netrc X`: Path to .netrc file for authentication
- `-Y X`, `--year X`: Years to download from input links file
- `-I`, `--isentropic`: Input data is over isentropic levels
- `-G`, `--gzip`: Input data is compressed
- `-t X`, `--timeout X`: Timeout in seconds for blocking operations
- `-l`, `--log`: Output log of files downloaded
- `-M X`, `--mode=X`: Permission mode of directories and files downloaded
