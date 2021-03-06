ucar_rda_jra55_surface.py
=========================

- Downloads JRA-55 products using a links list provided by the [NCAR/UCAR Research Data Archive (RDA)](https://rda.ucar.edu/)
- JRA-55 6-hour data is more regularly updated compared with the monthly means
- Combines 6-hour model outputs into monthly averages
- Will extract data files if compressed into a single tar file

#### Calling Sequence
```bash
python ucar_rda_jra55_surface.py --directory <path_to_directory> links_list_file
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/ucar_rda_jra55_surface.py)

#### Inputs
- `links_list_file`: UCAR links list file (`csh`)

#### Command Line Options
- `-U X`, `--user X`: Username for UCAR/NCAR RDA login
- `-P X`, `--password X`: Password for UCAR/NCAR RDA login
- `-N X`, `--netrc X`: Path to .netrc file for authentication
- `-D X`, `--directory X:` Full path to output directory
- `-Y X`, `--year X`: Years to download from input links file
- `-I`, `--interpolated`: Input data is interpolated to 1.25 degrees
- `-G`, `--gzip`: Input data is compressed
- `-l`, `--log`: Output log of files downloaded
- `-M X`, `--mode=X`: Permission mode of directories and files downloaded
