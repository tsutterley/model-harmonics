gesdisc_merra_monthly.py
========================

- Downloads MERRA-2 products using a links list provided by the Goddard Earth Sciences Data and Information Server Center (GES DISC)
- Combines daily model level outputs into monthly averages

#### Calling Sequence
```bash
python gesdisc_merra_monthly.py --directory <path_to_directory> links_list_file
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/gesdisc_merra_monthly.py)

#### Inputs
- `links_list_file`: GES DISC generated file listing files to download

#### Command Line Options
- `-U X`, `--user X`: username for NASA Earthdata Login
- `-P X`, `--password X`: Password for NASA Earthdata Login
- `-N X`, `--netrc X`: path to .netrc file for authentication
- `-D X`, `--directory X`: working data directory
- `-t X`, `--timeout X`: Timeout in seconds for blocking operations
- `-l`, `--log`: output log of files downloaded
- `-V`, `--verbose`: Output information for each output file
- `-M X`, `--mode X`: Local permissions mode of the files created
