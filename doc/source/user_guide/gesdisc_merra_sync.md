gesdisc_merra_sync.py
=====================

- Syncs MERRA-2 surface mass balance (SMB) related products from the Goddard Earth Sciences Data and Information Server Center (GES DISC)
- `tavgM_2d_int` (Vertically Integrated Diagnostics) collection:
    * `PRECCU` (convective rain)
    * `PRECLS` (large-scale rain)
    * `PRECSN` (snow)
    * `EVAP` (evaporation)
- `tavgM_2d_glc` (Land Ice Surface Diagnostics) collection:
    * `RUNOFF` (runoff over glaciated land)

#### Calling Sequence
```bash
python gesdisc_merra_sync.py --directory <path_to_directory>
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/SMB/gesdisc_merra_sync.py)

#### Command Line Options
- `-U X`, `--user X`: username for NASA Earthdata Login
- `-N X`, `--netrc X`: path to .netrc file for authentication
- `-D X`, `--directory X`: working data directory
- `-Y X`, `--year X`: years to sync
- `--log`: output log of files downloaded
- `--list`: print files to be transferred, but do not execute transfer
- `--clobber`: Overwrite existing data in transfer
- `-M X`, `--mode X`: permissions mode of the directories and files synced
