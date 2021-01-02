reanalysis_geopotential_heights.py
==================================

- Reads temperature and specific humidity data to calculate geopotential height and pressure difference fields at half levels from reanalysis

#### Calling Sequence
```bash
python reanalysis_geopotential_heights.py --directory <path_to_directory> ERA5 MERRA-2
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/reanalysis_geopotential_heights.py)

#### Inputs
- [ERA-Interim](http://apps.ecmwf.int/datasets/data/interim-full-moda)
- [ERA5](http://apps.ecmwf.int/data-catalogues/era5/?class=ea)
- [MERRA-2](https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/)

#### Command Line Options
- `-D X`, `--directory X`: working data directory
- `-Y X`, `--year X`: years to run
- `-V`, `--verbose`:  Output information for each output file
- `-M X`, `--mode X`: Permissions mode of the files created
