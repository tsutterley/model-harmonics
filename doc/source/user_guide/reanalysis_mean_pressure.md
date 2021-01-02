reanalysis_mean_pressure.py
===========================

- Calculates the mean surface pressure fields from reanalysis

#### Calling Sequence
```bash
python reanalysis_mean_pressure.py --directory <path_to_directory> ERA5 MERRA-2
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/reanalysis_mean_pressure.py)

#### Inputs
- [ERA-Interim](http://apps.ecmwf.int/datasets/data/interim-full-moda)
- [ERA5](http://apps.ecmwf.int/data-catalogues/era5/?class=ea)
- [MERRA-2](https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/)
- [NCEP-DOE-2](https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis2.html)
- [NCEP-CFSR](https://rda.ucar.edu/datasets/ds093.1/)
- [JRA-55](http://jra.kishou.go.jp/JRA-55/index_en.html)

#### Command Line Options
- `-D X`, `--directory X`: working data directory
- `--mean X`: Start and end year for mean
- `-V`, `--verbose`:  Output information for each output file
- `-M X`, `--mode X`: Permissions mode of the files created

