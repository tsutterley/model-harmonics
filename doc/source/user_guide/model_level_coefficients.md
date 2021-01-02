model_level_coefficients.py
===========================

- Creates a netCDF4 file of reanalysis A and B coefficients for model levels
- Model level coefficients are obtained using equation 3.17 of [Simmons and Burridge (1981)](https://doi.org/10.1175/1520-0493(1981)109<0758:AEAAMC>2.0.CO;2) and the methodology of [Trenberth et al (1993)](https://doi.org/10.5065/D6HX19NH)

#### Calling Sequence
```bash
python model_level_coefficients.py --directory <path_to_directory> ERA5 MERRA-2
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/model_level_coefficients.py)

#### Inputs
- [ERA-Interim](http://apps.ecmwf.int/datasets/data/interim-full-moda)
- [ERA5](http://apps.ecmwf.int/data-catalogues/era5/?class=ea)
- [MERRA-2](https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/)

#### Command Line Options
- `-D X`, `--directory X`: Working data directory
- `-M X`, `--mode X`: Permission mode of directories and files
