ecco_depth_version4.py
======================

- Interpolates [GEBCO bathymetry](https://www.bodc.ac.uk/data/hosted_data_systems/gebco_gridded_bathymetry_data/
) to ECCO Version 4 interpolated ocean model grids

#### Calling Sequence
```bash
python ecco_depth_version4.py --directory <path_to_directory> model_file
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/ECCO/ecco_depth_version4.py)

#### Inputs
- `model_file`: ECCO Version 4 Model File

#### Command Line Options
- `-D X`, `--directory X`: working data directory
- `-v X`, `--version X`: GEBCO bathymetry version
- `-M X`, `--mode X`: Permission mode of directories and files
