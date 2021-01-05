ecco_cube92_ocean_depth.py
==========================

- Interpolates [GEBCO bathymetry](https://www.bodc.ac.uk/data/hosted_data_systems/gebco_gridded_bathymetry_data/
) to ECCO2 Cube92 ocean model grids

#### Calling Sequence
```bash
python ecco_cube92_ocean_depth.py --directory <path_to_directory> model_file
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/ECCO/ecco_cube92_ocean_depth.py)

#### Inputs
- `model_file`: ECCO2 Cube92 Model File

#### Command Line Options
- `-D X`, `--directory X`: working data directory
- `-v X`, `--version X`: GEBCO bathymetry version
- `-M X`, `--mode X`: Permission mode of directories and files
