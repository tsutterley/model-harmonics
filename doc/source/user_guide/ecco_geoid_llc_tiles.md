ecco_geoid_llc_tiles.py
=======================

- Calculates geoid heights for ECCO ocean model LLC tiles using model coefficients from the [GFZ International Centre for Global Earth Models (ICGEM)](http://icgem.gfz-potsdam.de/home)

#### Calling Sequence
```bash
python ecco_geoid_llc_tiles.py --geoid <path_to_geoid_file> ECCO-Grid.nc ECCO-EGM2008.nc
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/ECCO/ecco_geoid_llc_tiles.py)

#### Inputs
- `input_file`: ECCO LLC tile grid file
- `output_file`: output geoid height file

#### Command Line Options
- `-G X`, `--geoid X`: gfc file from the GFZ ICGEM
- `-L X`, `--lmax X`: maximum spherical harmonic degree
- `-M X`, `--mode X`: Permission mode of directories and files
- `-V`, `--verbose`: Output information for each output file
