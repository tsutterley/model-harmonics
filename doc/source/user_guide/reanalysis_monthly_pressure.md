reanalysis_monthly_pressure.py
==============================

Reads daily atmospheric pressure fields from reanalysis and outputs monthly averages

#### Calling Sequence
```bash
python reanalysis_monthly_pressure.py --directory <path_to_directory> NCEP-DOE-2
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/reanalysis_monthly_pressure.py)

#### Inputs
- [NCEP-DOE-2](https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis2.html)

#### Command Line Options
- `-D X`, `--directory X`: working data directory
- `-Y X`, `--year X`: years to run
- `-V`, `--verbose`:  Output information for each output file
- `-M X`, `--mode X`: Permissions mode of the files created
