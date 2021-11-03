racmo_smb_cumulative.py
=======================

- Reads RACMO datafiles to calculate cumulative anomalies in derived surface mass balance products

#### Calling Sequence
```bash
python racmo_smb_cumulative.py --product smb --mean 1980 1995 model_file
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/SMB/racmo_smb_cumulative.py)

Inputs
######

- `model_file`: full path to input RACMO netCDF4 file

#### Command Line Options
- `-P X`, `--product X`: RACMO SMB product to calculate
- `--mean`: Start and end year of mean
- `-G`, `--gzip`: netCDF4 file is locally gzip compressed
- `-V`, `--verbose`: Output information for each output file
- `-M X`, `--mode X`: Local permissions mode of the directories and files
