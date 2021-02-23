harmonic_operators.py
====================

- Performs basic operations on spherical harmonic files

#### Calling Sequence
```bash
python harmonic_operators.py --operation add infile1 infile2 outfile
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/scripts/harmonic_operators.py)

#### Inputs
- path to input spherical harmonic files
- path to output spherical harmonic file

#### Command Line Options
- `--help`: list the command line options
- `-O X`, `--operation X`: Operation to run
    * `'add'`
    * `'subtract'`
    * `'multiply'`
    * `'divide'`
    * `'mean'`
    * `'destripe'`
- `-l X`, -`-lmax X`: maximum spherical harmonic degree
- `-m X`, `--mmax X`: maximum spherical harmonic order
- `-F X`, `--format X`: input and output data format
    * `'ascii'`
    * `'netCDF4'`
    * `'HDF5'`
- `-D`, `--date`: input and output files have date information
- `-V`, `--verbose`: Output information for each output file
- `-M X`, `--mode X`: Permission mode of directories and files
