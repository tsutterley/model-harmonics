merra_smb_mask.py
=================

- Creates a mask for [MERRA-2 land ice data](https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/M2C0NXASM.5.12.4/1980/MERRA2_101.const_2d_asm_Nx.00000000.nc4) using a set of shapefiles

#### Calling Sequence
```bash
python merra_smb_mask.py --shapefile <path_to_shapefiles> input_file output_file
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/SMB/merra_smb_mask.py)

#### Command Line Options
- `-v X`, `--variable X`: Variable from input netCDF4 file to extract
- `-F X`, `--shapefile X`: Shapefiles to use
- `-A X`, `--area X`: Minimum area threshold for polygons
- `-B X`, `--buffer X`: Distance to buffer polygons
- `-V`, `--verbose`: verbose output of processing run
- `-M X`, `--mode X`: Permissions mode of the files created
