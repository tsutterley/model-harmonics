merra_hybrid_regrid.py
======================

- Read and regrid MERRA-2 hybrid derived surface mass balance products

#### Calling Sequence
```bash
python merra_hybrid_regrid.py --directory <path_to_directory> --region gris
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/SMB/merra_hybrid_regrid.py)

#### Command Line Options
- `-D X`, `--directory X`: Working data directory
- `-R X`, `--region X`: Region to calculate
    * `gris`
    * `ais`
- `-v X`, `--version X`: Version of firn model to calculate
    * `v0`
    * `v1`
    * `v1.0`
    * `v1.1`
- `-P X`, `--product X`: MERRA-2 hybrid product to calculate
- `-Y X`, `--year X`: Years to run
- `--mask X`: netCDF4 mask files for reducing to regions
- `-S X`, `--spacing X`: spatial resolution of input data (dlon,dlat)
- `-I X`, `--interval X`: output grid interval
    * `1`: (0:360, 90:-90)
    * `2`: (degree spacing/2)
    * `3`: non-global grid (set with defined bounds)
- `-B X`, `--bounds X`: non-global grid bounding box (minlon,maxlon,minlat,maxlat)
- `-G`, `--gzip`: input netCDF4 file is gzip compressed
- `-V`, `--verbose`: Output information for each output file
- `-M X`, `--mode X`: Permissions mode of the files created