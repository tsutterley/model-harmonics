gldas_mask_arctic.py
====================

- Creates a mask for GLDAS data for Greenland, Svalbard, Iceland and the Russian High Arctic defined by a set of shapefiles

#### Calling Sequence
```bash
python gldas_mask_arctic.py --directory <path_to_directory> --shapefile <path_to_shapefiles>
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/gldas_mask_arctic.py)

#### Command Line Options
- `-D X`, `--directory X`: Working data directory
- `-S X`, `--spacing X`: Spatial resolution of models to run
    * `'10'`: 1.0 degrees latitude/longitude
    * `'025'`: 0.25 degrees latitude/longitude
- `-F X`, `--shapefile X`: Shapefiles to use
- `-V`, `--verbose`: verbose output of processing run
- `-M X`, `--mode X`: Permissions mode of the files created
