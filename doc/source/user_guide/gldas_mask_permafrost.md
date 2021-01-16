gldas_mask_permafrost.py
========================

- Creates a mask for GLDAS data based on the permafrost/surface classification from the [NSIDC Circum-Arctic Map of Permafrost and Ground-Ice Conditions](http://nsidc.org/data/ggd318.html)
    1. Continuous Permafrost
    2. Discontinuous Permafrost
    3. Isolated Permafrost
    4. Sporadic Permafrost
    5. Glaciated Area

#### Calling Sequence
```bash
python gldas_mask_permafrost.py --directory <path_to_directory> --shapefile <path_to_shapefile>
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/gldas_mask_permafrost.py)

#### Command Line Options
- `-D X`, `--directory X`: Working data directory
- `-S X`, `--spacing X`: Spatial resolution of models to run
    * `'10'`: 1.0 degrees latitude/longitude
    * `'025'`: 0.25 degrees latitude/longitude
- `-F X`, `--shapefile X`: Shapefile to use
- `-V`, `--verbose`: verbose output of processing run
- `-M X`, `--mode X`: Permissions mode of the files created
