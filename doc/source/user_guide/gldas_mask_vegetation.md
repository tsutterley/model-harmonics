gldas_mask_vegetation.py
========================

- Creates a mask for GLDAS data using the [GLDAS vegetation type binary files](https://ldas.gsfc.nasa.gov/gldas/vegetation-class-mask)
    1. Evergreen Needleleaf Forest
    2. Evergreen Broadleaf Forest
    3. Deciduous Needleleaf Forest
    4. Deciduous Broadleaf Forest
    5. Mixed Forest
    6. Closed Shrublands
    7. Open Shrublands
    8. Woody Savannas
    9. Savannas
    10. Grassland
    11. Permanent Wetland
    12. Cropland
    13. Urban and Built-Up
    14. Cropland/Natural Vegetation Mosaic
    15. Snow and Ice
    16. Barren or Sparsely Vegetated
    17. Ocean
    18. Wooded Tundra
    19. Mixed Tundra
    20. Bare Ground Tundra

#### Calling Sequence
```bash
python gldas_mask_vegetation.py --directory <path_to_directory>
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/gldas_mask_vegetation.py)

#### Command Line Options
- `-D X`, `--directory X`: Working data directory
- `-S X`, `--spacing X`: Spatial resolution of models to run
    * `'10'`: 1.0 degrees latitude/longitude
    * `'025'`: 0.25 degrees latitude/longitude
- `-V`, `--verbose`: verbose output of processing run
- `-M X`, `--mode X`: Permissions mode of the files created
