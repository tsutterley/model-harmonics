model-harmonics
===============

Python tools for obtaining and working with model synthetic spherical
harmonic coefficients for comparing with data from the the NASA/DLR
Gravity Recovery and Climate Experiment (GRACE) and the NASA/GFZ
Gravity Recovery and Climate Experiment Follow-On (GRACE-FO) missions

.. toctree::
    :maxdepth: 2
    :caption: Getting Started

    getting_started/Install.rst
    getting_started/Overview.rst
    getting_started/NASA-Earthdata.rst
    getting_started/Contributing.rst
    getting_started/Citations.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: API Reference

    api_reference/datum.rst
    api_reference/gen_atmosphere_stokes.rst
    api_reference/gen_point_pressure.rst
    api_reference/gen_pressure_stokes.rst
    api_reference/greens_kernel.rst
    api_reference/spatial.rst
    api_reference/utilities.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Utilities

    api_reference/harmonic_operators.rst
    api_reference/least_squares_mascons.rst
    api_reference/least_squares_mascon_timeseries.rst
    api_reference/regress_maps.rst
    api_reference/spatial_operators.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: GIA

    api_reference/GIA/GIA_GSFC_mascons.rst
    api_reference/GIA/GIA_spatial_maps.rst
    api_reference/GIA/GIA_uplift_ICESat2_ATL15.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: OBP

    api_reference/OBP/jpl_ecco_webdav.rst
    api_reference/OBP/jpl_ecco_sync.rst
    api_reference/OBP/jpl_ecco_cube92_sync.rst
    api_reference/OBP/jpl_ecco_llc_sync.rst
    api_reference/OBP/jpl_ecco_v4_sync.rst
    api_reference/OBP/ecco_mean_realtime.rst
    api_reference/OBP/ecco_mean_llc_tiles.rst
    api_reference/OBP/ecco_mean_version4.rst
    api_reference/OBP/ecco_read_realtime.rst
    api_reference/OBP/ecco_read_llc_tiles.rst
    api_reference/OBP/ecco_read_version4.rst
    api_reference/OBP/ecco_depth_version4.rst
    api_reference/OBP/ecco_geoid_llc_tiles.rst
    api_reference/OBP/ecco_llc_tile_harmonics.rst
    api_reference/OBP/ecco_monthly_harmonics.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Reanalysis

    api_reference/reanalysis/cds_reanalysis_retrieve.rst
    api_reference/reanalysis/ecmwf_reanalysis_retrieve.rst
    api_reference/reanalysis/gesdisc_merra_download.rst
    api_reference/reanalysis/gesdisc_merra_monthly.rst
    api_reference/reanalysis/gesdisc_merra_subset.rst
    api_reference/reanalysis/noaa_cdc_ncep_ftp.rst
    api_reference/reanalysis/ucar_rda_cfsr_surface.rst
    api_reference/reanalysis/ucar_rda_jra55_surface.rst
    api_reference/reanalysis/model_level_coefficients.rst
    api_reference/reanalysis/reanalysis_atmospheric_harmonics.rst
    api_reference/reanalysis/reanalysis_geopotential_heights.rst
    api_reference/reanalysis/reanalysis_inverse_barometer.rst
    api_reference/reanalysis/reanalysis_mean_harmonics.rst
    api_reference/reanalysis/reanalysis_mean_pressure.rst
    api_reference/reanalysis/reanalysis_monthly_harmonics.rst
    api_reference/reanalysis/reanalysis_monthly_pressure.rst
    api_reference/reanalysis/reanalysis_pressure_harmonics.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Seismic

    api_reference/seismic/seismic_monthly_harmonics.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: SMB

    api_reference/SMB/era5_smb_mean.rst
    api_reference/SMB/era5_smb_cumulative.rst
    api_reference/SMB/era5_smb_harmonics.rst
    api_reference/SMB/gemb_smb_cumulative.rst
    api_reference/SMB/gemb_smb_harmonics.rst
    api_reference/SMB/gesdisc_merra_sync.rst
    api_reference/SMB/merra_hybrid_cumulative.rst
    api_reference/SMB/merra_hybrid_harmonics.rst
    api_reference/SMB/merra_smb_mean.rst
    api_reference/SMB/merra_smb_cumulative.rst
    api_reference/SMB/merra_smb_mask.rst
    api_reference/SMB/merra_smb_harmonics.rst
    api_reference/SMB/racmo_downscaled_cumulative.rst
    api_reference/SMB/racmo_downscaled_mean.rst
    api_reference/SMB/racmo_smb_cumulative.rst
    api_reference/SMB/racmo_smb_harmonics.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: TWS

    api_reference/TWS/cds_land_retrieve.rst
    api_reference/TWS/era5_land_mean_monthly.rst
    api_reference/TWS/era5_land_read_monthly.rst
    api_reference/TWS/era5_land_mask_arctic.rst
    api_reference/TWS/era5_land_mask_permafrost.rst
    api_reference/TWS/era5_land_monthly_harmonics.rst
    api_reference/TWS/gesdisc_gldas_sync.rst
    api_reference/TWS/gldas_mean_monthly.rst
    api_reference/TWS/gldas_read_monthly.rst
    api_reference/TWS/gldas_mask_arctic.rst
    api_reference/TWS/gldas_mask_permafrost.rst
    api_reference/TWS/gldas_mask_vegetation.rst
    api_reference/TWS/gldas_monthly_harmonics.rst
    api_reference/TWS/gldas_scaling_factors.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Bibliography

    getting_started/Bibliography.rst
