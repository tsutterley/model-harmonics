========
Overview
========

This documentation is intended to explain how to compute spherical harmonics from model
outputs for comparing or correcting GRACE/GRACE-FO time-variable gravity measurements.
This software was developed with the goal of supporting science applications for
time-variable gravity.
`model-harmonics <https://github.com/tsutterley/model-harmonics>`__ consists of
extension routines for the set of
`read-GRACE-harmonics <https://github.com/tsutterley/read-GRACE-harmonics>`__ tools.

ECCO
====

Uses outputs from the NASA-JPL `Estimating the Circulation and Climate of the Ocean (ECCO) <https://ecco-group.org/>`_ model.
For ECCO near-real time Kalman-filtered (kf080i) and Rauch-Tung-Striebel (RTS) smoother (dr080i) models, reads 12-hour ocean bottom pressure data (OBP) and calculates monthly averages.
For ECCO version 4 models, reads monthly ocean bottom pressure potential anomalies and converts to estimates of absolute ocean bottom pressure (OBP).
Near-real time models are downloaded using the `jpl_ecco_sync.py <https://github.com/tsutterley/model-harmonics/blob/main/ECCO/jpl_ecco_sync.py>`_ program and
Version 4 models are downloaded using the `jpl_ecco_v4_sync.py <https://github.com/tsutterley/model-harmonics/blob/main/ECCO/jpl_ecco_v4_sync.py>`_ program.
Because Boussinesq-type models conserve volume rather than mass, the global area average of each monthly map is removed `(Greatbatch, 1994) <https://doi.org/10.1029/94JC00847>`_.
Monthly anomalies in ocean bottom pressure are calculated by removing a multi-annual mean (typically 2003 |ndash| 2007).
Ocean bottom pressure anomalies are converted to spherical harmonics following `Boy and Chao. (2005) <https://doi.org/10.1029/2002JB002333>`_ (Equations :eq:`1` and :eq:`2`).

.. math::
    :label: 1

	\left\{\begin{matrix}\tilde{C}_{lm}(t) \\ \tilde{S}_{lm}(t) \end{matrix} \right\} =
	\frac{3}{4\pi a\rho_{e}}\frac{1+k_l}{(2l+1)}\int \xi_l(\theta,\phi,t)~P_{lm}(\cos\theta)
	\left\{\begin{matrix}\cos{m\phi} \\ \sin{m\phi} \end{matrix} \right\}\,d\Omega

.. math::
    :label: 2

	\xi_l(\theta,\phi,t) = \left(\frac{a +  N(\theta,\phi) - d(\theta,\phi)}{a}\right)^{l+2}\frac{p_{bot}(\theta,\phi,t)}{g(\theta,\phi)}

.. graphviz::
    :caption: ECCO Spherical Harmonics Framework
    :align: center

    digraph {
        E [label="ECCO Ocean Bottom Pressure" shape=box style="filled" color="darkorchid"]
        N [label="Geoid Height" shape=box style="filled" color="darkorchid"]
        B [label="Ocean Bathymetry" shape=box style="filled" color="darkorchid"]
        M [URL="https://github.com/tsutterley/model-harmonics/blob/main/ECCO/ecco_mean_realtime.py"
            label="Calculate Temporal Mean" shape=box style="filled" color="gray"]
        R [URL="https://github.com/tsutterley/model-harmonics/blob/main/ECCO/ecco_read_realtime.py"
            label="Calculate Monthly Anomalies" shape=box style="filled" color="gray"]
        H [URL="https://github.com/tsutterley/model-harmonics/blob/main/ECCO/ecco_monthly_harmonics.py"
            label="Calculate Spherical Harmonics" shape=box style="filled" color="gray"]
        S [label="Spatial Maps" shape=box style="filled" color="mediumseagreen"]
        T [label="Time Series" shape=box style="filled" color="mediumseagreen"]
        E -> M
        M -> R
        E -> R
        R -> H
        N -> H
        B -> H
        H -> S
        H -> T
    }

GLDAS
=====

Uses `GLDAS model outputs <https://ldas.gsfc.nasa.gov/gldas>`_ from the NASA Goddard Space Flight Center (GSFC) Hydrological Sciences Laboratory (HSL)
`Global Land Data Assimilation System Version 2 (GLDAS-2) <https://disc.gsfc.nasa.gov/information/data-release?title=New%20and%20Reprocessed%20GLDAS%20Version%202%20Data%20Products%20Released>`_.
GLDAS outputs are downloaded using the `gesdisc_gldas_sync.py <https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/gesdisc_gldas_sync.py>`_ program.
GLDAS version 2.1 is forced with a combination of model and observation data.
Additionally, the GLDAS project produces two months of "early production stream" products that are run without the forcing data.
Here, monthly terrestrial water storage (TWS) estimates are calculated by combining the GLDAS soil moisture (SM), snow water equivalent (SWE) and total canopy storage outputs.
Monthly anomalies in terrestrial water storage are calculated by removing a multi-annual mean (typically 2003 |ndash| 2007).
Before converting to spherical harmonics, the GLDAS terrestrial water storage estimates are masked to remove glaciated and permafrost regions.
Terrestrial water storage anomalies are converted to spherical harmonics following `Wahr et al. (1998) <https://doi.org/10.1029/98JB02844>`_ (Equation :eq:`3`).

.. math::
    :label: 3

	\left\{\begin{matrix}\tilde{C}_{lm}(t) \\[-4pt] \tilde{S}_{lm}(t) \end{matrix} \right\} =
	\frac{3}{4\pi a\rho_{e}}\frac{1+k_l}{2l+1}\int\sigma(\theta,\phi,t)~P_{lm}(\cos\theta)
	\left\{\begin{matrix}\cos{m\phi} \\[-4pt] \sin{m\phi} \end{matrix} \right\}~d\Omega

.. graphviz::
    :caption: GLDAS Spherical Harmonics Framework
    :align: center

    digraph {
        E [label="GLDAS Land Surface\nModel Outputs" shape=box style="filled" color="darkorchid"]
        L [label="Vegetation and\nLand Surface Masks" shape=box style="filled" color="darkorchid"]
        M [URL="https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/gldas_mean_monthly.py"
            label="Calculate Temporal Mean" shape=box style="filled" color="gray"]
        R [URL="https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/gldas_read_monthly.py"
            label="Calculate Monthly Anomalies" shape=box style="filled" color="gray"]
        H [URL="https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/gldas_monthly_harmonics.py"
            label="Calculate Spherical Harmonics" shape=box style="filled" color="gray"]
        S [URL="https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/combine_harmonics.py"
            label="Spatial Maps" shape=box style="filled" color="mediumseagreen"]
        T [URL="https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/least_squares_mascon_timeseries.py"
            label="Time Series" shape=box style="filled" color="mediumseagreen"]
        E -> M
        E -> R
        M -> R
        R -> H
        L -> H
        H -> S
        H -> T
    }


Reanalysis
==========

`ERA-Interim <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era-interim>`_ is computed by ECMWF and is available starting from 1979.
`ERA5 <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>`_  is the latest reanalysis computed by ECMWF offering much higher spatial and temporal resolution and is available starting from 1950.
Differences between ERA-Interim and ERA5 are outlined `here <https://confluence.ecmwf.int/pages/viewpage.action?pageId=74764925>`_.
ERA-Interim outputs are downloaded using the `ecmwf_reanalysis_retrieve.py <https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/ecmwf_reanalysis_retrieve.py>`_ program following using the `ecmwf-api-client <https://confluence.ecmwf.int/display/WEBAPI/Access+ECMWF+Public+Datasets>`_ documentation.
ERA5 outputs are downloaded using the `cds_reanalysis_retrieve.py <https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/cds_reanalysis_retrieve.py>`_ program.
`MERRA-2 <https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/>`_ is computed by the NASA Global Modeling and Assimilation Office (GMAO) and is available starting from 1980.
MERRA-2 outputs are downloaded using the `gesdisc_merra_download.py <https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/gesdisc_merra_download.py>`_ or `gesdisc_merra_monthly.py <https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/gesdisc_merra_monthly.py>`_ programs.
`NCEP-DOE-2 <https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis2.html>`_ is computed by the National Centers for Environmental Prediction (NCEP) and is available starting from 1979.
NCEP-DOE-2 outputs are downloaded using the `noaa_cdc_ncep_ftp.py <https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/noaa_cdc_ncep_ftp.py>`_ program.
`NCEP-CFSR <https://cfs.ncep.noaa.gov/>`_ is computed by the National Centers for Environmental Prediction (NCEP) and is available starting from 1979 with Version 2 available from 2011 onward.
NCEP-CFSR outputs are downloaded using the `ucar_rda_cfsr_surface.py <https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/ucar_rda_cfsr_surface.py>`_ program.
`JRA-55 <http://jra.kishou.go.jp/JRA-55/index_en.html>`_ is computed by the Japan Meteorological Agency (JMA) and is available starting from 1958.
JRA-55 outputs are downloaded using the `ucar_rda_jra55_surface.py <https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/ucar_rda_jra55_surface.py>`_ program.

Spherical harmonics from reanalysis outputs are computed here using three different schemes of complexity following `Boy and Chao. (2005) <https://doi.org/10.1029/2002JB002333>`_:
1) a thin-layer 2D spherical geometry,
2) a thin-layer 2D geometry with realistic geometry incorporating model orography and estimates of geoid height (Equations :eq:`4` and :eq:`5`), and
3) a 3D atmospheric geometry integrating over the model layers (Equations :eq:`4` and :eq:`6`).
Anomalies for each reanalysis are calculated relative to a multi-annual mean (such as 2003 |ndash| 2014).

.. math::
    :label: 4

	\left\{\begin{matrix}\tilde{C}_{lm}(t) \\ \tilde{S}_{lm}(t) \end{matrix} \right\} =
	\frac{3}{4\pi a\rho_{e}}\frac{1+k_l}{(2l+1)}\int \xi_l(\theta,\phi,t)~P_{lm}(\cos\theta)
	\left\{\begin{matrix}\cos{m\phi} \\ \sin{m\phi} \end{matrix} \right\}\,d\Omega

.. math::
    :label: 5

	\xi_l(\theta,\phi,t) = \left(\frac{a + h(\theta,\phi) + N(\theta,\phi)}{a}\right)^{l+2}\frac{p_0(\theta,\phi,t)}{g(\theta,\phi)}

.. math::
    :label: 6

	\xi_l(\theta,\phi,t) = -\int_{p_0}^{0}\left(\frac{a + z(\theta,\phi) + N(\theta,\phi)}{a}\right)^{l+2}\frac{dp}{g(\theta,\phi,z)}

.. graphviz::
    :caption: Reanalysis Spherical Harmonics with Two-Dimensional Geometry Framework
    :align: center

    digraph {
        E [label="Reanalysis Surface Pressure" shape=box style="filled" color="darkorchid"]
        N [label="Geoid Height" shape=box style="filled" color="darkorchid"]
        O [label="Model Orography" shape=box style="filled" color="darkorchid"]
        M [URL="https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/reanalysis_mean_pressure.py"
            label="Calculate Temporal Mean" shape=box style="filled" color="gray"]
        H [URL="https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/reanalysis_pressure_harmonics.py"
            label="Calculate Spherical Harmonics" shape=box style="filled" color="gray"]
        S [URL="https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/combine_harmonics.py"
            label="Spatial Maps" shape=box style="filled" color="mediumseagreen"]
        T [URL="https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/least_squares_mascon_timeseries.py"
            label="Time Series" shape=box style="filled" color="mediumseagreen"]
        E -> M
        M -> H
        E -> H
        N -> H
        O -> H
        H -> S
        H -> T
    }

.. graphviz::
    :caption: Reanalysis Spherical Harmonics with Three-Dimensional Geometry Framework
    :align: center

    digraph {
        E [label="Reanalysis Temperature\nand Specific Humidity" shape=box style="filled" color="darkorchid"]
        L [URL="https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/model_level_coefficients.py"
            label="Model Level\nCoefficients" shape=box style="filled" color="darkorchid"]
        N [label="Geoid Height" shape=box style="filled" color="darkorchid"]
        O [label="Model Orography" shape=box style="filled" color="darkorchid"]
        G [URL="https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/reanalysis_geopotential_heights.py"
            label="Calculate Geopotential Heights\nand Pressure Differences" shape=box style="filled" color="gray"]
        M [URL="https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/reanalysis_mean_harmonics.py"
            label="Calculate Temporal Mean\nSpherical Harmonics" shape=box style="filled" color="gray"]
        H [URL="https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/reanalysis_atmospheric_harmonics.py"
            label="Calculate Spherical Harmonics" shape=box style="filled" color="gray"]
        S [URL="https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/combine_harmonics.py"
            label="Spatial Maps" shape=box style="filled" color="mediumseagreen"]
        T [URL="https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/least_squares_mascon_timeseries.py"
            label="Time Series" shape=box style="filled" color="mediumseagreen"]
        E -> G
        L -> G
        O -> G
        G -> M
        M -> H
        G -> H
        N -> H
        H -> S
        H -> T
    }

.. |ndash|    unicode:: U+2013 .. EN DASH