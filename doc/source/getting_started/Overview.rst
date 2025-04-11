========
Overview
========

This documentation is intended to explain how to compute spherical harmonics from model
outputs for comparing with or correcting time-variable gravity measurements from the
`GRACE/GRACE-FO <https://github.com/tsutterley/gravity-toolkit>`_ missions.
This software was developed with the goal of supporting science applications for
time-variable gravity.
The ``model-harmonics`` projects consists of extension routines for the set of ``gravity-toolkit`` tools.

OBP
===

Uses outputs from the NASA-JPL `Estimating the Circulation and Climate of the Ocean (ECCO) <https://ecco-group.org/>`_ model.
For ECCO near real-time Kalman-filtered (kf080i) and Rauch-Tung-Striebel (RTS) smoother (dr080i) models, reads 12-hour ocean bottom pressure data (OBP) and calculates monthly averages.
For ECCO version 4 models, reads monthly ocean bottom pressure potential anomalies and converts to estimates of absolute ocean bottom pressure (OBP).
Near real-time models are downloaded using the ``jpl_ecco_sync.py`` program,
monthly Cube92 models are calculated using the ``jpl_ecco_cube92_sync.py`` program,
interpolated Version 4 models are downloaded using the ``jpl_ecco_v4_sync.py`` program, and
monthly Version 4 and 5 models in LLC tile format are downloaded using the ``jpl_ecco_llc_sync.py`` program.
Because Boussinesq-type models conserve volume rather than mass, the global area average of each monthly map is removed :cite:p:`Greatbatch:1994dx`.
Monthly anomalies in ocean bottom pressure are calculated by removing a multi-annual mean (typically 2003 |ndash| 2007).
Ocean bottom pressure anomalies are converted to spherical harmonics following :cite:p:`Boy:2005el` (Equations :eq:`1` and :eq:`2`).

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
        E [label="ECCO Ocean Bottom Pressure"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        N [label="Geoid Height"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        B [label="Ocean Bathymetry"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        M [URL="https://github.com/tsutterley/model-harmonics/blob/main/OBP/ecco_mean_realtime.py"
            label="Calculate Temporal Mean"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        R [URL="https://github.com/tsutterley/model-harmonics/blob/main/OBP/ecco_read_realtime.py"
            label="Calculate Monthly Anomalies"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        H [URL="https://github.com/tsutterley/model-harmonics/blob/main/OBP/ecco_monthly_harmonics.py"
            label="Calculate Spherical Harmonics"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        S [URL="https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/combine_harmonics.py"
            label="Spatial Maps"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"]
        T [URL="https://github.com/tsutterley/model-harmonics/blob/main/scripts/least_squares_mascon_timeseries.py"
            label="Time Series"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"]
        E -> M [arrowsize=0.8]
        M -> R [arrowsize=0.8]
        E -> R [arrowsize=0.8]
        R -> H [arrowsize=0.8]
        N -> H [arrowsize=0.8]
        B -> H [arrowsize=0.8]
        H -> S [arrowsize=0.8]
        H -> T [arrowsize=0.8]
    }

TWS
===

Uses `GLDAS model outputs <https://ldas.gsfc.nasa.gov/gldas>`_ from the NASA Goddard Space Flight Center (GSFC) Hydrological Sciences Laboratory (HSL)
`Global Land Data Assimilation System Version 2 (GLDAS-2) <https://disc.gsfc.nasa.gov/information/data-release?title=New%20and%20Reprocessed%20GLDAS%20Version%202%20Data%20Products%20Released>`_
:cite:p:`Rodell:2004ke`.
GLDAS outputs are downloaded using the ``gesdisc_gldas_sync.py`` program.
GLDAS version 2.1 is forced with a combination of model and observation data.
Additionally, the GLDAS project produces two months of "early production stream" products that are run without the forcing data.
Here, monthly terrestrial water storage (TWS) estimates are calculated by combining the GLDAS soil moisture (`SM`), snow water equivalent (`SWE`) and total canopy storage outputs.
Monthly anomalies in terrestrial water storage are calculated by removing a multi-annual mean (typically 2003 |ndash| 2007).
Before converting to spherical harmonics, the GLDAS terrestrial water storage estimates are masked to remove
`urbanized <https://github.com/tsutterley/model-harmonics/blob/main/TWS/gldas_mask_vegetation.py>`_,
`glaciated <https://github.com/tsutterley/model-harmonics/blob/main/TWS/gldas_mask_arctic.py>`_ and
`permafrost <https://github.com/tsutterley/model-harmonics/blob/main/TWS/gldas_mask_permafrost.py>`_ regions.
Terrestrial water storage anomalies are converted to spherical harmonics following :cite:p:`Wahr:1998hy` (Equation :eq:`3`).

.. math::
    :label: 3

	\left\{\begin{matrix}\tilde{C}_{lm}(t) \\[-4pt] \tilde{S}_{lm}(t) \end{matrix} \right\} =
	\frac{3}{4\pi a\rho_{e}}\frac{1+k_l}{2l+1}\int\sigma(\theta,\phi,t)~P_{lm}(\cos\theta)
	\left\{\begin{matrix}\cos{m\phi} \\[-4pt] \sin{m\phi} \end{matrix} \right\}~d\Omega

.. graphviz::
    :caption: GLDAS Spherical Harmonics Framework
    :align: center

    digraph {
        E [label="GLDAS Land Surface\nModel Outputs"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        L [label="Vegetation and\nLand Surface Masks"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        M [URL="https://github.com/tsutterley/model-harmonics/blob/main/TWS/gldas_mean_monthly.py"
            label="Calculate Temporal Mean"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        R [URL="https://github.com/tsutterley/model-harmonics/blob/main/TWS/gldas_read_monthly.py"
            label="Calculate Monthly Anomalies"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        H [URL="https://github.com/tsutterley/model-harmonics/blob/main/TWS/gldas_monthly_harmonics.py"
            label="Calculate Spherical Harmonics"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        S [URL="https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/combine_harmonics.py"
            label="Spatial Maps"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"]
        T [URL="https://github.com/tsutterley/model-harmonics/blob/main/scripts/least_squares_mascon_timeseries.py"
            label="Time Series"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"]
        E -> M [arrowsize=0.8]
        E -> R [arrowsize=0.8]
        M -> R [arrowsize=0.8]
        R -> H [arrowsize=0.8]
        L -> H [arrowsize=0.8]
        H -> S [arrowsize=0.8]
        H -> T [arrowsize=0.8]
    }


Reanalysis
==========

`ERA-Interim <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era-interim>`_ is computed by ECMWF and is available starting from 1979.
`ERA5 <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>`_  is the latest reanalysis computed by ECMWF offering much higher spatial and temporal resolution and is available starting from 1950.
Differences between ERA-Interim and ERA5 are outlined `here <https://confluence.ecmwf.int/pages/viewpage.action?pageId=74764925>`_.
ERA-Interim outputs are downloaded using the ``ecmwf_reanalysis_retrieve.py`` program following using the `ecmwf-api-client <https://confluence.ecmwf.int/display/WEBAPI/Access+ECMWF+Public+Datasets>`_ documentation.
ERA5 outputs are downloaded using the ``cds_reanalysis_retrieve.py`` program following using the `cdsapi <https://cds.climate.copernicus.eu/api-how-to>`_ documentation.
`MERRA-2 <https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/>`_ is computed by the NASA Global Modeling and Assimilation Office (GMAO) and is available starting from 1980.
MERRA-2 outputs are downloaded using the ``gesdisc_merra_download.py`` or ``gesdisc_merra_monthly.py`` programs.
`NCEP-DOE-2 <https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis2.html>`_ is computed by the National Centers for Environmental Prediction (NCEP) and is available starting from 1979.
NCEP-DOE-2 outputs are downloaded using the ``noaa_cdc_ncep_ftp.py`` program.
`NCEP-CFSR <https://cfs.ncep.noaa.gov/>`_ is computed by the National Centers for Environmental Prediction (NCEP) and is available starting from 1979 with Version 2 available from 2011 onward.
NCEP-CFSR outputs are downloaded using the ``ucar_rda_cfsr_surface.py`` program.
`JRA-55 <http://jra.kishou.go.jp/JRA-55/index_en.html>`_ is computed by the Japan Meteorological Agency (JMA) and is available starting from 1958.
JRA-55 outputs are downloaded using the ``ucar_rda_jra55_surface.py`` program.

Spherical harmonics from reanalysis outputs are computed here using three different schemes of complexity following :cite:p:`Boy:2005el` and :cite:p:`Swenson:2002kf`:
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
        E [label="Reanalysis Surface Pressure"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        N [label="Geoid Height"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        O [label="Model Orography"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        M [URL="https://github.com/tsutterley/model-harmonics/blob/main/TWS/reanalysis_mean_pressure.py"
            label="Calculate Temporal Mean"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        H [URL="https://github.com/tsutterley/model-harmonics/blob/main/TWS/reanalysis_pressure_harmonics.py"
            label="Calculate Spherical Harmonics"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        S [URL="https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/combine_harmonics.py"
            label="Spatial Maps"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"]
        T [URL="https://github.com/tsutterley/model-harmonics/blob/main/scripts/least_squares_mascon_timeseries.py"
            label="Time Series"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"]
        E -> M [arrowsize=0.8]
        M -> H [arrowsize=0.8]
        E -> H [arrowsize=0.8]
        N -> H [arrowsize=0.8]
        O -> H [arrowsize=0.8]
        H -> S [arrowsize=0.8]
        H -> T [arrowsize=0.8]
    }

.. graphviz::
    :caption: Reanalysis Spherical Harmonics with Three-Dimensional Geometry Framework
    :align: center

    digraph {
        E [label="Reanalysis Temperature\nand Specific Humidity"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        L [URL="https://github.com/tsutterley/model-harmonics/blob/main/TWS/model_level_coefficients.py"
            label="Model Level\nCoefficients"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        N [label="Geoid Height"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        O [label="Model Orography"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        G [URL="https://github.com/tsutterley/model-harmonics/blob/main/TWS/reanalysis_geopotential_heights.py"
            label="Calculate Geopotential Heights\nand Pressure Differences"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        M [URL="https://github.com/tsutterley/model-harmonics/blob/main/TWS/reanalysis_mean_harmonics.py"
            label="Calculate Temporal Mean\nSpherical Harmonics"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        H [URL="https://github.com/tsutterley/model-harmonics/blob/main/TWS/reanalysis_atmospheric_harmonics.py"
            label="Calculate Spherical Harmonics"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        S [URL="https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/combine_harmonics.py"
            label="Spatial Maps"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"]
        T [URL="https://github.com/tsutterley/model-harmonics/blob/main/scripts/least_squares_mascon_timeseries.py"
            label="Time Series"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"]
        E -> G [arrowsize=0.8]
        L -> G [arrowsize=0.8]
        O -> G [arrowsize=0.8]
        G -> M [arrowsize=0.8]
        M -> H [arrowsize=0.8]
        G -> H [arrowsize=0.8]
        N -> H [arrowsize=0.8]
        H -> S [arrowsize=0.8]
        H -> T [arrowsize=0.8]
    }

SMB
===

Uses `MERRA-2 model outputs <https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/s>`_ from the NASA `Global Modeling and Assimilation Office (GMAO) <https://gmao.gsfc.nasa.gov/>`_,
or `ERA5 model outputs <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>`_  computed by ECMWF.
MERRA-2 `Vertically Integrated Diagnostics (M2TMNXINT) <https://disc.gsfc.nasa.gov/datasets/M2TMNXINT_5.12.4/summary>`_ and
`Land Ice Surface Diagnostics (M2TMNXGLC) <https://disc.gsfc.nasa.gov/datasets/M2TMNXGLC_5.12.4/summary>`_ are downloaded using the ``gesdisc_merra_sync.py`` program.
ERA5 precipitation and evaporation outputs are downloaded using the ``cds_reanalysis_retrieve.py`` program following using the `cdsapi <https://cds.climate.copernicus.eu/api-how-to>`_ documentation.
For MERRA-2, monthly surface mass balance (SMB) estimates are calculated by combining the
convective rain (`PRECCU`), large-scale rain (`PRECLS`), snow (`PRECSN`), evaporation (`EVAP`), and runoff over glaciated land (`RUNOFF`) variables.
For ERA5,  monthly surface mass balance (SMB) estimates are calculated by combining the total precipitation (`tp`) and evaporation (`e`) variables.
ERA5 surface mass balance estimates are not including runoff as those variables are presently `inaccurate over glaciated surfaces <https://confluence.ecmwf.int/pages/viewpage.action?pageId=208488132>`_.
Monthly cumulative anomalies in surface mass balance are calculated by removing a multi-annual mean (typically 1980 |ndash| 1995).
Before converting to spherical harmonics, the surface mass balance estimates are masked to isolate regions of interest.
Surface mass balance anomalies are converted to spherical harmonics following :cite:p:`Wahr:1998hy` (Equation :eq:`7`).

.. math::
    :label: 7

	\left\{\begin{matrix}\tilde{C}_{lm}(t) \\[-4pt] \tilde{S}_{lm}(t) \end{matrix} \right\} =
	\frac{3}{4\pi a\rho_{e}}\frac{1+k_l}{2l+1}\int\sigma(\theta,\phi,t)~P_{lm}(\cos\theta)
	\left\{\begin{matrix}\cos{m\phi} \\[-4pt] \sin{m\phi} \end{matrix} \right\}~d\Omega

.. graphviz::
    :caption: Surface Mass Balance Spherical Harmonics Framework
    :align: center

    digraph {
        E [label="MERRA-2 Reanalysis\nModel Outputs"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        L [label="Region Masks"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        M [URL="https://github.com/tsutterley/model-harmonics/blob/main/SMB/merra_smb_mean.py"
            label="Calculate Temporal Mean"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        R [URL="https://github.com/tsutterley/model-harmonics/blob/main/SMB/merra_smb_cumulative.py"
            label="Calculate Cumulative Anomalies"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        H [URL="https://github.com/tsutterley/model-harmonics/blob/main/SMB/merra_smb_harmonics.py"
            label="Calculate Spherical Harmonics"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        S [URL="https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/combine_harmonics.py"
            label="Spatial Maps"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"]
        T [URL="https://github.com/tsutterley/model-harmonics/blob/main/scripts/least_squares_mascon_timeseries.py"
            label="Time Series"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"]
        E -> M [arrowsize=0.8]
        E -> R [arrowsize=0.8]
        M -> R [arrowsize=0.8]
        R -> H [arrowsize=0.8]
        L -> H [arrowsize=0.8]
        H -> S [arrowsize=0.8]
        H -> T [arrowsize=0.8]
    }

.. |ndash|    unicode:: U+2013 .. EN DASH
