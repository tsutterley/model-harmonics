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