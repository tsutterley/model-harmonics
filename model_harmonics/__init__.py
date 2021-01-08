"""
A spherical harmonics toolkit for Python
========================================

model_harmonics contains Python tools for obtaining and working with
model synthetic spherical harmonic coefficients for comparing with
data from the the NASA/DLR Gravity Recovery and Climate Experiment
(GRACE) and the NASA/GFZ Gravity Recovery and Climate Experiment
Follow-On (GRACE-FO) missions

Documentation is available at https://model-harmonics.readthedocs.io
"""
import model_harmonics.utilities
from model_harmonics.gen_atmosphere_stokes import gen_atmosphere_stokes
from model_harmonics.gen_pressure_stokes import gen_pressure_stokes
