#!/usr/bin/env python
u"""
test_ellispoid.py (12/2022)
Tests the that ellipsoid estimates are equivalent
"""
import pytest
import geoid_toolkit as geoidtk
import model_harmonics as mdlhmc

# PURPOSE: test ellipsoid outputs
@pytest.mark.parametrize("ellipsoid", ['WGS84'])
@pytest.mark.parametrize("units", ['MKS','CGS'])
def test_ellipsoid(ellipsoid, units):
    # get dictionary output from geoid-toolkit
    ellipsoid_params = geoidtk.ref_ellipsoid(ellipsoid, UNITS=units)
    # get class output from model-haromnics
    units = mdlhmc.datum(ellipsoid=ellipsoid, units=units)
    # assert values
    assert ellipsoid_params['a'] == units.a_axis
    assert ellipsoid_params['b'] == units.b_axis
    assert ellipsoid_params['f'] == units.flat
    assert ellipsoid_params['rad_e'] == units.rad_e
    assert ellipsoid_params['rad_p'] == units.rad_p
    assert ellipsoid_params['ratio'] == units.ratio
    assert ellipsoid_params['GM'] == units.GM
    assert ellipsoid_params['omega'] == units.omega
    assert ellipsoid_params['C20'] == units.C20
    assert ellipsoid_params['J2'] == units.J2
    assert ellipsoid_params['U0'] == units.U0
    assert ellipsoid_params['dk'] == units.dk
    assert ellipsoid_params['norm_a'] == units.ga
    assert ellipsoid_params['norm_b'] == units.gb
    assert ellipsoid_params['mp'] == units.m
    assert ellipsoid_params['q'] == units.q
    assert ellipsoid_params['q0'] == units.q0
    assert ellipsoid_params['ecc'] == units.ecc
    assert ellipsoid_params['ecc1'] == units.ecc1
    assert ellipsoid_params['ecc2'] == units.ecc2
    assert ellipsoid_params['area'] == units.area
    assert ellipsoid_params['volume'] == units.volume
    assert ellipsoid_params['rho_e'] == units.rho_e
