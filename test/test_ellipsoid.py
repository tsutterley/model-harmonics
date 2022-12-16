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
    refell = mdlhmc.constants(ellipsoid=ellipsoid, units=units)
    # assert values
    assert ellipsoid_params['a'] == refell.a_axis
    assert ellipsoid_params['b'] == refell.b_axis
    assert ellipsoid_params['f'] == refell.flat
    assert ellipsoid_params['rad_e'] == refell.rad_e
    assert ellipsoid_params['rad_p'] == refell.rad_p
    assert ellipsoid_params['ratio'] == refell.ratio
    assert ellipsoid_params['GM'] == refell.GM
    assert ellipsoid_params['omega'] == refell.omega
    assert ellipsoid_params['C20'] == refell.C20
    assert ellipsoid_params['J2'] == refell.J2
    assert ellipsoid_params['U0'] == refell.U0
    assert ellipsoid_params['dk'] == refell.dk
    assert ellipsoid_params['norm_a'] == refell.ga
    assert ellipsoid_params['norm_b'] == refell.gb
    assert ellipsoid_params['mp'] == refell.m
    assert ellipsoid_params['q'] == refell.q
    assert ellipsoid_params['q0'] == refell.q0
    assert ellipsoid_params['ecc'] == refell.ecc
    assert ellipsoid_params['ecc1'] == refell.ecc1
    assert ellipsoid_params['ecc2'] == refell.ecc2
    assert ellipsoid_params['area'] == refell.area
    assert ellipsoid_params['volume'] == refell.volume
    assert ellipsoid_params['rho_e'] == refell.rho_e
