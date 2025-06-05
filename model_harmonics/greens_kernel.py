#!/usr/bin/env python
"""
greens_kernel.py
Written by Tyler Sutterley (11/2024)

Calculate a Green's function kernel for a given set of Love Numbers

CALLING SEQUENCE:
    X, Y, G = greens_kernel(LMAX, WIDTH=[wx,wy],
        SPACING=[dx,dy], LOVE=(hl,kl,ll))

INPUTS:
    LMAX: Maximum spherical harmonic degree

OPTIONS:
    SPACING: Grid spacing in x and y directions (meters)
    WIDTH: Grid width in x and y directions (meters)
    LOVE: List of Love Numbers kl, hl, and ll

OUTPUTS:
    X: X-coordinates of the kernel (meters)
    Y: Y-coordinates of the kernel (meters)
    G: Green's function kernel (m/kg)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/

PROGRAM DEPENDENCIES:
    units.py: class for converting units
    legendre_polynomials.py: computes fully normalized Legendre polynomials

UPDATE HISTORY:
    Updated 11/2024: use bessel function for points within cutoff distance
    Written 11/2024
"""

import numpy as np
import scipy.special
import gravity_toolkit as gravtk

def greens_kernel(LMAX, SPACING=[], WIDTH=[], LOVE=None):
    """
    Calculate the Green's function for a given set of Love Numbers
    following :cite:t:`Farrell:1972cm,Farrell:1973ui,Longman:1962ev`

    Parameters
    ----------
    LMAX: int
        Maximum spherical harmonic degree
    SPACING: list, default []
        Grid spacing in x and y directions
    WIDTH: list, default []
        Grid width in x and y directions
    LOVE: list or None, default None
        List of Love Numbers kl, hl, and ll

    Returns
    -------
    X: numpy.ndarray
        X-coordinates of the kernel
    Y: numpy.ndarray
        Y-coordinates of the kernel
    G: numpy.ndarray
        Green's function kernel
    """
    # get Earth parameters
    # radius of the Earth in meters
    rad_e = gravtk.units().rad_e/100.0
    # average density of the Earth in kg/m^3
    rho_e = gravtk.units().rho_e*1000.0
    # scale factor to convert a mass load to uplift
    scale = 3.0/(4.0*np.pi*rho_e*rad_e**2)
    # verify values are close to expected (a/M_e)
    assert np.isclose(scale, 6.371e6/5.972e24)
    # extract arrays of kl, hl, and ll Love Numbers
    hl, kl, ll = LOVE
    # spherical harmonic degrees
    l = np.arange(LMAX+1)
    # scale used to originally normalize the Legendre polynomials
    norm = np.sqrt(2.0*l + 1)
    # grid spacing
    dx,dy = np.broadcast_to(np.atleast_1d(SPACING),(2,))
    cutoff = np.sqrt(dx*dy)
    # grid width
    W = np.broadcast_to(np.atleast_1d(WIDTH),(2,))
    # centered coordinates
    X = np.arange(0, W[0] + dx, dx) - W[0]/2.0
    Y = np.arange(0, W[1] + dy, dy) - W[1]/2.0
    # create a grid of coordinates
    gridx, gridy = np.meshgrid(X, Y)
    # calculate distance from central point
    D = np.sqrt(gridx**2 + gridy**2)
    # allocate for output Green's function
    nx = np.int64(W[0]//dx) + 1
    ny = np.int64(W[1]//dy) + 1
    G = np.zeros((ny, nx))
    # calculate Green's function for each point
    for i,x in enumerate(X):
        for j,y in enumerate(Y):
            # check if distance is within cutoff
            if (D[j,i] < cutoff):
                # adjustment to potentially avoid singularity
                # equivalent radius of a disc load
                radius = np.sqrt(dx*dy/np.pi)
                # calculate distance at one half grid spacing
                alpha = 0.5*radius/rad_e
                # Bessel function up to LMAX
                # multiply by 2.0 to account for the adjustment
                Pl = 2.0*scipy.special.jv(0, (l + 0.5)*alpha)
            else:
                # angular distance from central point
                alpha = np.cos(D[j,i]/rad_e)
                # Legendre polynomials up to LMAX
                Ptemp,_ = gravtk.legendre_polynomials(LMAX, alpha)
                # unnormalizing Legendre polynomials
                Pl = np.squeeze(Ptemp)/norm
            # calculate Green's function
            G[j,i] = scale*np.sum(hl*Pl)
    # return the Green's function and the coordinates
    return (X, Y, G)
