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
    CUTOFF: Distance from central point to use a disc load (meters)
    LOVE: List of Love Numbers kl, hl, and ll

OUTPUTS:
    X: X-coordinates of the kernel (meters)
    Y: Y-coordinates of the kernel (meters)
    G: Green's function kernel (m/kg)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    units.py: class for converting units
    legendre_polynomials.py: computes fully normalized Legendre polynomials

UPDATE HISTORY:
    Written 11/2024
"""

import numpy as np
import gravity_toolkit as gravtk

def greens_kernel(LMAX, SPACING=[], WIDTH=[], CUTOFF=0.0, LOVE=None):
    """
    Calculate the Green's function for a given set of Love Numbers
    following :cite:p:`Farrell:1972cm`, :cite:p:`Farrell:1973ui`
    and :cite:p:`Longman:1962ev`

    Parameters
    ----------
    LMAX: int
        Maximum spherical harmonic degree
    SPACING: list, default []
        Grid spacing in x and y directions
    WIDTH: list, default []
        Grid width in x and y directions
    CUTOFF: float, default 0.0
        Distance from central point to use a disc load
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
    # scale used to originally normalize the Legendre polynomials
    norm = np.sqrt(2.0*np.arange(LMAX+1) + 1)
    # grid spacing
    dx,dy = np.broadcast_to(np.atleast_1d(SPACING),(2,))
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
            # use disc load if within cutoff following Farrell (1973)
            if (D[j,i] < CUTOFF):
                # for interior points: use a disc load
                # to potentially avoid singularity
                Pl = np.zeros((LMAX+1))
                # divide pixel area by area of a half sphere
                alpha = (1.0 - (dx*dy)/(2.0*np.pi*rad_e**2))
                # l=0 is a special case
                Pl[0] = (1.0 - alpha)/2.0
                # Legendre polynomials up to LMAX+1
                Ptemp,_ = gravtk.legendre_polynomials(LMAX+1, alpha)
                for l in range(1, LMAX+1):
                    # unnormalizing Legendre polynomials
                    # sqrt(2*l - 1) == sqrt(2*(l-1) + 1)
                    # sqrt(2*l + 3) == sqrt(2*(l+1) + 1)
                    pll = Ptemp[l-1]/np.sqrt(2.0*l - 1.0)
                    plu = Ptemp[l+1]/np.sqrt(2.0*l + 3.0)
                    Pl[l] = (pll - plu)/2.0
            else:
                # Legendre polynomials up to LMAX
                alpha = np.cos(D[j,i]/rad_e)
                Ptemp,_ = gravtk.legendre_polynomials(LMAX, alpha)
                # unnormalizing Legendre polynomials
                Pl = np.squeeze(Ptemp)/norm
            # calculate Green's function
            G[j,i] = scale*np.sum(hl*Pl)
    # return the Green's function and the coordinates
    return (X, Y, G)
