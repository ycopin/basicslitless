# Geometry Module

import numpy as np

from shapely.geometry import Polygon
from shapely.geometry.point import Point


def width(mag, maglim=20, seeing=1):
    """
    Gives the approximate size of a star on a captor, based on its magnitude,
    the highest magnitude visible and the seeing of the captor.

    Parameters
    ----------
    mag : float
       magnitude of the star.
    maglim : float
       highest magnitude visible by the captor.
    seeing : float
       seeing, FWHM of the point spread function of the atmosphere.

    Returns
    -------
    out : the size of the point on the captor, in arcsec.
    """

    if mag >= maglim:
        w = 0
    else:
        w = 0.58 * seeing * np.sqrt(maglim - mag)
    return w  # arcsec


def order_shape(n: int, x, y, mag, config, maglim, seeing):
    """
    Generates shapes representing the visible n'th order of the point spread
    function of a star on the captor in the case of slitless spectroscopy.

    Parameters
    ----------
    n : int
       spectrum order
    x, y : floats
       position of the star on the captor (along 0x and 0y axes).
    mag : float
       apparent magnitude of the star.
    config : Configuration object
       made from a configuration file.
    maglim : float
       highest visible magnitude on the captor.
    seeing : float
       seeing, FWHM of the point spread function of the telescope.

    Returns
    -------
    out : shapely Polygon object
    """

    try:
        disperserate = False
        lmin, lmax = config.lambda_min, config.lambda_max
        gpm = config.grooves_per_mm
        d2ccd = config.distance2ccd
        p2m, p2a = config.pixel2mm, config.pixel2arcsec
    except BaseException:
        disperserate = True
        lmin, lmax = config.lambda_min, config.lambda_max
        dr, p2a = config.dispersion_ratio, config.pixel2arcsec

    w = width(mag, maglim, seeing)/p2a

    if n == 0:
        order = Point(x, y).buffer(w)
    elif isinstance(n, int):
        if disperserate:
            hstart, hstop = lmin/dr, lmax/dr
        else:
            hstart = n*np.tan(np.arcsin(lmin*gpm))*d2ccd/p2m
            hstop = n*np.tan(np.arcsin(lmax*gpm))*d2ccd/p2m

        xmin, xmax = x + hstart, x + hstop
        ymin, ymax = y - w/abs(n), y + w/abs(n)
        order = Polygon([[xmin, ymin], [xmin, ymax],
                         [xmax, ymax], [xmax, ymin]])

    return order


def rotate_around(matrix, centre, angle):
    """
    Rotates a collection of points (2xn matrix) by an angle around a centre.

    Parameters
    ----------
    matrix : 2xn array_like
       concatenation of n points on a 2D plain.
    centre : 2x1 array_like
       centre around which the points will rotate.
    angle : float
       rotation angle, in radian.

    Returns
    -------
    out : 2xn numpy array.
       the same concatenation of points but rotated by the given angle around
       the centre.
    """

    RotMat = np.array(((np.cos(angle), -np.sin(angle)),
                       (np.sin(angle), np.cos(angle))))
    n = matrix.shape[-1]
    X0, Y0 = centre
    CentreMat = np.repeat([[X0], [Y0]], n, axis=1)

    return np.dot(RotMat, matrix - CentreMat) + CentreMat
