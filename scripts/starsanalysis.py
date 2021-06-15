# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module for slitless spectroscopy basic simulations """

__author__ = "Sadek Al-Jibouri <sadek.aljibouri@ens-lyon.fr>"

# Importation des librairies

from classes import *

# Classes :


def total_overlap(cstar, residues, angle):
    """
    Gives the overlap on the PSF of the central star caused by its
    surroundings.

    Parameters
    ----------
    cstars : Star object
       The central star
    residues :list of stars
       The list of residual stars, the ones that could contaminate the central
       one
    angle : float
       angle of dispersion of the PSF on the captor in radian, an angle of zero
       corresponds to a dispersion along the X axis.

    Returns
    -------
    out : the percentage of overlapping on the central star spectrum caused by
    its surroundings.
    """
    try:
        residues[0]
    except BaseException:
        residues = [residues]
    Xl = np.array([star.x for star in residues])
    Yl = np.array([star.y for star in residues])
    orders = [star.all_orders for star in residues]
    n = len(residues)
    cx, cy, shape0 = cstar.x, cstar.y, cstar.all_orders
    mag0 = cstar.mag

    Pos = np.vstack([Xl, Yl])

    M = rotate_around(Pos, [cx, cy], -angle)

    xmin, ymin, xmax, ymax = shape0.bounds
    w = xmax - xmin
    h = ymax - ymin

    filtre = (np.abs(M[0]-cx) <= 2*w)*(np.abs(M[1]-cy) <= 2*h)
    L = []
    for x, x0, y, y0, i in zip(M[0][filtre], Xl[filtre], M[1][filtre],
                               Yl[filtre], np.arange(0, n)[filtre]):
        L.append(translate(orders[i], xoff=x - x0, yoff=y - y0))
    shape = unary_union(L)

    t = shape0.intersection(shape).area

    return 100*t/shape0.area


def best_angle(stars, config):
    """returns the optimal angle for the dispersion of a spectrum around a star along with how far one can deviate from this angle.

    Parameters
    ----------
    stars : list of star objects
       The first star will be considered as the important one
    config : Configuration
       made from a configuration file.

    Returns
    -------
    b_a : float
       the angle (in radian) with the less contamination on the central star. 0 corresponds to a dispersion along the X axis and the angle is then mesured counterclockwise.
    wid : float
       the maximal deviation from this angle one can take without having too high a dispersion
    Overlap_list : list
       the list of contamination for 100 angles evenly spaced in between 0 and pi
    """

    cstar = stars[0]
    residues = stars[1:]

    LTh = np.linspace(-np.pi/2, np.pi/2, 100)

    Overlap_list = np.array([total_overlap(cstar, residues, th) for th in LTh])

    Mi, Ma, n = min(Overlap_list), max(Overlap_list), len(Overlap_list)
    thres = (Ma - Mi)/5 + Mi
    if max(Overlap_list) == min(Overlap_list):
        return 0, LTh, Overlap_list

    imax = [i for i, j in enumerate(Overlap_list) if j == Ma][0]
    LMin = []
    List_index = [i % n for i in range(imax, imax + n)]

    i = 0
    while i < n:
        if Overlap_list[List_index[i]] <= thres and Overlap_list[List_index[i - 1]] >= thres:
            LMin.append([List_index[i], 0])
        elif Overlap_list[List_index[i]] <= thres:
            LMin[-1][-1] += 1
        i += 1
    LMin = np.array(LMin)
    longest = max(LMin[:, 1])
    best = [i for i, j in enumerate(LMin[:, 1]) if j == longest][0]

    angle, wid = LTh[LMin[best, 0]], longest*np.pi/200

    if (angle + wid) % (np.pi) > np.pi/2:
        b_a = - np.pi + (angle + wid) % (np.pi)
    else:
        b_a = (angle + wid) % (np.pi)
    return b_a, LTh, Overlap_list
