#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculations and plots for only a config file with no fits given
"""

# Importation des librairies
import numpy as np
import matplotlib.pyplot as plt

from astropy.visualization import ZScaleInterval
from astropy.io import fits
import astropy.units as u

from astroquery.skyview import SkyView
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad

from shapely.ops import unary_union
from polygon import PolygonCollection

from starsanalysis import best_angle
from classes import Configuration, list_stars

customSimbad = Simbad()
customSimbad.add_votable_fields('ra(d), dec(d), flux(V)')



def intensity2mag(I, refmag, refI):
    """Converts an intensity to a magnitude.

    Parameters
    ----------
    I : float
       intensity to convert.
    refmag : float
       reference magnitude.
    refI : float
       corresponding intensity

    Returns
    -------
    A magnitude corresponding to I
    """
    
    return refmag - 2.5*np.log10(I/refI)


def mag2intensity(mag, refI, refmag):
    """Converts a magnitude to an intensity.

    Parameters
    ----------
    mag : float
       magnitude to convert.
    refI : float
       reference intensity
    refmag : float
       corresponding magnitude.

    Returns
    -------
    An intensity corresponding to mag
    """
    
    return refI * 10**((refmag-mag)/2.5)


def genSimulation(configpath, starname, orders=[-1, 0, 1], seeing=1., magoffset=5.):
    """
    Generates the arguments for plot_all to simulate the spectrogram of star
    through a slitless spectrogram which properties are contained in the
    configuration file.

    Parameters
    ----------
    configpath : string
       path to the slitless spectrograph configutation file
    starname : string
       the star name, to be resolved
    orders : list of integers
       the list of orders to represent.
    seeing : float
       the wanted seeing for the spectrogram in arcseconds, 1 arcsec by default

    Returns
    -------
    stars : list
       a list of Star objects, containing the stars to take into account
       with their position on the detector,
    img : 2d ndarray
       a background image for the stars,
    angle : float
       the angle of dispersion used for the spectrogram [rad],
    config : Configuration
       the configuration used for the spectrogram,
    angles : list
       a list of angles [deg] for contamination,
    contaminations : list
       the list of the total overlap caused by surrounding stars on the
       studied star spectrogram
    """

    print(f"Reading configuration {configpath!r}")
    config = Configuration(configpath)
    
    imsize, pix2ars = config.ccd_imsize, config.pixel2arcsec
    print(f"CCD size: {imsize}, scale: {pix2ars}\"/px")

    print(f"Querying Simbad for {starname!r}...")
    query = customSimbad.query_object(starname)
    if query is None:
        raise RuntimeError(f"Unresolved object {starname!r}.")

    query = query.filled()  # Turn masked values to NaN
    ra0, dec0, mag0 = query["RA_d"][0], query["DEC_d"][0], query["FLUX_V"][0]

    print(f"{starname!r}: RA={ra0:.2f}, Dec={dec0:.2f}, V={mag0:.2f}")
    if np.isnan(mag0):
        print("No V-mag available, default to 17.")
        mag0 = 17.

    r = 2 * imsize * pix2ars / 3600
    maglim = mag0 + magoffset

    print(f"Querying Gaia for {starname!r}, mag < {maglim:.2f}...")
    job = Gaia.launch_job_async(
        f"SELECT ra, dec, DISTANCE(POINT('ICRS', ra, dec), "
        f"POINT('ICRS', {ra0}, {dec0})) AS dist, "
        f"phot_g_mean_mag AS flux "
        f"FROM gaiadr2.gaia_source "
        f"WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra0}, {dec0}, {r})) = 1 AND "
        f"phot_g_mean_mag < {maglim} ORDER BY dist")
    results = job.get_results()
    print(f"{len(results)} entries")

    ra0, dec0 = results['ra'][0], results['dec'][0]
    config.wcs.wcs.crval = [ra0, dec0]

    cra, cdec = config.wcs.wcs_pix2world([[imsize/2, imsize/2]], 0)[0]
    stars = list_stars(results, config, maglim, seeing, orders)
    angle, angles, contaminations = best_angle(stars, config)

    for star in stars:
        star.rotate_orders(angle, use_radians=True)

    angles = angles * 180 / np.pi # [deg]
    img = SkyView.get_images(position='{}, {}'.format(cra, cdec),
                             survey='DSS',
                             width=imsize*pix2ars*u.arcsec,
                             height=imsize*pix2ars*u.arcsec)[0][0].data

    return (stars, img, angle, config, angles, contaminations)


def plot_all(stars, img, angle, config, angles=None, contaminations=None):
    """
    plots a simplistic view of the spectra on the region of a star of interest,
    with dispersion of the spectra being tilted with an angle 'angle' along with
    how the contamination evolves for angles between 0 and pi. Angle counted
    counterclockwise from the X axis.

    Parameters
    ----------
    stars : list of star objects
       list of stars to plot
    img : ndarray
       image to place under the stars (the stars have to be aligned with the
       stars on this image)
    config : Configuration
       made from a configuration file.
    angle : float
       dispersion angle in radian
    angles : list
       List of angles [deg] for contamination.
    contaminations : list
       the list(s) of overlap for the angles.

    Returns
    -------
    out : None,
       Plots contamination as function of angle and a simulation of the spectrogram
       on the detector over img.
    """

    name, imsize = config.obs_name, config.ccd_imsize
    n, m = img.shape

    if (n, m) == (300, 300):
        n, m = imsize, imsize

    X0, Y0 = stars[0].x, stars[0].y

    fig = plt.figure(figsize=(14, 6))

    if angles is None:
        ax2 = fig.add_subplot(111)
    else:
        ax1, ax2 = fig.add_subplot(121), fig.add_subplot(122)
        ax1.plot(angles, contaminations)
        ax1.set_title("Evolution of the contamination for different tilts")
        ax1.set_xlabel("Tilt angle [°]")
        ax1.set_ylabel("Contamination [%]")

    orders_shape = unary_union([star.all_orders for star in stars[1:]])
    ax2.add_collection(PolygonCollection(orders_shape, alpha=1,
                                         color='white', ec='black'))

    ax2.add_collection(PolygonCollection(
        unary_union([star.order[0] for star in stars[1:]]),
        color='white'))

    ax2.add_collection(PolygonCollection(stars[0].all_orders, alpha=1, color='red'))

    ax2.imshow(img, extent=(0, m, n, 0), cmap='gray_r')

    ax2.axis('square')
    ax2.set_xlabel("X (pixels)")
    ax2.set_ylabel("Y (pixels)")
    ax2.set_xlim([0, m])
    ax2.set_ylim([0, n])
    ax2.set_title(f"{name}, dispersion direction: {angle*180/np.pi:.2f}°")

    return fig
