#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculations and plots for only a config file with no fits given
"""

# Importation des librairies
import numpy as np
import matplotlib.pyplot as plt

#from astropy.visualization import ZScaleInterval
#from astropy.io import fits
import astropy.units as u
from astropy.wcs import WCS

from astroquery.skyview import SkyView
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad

from shapely.ops import unary_union
from polygon import PolygonCollection

import starsanalysis as SA
from classes import Configuration, list_stars

customSimbad = Simbad()
customSimbad.add_votable_fields('ra(d), dec(d), flux(V)')


def read_config(configpath):

    print(f"Reading configuration {configpath!r}")
    config = Configuration(configpath)

    return config


def genSimulation(config, starname, orders=[-1, 0, 1], seeing=1., magoffset=5., survey="DSS"):
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
    magoffset : float
       limit magnitude offset, 5 mag by default
    survey : string
       SkyView survey for the background image
       (see https://astroquery.readthedocs.io/en/latest/skyview/skyview.html)

    Returns
    -------
    stars : list
       a list of Star objects, containing the stars to take into account
       with their position on the detector,
    img : 2d ndarray
       a background image for the stars,
    angle : float
       the angle of dispersion used for the spectrogram [rad],
    angles : list
       a list of angles [deg] for contamination,
    contaminations : list
       the list of the total overlap caused by surrounding stars on the
       studied star spectrogram
    """

    print(f"CCD size: {config.ccd_imsize}, scale: {config.pixel2arcsec}\"/px")

    print(f"Querying Simbad for {starname!r}...")
    query = customSimbad.query_object(starname)
    if query is None:
        raise RuntimeError(f"Unresolved object {starname!r}.")

    query = query.filled()  # Turn masked values to NaN
    ra0, dec0, mag0 = query["RA_d"][0], query["DEC_d"][0], query["FLUX_V"][0]

    print(f"{starname!r}: RA={ra0:.6f}, Dec={dec0:.6f}, V={mag0:.2f}")
    if np.isnan(mag0):
        print("No V-mag available, default to 17.")
        mag0 = 17.

    radius = config.ccd_imsize * config.pixel2arcsec / 3600  # [deg]
    maglim = mag0 + magoffset

    catalog = "gaiadr2.gaia_source"
    mag = "phot_g_mean_mag"
    print(f"Querying Gaia:{catalog} for {starname!r}, {mag} < {maglim:.2f}:")
    query = (f"SELECT ra, dec, "
             f"DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {ra0}, {dec0})) AS dist, "
             f"{mag} AS mag "
             f"FROM {catalog} "
             f"WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra0}, {dec0}, {radius})) = 1 "
             f"AND {mag} < {maglim} "
             f"ORDER BY dist")
    print(query)
    job = Gaia.launch_job_async(query)
    results = job.get_results()
    print(f"{len(results)} entries")

    ra0, dec0 = results['ra'][0], results['dec'][0]
    config.wcs.wcs.crval = [ra0, dec0]

    stars = list_stars(results, config, maglim, seeing, orders)
    best_angle, angles, contaminations = SA.best_angle(stars, config)

    for star in stars:
        star.rotate_orders(best_angle, use_radians=True)

    angles = angles * 180 / np.pi # [deg]

    cra, cdec = config.wcs.wcs_pix2world([[config.ccd_imsize/2] * 2], 0)[0]
    width  = config.ccd_imsize * config.pixel2arcsec * u.arcsec  # [arcsec]
    print(f"Querying SkyView.{survey} around RA={cra:.6f}, Dec={cdec:.6f} "
          f"[{width:.2f}\"×{width:.2f}\"]...")
    hdu = SkyView.get_images(position=f"{cra}, {cdec}",
                             survey=survey, width=width, height=width)[0][0]

    return (stars, hdu, best_angle, angles, contaminations)


def plot_all(stars, hdu, best_angle, config, angles=None, contaminations=None):
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
    best_angle : float
       best dispersion angle in radian
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

    n, m = hdu.data.shape
    if (n, m) == (300, 300):
        n = m = config.ccd_imsize

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

    # Plot spectrogram of reference star
    ax2.add_collection(PolygonCollection(stars[0].all_orders,
                                         alpha=1, lw=3, color='red'))

    # Plot spectrograms of all neighboring stars
    # ax2.add_collection(PolygonCollection(
    #     unary_union([star.all_orders for star in stars[1:]]),
    #     alpha=1, color='white', ec='black'))
    ax2.add_collection(PolygonCollection(
        unary_union([star.all_orders for star in stars[1:]]),
        alpha=1, lw=2, color='blue', ec='blue'))

    # Background image
    ax2.imshow(hdu.data, extent=(0, m, n, 0), cmap='gray_r')

    ax2.axis('square')
    ax2.set(xlabel="X [px]", ylabel="Y [px]",
            xlim=[0, m], ylim=[0, n],
            title=f"{config.obs_name}, dispersion direction: {best_angle*180/np.pi:.2f}°")

    return fig


def plot_angles(angles, contaminations, best_angle, ax=None):

    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(angles.to(u.degree), contaminations)
    ax.set(title=f"Best angle: {best_angle.to(u.deg):.1f}",
           xlabel="Dispersion angle [°]",
           ylabel="Contamination [%]")

    return ax


def show_scene(hdu, stars, config, ax=None):

    bkgimg = hdu.data
    n, m = bkgimg.shape
    if (n, m) == (300, 300):
        n = m = config.ccd_imsize

    # Background image
    print("Background image WCS")
    wcs = WCS(hdu.header)
    wcs.printwcs()

    print("Configuration WCS")
    config.wcs.printwcs()

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, projection=wcs)
    else:                       # Replace old ax by new one at same position
        ax2 = ax.figure.add_axes(ax.get_position().bounds, projection=wcs)
        ax.remove()
        ax = ax2

    ax.imshow(hdu.data, cmap='gray_r')
    ax.coords.grid(True, color='green', ls='solid')
    ax.coords[0].set_axislabel(wcs.wcs.ctype[0])
    ax.coords[1].set_axislabel(wcs.wcs.ctype[1])

    # ra, dec, mag = np.array([ [s.ra, s.dec, s.mag] for s in stars ]).T
    # ax.scatter(ra, dec, s=20*10**(0.4*(mag[0] - mag)), marker='+',
    #            transform=ax.get_transform('world'))

    # Plot spectrogram of reference star
    ax.add_collection(PolygonCollection(stars[0].all_orders,
                                        alpha=1, lw=3, color='red',
                                        transform=ax.get_transform(config.wcs)))

    # Plot spectrograms of all neighboring stars
    ax.add_collection(PolygonCollection(
        unary_union([star.all_orders for star in stars[1:]]),
        alpha=1, lw=2, color='blue', ec='blue',
        transform=ax.get_transform(config.wcs)))

    return ax
