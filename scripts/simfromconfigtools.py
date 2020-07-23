#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculations and plots for only a config file with no fits given
"""

# Importation des librairies
import numpy as np
import matplotlib.pyplot as plt

from cv2 import resize

from astroquery.skyview import SkyView

from astroquery.simbad import Simbad
customSimbad = Simbad()
customSimbad.add_votable_fields('ra(d), dec(d), flux(V)')
from astroquery.gaia import Gaia


from starsanalysis import *

# Definition des fonctions


def genSimulation(configpath, star, seeing = 1, orders = [-1,0,1]):
   config = configuration(configpath)
   imsize, pix2ars = config.ccd_imsize, config.pixel2arcsec

   query = customSimbad.query_object(star)
   ra0, dec0, mag0 = query["RA_d"][0], query["DEC_d"][0], query["FLUX_V"][0]
   r = 2*imsize*config.pixel2arcsec/3600
   
   maglim = mag0 + 5
   job = Gaia.launch_job_async("SELECT ra, dec, DISTANCE(POINT('ICRS', ra, dec),POINT('ICRS', {0}, {1})) AS dist, phot_g_mean_mag AS flux FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS', {0}, {1}, {2})) = 1 AND phot_g_mean_mag < {3} ORDER BY dist".format(ra0, dec0, r, maglim))
   result = job.get_results()
   
   ra0, dec0 = result['ra'][0], result['dec'][0]
   config.set_wcs(ra0,dec0)
   
   cra, cdec = config.wcs.wcs_pix2world([[imsize/2, imsize/2]], 0)[0]
   
   stars = list_stars(result,config,maglim,seeing, orders)
   angle, LTh, Overlap_list = best_angle(stars, config)
   
   for star in stars :
      star.rotate_orders(angle, use_radians = True)
   
   LTh = LTh*180/np.pi
   img = SkyView.get_images(position = '{}, {}'.format(cra, cdec), survey = 'DSS', width = imsize*pix2ars*u.arcsec, height = imsize*pix2ars*u.arcsec)[0][0].data
   
   
   return (stars, img, config, angle, LTh, Overlap_list)



def plot_all(stars, img, config, angle, LTh = [], Overlap_list = []):
   """
   plots a simplistic view of the spectra on the region of a star of interest, with dispersion of the spectra being tilted with an angle 'angle' along with how the contamination evolves for angles between 0 and pi. Angle counted counterclockwise from the X axis.
   
   Parameters
   ----------
   stars : list of star objects
      list of stars to plot
   img : ndarray
      image to place under the stars (the stars have to be aligned with the stars on this image)
   config : configuration object
      made from a configuration file.
   angle : float
      dispersion angle in radian
   LTh : list
      List of angles for the Overlap list
   Overlap_list : list
      the list(s) of overlap for the angles in LTh
   
   Returns
   -------
   out : None, 
      Plots the graph of Overlap_list(LTh) and the simulation of the spectra on the captor over img."""
   angle = 0
   name, imsize = config.obs_name, config.ccd_imsize
   n,m = img.shape
   
   if (n,m) == (300,300):
      n,m = imsize, imsize
   
   X0,Y0 = stars[0].x, stars[0].y
   
   fig = plt.figure(figsize = (14, 6))
   
   if len(LTh) == 0:
      ax2 = fig.add_subplot(111)
   else : 
      ax1, ax2 = fig.add_subplot(121), fig.add_subplot(122)
      ax1.plot(LTh, Overlap_list)
      ax1.set_title("Evolution of the contamination for different tilts")
      ax1.set_xlabel("Tilt angle (degree)")
      ax1.set_ylabel("Contamination (% of central spectrum surface area)")

   orders_shape = unary_union([star.all_orders for star in stars[1:]])
   ax2.add_patch(PolygonPatch(orders_shape, alpha = 0.7))
   
   ax2.add_patch(PolygonPatch(stars[0].all_orders, alpha = 0.7, color = 'red'))

   ax2.imshow(img, extent = (0, m, n, 0), cmap = 'gray')
   
   
   """for th, s in zip([0, angle], [ -1, 1]):
      Angle_0 = Rotation_Around(np.vstack([[X0, X0], [ - n+m, n+m]]), [X0, Y0], th - np.pi/2)
      ax2.plot(*Angle_0, 'white')
      M0 = np.abs(np.cos(th))*(Angle_0 - np.array([X0, Y0]))/5 + np.array([X0, Y0])
      Xt0, Yt0 = M0[0, 0] + s*200, M0[1, 0] + s*200 
      plt.text(Xt0, Yt0, 'Tilt : {:.0f}°'.format(180*(th)/np.pi), rotation = th*180/np.pi, rotation_mode = 'anchor', bbox = dict(boxstyle = "square", ec = (1., 0.5, 0.5), fc = (1., 0.8, 0.8), alpha = 0.5, color = 'white'))
   """

   ax2.axis('square')
   ax2.set_xlabel("X (pixels)")
   ax2.set_ylabel("Y (pixels)")
   ax2.set_xlim([0, m])
   ax2.set_ylim([0, n])
   ax2.set_title("{}, Deviation used : {:.2f}°".format(name, angle*180/np.pi))