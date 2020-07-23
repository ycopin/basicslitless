#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Descriptif du fichier
"""

# Importation des librairies
import numpy as np
import matplotlib.pyplot as plt
import cv2,time

from astropy.io import fits
from astropy.visualization import ZScaleInterval
zscale = ZScaleInterval()

from astroquery.simbad import Simbad
customSimbad = Simbad()
customSimbad.add_votable_fields('ra(d), dec(d)')
from astroquery.gaia import Gaia

from starsanalysis import *

# Definition des fonctions

def align_order1(img,stars):
   """Returns the angle of dispersion for the first order of a a given star
   The first star must be aligned with its representation on the image (see the wcs to do that).
   
   The angle is determined by finding the maximum of S(angle) with S the sum of the pixels' values under the first order shape of the star. Since the first order can be approximated as a long rectangle of bright pixels, this method works most of the time.
   
   Parameters
   ----------
   img : ndarray-like
      Gray scale numpy image of the stars on a ccd.
   stars : list of star_class object
      The stars containing the first orders to align. stars[0].x and stars[0].y must correspond to a bright star position on the ccd.
   """
   n,m = img.shape
   
   order1 = stars[0].order[1]
   
   X0, Y0 = stars[0].x, stars[0].y

   angles = np.linspace(-np.pi,np.pi,100)
   L = [0]*len(angles)

   x, y = np.meshgrid(np.arange(m), np.arange(n)) # make a canvas with coordinates
   x, y = x.flatten(), y.flatten()
   points = np.vstack((x,y)).T
   grid = PolygonPatch(order1).contains_points(points)
   in_points = points[grid]
   
   tot = np.sum(img)
   
   img_copy = np.copy(img).T
   
   for i in range(len(angles)) :
      angle = angles[i]
      newpoints = np.floor(Rotation_Around(in_points.T, [X0,Y0], angle).T).astype('int')
      newpoints = newpoints[(newpoints[:,0]>=0) * (newpoints[:,1]>=0) * (newpoints[:,0]<m) * (newpoints[:,1]<n)]
      L[i] = img_copy[newpoints[:,0],newpoints[:,1]].sum()
   i_max = [i for i, j in enumerate(L) if j == max(L)][0]
   
   dispersion_angle = angles[i_max]
   for star in stars :
      star.rotate_orders(dispersion_angle, use_radians = True)
   return(L, dispersion_angle)



def align_stars(img = None, stars = None):
   """Aligns a set of stars to a picture, with the first star being on the right place
   Parameters 
   ----------
   img : 2d ndarray
      image to align the stars to, the projection must be the same as the one used for the stars. Only the angle must be off.
   stars : list of star objects to align
      The first star has to be on the right positon, with the right dispersion.
   """
   print("Work in progress, still doesn't work")
   return None
   cstar = stars[0]
   order = cstar.all_orders
   X0, Y0 = cstar.x, cstar.y
   
   
   img_copy = ((img - np.min(img))/(np.max(img)-np.min(img))*255).astype('uint8').T
   
   
   n,m = img_copy.shape
   x, y = np.meshgrid(np.arange(n), np.arange(m)) # make a canvas with coordinates
   x, y = x.flatten(), y.flatten()
   points = np.vstack((x,y)).T
   grid = PolygonPatch(order).contains_points(points)
   in_points = points[grid]


   med = np.median(img_copy)


   img_copy[in_points[:,0], in_points[:,1]] = med
   kernel = np.ones((5,5),np.float32)/25
   dst = cv2.filter2D(img_copy,-1,kernel)
   plt.imshow(dst) , plt.show()
   
   
   order0s = [star.order[0] for star in stars]
   
   
   angles = np.linspace(-np.pi,np.pi,100)
   L = [0]*len(angles)
   
   
   tot = np.sum(dst)
   
   for i in range(len(angles)) :
      angle = angles[i]
      for order0 in order0s:
         shape = rotate(order0, angle, [X0,Y0],True)
         a,b,c,d = np.floor(shape.bounds)
         if a>0 and c<n and b>0 and d<m :
            L[i] += dst[int(a):int(c),int(b):int(d)].sum()
   
   i_max = [i for i, j in enumerate(L) if j == max(L)][0]
   
   ccd_angle = angles[i_max]
   
   rotate_stars(stars, [X0,Y0], angle)
   return(L, ccd_angle)
   
   
def max_in_area(image,xmin,ymin,xmax,ymax):
   """Gives the maximum value inside a delimited rectangle in a gray-scale image
   image : 2d ndarray
      image to analyse
   xmin, ymin : integers
      bottom left coordinates of the rectangle
   xmax, ymax : integers
      top right coordinates of the rectangle
   
   Returns
   -------
   M : maximum values inside of the rectangle.
   """
   M = np.max(image[int(ymin):int(ymax),int(xmin):int(xmax)])
   L = []
   n,m = np.shape(image)
   for i in range(n):
      if M in image[i]:
         for j in [k for k, j in enumerate(image[i]) if j == M]:
            L.append([i,j])
   return M


def intensity2mag(I,refmag,refI):
   """Converts an intensity to a magnitude
   
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
   """Converts a magnitude to an intensity
   
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

def genSuperposition(configpath, fitspath, star, pos):
   """Generates a superposition of the geometric representation of  the arguments for plotall from these parameters
   Parameters
   ----------
   configpath : string
      Path to the configuration file (.ini) for the telescope used.
   fitspath : string
      Path to the fits file to be covered
   star : string
      Identifier of the most important star on the fits (brightest, or just the most interesting one)
   pos : [int, int]
      coordinates of the star on the fits
   """
   config = configuration(configpath)
   imsize, pix2ars = config.ccd_imsize, config.pixel2arcsec
   X0, Y0 = pos

   # Init picture: 
   head = fits.getheader(fitspath)
   img = fits.getdata(fitspath)
   n,m = img.shape
   # Init Geometry: 
   query = customSimbad.query_object(star)
   ra0, dec0 = query["RA_d"][0], query["DEC_d"][0]

   mag0 = Gaia.launch_job_async("SELECT phot_g_mean_mag AS flux, POWER(ra - {0},2) +  POWER(dec - {1},2) AS dist FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS', {0}, {1}, {2})) = 1 ORDER BY dist".format(ra0, dec0, 1/360)).get_results()['flux'][0]
   
   I0 = max_in_area(img,X0-n/10, Y0 - m/10, X0 + n/10, Y0 + m/10)
   
   maglim = intensity2mag(np.mean(img), mag0, I0)
   
   try :
      seeing = head['SEEING']
   except :
      seeing = 1
   # Finding stars: 
   r = 2*max(n,m)*config.pixel2arcsec/3600

   job = Gaia.launch_job_async("SELECT ra, dec, POWER(ra - ({0}),2) +  POWER(dec - ({1}),2) AS dist, phot_g_mean_mag AS flux FROM gaiadr1.gaia_source WHERE CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS', {0}, {1}, {2})) = 1 AND phot_g_mean_mag < {3} ORDER BY dist".format(ra0, dec0, r, maglim))
   
   result = job.get_results()

   # Init WCS: 
   ra0, dec0 = result['ra'][0], result['dec'][0]
   config.set_wcs(ra0,dec0)
   config.wcs.wcs.crpix = [X0, Y0]

   # Applying WCS
   
   stars = list_stars(result,config,maglim,seeing,[-1,0,1])
   
   # Angles Determination:
   L1,A1 = align_order1(img,stars)
   print('Aligned Orders')
   L2,A2 = align_stars(img, stars)
   print("Aligned Stars")
   
   LTh = np.linspace(-180,180,len(L2))

   return (stars, zscale(img), config, A1, LTh, L2)