import numpy as np
import matplotlib.pyplot as plt
from configparser import ConfigParser
from astropy.wcs import WCS


from geometrysl import *

#Classes definition

class configuration :
   """Regroups useful informations on a disperser"""
   def __init__(self, fichier = "default_conf"):
      """Gathers informations from a configuration file. \nfichier must be the path to a .ini file, containing the requiered informations (listed in the config/default.ini file).
      use "default_conf" as fichier to have a default configuration object"""
      parser = ConfigParser()
      if fichier == "default_conf":
         self.obs_name = "DEFAULT" # string
         self.ccd_imsize = 2048 # px
         self.pixel2mm = 24e-3 # mm/px
         self.pixel2arcsec = 0.401 # arcsec/px
         self.distance2ccd = 55 # mm
         self.grooves_per_mm = 350 # mm^ - 1
         self.lambda_min = 300*1e-6 # mm
         self.lambda_max = 1100*1e-6 # mm
      else : 
         parser.read(fichier)
         self.obs_name = parser.get('instrument', 'OBS_NAME') # string
         self.ccd_imsize = int(parser.get('CCD', 'CCD_IMSIZE')) # px
         self.pixel2arcsec = float(parser.get('CCD', 'CCD_PIXEL2ARCSEC')) # arcsec/px
         self.lambda_min = float(parser.get('spectrum range', 'LAMBDA_MIN'))*1e-6 # mm
         self.lambda_max = float(parser.get('spectrum range', 'LAMBDA_MAX'))*1e-6 # mm
         try :
            self.pixel2mm = float(parser.get('CCD', 'CCD_PIXEL2MM')) # mm/px
            self.distance2ccd = float(parser.get('dispersers', 'DISTANCE2CCD')) # mm
            self.grooves_per_mm = float(parser.get('dispersers', 'GROOVES_PER_MM')) # mm^ - 1
         except :
            self.dispersion_ratio = float(parser.get('dispersers', 'DISPERSION_RATIO'))*1e-7 # mm/px
      self.wcs = WCS(naxis = 2)
      self.wcs.wcs.cdelt = [ - self.pixel2arcsec/3600, self.pixel2arcsec/3600]
      self.wcs.wcs.crpix = [self.ccd_imsize//2 + 1, self.ccd_imsize//2 + 1]
      self.wcs.wcs.crval = [0, 0]
      self.wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]


class star_class :
   """Contains useful caracteristics of a star observed through a slitless spectroscope"""
   def __init__(self, ra, dec, mag,config):
      """Sets a star object from the given informations. Set all parameters to 'd' for a default star (HD107696 here)
      
      Parameters
      ----------
      ra : float
         Right ascension of the star in degree.
      dec : float
         Declination of the star in degree.
      mag : float
         Visible magnitude of the star.
      config : configuration
         configuration object containing useful informations for the analysis.
      """
      if ra == 'd' and dec == 'd' and mag == 'd':
         self.ra = 0
         self.dec = 0
         self.mag = 0
      else : 
         self.ra = ra
         self.dec = dec
         self.mag = mag
      self.order = {}
      if config == 'd':
         config = configuration()
      self.x, self.y = config.wcs.wcs_world2pix([[self.ra,self.dec]],0)[0]
      self.all_orders = Polygon()
   
   def set_orders(self, Ln: list, config, maglim = 15, seeing = 1):
      """Sets an group of polygons representing the orders listed in Ln.
      
      Parameters
      ----------
      Ln : list
         List of integers, each integers represents an order to set.
      config : configuration
         configuration object containing useful informations for the analysis.
      maglim : float 
         highest magnitude visible by the captor.
      seeing : float
         seeing, FWHM of the point spread function of the telescope.
      """
      for n in Ln:
         self.order[n] = Order(n,self.x,self.y,self.mag,config, maglim, seeing)
         self.all_orders = unary_union([self.order[n] for n in self.order])
   def rotate_orders(self, angle, use_radians = False):
      for n in self.order:
         self.order[n] = rotate(self.order[n], angle, [self.x, self.y], use_radians)
      self.all_orders = rotate(self.all_orders, angle, [self.x, self.y], use_radians)
      
   def translate(self, x_off = 0,y_off = 0):
      self.all_orders = translate(self.all_orders,x_off, y_off)
      for i in self.order:
         self.order[i] = translate(self.order[i],x_off, y_off)
      self.x, self.y = self.x + x_off, self.y + y_off


#Functions definition

def gen_orders(Ln, stars, config, maglim, seeing):
   """Sets the orders contained in Ln for all the stars in stars.
   
   Parameters
   ----------
   Ln : list
      List of integers, each integers represents an order to set.
   stars : list
      List of stars to set orders to.
   config : configuration
      configuration object containing useful informations for the analysis.
   maglim : float 
      highest magnitude visible by the captor.
   seeing : float
      seeing, FWHM of the point spread function of the telescope.
   """
   for i in stars :
      i.set_orders(Ln, config, maglim, seeing)

def list_stars(table,config, maglim, seeing, orders):
   """Sets and lists stars from a table in star objects
   
   Parameters
   ----------
   table : astropy table or equivalent
      Contains informations on the visible stars on the ccd, the column containing the right ascensions of the stars must contain "ra" in its name, "dec" for the declination and "flux" for the visible magnitude. These three informations must be present for all of the stars. The angles must be given in degree
   config : configuration
      configuration object containing useful informations for the analysis.
   maglim : float 
      highest magnitude visible by the captor.
   seeing : float
      seeing, FWHM of the point spread function of the telescope.
   orders : list
      orders to take into account
      
   Returns
   -------
   L : list
      list of star objects with their orders from the orders list set.
   """
   L = []
   
   ra_colname = [table.colnames[i] for i, j in enumerate(table.colnames) if 'ra' in j.lower()][0]
   dec_colname = [table.colnames[i] for i, j in enumerate(table.colnames) if 'dec' in j.lower()][0]
   mag_colname = [table.colnames[i] for i, j in enumerate(table.colnames) if 'flux' in j.lower()][0]
   
   
   for ra, dec, mag in zip(table[ra_colname], table[dec_colname],table[mag_colname]):
      if np.isreal(mag) :
         L.append(star_class(ra,dec,mag,config))
   gen_orders(orders,L,config,maglim,seeing)
   return L

def rotate_stars(stars, centre, angle):
   """Rotates the positions of a set of stars by an angle "angle" 
   
   Parameters
   ----------
   stars : list of star objects
   centre : [float, float]
      Position of the center of rotation in pixel
   angle : float
      the angle of rotation in radian
   """
   n = len(stars)
   M = np.zeros([2,n])
   for i in range(n):
      M[:,i] = [stars[i].x,stars[i].y]
   Mr = Rotation_Around(M, centre, angle)
   for i in range(n):
      x, y = stars[i].x, stars[i].y
      new_x, new_y = Mr[:,i]
      stars[i].x, stars[i].y = Mr[:,i]
      for n in stars[i].order :
         stars[i].order[n] = translate(stars[i].order[n], new_x - x, new_y - y)
      stars[i].all_orders = translate(stars[i].all_orders, new_x - x, new_y - y)
   
