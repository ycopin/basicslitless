from simfromconfigtools import plot_all
from simfromfitstools import genSuperposition
from classes import *

import numpy as np
import matplotlib.pyplot as plt

import time




if __name__ == "__main__":
   from argparse import ArgumentParser
   
   parser = ArgumentParser()
   
   parser.add_argument("config", help="Configuration file path", default = "./config/auxtel.ini")
   parser.add_argument("fits", help="FITS file path", default = "./ressources/auxtel_first_light.fits")
   parser.add_argument("star", help="Star identifier", default = "HD107696")
   parser.add_argument("pos", help="Observed star position on the CCD (format X,Y)", default = "0,0")


   args = parser.parse_args()
   
   xs, ys = args.pos.split(',')
   X0, Y0 = int(xs), int(ys)

   stars, img, config, angles, LTh, list = genSuperposition(args.config, args.fits, args.star, [X0, Y0])
   

   plot_all(stars, img, config, angles, LTh, list)
   plt.savefig("./Examples/Sim_from_fits_{}.png".format(time.strftime("%m_%d_%H%M%S")))
   plt.show()