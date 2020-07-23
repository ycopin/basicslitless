from simfromconfigtools import plot_all, genSimulation
from classes import *

import numpy as np
import matplotlib.pyplot as plt

import time




if __name__ == "__main__":
   from argparse import ArgumentParser
   
   parser = ArgumentParser()
   
   parser.add_argument("config", help="Configuration file path", default = "./Ressources/auxtel.ini")
   parser.add_argument("star", help="Star identifier", default = "HD107696")
   parser.add_argument("-m", "--maglim", help="Highest visible magnitude", type = float, default=15)
   parser.add_argument("-s", "--seeing", help="Seeing in arcsec", default = 1, type = float)
   parser.add_argument("-o", "--orders", help ="Orders to take into account as a comma separated list (default : (-1,0,1))", default = "(-1,0,1)")


   args = parser.parse_args()
   
   order_list =[int(i) for i in args.orders[1:-1].split(',')]
   

   stars, img, angles, config, LTh, Overlap_list = genSimulation(args.config, args.star, args.maglim, args.seeing, order_list)
   

   plot_all(stars, img, angles, config, LTh, Overlap_list)
   plt.savefig("./Examples/Sim_from_config_{}.png".format(time.strftime("%m_%d_%H%M%S")))
   plt.show()