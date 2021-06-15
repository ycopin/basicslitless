# Time-stamp: <2021-06-08 12:57:21 ycopin>

from simulationtools import plot_all, genSimulation
from classes import *

import matplotlib.pyplot as plt

import time

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

parser.add_argument("config", help="Configuration file.ini")
parser.add_argument("star", help="Star identifier")
parser.add_argument("-s", "--seeing",
                    help="Seeing in arcsec",
                    default=1., type=float)
parser.add_argument("-O", "--orders",
                    help="Dispersion orders",
                    default="-1,0,+1")
parser.add_argument("-o", "--output",
                    help="Output figure")

args = parser.parse_args()

orders = [ int(i) for i in args.orders.split(',') ]

stars, img, angles, config, LTh, Overlap_list = genSimulation(
    args.config, args.star, orders, args.seeing)

fig = plot_all(stars, img, angles, config, LTh, Overlap_list)
if args.output:
    fig.savefig(args.output)
else:
    plt.show()
