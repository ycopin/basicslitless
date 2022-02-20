# Time-stamp: <2021-06-08 12:57:21 ycopin>

"""
Main graphical script.
"""

import argparse

from simulationtools import plot_all, genSimulation
#from classes import *

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument("config", help="Configuration file.ini")
parser.add_argument("starname", help="Star name (to be resolved)")
parser.add_argument("-s", "--seeing", help="Seeing in arcsec [%(default)d\"]",
                    type=float, default=1.)
parser.add_argument("-O", "--orders", help="Dispersion orders [%(default)s]",
                    default="-1,0,+1")
parser.add_argument("-m", "--magoffset", help="Magnitude offset [%(default)s]",
                    type=float, default=5.)
parser.add_argument("-o", "--output", help="Output figure")

args = parser.parse_args()

orders = [ int(i) for i in args.orders.split(',') ]
print(f"Dispersion orders:", orders)

stars, img, best_angle, config, angles, contaminations = genSimulation(
    args.config, args.starname, orders=orders, seeing=args.seeing, magoffset=args.magoffset)
print(f"{len(stars)} stars")
print(f"Best angle: {best_angle*57.29578:.1f}°")

fig = plot_all(stars, img, best_angle, config, angles, contaminations)
if args.output:
    fig.savefig(args.output)
else:
    plt.show()
