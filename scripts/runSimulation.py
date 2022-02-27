# Time-stamp: <2022-02-25 20:28:39 ycopin>

"""
Main graphical script.
"""

import argparse

import matplotlib.pyplot as plt
import astropy.units as u

import simulationtools as ST

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument("config", help="Configuration file.ini")
parser.add_argument("starname", help="Star name (to be resolved)")
parser.add_argument("-s", "--seeing", help="Seeing in arcsec [%(default)d\"]",
                    type=float, default=1.)
parser.add_argument("-O", "--orders", help="Dispersion orders [%(default)s]",
                    default="-1,0,+1")
parser.add_argument("-m", "--magoffset", help="Magnitude offset [%(default)s]",
                    type=float, default=5.)
parser.add_argument("-S", "--survey", help="Background image [%(default)s]",
                    default="DSS")
parser.add_argument("-o", "--output", help="Output figure")

args = parser.parse_args()

orders = [ int(i) for i in args.orders.split(',') ]
print("Dispersion orders:", orders)

config = ST.read_config(args.config)
stars, hdu, best_angle, angles, contaminations = ST.genSimulation(
    config, args.starname, orders=orders,
    seeing=args.seeing, magoffset=args.magoffset, survey=args.survey)
print(f"{len(stars)} stars")
print(f"Best angle: {best_angle*57.29578:.1f}°")  # [rad] to [deg]

fig = ST.plot_all(stars, hdu, best_angle, config, angles, contaminations)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
ST.plot_angles(angles * u.degree, contaminations * u.percent, best_angle * u.radian, ax=ax1)
ST.show_scene(hdu, stars, config, ax=ax2)
fig.tight_layout()

if args.output:
    fig.savefig(args.output)
else:
    plt.show()
