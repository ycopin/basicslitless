## Overlap And Query Tests :

from ..scripts.stars_analysis import *
from runSimulation import *

query = customSimbad.query_object('HD107696')
config = configuration('./Ressources/auxtel.ini')
r_pix = 2*config.ccd_imsize
arcsperpix = config.pixel2arcsec
r = r_pix*arcsperpix/3600
ra0, dec0 = ra0, dec0 = customSimbad.query_object("HD107696")["RA_d"][0], customSimbad.query_object('HD107696')["DEC_d"][0]

job = Gaia.launch_job_async("SELECT ra, dec, POWER(ra - {0},2) +  POWER(dec - {1},2) AS dist, phot_g_mean_mag AS flux FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS', {0}, {1}, {2})) = 1 AND phot_g_mean_mag < {3} ORDER BY dist".format(ra0, dec0, r, 12))

result = job.get_results()
ra0, dec0 = result['ra'][0], result['dec'][0]

config.set_wcs(ra0,dec0)
stars = list_stars(result, config, 15,1)

cstar = stars[0]

residues = stars[1:]

t = best_angle(stars, config)