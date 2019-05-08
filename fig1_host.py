""" Plot a host image from SDSS
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.wcs
from astropy import coordinates as coords
from astropy.visualization import make_lupton_rgb

ra = 178.1817541
dec = 25.6750332

ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18aaqjovh/data/host_sdss"

gim = fits.open(ddir + "/frame-g-005116-5-0061.fits.bz2")
g = gim[0].data # 0 to 0.1 looks good
u = fits.open( 
        ddir + "/frame-u-005116-5-0061.fits.bz2")[0].data
r = fits.open( # 0 to 0.2 looks good
        ddir + "/frame-r-005116-5-0061.fits.bz2")[0].data
z = fits.open( # 0 to 0.6 looks good
        ddir + "/frame-z-005116-5-0061.fits.bz2")[0].data

# Figure out pos from header
head = gim[0].header
wcs = astropy.wcs.WCS(head)
target_pix = wcs.all_world2pix([(np.array([ra,dec], np.float_))], 1)[0]
xpos = target_pix[0]
ypos = target_pix[1]
 
# Plot cutout
im = g[int(ypos-100):int(ypos+100),int(xpos-100):int(xpos+100)]
rgb = make_lupton_rgb(
        z/6, r/2, g, Q=10, stretch=0.5)
plt.imshow(
        rgb[int(ypos-100):int(ypos+100),int(xpos-100):int(xpos+100)], 
        origin='lower')

plt.show()
