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

gim = fits.open(ddir + "/frame-g-005116-5-0061.resamp.fits")
g = gim[0].data # 0 to 0.1 looks good
r = fits.open( # 0 to 0.2 looks good
        ddir + "/frame-r-005116-5-0061.resamp.fits")[0].data
z = fits.open( # 0 to 0.6 looks good
        ddir + "/frame-z-005116-5-0061.resamp.fits")[0].data
u = fits.open( 
        ddir + "/frame-u-005116-5-0061.resamp.fits")[0].data
i = fits.open( 
        ddir + "/frame-i-005116-5-0061.resamp.fits")[0].data

# Figure out pos from header
head = gim[0].header
wcs = astropy.wcs.WCS(head)
target_pix = wcs.all_world2pix([(np.array([ra,dec], np.float_))], 1)[0]
xpos = target_pix[0]
ypos = target_pix[1]
 
imsize = 120
# Plot cutout
gcut = g[int(ypos-imsize):int(ypos+imsize),int(xpos-imsize):int(xpos+imsize)]
rcut = r[int(ypos-imsize):int(ypos+imsize),int(xpos-imsize):int(xpos+imsize)]
zcut = z[int(ypos-imsize):int(ypos+imsize),int(xpos-imsize):int(xpos+imsize)]
ucut = u[int(ypos-imsize):int(ypos+imsize),int(xpos-imsize):int(xpos+imsize)]
icut = i[int(ypos-imsize):int(ypos+imsize),int(xpos-imsize):int(xpos+imsize)]

# icut: vmin = 0, vmax = 0.5
# rcut: vmin = 0, vmax = 0.3
# gcut: vmin = 0, vmax = 0.2

rgb = make_lupton_rgb(
        icut/2, rcut/1.5, gcut, Q=4, stretch=0.2)
plt.imshow(rgb, origin='lower')

#plt.show()
plt.savefig("host.eps", format='eps', dpi=1000)

