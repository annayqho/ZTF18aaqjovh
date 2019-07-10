""" Constrain the fraction of Ic-BL SNe with emission comparable
to SN2006aj and SN1998bw """


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.cosmology import Planck15
from astropy.table import Table


ddir = "/Users/annaho/Dropbox/Projects/Research/IcBL/data/radio_compilations"
dat = Table.read(
    ddir + "/1998bw.dat",
    delimiter="&", format='ascii.fixed_width')
freq = dat['freq']
choose = freq == 2.49 # closest to 3 GHz
t = dat['dt'][choose]
flux = dat['flux'][choose] * 1e-3 * 10**(-23)
d = Planck15.luminosity_distance(z=0.0085).cgs.value
lum = flux * 4 * np.pi * d**2
ax.plot(t, lum, c=col, label="_nolegend_")
