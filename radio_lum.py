""" Plot of the radio luminosity over time

Include:
LLGRBs
1998bw
060218
100316D

relativistic SNe
2009bb
2012ap

these are around 5 GHz
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from astropy.io import ascii
rc("font", family="serif")
from astropy.cosmology import Planck15
import astropy.units as u


# current values from the data, uJy, z=0.054
dt = np.array([6, 11, 23])
f = np.array([32.5, 29.6, 26.6])
flux = f * 10**-6 * 10**-23
z = 0.054
dL = Planck15.luminosity_distance(z=z)
lum = 4 * np.pi * dL.cgs.value**2 * flux / (1+z)
plt.scatter(dt, lum, c='k', marker='s')
#plt.scatter(dt+10, lum, marker='s', facecolors='white', edgecolors='k')


# mJy, z=0.105
dat = ascii.read("data/031203_obs.dat", delimiter=',')
dt = dat['dt']
flux = dat['flux'] * 10**-3 * 10**-23
z = 0.105
dL = Planck15.luminosity_distance(z=0.105)
lum = 4 * np.pi * dL.cgs.value**2 * flux / (1+z)
plt.plot(dt, lum, c='grey', alpha=0.5)
plt.text(dt[0], lum[0]*1.5, s='031203')


# uJy, z = 0.033
dat = ascii.read("data/060218_obs.dat", delimiter=',')
dt = dat['dt']
flux = dat['flux'] * 10**-6 * 10**-23
z = 0.033
dL = Planck15.luminosity_distance(z=z)
lum = 4 * np.pi * dL.cgs.value**2 * flux / (1+z)
plt.plot(dt, lum, c='grey', alpha=0.5)
plt.text(dt[2], lum[2]*1.2, s='060218', horizontalalignment='center')

# uJy, z = 0.0591
dat = ascii.read("data/100316D_obs.dat", delimiter=',')
dt = dat['dt']
flux = dat['flux'] * 10**-6 * 10**-23
z = 0.0591
dL = Planck15.luminosity_distance(z=z)
lum = 4 * np.pi * dL.cgs.value**2 * flux / (1+z)
choose = dat['fluxerr'] > 0
plt.plot(dt[choose], lum[choose], c='grey', alpha=0.5)
plt.text(dt[2], lum[2]*1.2, s='100316D', horizontalalignment='center')

# mJy, z = 0.0085
dat = ascii.read("data/1998bw_obs.dat", delimiter=',')
dt = dat['dt']
flux = dat['flux'] * 10**-3 * 10**-23
z = 0.0085
dL = Planck15.luminosity_distance(z=z)
lum = 4 * np.pi * dL.cgs.value**2 * flux / (1+z)
plt.plot(dt, lum, c='grey', alpha=0.5)
plt.text(dt[2], lum[2]*1.2, s='1998bw', horizontalalignment='center')

# mJy, z = 0.009
dat = ascii.read("data/2009bb_obs.dat", delimiter=',')
dt = dat['dt']
flux = dat['flux'] * 10**-3 * 10**-23
z = 0.009
dL = Planck15.luminosity_distance(z=z)
lum = 4 * np.pi * dL.cgs.value**2 * flux / (1+z)
plt.plot(dt, lum, c='grey', alpha=0.5)
plt.text(dt[1], lum[1]*1.2, s='2009bb', horizontalalignment='center')


# mJy, z = 0.009
dat = ascii.read("data/2012ap_obs.dat", delimiter=',')
dt = dat['dt']
flux = dat['flux'] * 10**-3 * 10**-23
z = 0.009
dL = Planck15.luminosity_distance(z=z)
lum = 4 * np.pi * dL.cgs.value**2 * flux / (1+z)
plt.scatter(dt, lum, c='grey', alpha=0.5)
plt.text(dt[0], lum[0]*1.2, s='2012ap', horizontalalignment='center')


plt.xlim(0,40)
plt.xlabel("Days since discovery", fontsize=16)
#plt.ylabel("$4 \pi d_L^2 F_\mathrm{radio} / (1+z)$ erg/Hz/s", fontsize=16)
plt.ylabel(r"$L_{\nu=5\,\mathrm{GHz}}$ [erg/Hz/s]", fontsize=16)
plt.yscale('log')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.tight_layout()

#plt.show()
plt.savefig("radio_lum.png")
