""" Constrain the fraction of Ic-BL SNe with emission comparable
to SN2006aj and SN1998bw """


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.cosmology import Planck15
from astropy.table import Table
from radio_lc import ujy_to_flux


def plot_98bw(ax):
    dat = Table.read(
        ddir + "/1998bw.dat",
        delimiter="&", format='ascii.fixed_width')
    freq = dat['freq']
    choose = freq == 4.9 # closest to 3 GHz
    t = dat['dt'][choose]
    flux = dat['flux'][choose] * 1e-3 * 10**(-23)
    d = Planck15.luminosity_distance(z=0.0085).cgs.value
    lum = flux * 4 * np.pi * d**2
    col = orag
    ax.plot(t, lum, c=col, label="_nolegend_", lw=2)


def plot_06aj(ax):
    z = 0.03342
    dat = Table.read(ddir + "/060218_obs.dat", format='ascii.csv')
    col = orag
    freq = dat['freq']
    choose = freq== 4.86
    t = dat['dt'][choose]
    lum = ujy_to_flux(dat['flux'][choose], z)
    toplot = dat['fluxerr'][choose] > 0
    print(lum[toplot])
    ax.plot(t[toplot], lum[toplot], c=col, label="_nolegend_")


fig,ax = plt.subplots(1,1,figsize=(8,6))

dark = '#140b34'
purp = '#84206b'
orag = '#e55c30'
yell = '#f6d746'

ddir = "/Users/annaho/Dropbox/Projects/Research/IcBL/data/radio_compilations"

ax.set_yscale('log')
ax.tick_params(axis='both', labelsize=16)
ax.set_xlabel(r"$\Delta t$ (days)", fontsize=16)
ax.set_ylabel(r'4.9 GHz Radio Luminosity ($10^{27}$ erg/s)', fontsize=16)

plt.show()
