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
    ax.plot(t[toplot], lum[toplot], c=col, label="_nolegend_", lw=2)


def plot_ztf_lims(ax):
    dt = [5, 19, 46, 100.4, 68.2, 37.95]
    lums = [7.599669433717404E26, 8.327790693956049E26, 
        3.0907667047489497E26, 1.237984570092877E27,
        2.3488805252832247E27, 3.711388552298804E27]
    ax.scatter(dt, lums, marker='v', c='k')


def plot_ptf_lims(ax):
    dat = np.loadtxt("../../data/ptf_lims.txt", delimiter=',')
    ax.scatter(dat[:,0], dat[:,1], marker='v', 
            facecolor='white', edgecolor='k')


def plot_ztf_det(ax):
    ax.plot(
            [16, 22, 34], 
            ujy_to_flux(np.array([32.5, 29.6, 26.6]), 0.05403),
            c=dark, lw=2, label="_nolegend_")
    ax.errorbar(
        [16, 22, 34],
        ujy_to_flux(np.array([32.5, 29.6, 26.6]), 0.05403),
        yerr=ujy_to_flux(np.array([7.1, 5.3, 5.4]), 0.05403),
        c=dark, fmt='o')


def plot_iptf_det(ax):
    # iPTF14dby
    dt = [20.42, 44.18, 65.44, 75.41, 106.30, 154.16]
    lum = [4.3E28, 6.8E28, 5.1E28, 5.2E28, 5.4E28, 4.7E28]
    ax.plot(dt, lum, c=dark, lw=2, ls='--', label="_nolegend_")


fig,ax = plt.subplots(1,1,figsize=(8,6))

dark = '#140b34'
purp = '#84206b'
orag = '#e55c30'
yell = '#f6d746'

ddir = "/Users/annaho/Dropbox/Projects/Research/IcBL/data/radio_compilations"

plot_98bw(ax)
plot_06aj(ax)
plot_ztf_lims(ax)
plot_ptf_lims(ax)
plot_ztf_det(ax)
plot_iptf_det(ax)

ax.set_yscale('log')
ax.set_xscale('log')
ax.tick_params(axis='both', labelsize=16)
ax.set_xlabel(r"Time Since GRB/Discovery (days)", fontsize=16)
ax.set_ylabel(r'4--6 GHz Radio Luminosity (erg/s/Hz)', fontsize=16)

plt.show()
