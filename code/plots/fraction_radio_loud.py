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


def plot_09bb(ax):
    z = 0.009937
    dat = np.loadtxt(ddir + "/2009bb_49ghz.txt", delimiter=',')
    col = orag
    lum = ujy_to_flux(dat[:,1]*1000, z)
    ax.plot(dat[:,0], lum, c=col, label="_nolegend_", lw=2)


def plot_100316D(ax):
    z = 0.0593
    dat = Table.read(ddir + "/100316D_obs.dat", format='ascii.csv')

    col = orag
    freq = dat['freq']
    choose = freq== 5.4
    t = dat['dt'][choose]
    lum = ujy_to_flux(dat['flux'][choose], z)
    toplot = dat['fluxerr'][choose] > 0
    ax.plot(t[toplot], lum[toplot], c=col, label="_nolegend_")


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
    dat = np.loadtxt("../../data/ptfdby.txt", delimiter=',')
    ax.plot(dat[:,0], dat[:,1], c=dark, lw=0.5, ls='--', label="_nolegend_")

    # iPTF11qcj
    dat = np.loadtxt("../../data/ptf11qcj.txt", delimiter=',')
    ax.plot(dat[:,0], dat[:,1], c=dark, lw=0.5, ls='--', label="_nolegend_")

    # 11cmh
    ax.plot(
            [171,991], [4E28, 4.5E27], lw=0.5, 
            c=dark, ls='--', label="_nolegend_")


fig,ax = plt.subplots(1,1,figsize=(6,5))

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
plot_100316D(ax)

ax.set_yscale('log')
ax.set_xscale('log')
ax.tick_params(axis='both', labelsize=16)
ax.set_xlim(2,160)
ax.set_xlabel(r"Time Since GRB/Discovery (days)", fontsize=16)
ax.set_ylabel(r'4--6 GHz Radio Luminosity (erg/s/Hz)', fontsize=16)
plt.tight_layout()

plt.show()
#plt.savefig("frac_det.eps", format='eps', dpi=500)
