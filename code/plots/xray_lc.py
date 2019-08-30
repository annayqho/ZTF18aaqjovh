""" Plot a grid of radio light curves """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15
import glob

dark = '#140b34'
purp = '#84206b'
orag = '#e55c30'
yell = '#f6d746'

ddir = "/Users/annaho/Dropbox/Projects/Research/IcBL/data/radio_compilations"


def plot_source(ax):
    # Chandra points
    ax.scatter(33, 1.9E40, c=dark, marker='v')


def plot_98bw(ax, background=False):
    """ Plot SN 1998bw radio light curves """
    col = 'lightgrey'
    if background is False:
        col=purp
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
 
    if background is False:
        col=dark
    choose = freq == 4.9 # closest to 3 GHz
    t = dat['dt'][choose]
    flux = dat['flux'][choose] * 1e-3 * 10**(-23)
    d = Planck15.luminosity_distance(z=0.0085).cgs.value
    lum = flux * 4 * np.pi * d**2
    ax.plot(t, lum, c=col, label="_nolegend_")

    if background is False:
        col=yell
    choose = freq == 8.64 # X-band
    t = dat['dt'][choose]
    flux = dat['flux'][choose] * 1e-3 * 10**(-23)
    d = Planck15.luminosity_distance(z=0.0085).cgs.value
    lum = flux * 4 * np.pi * d**2
    ax.plot(t, lum, c=col, label="_nolegend_")

    if background==False:
        ax.text(0.1, 0.1, "SN1998bw", fontsize=12, transform=ax.transAxes)


def plot_09bb(ax, background=False):
    """ Plot SN2009bb radio light curves """
    dat = Table.read(
        ddir + "/2009bb.dat",
        delimiter="&", format='ascii.fixed_width')

    col = 'lightgrey'
    if background is False:
        col = dark
    freq = dat['freq']
    choose = freq== 4.86
    t = dat['dt'][choose]
    flux = dat['flux'][choose] * 1e-3 * 10**(-23) 
    d = Planck15.luminosity_distance(z=0.009937).cgs.value
    lum = flux * 4 * np.pi * d**2 
    ax.plot(t, lum, c=col, label='_nolegend_')

    if background is False:
        col = yell
    freq = dat['freq']
    choose = freq== 8.46
    t = dat['dt'][choose]
    flux = dat['flux'][choose] * 1e-3 * 10**(-23) 
    d = Planck15.luminosity_distance(z=0.009937).cgs.value
    lum = flux * 4 * np.pi * d**2 
    ax.plot(t, lum, c=col, label='_nolegend_')

    if background==False:
        ax.text(0.1, 0.1, "SN2009bb", fontsize=12, transform=ax.transAxes)



def plot_0316d(ax, background=False):
    """ Plot GRB 100316D radio light curves """
    z = 0.0593
    dat = Table.read(ddir + "/100316D_obs.dat", format='ascii.csv')

    col = 'lightgrey'
    if background is False:
        col = dark
    freq = dat['freq']
    choose = freq== 5.4
    t = dat['dt'][choose]
    lum = dat['flux'][choose]
    toplot = dat['fluxerr'][choose] > 0
    ax.plot(t[toplot], lum[toplot], c=col, label="_nolegend_")

    if background is False:
        col = yell
    freq = dat['freq']
    choose = freq== 9
    t = dat['dt'][choose]
    lum = dat['flux'][choose]
    toplot = dat['fluxerr'][choose] > 0
    ax.scatter(
            t[toplot], lum[toplot], c=col, marker='.', label="_nolegend_") 

    if background==False:
        ax.text(0.1, 0.1, "SN2010bh", fontsize=12, transform=ax.transAxes)


def plot_06aj(ax, background=False):
    """ Plot GRB 060218 radio light curves """
    z = 0.03342
    dat = Table.read(ddir + "/060218_obs.dat", format='ascii.csv')

    col = 'lightgrey'
    if background is False:
        col = dark
    freq = dat['freq']
    choose = freq== 4.86
    t = dat['dt'][choose]
    lum = dat['flux'][choose]
    toplot = dat['fluxerr'][choose] > 0
    ax.plot(t[toplot], lum[toplot], c=col, label="_nolegend_")

    if background is False:
        col = yell
    freq = dat['freq']
    choose = freq== 8.46
    t = dat['dt'][choose]
    lum = dat['flux'][choose]
    toplot = dat['fluxerr'][choose] > 0
    ax.plot(t[toplot], lum[toplot], c=col, label="_nolegend_")

    col = 'lightgrey'
    if background is False:
        col = orag
    freq = dat['freq']
    choose = freq== 15
    t = dat['dt'][choose]
    lum = dat['flux'][choose]
    toplot = dat['fluxerr'][choose] > 0
    ax.plot(t[toplot], lum[toplot], c=col, label="_nolegend_")

    if background==False:
        ax.text(0.1, 0.1, "SN2006aj", fontsize=12, transform=ax.transAxes)


def plot_17cw(ax, background=False):
    """ Plot iPTF17cw radio light curves """
    z = 0.093
    col = 'lightgrey'
    if background is False:
        col = dark
    t = np.array([12.6, 15.7, 21.6, 30.7, 41.6])
    flux = np.array([38.1, 30.4, 19.9, 22.4, 19])
    lum = flux
    ax.plot(t, lum, c=col, label='_nolegend_')

    if background is False:
        col = purp
    t = np.array([15.9, 24.9, 41.6])
    flux = np.array([50, 41.4, 44])
    lum = flux
    ax.plot(t, lum, c=col, label='_nolegend_')

    if background is False:
        col = orag
    t = np.array([15.9, 24.9, 65.6, 105])
    flux = np.array([25.4, 21.1, 20.2, 10.7])
    lum = flux
    ax.plot(t, lum, c=col, label='_nolegend_')

    if background==False:
        ax.text(0.1, 0.1, "iPTF17cw", fontsize=12, transform=ax.transAxes)


def plot_12ap(ax, background=False):
    """ Plot SN2012 radio light curves. I'm using the peaks reported
    in Chakraborti+ 2012, because they don't give the actual fluxes """
    z = 0.009
    col = 'lightgrey'
    if background is False:
        col = yell
    t = np.array([12, 18.2])
    flux = np.array([5.85, 5.69]) * 1E3
    lum = flux
    ax.plot(t, lum, c=col, label='_nolegend_')

    if background is False:
        col = dark
    t = np.array([27, 37.9])
    flux = np.array([4.71, 4.20]) * 1E3
    lum = flux
    ax.plot(t, lum, c=col, label='_nolegend_')

    if background==False:
        ax.text(0.1, 0.1, "SN2012ap", fontsize=12, transform=ax.transAxes)


if __name__=="__main__":
    fig,axarr = plt.subplots(
            3, 2, figsize=(6,6), sharex=True, sharey=True)

    for ii,ax in enumerate(axarr.reshape(-1)):
        # LC of ZTF18aaqjovh
        plot_source(ax)
        ax.yaxis.set_tick_params(labelsize=14)
        ax.xaxis.set_tick_params(labelsize=14)
        ax.set_yscale('log')
        ax.set_xlim(0, 40)
        ax.set_ylim(1E39, 1E45)

    fig.text(0.5, 0.04, r"$\Delta t$ (days)", ha='center', fontsize=16) 
    fig.text(
            0.04, 0.5, r'0.3--10 keV X-ray Luminosity (erg/s)', 
            fontsize=16, rotation='vertical', horizontalalignment='center',
            verticalalignment='center')

    fig.subplots_adjust(wspace=0.1, hspace=0.1)

    plt.show()
    #plt.savefig(
    #    "xray_lc.eps", format='eps', dpi=500, bbox_inches='tight',
    #    pad_inches=0.1)
