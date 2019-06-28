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


def ujy_to_flux(ujy, z):
    d = Planck15.luminosity_distance(z=z).cgs.value
    return ujy*1E-6*1E-23*4*np.pi*d**2


def plot_source(ax):
    lw=3
    # 6 GHz LC
    ax.plot(
            [16, 22, 34], 
            ujy_to_flux(np.array([32.5, 29.6, 26.6]), 0.05403)/1E27, 
            c=dark, lw=lw, label="6 GHz")
    ax.errorbar(
            [16, 22, 34], 
            ujy_to_flux(np.array([32.5, 29.6, 26.6]), 0.05403)/1E27, 
            yerr=ujy_to_flux(np.array([7.1, 5.3, 5.4]), 0.05403)/1E27,
            c=dark, fmt='o')

    # 3 GHz LC
    ax.plot(
            [21, 36], 
            ujy_to_flux(np.array([26, 34.6]), 0.05403)/1E27, 
            c=purp, lw=lw, label="3 GHz")
    ax.errorbar(
            [21, 36], 
            ujy_to_flux(np.array([26, 34.6]), 0.05403)/1E27, 
            yerr=ujy_to_flux(np.array([6.9, 4.8]), 0.05403)/1E27,
            c=purp, fmt='s')

    # 15 GHz LC
    ax.errorbar(21, ujy_to_flux(15.2, 0.05403)/1E27, 
            yerr=ujy_to_flux(5.2, 0.05403)/1E27,
            c=orag, fmt='D')


fig,axarr = plt.subplots(4, 3, figsize=(8,8), sharex=True, sharey=True)

for ii,ax in enumerate(axarr.reshape(-1)):
    # LC of ZTF18aaqjovh
    plot_source(ax)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.xaxis.set_tick_params(labelsize=14)
    #ax.set_yscale('log')
    ax.set_xlim(0,50)
    ax.set_ylim(0, 5)


fig.text(0.5, 0.04, r"$\Delta t$ (days)", ha='center', fontsize=16) 
fig.text(
        0.04, 0.5, r'Radio Luminosity ($10^{27}$ erg/s)', 
        fontsize=16, rotation='vertical', horizontalalignment='center',
        verticalalignment='center')

plt.subplots_adjust(hspace=0, wspace=0)

plt.show()
