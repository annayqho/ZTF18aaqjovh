""" Plot the light curve of ZTF18aaqjovh
This is Fig 2 of our paper.

Light curve was retrieved in ZTF_Tools/query_lc.py,
and saved to a text file lc.dat
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.io import ascii
from astropy.cosmology import Planck15
#from ztfquery import marshal
#import extinction


data_dir = "/Users/annaho/Dropbox/Projects/Research/ZTF18aaqjovh/data/phot"
# time of the last non-detection
t0 = 58233.17615 # in MJD
# redshift
z = 0.05403


def plot_98bw(ax):
    offset = 0.5
    dat = ascii.read(data_dir + "/sn1998bw.dat", delimiter=';')
    jd = dat['JD']
    rband = dat['Rcmag']
    erband = dat['e_Rcmag']
    # Extinction is 0.127 in R-band in this direction
    dm = Planck15.distmod(z=z).value-Planck15.distmod(z=0.0085).value
    ax.plot(jd-jd[0], rband+dm-0.127, color='#e55c30', lw=1.5, 
            label="98bw $Rc$")
    ax.plot(
            (jd-jd[0])/(1.0085), (rband+dm)+offset-0.127, color='#e55c30', 
            lw=0.5, ls='--', label="98bw $Rc$+%s mag" %offset)
    gband = dat['Bmag']
    egband = dat['e_Bmag']
    # Extinction is 0.212 in B-band in this direction
    ax.plot(jd-jd[0], gband+dm-0.212, color='#140b34', lw=1.5, 
            label="98bw $B$")
    ax.plot(
            (jd-jd[0])/(1.0085), (gband+dm)+offset-0.212, color='#140b34', 
            lw=0.5, ls='--', label="98bw $B$+%s mag" %offset)
    ax.axvline(x=-0.1, c='k', lw=0.5)
    ax.text(0,18.7,'GRB 980425', fontsize=10, rotation=270)
    # ax.fill_between(
    #         jd-jd[0], rband+dm-erband, 
    #         rband+dm+erband, color='lightgrey')



def load_lc(ax):

    # Plot r-band light curve
    f = data_dir + "/ZTF18aaqjovh_force_phot_lc_r.txt"
    lc = np.loadtxt(f)
    dt = lc[:,0]-2400000.5-t0
    mag = lc[:,1]
    emag = lc[:,2]
    choose = emag < 1
    ax.errorbar(
            dt[choose], mag[choose], yerr=emag[choose], 
            fmt='o', c='#e55c30', label="ZTF18aaqjovh, $r$")

    # Plot g-band light curve
    f = data_dir + "/ZTF18aaqjovh_force_phot_lc_g.txt"
    lc = np.loadtxt(f)
    dt = lc[:,0]-2400000.5-t0
    mag = lc[:,1]
    emag = lc[:,2]
    choose = emag < 1
    ax.errorbar(
            dt[choose], mag[choose], yerr=emag[choose], 
            fmt='s', c='#140b34', label="ZTF18aaqjovh, $g$")

    # Show some spectral epochs with an S
    sp = [14, 19, 20, 44]
    for s in sp:
        ax.text(s, 17.7, "S", fontsize=14)

    # # Plot the last non-detection
    ax.scatter(0, 20.6, marker='v', c='k')


if __name__=="__main__":
    fig,ax = plt.subplots(1,1,figsize=(7, 6))

    load_lc(ax)
    plot_98bw(ax)


    ax.set_xlabel(
            "Days (Rest-frame) Since Last Non-Detection of ZTF18aaqjovh", 
            fontsize=16)
    ax.set_ylabel("Apparent Mag (AB)", fontsize=16)
    ax.set_ylim(17.5,21.5)
    ax.set_xlim(-1,52)
    ax.tick_params(axis='both', labelsize=14)
    ax.invert_yaxis()
    #ax2.invert_yaxis()

    # Put a twin axis with the absolute magnitude
    ax2 = ax.twinx()
    y_f = lambda y_i: y_i-Planck15.distmod(z=z).value
    ymin, ymax = ax.get_ylim()
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.plot([],[])
    ax2.tick_params(axis='both', labelsize=14)
    ax2.set_ylabel(
            "Absolute Mag (AB)", fontsize=16, rotation=270, va='bottom')
    ax2.set_xlim(-1,52)

    ax.legend(loc='lower center', fontsize=12, ncol=2)
    fig.tight_layout()
    #plt.show()
    plt.savefig('lc.eps', dpi=300)
