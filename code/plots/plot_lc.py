""" Plot the light curve of ZTF18aaqjovh

Light curve was retrieved in ZTF_Tools/query_lc.py,
and saved to a text file lc.dat
"""

import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.io import ascii
from astropy.cosmology import Planck15
#from ztfquery import marshal
#import extinction


# time of the last non-detection
t0 = 58233.17615 # in MJD
# redshift
z = 0.05403


def plot_98bw(ax):
    datadir = "/Users/annaho/Dropbox/Projects/Research/IcBL/data/optical_compilations"
    dat = ascii.read(datadir + "/sn1998bw.dat", delimiter=';')
    jd = dat['JD']
    rband = dat['Rcmag']
    erband = dat['e_Rcmag']
    dm = Planck15.distmod(z=z).value-Planck15.distmod(z=0.0085).value
    ax.plot(jd-jd[0], rband+dm, color='k', lw=0.5, label="98bw")
    ax.plot(
            (jd-jd[0])/(1.0085), (rband+dm)+0.5, color='k', 
            lw=0.5, ls='--', label="98bw+0.5 mag")
    print(jd[0])
    # ax.fill_between(
    #         jd-jd[0], rband+dm-erband, 
    #         rband+dm+erband, color='lightgrey')


#def download_lc():
#     """ Download the light curve """
#     # connect to databases
#     m = marshal.MarshalAccess()
#     print("Connected")
# 
#     # download light curves
#     marshal.download_lightcurve('ZTF18aaqjovh')


def load_marshal_lc():
    f = "Data/marshal/lightcurves/ZTF18aaqjovh/marshal_plot_lc_lightcurve_ZTF18aaqjovh.csv"
    lc = ascii.read(f)
    det = lc['mag'] < 99
    dt = lc['jdobs'][det]-lc['jdobs'][det][0]
    mag = lc['mag'][det]
    emag = lc['emag'][det]
    dt_nondet = lc['jdobs'][~det]-lc['jdobs'][det][0]
    lm = lc['limmag'][~det]


def load_danny_lc(ax):
    data_dir = "/Users/annaho/Dropbox/Projects/Research/ZTF18aaqjovh/data"
    f = data_dir + "/ZTF18aaqjovh.csv"
    lc = ascii.read(f)
    det = lc['mag'] < 99
    dt = lc['mjd'][det]-t0
    mag = lc['mag'][det]
    emag = lc['magerr'][det]
    dt_nondet = lc['mjd'][~det]-t0
    lm = lc['lim_mag'][~det]

    # Plot r-band light curve
    choose = lc['filter'][det] == 'ztfr'
    ax.errorbar(
            dt[choose], mag[choose], yerr=emag[choose], 
            fmt='o', c='k', label="ZTF18aaqjovh")

    # Plot r-band upper limits
    choose = lc['filter'][~det] == 'ztfr'
    ax.errorbar(
            (dt_nondet[choose])/1.05403, lm[choose], fmt='v', c='k', label=None)

    # Show some spectral epochs with an S
    sp = [14, 19, 20, 44, 105]
    for s in sp:
        ax.text(s, 17.7, "S", fontsize=12)



if __name__=="__main__":
    fig,ax = plt.subplots(1,1,figsize=(9, 7))

    load_danny_lc(ax)
    plot_98bw(ax)

    ax.set_xlabel("Days (Rest-frame) Since Last Non-Detection of ZTF18aaqjovh", fontsize=16)
    ax.set_ylabel("Apparent Mag (AB)", fontsize=16)
    ax.set_ylim(17.5,21)
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

    ax.legend(loc='lower right', fontsize=14)
    #plt.show()
    plt.savefig('lc.png', dpi=100)
