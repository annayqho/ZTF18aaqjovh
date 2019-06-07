""" Plot the light curve of ZTF18aaqjovh

Light curve was retrieved in ZTF_Tools/query_lc.py,
and saved to a text file lc.dat
"""

import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.io import ascii
from ztfquery import query
from ztfquery import marshal
import extinction


# time of the last non-detection
t0 = 58233.17615 # in MJD
# redshift
z = 0.05403


def download_lc():
    """ Download the light curve """
    # connect to databases
    m = marshal.MarshalAccess()
    print("Connected")

    # download light curves
    marshal.download_lightcurve('ZTF18aaqjovh')


def load_marshal_lc():
    f = "Data/marshal/lightcurves/ZTF18aaqjovh/marshal_plot_lc_lightcurve_ZTF18aaqjovh.csv"
    lc = ascii.read(f)
    det = lc['mag'] < 99
    dt = lc['jdobs'][det]-lc['jdobs'][det][0]
    mag = lc['mag'][det]
    emag = lc['emag'][det]
    dt_nondet = lc['jdobs'][~det]-lc['jdobs'][det][0]
    lm = lc['limmag'][~det]


def load_danny_lc():
    data_dir = "/Users/annaho/Dropbox/Projects/Research/ZTF18aaqjovh/data"
    f = data_dir + "/ZTF18aaqjovh.csv"
    lc = ascii.read(f)
    det = lc['mag'] < 99
    dt = lc['mjd'][det]-t0
    mag = lc['mag'][det]
    emag = lc['magerr'][det]
    dt_nondet = lc['mjd'][~det]-t0
    lm = lc['lim_mag'][~det]
    fig,ax = plt.subplots(figsize=(7, 4))

    # Plot r-band light curve
    choose = lc['filter'][det] == 'ztfr'
    plt.errorbar(
            dt[choose], mag[choose], yerr=emag[choose], 
            fmt='o', c='k', label="$r$")

    # Plot r-band upper limits
    choose = lc['filter'][~det] == 'ztfr'
    plt.errorbar(
            dt_nondet[choose], lm[choose], fmt='v', c='k', label=None)

    # Show some spectral epochs with an S
    sp = [14, 19, 20, 44, 105]
    for s in sp:
        ax.text(s, 18.1, "S", fontsize=12)

    plt.xlabel("Days Since Last Non-Detection", fontsize=16)
    plt.ylabel("Apparent Mag (AB)", fontsize=16)
    plt.ylim(17.9,21)
    plt.xlim(-1,52)
    plt.gca().invert_yaxis()
    plt.tick_params(axis='both', labelsize=14)
    # fig.text(0.94, 0.5, "Absolute Magnitude (AB)",
    #          ha='center', va='center', fontsize=16, rotation='vertical')
    plt.tight_layout()
    plt.legend(loc='lower right', fontsize=14)
    #plt.show()
    plt.savefig('lc.eps', format='eps', dpi=1000)


if __name__=="__main__":
    load_danny_lc()
