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
    f = "/Users/annaho/Dropbox/Projects/Research/ZTF18aaqjovh/data/ZTF18aaqjovh.csv"

    fig = plt.figure(figsize=(7, 4))

    # Plot r-band light curve
    choose = lc['filter'][det] == 'r'
    plt.errorbar(
            dt[choose], mag[choose], yerr=emag[choose], 
            fmt='o', c='k', label="$r$")

    # Plot r-band upper limits
    choose = lc['filter'][~det] == 'r'
    plt.errorbar(
            dt_nondet[choose], lm[choose], 
            fmt='v', c='k', label="$r$")


    plt.xlabel("Days Since First ZTF Detection", fontsize=16)
    plt.ylabel("Apparent Mag (AB)", fontsize=16)
    plt.ylim(17.9,21)
    plt.xlim(-12,63)
    plt.gca().invert_yaxis()
    plt.tick_params(axis='both', labelsize=14)
    plt.tight_layout()
    #plt.show()
    plt.savefig('lc.eps', format='eps', dpi=1000)
