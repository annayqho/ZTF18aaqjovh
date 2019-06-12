""" Plot the spectral sequence of ZTF18aaqjovh
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.io import ascii
from ztfquery import query
from ztfquery import marshal
import extinction
import glob
from astropy.time import Time
sys.path.append("/Users/annaho/Github/Spectra")
from normalize import smooth_spec
from measure_snr import get_snr


z = 0.03154
SPEC_DIR = "Data/marshal/spectra/ZTF18aaqjovh"


def get_res(tel):
    """ Here, this means the width of a line in Angstroms """
    if tel == 'LT':
        res = 18 # Angstrom, res at central wavelength
        res = 30 # add a couple of Ang?
    elif tel == 'P200':
        res = 10 # determined by eye from the spectrum
        # basically, width of a galaxy emission line is 10 AA
        # and each pixel is 1 AA
    elif tel == 'Keck1':
        res = 7*2 # determined by eye from spectrum
        # width of a line is around 7 pixels
        # and each pixel is 2 Angstroms
    elif tel == 'NOT':
        # width of a line is around 8 pixels
        # and each pixel is around 2.63 Ang
        res = 8*2.63
    elif tel == 'DCT':
        # width of a line is around 7 pixels
        # and each pixel is 2.2 Ang
        res = 7*2.2
    elif tel == 'P60':
        res = 20
    else:   
        print("I don't have this telescope")
    return res


def download_spec():
    """ Download the spectra """
    # connect to databases
    m = marshal.MarshalAccess()
    print("Connected")

    # download light curves
    marshal.download_spectra('ZTF18aaqjovh')

    # return filenames
    ddir = "Data/marshal/spectra/ZTF18aaqjovh"
    f = glob.glob(ddir + "/*.ascii")
    print(f)
    return f


def get_files(sind, eind):
    """ 
    start_ind: starting index
    end_ind: end index
    """
    files = np.array(glob.glob(SPEC_DIR + "/*.ascii"))
    dt = np.zeros(len(files))
    tels = []
    cols = np.array([""]*len(dt), dtype='U10')

    # Read in all of the files, pull out the corresponding dates, and sort by date
    t0 = 58233.17615 # in MJD
    for ii,f in enumerate(files):
        tel = f.split("_")[2]
        tels.append(tel)
        alldat = open(f).readlines()
        if tel == 'LT': 
            for line in alldat:
                if 'DATE-OBS' in line:
                    obsdate = line[13:36]
                    t = Time(obsdate, format='isot').mjd
                    dt[ii] = t-t0
            cols[ii] = 'magenta'
        elif tel == 'P200':
            for line in alldat:
                if 'UTSHUT' in line:
                    obsdate = line[11:]
                    t = Time(obsdate, format='isot').mjd
                    dt[ii] = t-t0
            cols[ii] = 'lightblue'
        elif tel == 'Keck1':
            for line in alldat:
                if 'DATE_BEG' in line:
                    obsdate = line[13:32]
                    t = Time(obsdate, format='isot').mjd
                    dt[ii] = t-t0
            cols[ii] = 'red'
        elif tel == 'DCT':
            obsdate = '2018-09-14T00:00:00' # temporary
            t = Time(obsdate, format='isot').mjd
            dt[ii] = t-t0
            cols[ii] = 'yellow'
        elif tel == 'NOT':
            obsdate = '2018-05-15T00:00:00' # temporary
            t = Time(obsdate, format='isot').mjd
            dt[ii] = t-t0
            cols[ii] = 'green'
        elif tel == 'P60':
            for line in alldat:
                if 'OBSUTC' in line:
                    temp = line[10:].split(" ")
                    obsdate = temp[0] + "T" + temp[1]
                    t = Time(obsdate, format='isot').mjd
                    dt[ii] = t-t0
            cols[ii] = 'black'
        else:
            print("couldn't find telescope")
            print(tel)
    order = np.argsort(dt)
    files_sorted = files[order]
    dt_sorted = dt[order]
    tel_sorted = np.array(tels)[order]
    cols = cols[order]
    return files_sorted[sind:eind], dt_sorted[sind:eind], tel_sorted[sind:eind]


def load_spec(f):
    lc = np.loadtxt(f)
    wl = lc[:,0]
    f = lc[:,1]
    if tel == 'Keck':
        eflux = lc[:,3]
    else:
        # estimate uncertainty from scatter in a continuum region
        eflux = np.array([get_snr(wl, f, 6300, 6500)]*len(wl))
    ivar = 1/eflux**2
    return wl, f, ivar


def plot_smoothed_spec(ax, x, y, ivar, tel, epoch, ls='-', lw=0.5, c='black', label=None, text=True):
    """ plot the smoothed spectrum """
    res = get_res(tel)
    smoothed = smooth_spec(x, y, ivar, res*3)
    ax.plot(
            x, smoothed, c=c,
            drawstyle='steps-mid', lw=lw, ls=ls, alpha=1.0, label=label,
            zorder=10)
    dt_str = r"+%s\,d" %(
            str(np.round(epoch, 1)))
    if text:
        ax.text(
                x[-1]+100, smoothed[-1],  s=dt_str,
                horizontalalignment='left', verticalalignment='center',
                fontsize=12)
    return smoothed


def plot_spec(ax, x, y, tel, epoch):
    """ plot the spectrum """
    ax.plot(
            x, y, c='lightgrey',
            drawstyle='steps-mid', lw=0.5, alpha=0.4)
    return ax


if __name__=="__main__":
    fig,ax = plt.subplots()
    files, epochs, tels = get_files(0, 6)
    nfiles = len(files)
    shift = [0, 1, 1.5, 2.5, 3.5, 5, 6, 7, 8]
    for ii,f in enumerate(files):
        tel = tels[ii]
        dt = epochs[ii]
        wl, flux, ivar = load_spec(f)
        choose = np.logical_and(wl>3660, wl < 9000)
        scale = flux[wl>4100][0]
        shifted = flux/scale-shift[ii]
        plot_spec(ax, wl[choose], (shifted-shift[ii])[choose], tel, dt)
        plot_smoothed_spec(
                ax, wl[choose], (shifted-shift[ii])[choose], 
                ivar[choose], tel, dt)

    plt.tick_params(axis='both', labelsize=14)
    plt.xlim(3660, 10140)
    plt.ylim(-10, 0)
    plt.xlabel(r"Observed Wavelength (\AA)", fontsize=16)
    plt.ylabel(r"Scaled $F_{\lambda}$ + const.", fontsize=16)
    ax.get_yaxis().set_ticks([])
    plt.tight_layout()
    #plt.show()
    plt.savefig("spec_sequence.eps", format='eps', dpi=500)
