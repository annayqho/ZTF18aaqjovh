""" Plot the spectral sequence of ZTF18aaqjovh
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.io import ascii
from ztfquery import query
from ztfquery import marshal
import extinction
import glob
from astropy.time import Time


z = 0.03154
SPEC_DIR = "Data/marshal/spectra/ZTF18aaqjovh"


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
            obsdate = '2018-09-17T00:00:00' # temporary
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
    return wl, f


if __name__=="__main__":
    fig,ax = plt.subplots()
    files, epochs, tels = get_files(0, 6)
    nfiles = len(files)
    shift = [0, 1, 1.5, 2.5, 3.5, 5, 6, 7, 8]
    for ii,f in enumerate(files):
        tel = tels[ii]
        dt = epochs[ii]
        wl, flux = load_spec(f)
        scale = flux[wl>4100][0]
        shifted = flux/scale-shift[ii]
        plt.plot(
            wl, shifted-shift[ii], c='lightgrey', 
            drawstyle='steps-mid', lw=0.4, alpha=1.0)

    plt.tick_params(axis='both', labelsize=14)
    plt.xlim(3600, 9500)
    plt.ylim(-10, 0)
    plt.xlabel(r"Observed Wavelength (\AA)", fontsize=16)
    plt.ylabel(r"Scaled $F_{\lambda}$ + const.", fontsize=16)
    plt.show()
