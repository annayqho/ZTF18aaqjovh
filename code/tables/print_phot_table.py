""" Print table of photometry, ground-based and UVOT """

import numpy as np
from math import floor, log10
from astropy.time import Time
from astropy.cosmology import Planck15
from astropy.io import ascii
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from load_lc import get_uv_lc


def round_sig(x, sig=2):
    print(x)
    if x < 0:
        return -round(-x, sig-int(floor(log10(-x)))-1)
    return round(x, sig-int(floor(log10(x)))-1)


def ndec(num):
    dec = str(num).split('.')[-1]
    return len(dec)


d = Planck15.luminosity_distance(z=0.05403).cgs.value

headings = np.array(
        ['Date (JD)', '$\Delta t$', 'Instrument', 'Filter', 
         'AB Mag', 'Error in AB Mag'])
label = "opt-phot"
caption = "Optical photometry for ZTF18aaqjovh"

# Print the table headers
ncol = len(headings)
colstr = ""
colstr += 'l'
for col in np.arange(ncol-1): colstr+="r"
print(colstr)

colheadstr = ""
for col in np.arange(ncol-1):
    colheadstr += "\colhead{%s} & " %headings[col]
colheadstr += "\colhead{%s}" %headings[-1]

rowstr = ""
for col in np.arange(ncol-1):
    rowstr += "%s & "
rowstr += "%s \\\ \n"

outputf = open("table_%s.txt" %label, "w")
outputf.write("\\startlongtable \n")
outputf.write("\\begin{deluxetable*}{%s} \n" %colstr)
outputf.write("\\tablecaption{%s\label{tab:%s}} \n" %(caption,label))
outputf.write("\\tablewidth{0pt} \n")
outputf.write("\\tablehead{ %s } \n" %colheadstr)
#outputf.write("\\rotate \n")
outputf.write("\\tabletypesize{\scriptsize} \n")
outputf.write("\startdata \n")

data_dir = "/Users/annaho/Dropbox/Projects/Research/ZTF18aaqjovh/data"
t0 = 58233.17615 # time of the last non-detection, in MJD
f = data_dir + "/ZTF18aaqjovh.csv"
lc = ascii.read(f)
det = lc['mag'] < 99
mjd = lc['mjd'][det]
tel = np.array(['P48']*len(mjd))
mag = lc['mag'][det]
emag = lc['magerr'][det]

# add P60 photometry
mjd = np.hstack((mjd, [2458247.8588-2400000.5, 2458248.8353-2400000.5]))
tel = np.hstack((tel, ['P60', 'P60']))
mag = np.hstack((mag, [18.30, 18.21]))
emag = np.hstack((emag, [0.04, 0.03]))

# generate dt
dt = mjd-t0

# sort
order = np.argsort(dt)
mjd = mjd[order]
dt = dt[order]
mag = mag[order]
emag = emag[order]
tel = tel[order]

for ii in np.arange(len(dt)):
    mjd_str = round_sig(mjd[ii], 11)
    dt_str = np.round(dt[ii], 2)
    mag_str = str(round_sig(mag[ii], 4)).zfill(5)
    emag_str = str(np.round(emag[ii], ndec(mag_str))).zfill(4)
    row = rowstr %(
            mjd_str, dt_str, tel[ii], r"$r$", 
            mag_str, emag_str)
    outputf.write(row)

outputf.write("\enddata \n")
outputf.write("\end{deluxetable*} \n")
