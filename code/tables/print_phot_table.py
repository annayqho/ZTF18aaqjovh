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

data_dir = "/Users/annaho/Dropbox/Projects/Research/ZTF18aaqjovh/code/forced_phot"
t0 = 58233.17615 # time of the last non-detection, in MJD

# g-band P48 photometry
f = data_dir + "/ZTF18aaqjovh_force_phot_lc_g.txt"
lc = np.loadtxt(f)
jd = lc[:,0]
mag = lc[:,1]
det = ~np.isnan(mag) 
mjd = jd[det]-2400000.5
mag = mag[det]
emag = lc[:,2][det]
ndet = sum(det)
tel = np.array(['P48']*ndet)
filt = np.array(['g']*ndet)

# add r-band P48 photometry
f = data_dir + "/ZTF18aaqjovh_force_phot_lc_r.txt"
lc = np.loadtxt(f)
jd = lc[:,0]
mag_r = lc[:,1]
det = ~np.isnan(mag_r) 
ndet = sum(det)
mjd = np.hstack((mjd, jd[det]-2400000.5))
mag = np.hstack((mag, mag_r[det]))
emag = np.hstack((emag, lc[:,2][det]))
tel = np.hstack((tel, np.array(['P48']*ndet)))
filt = np.hstack((filt, np.array(['r']*ndet)))

# add P60 photometry
mjd = np.hstack((mjd, [2458247.8588-2400000.5, 2458248.8353-2400000.5]))
tel = np.hstack((tel, ['P60', 'P60']))
mag = np.hstack((mag, [18.30, 18.21]))
emag = np.hstack((emag, [0.04, 0.03]))
filt = np.hstack((filt, np.array(['r']*ndet)))

# generate dt
dt = mjd-t0

# sort
order = np.argsort(dt)
mjd = mjd[order]
dt = dt[order]
mag = mag[order]
emag = emag[order]
tel = tel[order]
filt = filt[order]

# remove duplicate entries
mjd,ind = np.unique(mjd, return_index=True)
dt = dt[ind]
mag = mag[ind]
emag = emag[ind]
tel = tel[ind]
filt = filt[ind]

for ii in np.arange(len(dt)):
    if(dt[ii] < 55):
        mjd_str = '{:<08f}'.format(round_sig(mjd[ii], 11)) # pad with zeros
        dt_str = '{:.2f}'.format(np.round(dt[ii], 2))
        mag_str = '{:.2f}'.format(round_sig(mag[ii], 4))
        emag_str = '{:.2f}'.format(np.round(emag[ii], ndec(mag_str)))
        row = rowstr %(
                mjd_str, dt_str, tel[ii], filt[ii], 
                mag_str, emag_str)
        outputf.write(row)

outputf.write("\enddata \n")
outputf.write("\end{deluxetable*} \n")
