""" Get forced photometry for ZTF18aaqjovh,
taking away the contribution of SN light in the reference. """

import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF_fast_transient_search/code/forced_phot")
from run_forced_phot import *
from ztfquery import query
import glob
import multiprocessing

# Run forced photometry
name = 'ZTF18aaqjovh'
ra = 178.181754
dec = 25.675033
jdobs = Time("2018-04-25", format='isot').jd
filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs)
choose_r = filt == 'r'
choose_g = filt == 'g'

# Now, calculate the offsets
choose = np.logical_and(filt=='g', jd>2.45855E6)
mean,wsum = np.average(flux[choose], weights=1/eflux[choose]**2, returned=True)
efmean = np.sqrt(1/wsum)
print(max(jd[choose]))
# mean = -270.81 +/- 0.88 for g-band

choose = np.logical_and(filt=='r', jd>2.45855E6)
print(max(jd[choose]))
mean,wsum = np.average(flux[choose], weights=1/eflux[choose]**2, returned=True)
efmean = np.sqrt(1/wsum)
# mean = 1.67 +/- 1.00 for r-band

# Apply the offsets to create an r-band light curve.
# (In the forced phot code, I shifted everything to a ZP of 25)
mag_g = -2.5*np.log10(flux[choose_g]+270.81)+25
mag_r = -2.5*np.log10(flux[choose_r]-1.67)+25

# Calculate error bars
fluxerr = np.sqrt(eflux[choose_r]**2 + 1**2)
flux_r = flux[choose_r]-1.67
emag_r =np.abs(-1.0855*fluxerr/flux_r)

fluxerr = np.sqrt(eflux[choose_g]**2 + 0.88**2)
flux_g = flux[choose_g]+270.81
emag_g =np.abs(-1.0855*fluxerr/flux_g)

# It's identical to the previous r-band LC, which is a relief.

# Now, save the photometry and subtract Galactic extinction
ext_g = 0.053
ext_r = 0.037
np.savetxt(
    "ZTF18aaqjovh_force_phot_lc_r.txt", 
    np.array([jd[choose_r],mag_r-0.037,emag_r]).T)
np.savetxt(
    "ZTF18aaqjovh_force_phot_lc_g.txt", 
    np.array([jd[choose_g],mag_g-0.053,emag_g]).T)
