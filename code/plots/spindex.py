""" Constrain the spectral index from the radio to the X-ray observations """

from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import Planck15

bigsize=16
smallsize=14

d = Planck15.luminosity_distance(z=0.05403).cgs.value

fig = plt.figure(figsize=(7,4.7))

# May 28 / Day 33: Chandra limit < 2.8E-15 erg/s/cm2
# May 29: VLA limit 26.6 \pm 5.4 at 6 GHz

# Xray: 2.4E18 Hz is 10 keV, 7.2E16 Hz is 0.3 keV
xray_nu = np.sqrt(2.4E18*7.2E16) # geometric mean
xray_nulnu = 2.8E-15*4*np.pi*d**2
plt.scatter(
        xray_nu, xray_nulnu, marker='o', c='k')
plt.arrow(
        xray_nu, xray_nulnu, 0, -xray_nulnu/1.2, length_includes_head=True, 
        head_width=xray_nu/2, head_length=xray_nulnu/5, fc='k')

# Radio
radio_nu = 6E9
radio_nulnu = radio_nu * 26.2E-6 * 1E-23 * 4 * np.pi * d**2
enulnu = radio_nu * 5.4E-6 * 1E-23 * 4 * np.pi * d**2 
plt.errorbar(
        radio_nu, radio_nulnu, yerr=enulnu, fmt='o', c='k', label="Day 33",
        zorder=1)

# Calculate spectral index
alpha = np.log10(xray_nulnu/radio_nulnu) / np.log10(xray_nu/radio_nu)
plt.plot([radio_nu, xray_nu], [radio_nulnu, xray_nulnu], c='k', ls='--', lw=0.5)
plt.text(
    1E13, 1E39, 
    r'$\nu L_\nu \propto \nu^{%s}$' %np.round(alpha,2), fontsize=smallsize)

# We have no comparable radio measurements
# July 24 / Day 90: Chandra limit < 2.9E-15 erg/s/cm2

# Instead maybe plot the data from May 17, which is Day 22
nus = np.array([3E9, 6E9, 15E9])
fluxes = np.array([26, 29.6, 15.1]) * 1E-6 * 1E-23 * 4 * np.pi * d**2
efluxes = np.array([6.9, 5.3, 5.2]) * 1E-6 * 1E-23 * 4 * np.pi * d**2
plt.errorbar(
        nus, nus*fluxes, yerr=nus*efluxes, 
        c='#e55c30', marker='s', label="Day 22", zorder=0)

plt.xscale('log')
plt.yscale('log')
plt.xlabel("Frequency [Hz]", fontsize=bigsize)
plt.ylabel(r"$\nu L_\nu$ (erg/s)", fontsize=bigsize)
plt.xticks(fontsize=bigsize)
plt.yticks(fontsize=bigsize)
plt.legend(fontsize=bigsize, loc='upper left')
plt.xlim(1E9, 1E18)
plt.ylim(1E36, 1E41)
plt.tight_layout()
plt.savefig("sed.eps", format='eps', dpi=500)
plt.close()
#plt.show()
