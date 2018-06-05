""" Plot the LRIS May 14 spectrum of the source """

from astropy.io import ascii
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")


dat = ascii.read("data/ZTF18aaqjovh_20180514_Keck1_v1.ascii")

fig = plt.figure(figsize=(8,4))
plt.step(dat['col1'], dat['col2']/10**-16, c='k', where='mid', lw=0.5)
plt.ylabel("Flux [$\\times 10^{-16}$ erg/s/cm${}^2/\AA$]", fontsize=16)
plt.xlabel("Wavelength $\lambda (\AA)$", fontsize=16)
plt.tick_params(axis='both', labelsize=14)

plt.tight_layout()
plt.savefig("spec.png")
