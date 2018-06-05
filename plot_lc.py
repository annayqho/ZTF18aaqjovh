""" Plot the light curve of ZTF18aaqjovh

Light curve was retrieved in ZTF_Tools/query_lc.py,
and saved to a text file lc.dat
"""

import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
from astropy.io import ascii


fig = plt.figure(figsize=(7, 4))
lc = ascii.read('lc.dat', delimiter=" ")
plt.errorbar(lc['t']-lc['t'][0], lc['mag'], yerr=lc['emag'], fmt='.', c='k')
plt.xlabel("Days since first detection", fontsize=16)
plt.ylabel("ZTF $r$-mag", fontsize=16)
plt.ylim(17.9,18.8)
plt.gca().invert_yaxis()
plt.tick_params(axis='both', labelsize=14)
plt.tight_layout()
plt.savefig('lc.png')
