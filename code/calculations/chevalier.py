""" perform the Chevalier calculations """

import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18aaqjovh/code/plots")
from astropy.cosmology import Planck15
from radio_lc import ujy_to_flux

eps_B = 1/3
eps_e = 1/3
f = 1/2

def v(x, y):
    return 10**((9/19) * (np.log10(y) - 26 - (19/9)*np.log10(x)))

def mdot(x,y):
    return 10**((4/19) * 
            ((26-np.log10(y)) + 19/4 * np.log10(0.0005/eps_B) + \
            (2*19/4)*np.log10(x)))


def get_B(Fp, nup, d_mpc):
    """ This is Equation 14 in C98
    Fp in Jy
    nup in GHz
    
    Returns Bp in Gauss
    """
    alpha = eps_e / eps_B
    Bp = 0.58 * alpha**(-4/19) * (f/0.5)**(-4/19) * (Fp)**(-2/19) * \
         d_mpc**(-4/19) * (nup/5)
    return Bp

# The highest velocity will come from the minimum value of x
# together with the maximum value of y

x = 20*(3/5)
y = ujy_to_flux(30, 0.05403)
print(v(x,y))

# The lowest velocity will come from the highest value of x
# together with the minimum value of y

x = 20*(15/5)
y = ujy_to_flux(20, 0.05403)
print(v(x,y))

# The max Mdot will come from the max value of x
# together with the max value of y
x = 20*(15/5)
y = ujy_to_flux(30, 0.05403)
print(mdot(x,y))

# The min Mdot will come from the min value of x
# together with the min value of y
x = 20*(3/5)
y = ujy_to_flux(20, 0.05403)
print(mdot(x,y))

# The max B is when Fp is min and nu_p is max
x = 20*(15/5)
y = ujy_to_flux(20, 0.05403)
dmpc = Planck15.luminosity_distance(z=0.05402).value
print(get_B(y, x, dmpc))
 
# The min B is when Fp is max and nu_p is min
x = 20*(3/5)
y = ujy_to_flux(30, 0.05403)
dmpc = Planck15.luminosity_distance(z=0.05402).value
print(get_B(y, x, dmpc))
