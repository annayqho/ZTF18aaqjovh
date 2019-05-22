""" Align images using SWARP """

import os

swarp_path = '/usr/local/optical/swarp/bin/swarp'

gfits = "../data/host_sdss/frame-g-005116-5-0061.fits"
rfits = "../data/host_sdss/frame-r-005116-5-0061.fits"
ifits = "../data/host_sdss/frame-i-005116-5-0061.fits"
zfits = "../data/host_sdss/frame-z-005116-5-0061.fits"
ufits = "../data/host_sdss/frame-u-005116-5-0061.fits"

# centering
ra = 178.181754 
dec = 25.675033
size = 600

swarp_command = swarp_path \
        + " %s %s %s %s %s -c config.swarp -CENTER '" %(gfits,rfits,ifits,zfits,ufits)\
        + str(ra) + " " + str(dec) \
        + "' -SUBTRACT_BACK Y -RESAMPLE Y -COMBINE N -IMAGE_SIZE '" \
        + str(size) + "," + str(size) + "'"
os.system(swarp_command)
