""" Align images using SWARP """

import os

swarp_path = '/usr/local/optical/swarp/bin/swarp'

gfits = "../data/host_sdss/frame-g-005116-5-0061.fits"
rfits = "../data/host_sdss/frame-r-005116-5-0061.fits"
ifits = "../data/host_sdss/frame-i-005116-5-0061.fits"

swarp_command = swarp_path \
        + " g.fits r.fits i.fits -c config.swarp -CENTER'" \
        + ra + " " + dec \
        + "; -SUBTRACT_BACK Y -RESAMPLE Y -COMBINE Y -IMAGE_SIZE'" \
        + '1700' + "," + '1700' + "'"
os.system(swarp_command)
