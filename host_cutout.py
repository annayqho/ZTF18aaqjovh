""" cutout of the host """

from matplotlib.colors import LogNorm
import glob
import requests
from requests.exceptions import ConnectionError
from astroquery.mast import Observations
from astropy import coordinates
import astropy.units as u
from astropy.io import fits
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")


def mast():
    """ access PS1 images through astroquery/MAST 
    the problem is that I don't think I can get image cutouts this way,
    only the full stacked image in the region.
    I don't want to have to do the cutting myself,
    so for now I'm going to use the PS1 image access webpage.
    """
    # radius should be 25 arcseconds, which is 0.00694 degrees
    obsTable = Observations.query_criteria(
            dataproduct_type=["image"],
            obs_collection=["PS1"],
            coordinates="331.675712 +36.208096", 
            radius="0.00694 deg")
    obsids = obsTable['obsid']

    dataProductsByID = Observations.get_product_list(obsids)
    Observations.download_products(
            obsids, curl_flag=True, # only the code to download it later, 
            mrp_only=True) # minimum recommended products
    # you get a .sh script that you can then run to download the products


def ps1(ra, dec):
    """ 
    use the PS1 image cutout server.
    
    documentation:
    http://hla.stsci.edu/fitscutcgi_interface.html    
    """

    # step 1: image list service
    s = requests.Session()
    www = "http://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    payload = {}
    payload["ra"] = ra
    payload["dec"] = dec
    payload["filters"] = 'i'
    payload["type"] = 'stack'
    payload["sep"] = 'comma'
    try:
        r = s.get(www, params=payload, timeout=(3.05,10))
    except requests.exceptions.RequestException as e: 
        print(e)
        sys.exit(1)
    print(r.url)
    if r.status_code != requests.codes.ok:
        print("status code" + str(r.status_code))
    s.close()
    raw = r.text.split("\n")
    header = np.array(raw[0].split(","))
    content = np.array(raw[1].split(","))
    sname = content[header=='shortname'][0]
    fname = content[header=='filename'][0]
    print(fname)
    
    # step 2: download the FITS file of the full image
    # server_name = "http://ps1images.stsci.edu"
    # url = server_name + fname
    # wget url

    # step 2: download the FITS file of the image cutout
    # 30 arcseconds by 30 arcseconds
    if glob.glob(sname):
        print("already saved")
    else:
        www = "http://ps1images.stsci.edu/cgi-bin/fitscut.cgi" 
        payload = {}
        payload["red"] = fname
        payload["format"] = 'fits'
        payload["size"] = '120'
        payload["ra"] = ra
        payload["dec"] = dec
        try:
            r = s.get(www, params=payload, timeout=(3.05,10), allow_redirects=True)
            open(sname, 'wb').write(r.content)
        except requests.exceptions.RequestException as e: 
            print(e)
            sys.exit(1)
        print(r.url)
        if r.status_code != requests.codes.ok:
            print("status code" + str(r.status_code))
        s.close()

    return sname # saved file


def plot(ax, fname, vmin_val, vmax_val):
    """ plot image cutout in Python """
    hdu_list = fits.open(fname)
    image_data = hdu_list[0].data
    hdu_list.close()
    values = np.ma.masked_invalid(image_data.flatten())
    counts = np.ma.masked_invalid(image_data[::-1])
    ax.imshow(
            np.arcsinh(counts), cmap='Greys_r', 
            vmin=np.percentile(np.arcsinh(values), vmin_val),
            vmax=np.percentile(np.arcsinh(values), vmax_val))
    # N up, E left
    ax.axis('off')

    ax.scatter(60, 60, marker='x', color='white', s=60, linewidths=1.0)

    # compass markings
    # ax.annotate("", xy=(10, 9), xytext=(35,9),
    #         arrowprops=dict(facecolor='black', width=1))
    # ax.annotate("", xy=(35,34), xytext=(35,9),
    #         arrowprops=dict(facecolor='black', width=1))
    # ax.annotate("E", xy=(7,9), fontsize=11,
    #         horizontalalignment='right',
    #         verticalalignment='center')
    # ax.annotate(
    #         "S", xy=(35,38), fontsize=11, 
    #         horizontalalignment='center',
    #         verticalalignment='top')

    # line to show scale
    ax.annotate("", xy=(70,110), xytext=(110,110),
            arrowprops=dict(arrowstyle="-", color='white'))
    ax.annotate(
            "10''", xy=(90,110), fontsize=11, 
            horizontalalignment='center',
            verticalalignment='bottom', color='white')


if __name__=="__main__":
    ra = 178.181754 
    dec = 25.675033

    vmin = 90
    vmax = 100

    fig,ax= plt.subplots(1,1, figsize=(6,5))
    fname = ps1(ra, dec) # file saved by this function
    plot(ax, fname, vmin, vmax)
    plt.subplots_adjust(wspace=0.1, hspace=0.01)
    plt.savefig("host_cutout.png")
    #plt.show()

