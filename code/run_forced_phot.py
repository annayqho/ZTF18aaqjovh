import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF_fast_transient_search/code")
import ForcePhotZTF
from ForcePhotZTF.phot_class import ZTFphot
from ztfquery import query
import glob


def get_forced_phot(name,ra,dec,jdobs):
    """ 
    Uses Yuhan's code to perform forced photometry 
    for a given RA and Dec and central JD (day of the transient)
    Over a window of +/- 10 days
    
    Parameters
    ----------
    name: name of the source 
    ra: ra in decimal degrees
    dec: dec in decimal degrees
    """
    start_jd = jdobs-10
    end_jd = jdobs+100

    zquery = query.ZTFQuery()
    zquery.load_metadata(
            radec=[ra, dec], size=0.0001, 
            sql_query="obsjd>%s and obsjd<%s" %(start_jd,end_jd))
    zquery.download_data("scimrefdiffimg.fits.fz")
    zquery.download_data("diffimgpsf.fits")

    print("finished downloading data")

    filefracday = zquery.metatable['filefracday'].values

    jd = []
    flux = []
    eflux = []
    mag = []
    emag = []
    filt = []

    # Normalize all fluxes to the same zp
    ZP = 25

    for ffd in filefracday:
        ffd = str(ffd)
        fpath = 'Data/sci/' + ffd[0:4] + '/' + ffd[4:8] + '/' + ffd[8:]

        # Check if directory is empty
        if len(glob.glob(fpath + "/*")) > 0:
            imgpath = glob.glob(fpath + "/*.fits.fz")[0]
            psfpath = glob.glob(fpath + "/*diffimgpsf.fits")[0]

            pobj = ZTFphot(name, ra, dec, imgpath, psfpath)
            pobj.load_source_cutout() 
            pobj.load_bkg_cutout()
            pobj.get_scr_cor_fn()  
            pobj.fit_psf()

            jd.append(pobj.obsjd)
            zp = pobj.zp
            new_flux = pobj.Fpsf * 10**(-0.4*(zp-ZP))
            flux.append(new_flux)
            new_eflux = pobj.eFpsf * 10**(-0.4*(zp-ZP))
            eflux.append(new_eflux)
            mag.append(pobj.mag)
            emag.append(pobj.mag_unc)
            filt.append(pobj.filter)

    jd = np.array(jd)
    flux = np.array(flux)
    eflux = np.array(eflux)
    mag = np.array(mag)
    emag =np.array(emag)
    filt = np.array(filt)
    
    return filt,jd,flux,eflux,mag,emag


if __name__=="__main__":
    name = 'ZTF18aaqjovh'
    ra = 178.181754
    dec = 25.675033
    jdobs = 2458256.7235
    filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs)
