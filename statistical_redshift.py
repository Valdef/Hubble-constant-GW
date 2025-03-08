#!/usr/bin/env python
# coding: utf-8

# In[1]:


# convert dec from degree to ° ' " 
def dec_to_degminsec(degree):
    deg = int(degree)
    min = int((degree - deg)*60)
    sec = (((degree - deg)*60) - min) * 3600

    return deg, min, sec


# In[2]:


# convert ra from degree to hours min sec
def ra_to_hhmmss(deg):
    hours = int(deg/15)
    min = int(((deg/15) - hours) * 60)
    sec = ((((deg/15) - hours) * 60) - min) * 60

    return hours, min, sec


# In[3]:


import numpy as np
import pandas as pd
import healpy as hp
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u

name_event = 'GW170817'

# open file
with fits.open(f"./skymap/{name_event}_skymap.fits.gz") as hdul:
    data = hdul[1].data
    prob = data['PROB'] # get proba per pixels
    header = hdul[1].header
    nside = header['NSIDE'] # resolution of map

# Define acceptable distance range
min_dist = 21
max_dist = 88

# Find pixel with maximum probability
max_pix = np.argmax(prob)

# Convert pixel index to angular coordinates
theta, phi = hp.pix2ang(nside, max_pix, nest=True) 
ra  = np.degrees(phi) 
dec = 90 - np.degrees(theta)

# convert ra and dec and hhmmss and ° ' "
deg, m, s = dec_to_degminsec(dec)
h, min, sec = ra_to_hhmmss(ra)

print(f"Most probable {name_event} event location: RA={ra:.2f}, Dec={dec:.2f}")
print('')
print(f"Most probable {name_event} event in hhmmss and ° ': "f"RA={h}h {min}min {sec:.2f}sec, Dec={deg}° {m}min, {s:.2f}sec")
print('')
print(f"Distance range: {min_dist:.2f} - {max_dist:.2f} Mpc")


# In[4]:


# Define column names for the GLADE+ catalog
col_names = [
    "GLADE_no", "PGC_no", "GWGC_name", "HyperLEDA_name", "2MASS_name", "WISExSCOS_name",
    "SDSS-DR16Q_name", "Object_type_flag", "RA", "Dec", "B_mag", "B_err", "B_flag", "B_Abs",
    "J_mag", "J_err", "H_mag", "H_err", "K_mag", "K_err", "W1_mag", "W1_err", "W2_mag", "W2_err",
    "W1_flag", "BJ_mag", "BJ_err", "z_helio", "z_cmb", "z_flag", "v_err", "z_err", "d_L",
    "d_L_err", "dist_flag", "M_star", "M_star_err", "M_star_flag", "Merger_rate", "Merger_rate_err"
]


# In[5]:


# Load high-probability GW sky localization
def load_gw_skymap(fits_file, prob_threshold=0.5):
    with fits.open(fits_file) as hdul:
        data = hdul[1].data
        prob = data["PROB"]
        nside = hdul[1].header["NSIDE"]

        prob_thresh = np.percentile(prob, 100 * (1 - prob_threshold)) 
        significant_pixels = np.where(prob >= prob_thresh)[0] # get pixels above threshold

        theta, phi = hp.pix2ang(nside, significant_pixels, nest=True)
        ra = np.rad2deg(phi)
        dec = 90 - np.rad2deg(theta)

    return SkyCoord(ra=ra * u.degree, dec=dec * u.degree)


# In[6]:


gw_coords = load_gw_skymap(f"./skymap/{name_event}_skymap.fits.gz", 0.683)


# In[7]:


def match_galaxies_chunk(chunk, gw_coords, min_dist, max_dist):
    chunk = chunk.dropna(subset=["RA", "Dec", "d_L", "z_helio", "z_err"])
    
    # Filter by distance
    chunk = chunk[(chunk["d_L"] >= min_dist) & (chunk["d_L"] <= max_dist)]
    if chunk.empty:
        return None

    # Match galaxies to GW sky region
    galaxy_coords = SkyCoord(ra=chunk["RA"].values * u.degree, dec=chunk["Dec"].values * u.degree)
    idx, d2d, _ = gw_coords.match_to_catalog_sky(galaxy_coords)
    
    # Keep galaxies within 1 degree of GW region
    matched = chunk.iloc[idx[d2d < 1 * u.degree]]
    return matched


# In[44]:


from tqdm import tqdm
redshifts = []
redshifts_err = []
# find redshift of galaxies by chunks
for chunk in tqdm(pd.read_csv("./GLADE+.txt", sep=" ", header=None, names=col_names, chunksize=10**5)):
    matched = match_galaxies_chunk(chunk, gw_coords, min_dist, max_dist)
    if matched is not None:
        redshifts.extend(matched["z_helio"].values)
        redshifts_err.extend(matched["z_err"].values)


# In[53]:


redshifts_err = np.array(redshifts_err)
redshifts = np.array(redshifts)
redshifts_err1 = redshifts_err[redshifts_err<=redshifts*0.5]


# In[ ]:


# get statistical redshift
from uncertainties import unumpy
redshift_w_error = unumpy.uarray(redshifts, redshifts_err)
redshift_w_error = np.mean(redshift_w_error)
print(f"Mean redshift: {redshift_w_error}")


# In[ ]:


c = 3*10**5 # km/s
v = c*redshift_w_error # recessional velocity


# In[ ]:


# distances from other jupyter notebook
distance = 48 # Mpc
distance_err1 = 40
distance_err2 = -27


# In[ ]:


from uncertainties import ufloat, unumpy
distance_w_err1 = ufloat(distance, distance_err1)
distance_w_err2 = ufloat(distance, abs(distance_err2))


# In[ ]:


# calculate H_0 with uncertainties
H_0 = unumpy.nominal_values(v)/unumpy.nominal_values(distance_w_err1) # kmsMpc
H_err1 = v/distance_w_err1
H_err2 = v/distance_w_err2

print(f"H_0 = {unumpy.nominal_values(H_0):.3f} (+{unumpy.std_devs(H_err1):.3f}, -{unumpy.std_devs(H_err2):.3f}) Km/s/Mpc")


# In[ ]:


z=ufloat(0.009877, 1.67*10**-5)
v = c*z # recessional velocity
H_0 = unumpy.nominal_values(v)/unumpy.nominal_values(distance_w_err1)
H_err1 = v/distance_w_err1
H_err2 = v/distance_w_err2

print(f"H_0 = {unumpy.nominal_values(H_0):.3f} (+{unumpy.std_devs(H_err1):.3f}, -{unumpy.std_devs(H_err2):.3f}) Km/s/Mpc")


# In[ ]:


import numpy as np
x = np.array([0, 1, 2, 3,4,5,6,7])
redshift_w_error = x[x<=4]
print(redshift_w_error)


# In[ ]:




