
import pandas as pd
import numpy as np
import healpy as hp
from astropy import units as u
from astropy.coordinates import SkyCoord
import random

valid_status = True #If True, load cleaned catalog
indexed_status = True #If True, load indexed catalog
sampled_status = True #If True, load sampled catalog
base_status = True #If True, load local bases
catalog_directory = '/Volumes/NTU_Astro/Catalog'

#%%
class catalog_process:
        
    def validity_selection(catalog):
        catalog_cleaned = catalog.copy()
        good_parallax = np.where(abs(catalog_cleaned['parallax']/catalog['parallax_err']) >= 3)[0]
        catalog_cleaned = catalog_cleaned.iloc[good_parallax].reset_index(drop = True)
        no_feh = np.where(catalog_cleaned['feh'] < -10)[0]
        catalog_cleaned = catalog_cleaned.drop(no_feh).reset_index(drop = True)
        low_snr = np.where(catalog_cleaned['SNR'] < 10)[0]
        catalog_cleaned = catalog_cleaned.drop(low_snr).reset_index(drop = True)
        catalog = catalog_cleaned
        return catalog_cleaned
    
    def indexing(catalog_cleaned, teff_bins, feh_bins, logg_bins, dist_bins, nside):
        catalog_indexed = catalog_cleaned.copy()
        idx = index(catalog_cleaned)
        grid_indices = idx.grid_idx(teff_bins, feh_bins, logg_bins)
        dist_indices = idx.dist_idx(dist_bins)
        pix_indices = idx.pix_idx(nside)
        catalog_indexed['grid_idx'] = grid_indices
        catalog_indexed['dist_idx'] = dist_indices
        catalog_indexed['pix_idx'] = pix_indices
        return catalog_indexed
        
    def sampling(catalog_indexed, size):
        sampled = random.sample(list(catalog_indexed.index), size)
        catalog_sampled = catalog_indexed.iloc[sampled].reset_index(drop = True)
        return catalog_sampled
    
    def dust_mapping(catalog_indexed, sky_map, nside):
        catalog_dust = dust_map(nside, catalog_indexed, sky_map)
        return catalog_dust

        
#%%

class index:
    def __init__(self, catalog):
        self.catalog = catalog
        
    def grid_idx(self, teff_bins, feh_bins, logg_bins):
        self.grid_indices = np.empty((len(self.catalog),), dtype=np.int64)

        # Find the bin indices for each value
        self.teff_indices = np.digitize(self.catalog['teff'], teff_bins) 
        self.logg_indices = np.digitize(self.catalog['logg'], logg_bins) 
        self.feh_indices = np.digitize(self.catalog['feh'], feh_bins) 
        self.num_logg_bins = len(logg_bins) 
        self.num_feh_bins = len(feh_bins) 
        self.num_teff_bins = len(teff_bins) 
        # Calculate the corresponding 1D index in the 3D grid space
        a, b, c  = self.teff_indices, self.logg_indices, self.feh_indices
        self.grid_indices = (a * (self.num_logg_bins * self.num_feh_bins)) + (b * self.num_feh_bins) + c
        return self.grid_indices

    def dist_idx(self, dist_bins):
        self.dist_indices = np.digitize(1/self.catalog['parallax'], dist_bins) - 1
        return self.dist_indices
    
    def galactic_coord(self, ra, dec):
        eq = SkyCoord(ra, dec, unit=u.deg)
        gal = eq.galactic
        l, b = gal.l.radian, gal.b.radian
        hp_b = -(b - np.pi/2)
        return l, b, hp_b   
    
    def pix_idx(self, nside):
        l, _, hp_b = self.galactic_coord(self.catalog['ra'], self.catalog['dec'])
        self.pix = hp.ang2pix(nside, hp_b, l, nest = True)
        return self.pix

    def index_to_digit(self, unique_index):
        c = unique_index % self.num_feh_bins
        temp = (unique_index - c) // self.num_feh_bins
        b = temp % self.num_logg_bins
        a = (temp - b) // self.num_logg_bins
        return a, b, c

def dust_map(nside, catalog, sky_map):
    dust = sky_map['Planck_EBV']
    dust_dg = hp.ud_grade(dust, nside_out=nside)
    HI = sky_map['HI']
    HI_dg = hp.ud_grade(HI, nside_out=nside)
    dust_panstar = sky_map['Panstarrs_EBV']
    dust_panstar_dg = hp.ud_grade(dust_panstar, nside_out=nside)
    # for catalog, use catalog['pix_idx'] to match dust_dg with catalog
    ebv = []
    ebv_panstar = []
    HI = []
    for i in range(0, len(catalog)):
        ebv.append(dust_dg[catalog['pix_idx'][i]])
        ebv_panstar.append(dust_panstar_dg[catalog['pix_idx'][i]])
        HI.append(HI_dg[catalog['pix_idx'][i]])
    catalog['Planck_EBV'] = ebv
    catalog['Panstarrs_EBV'] = ebv_panstar
    catalog['HI'] = HI
    return catalog

def catalog_by_criteria(criteria, catalog):
    selected_catalog = catalog
    for key in criteria.keys():
        # if the value for the key is a list, then it is a range. Else, the selection criteria should be equal to the value
        if type(criteria[key]) == list:
            selected_catalog = selected_catalog[(selected_catalog[key] >= criteria[key][0]) & (selected_catalog[key] <= criteria[key][1])].reset_index(drop=True)
        else:
            selected_catalog = selected_catalog[selected_catalog[key] == criteria[key]].reset_index(drop=True)
    return selected_catalog
