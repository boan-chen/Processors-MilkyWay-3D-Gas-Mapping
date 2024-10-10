#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from SpecReader import SpecReader
from rv_shift import resampling
from SpecDownload import download
from continuum import continuum
import healpy as hp

from SpecTools import *
from tqdm import tqdm
dir_dic = {'catalog': '/Volumes/NTU_Astro/Catalog/',
           'stellar_spectra' : '/Volumes/NTU_Astro/Spectra/Stellar_Spectra/',
           'wave': '/Volumes/NTU_Astro/Spectra/segue_wave/'
           }
st_wave = np.load(dir_dic['wave'] + 'segue_standard_wave.npy')

#%%
def extract_local(catalog, dist_idx, grid_size):
    local_bases = {}
    local_bases_counter = {}

    for i in tqdm(range(grid_size)):
        catalog_loaded = catalog[catalog['grid_idx'] == i]
        counter = 0

        if catalog_loaded.empty:
            local_bases[i] = None
            local_bases_counter[i] = counter
            continue

        local_catalog = catalog_loaded[catalog_loaded['dist_idx'] <= dist_idx].reset_index(drop=True)

        if local_catalog.empty:
            local_bases[i] = None
            local_bases_counter[i] = counter
            continue

        local = SpecReader(local_catalog, dir_dic)
        local_catalog, local_flux, local_ivar, local_model = local.load_data()

        if local_flux is None:
            local_bases[i] = None
            local_bases_counter[i] = 0
            continue

        resampled_flux = resampling(st_wave, local_flux, local_catalog['RV'], reverse=False)
        local_bases_counter[i] = len(local_flux)

        normalizer = continuum(st_wave, resampled_flux)
        local_norm_flux = normalizer.norm_all_flux()

        if not local_norm_flux:
            local_bases[i] = None
            local_bases_counter[i] = 0
            continue

        local_bases_counter[i] = len(local_norm_flux)
        local_bases[i] = composite(local_norm_flux)

    print('Done: local bases extracted')
    return local_bases, local_bases_counter


def extract_residual(catalog, nside, dist_idx, local_bases, local_bases_counter, dir_dic, threshold=1):
    npix = hp.nside2npix(nside)
    residuals = {}
    residuals_counter = {}
    for i in tqdm(range(npix)):
        subcatalog_init = catalog[catalog['pix_idx'] == i].reset_index(drop=True)
        subcatalog = subcatalog_init[subcatalog_init['dist_idx'] >= dist_idx].reset_index(drop=True)  # new_added
        counter = 0
        if subcatalog.empty:
            residuals_counter[i] = counter
            continue
        pix_loader = SpecReader(subcatalog, dir_dic)
        subcatalog, subflux, _ = pix_loader.load_data()
        if not subcatalog:
            residuals_counter[i] = counter
            continue
        normalizer = continuum(st_wave, subflux)
        pix_flux_norm = normalizer.norm_all_flux()
        pix_flux_shifted = resampling(st_wave, pix_flux_norm, subcatalog['RV'], reverse=False)

        pix_flux_residual_r = np.zeros_like(pix_flux_shifted)

        for j in range(len(subcatalog)):
            local_base = local_bases[subcatalog['grid_idx'].iloc[j]]

            if local_base is None or local_bases_counter[subcatalog['grid_idx'].iloc[j]] <= threshold:
                pix_flux_residual_r[j] = np.array([np.nan] * len(pix_flux_shifted[j]))
                continue
            else:
                pix_flux_residual_r[j] = pix_flux_shifted[j] / local_base

        nan_rows = np.all(np.isnan(pix_flux_residual_r), axis=1)

        if len(nan_rows) == len(pix_flux_residual_r):
            residuals[i] = None
            residuals_counter[i] = 0
            continue

        pix_flux_residual_r = pix_flux_residual_r[~nan_rows]
        pix_flux_residual_o = resampling(st_wave, pix_flux_residual_r, subcatalog['RV'], reverse=True)
        pix_flux_residual = composite(pix_flux_residual_o)

        residuals_counter[i] = len(pix_flux_residual_r)
        residuals[i] = pix_flux_residual

    print('Done: residuals extracted')
    return residuals, residuals_counter


# %%
