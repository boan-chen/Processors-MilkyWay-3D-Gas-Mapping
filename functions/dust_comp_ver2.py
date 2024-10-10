#%%
import numpy as np
import pandas as pd
import healpy as hp
import matplotlib.pyplot as plt
import os
# go to parent directory and import continuum
import sys
sys.path.append('..')

from SpecReader import SpecReader
from catalog_process import catalog_process, catalog_by_criteria
import SpecTools as st
from rv_shift import resampling
from SpecNormalization import poly_normalization as pn
from Visualization import spectrum_visualization as sv
from Visualization import spectra_visualization as ssv
from tqdm import tqdm
import concurrent.futures
from functools import partial
from server_shifting import server

nside = 16
serv, servers, dir_dic = server()
st_wave = np.load(dir_dic['wave'] + 'segue_standard_wave.npy')
catalog = pd.read_csv(dir_dic['catalog'] + 'boan_standard_catalog_ver2.csv')


def distribution_matching(data_a, data_b):
    hist_a, bin_edges = np.histogram(data_a, bins='auto')
    indecies = np.searchsorted(bin_edges, data_b, side = 'right')
    indecies_prob = np.zeros(len(bin_edges) + 1)
    indecies_prob[1:-1] = hist_a
    prob = np.empty_like(data_b)
    for i in range(0, len(data_b)):
        prob[i] = indecies_prob[indecies[i]]
    prob = prob / prob.sum()
    return prob

def residual_strapped(local_catalog, distant_catalog, shift = True, match = True):
    readr_local = SpecReader(local_catalog, st_wave, path = dir_dic)
    c_local, f_local, i_local, m_local = readr_local.load_data(tqdm_disable=False)
    norm_local, _ = pn(st_wave, f_local, i_local, tqdm_disable=True)
    readr_dist = SpecReader(distant_catalog, st_wave, path = dir_dic)
    c_dist, f_dist, i_dist, m_dist = readr_dist.load_data(tqdm_disable=False)
    norm_dist, _ = pn(st_wave, f_dist, i_dist, tqdm_disable=True)
    if shift is True:
        norm_dist = resampling(st_wave, norm_dist, distant_catalog['RV'], reverse=False)
    if match is True:
        prob = distribution_matching(distant_catalog['feh'], local_catalog['feh'])
        idx = np.random.choice(len(norm_local), size=len(norm_local), p=prob)
        norm_local = norm_local[idx]
        local_feh_m = local_catalog['feh'][idx]
        # plt.hist(local_feh_m, bins='auto', alpha = 0.5,  label = 'matched')
        # plt.hist(distant_catalog['feh'], bins='auto', alpha = 0.5, label = 'distant')
        # plt.hist(local_catalog['feh'], bins='auto', alpha = 0.5, label = 'local')
        # plt.legend()
    model = st.composite(norm_local)
    if shift is True:
        dist_flux_shifted = resampling(st_wave, norm_dist, distant_catalog['RV'], reverse=False)
        residual_unstacked = dist_flux_shifted/model
        residual_obs = resampling(st_wave, residual_unstacked, distant_catalog['RV'], reverse=True)
    else:
        residual_obs = norm_dist/model
    if len(residual_obs) == 0:
        print('No residual')
        return None, None, None, None, None
    return model, residual_obs
    
#%%
# grids = [4984,  4985,  5666,  4983,  6348,  5667,  4271,  4270,  4952,  5665,\
#     2134,  1260,  2133,  1261,  1259,  2135,  2103,  1292,  1293,  2132]

def process_single_grid_idx(j, catalog, distant_catalog_all):
    local_catalog = catalog_by_criteria({'grid_idx': j, 'dist_idx': [0, 1.5]}, catalog)
    if len(local_catalog) <= 30:
        return None
    distant_catalog = catalog_by_criteria({'grid_idx': j}, distant_catalog_all)
    if len(distant_catalog) <= 10:
        return None
    model, residual_obs = residual_strapped(local_catalog, distant_catalog, shift=True, match=True)
    return residual_obs


for i in range(0, 5):
    b = [i * 20, (i + 1) * 20]
    criteria_distant = {'b': b, 'dist_idx': [2, np.inf], 'teff': [5000, 6000]}
    distant_catalog_all = catalog_by_criteria(criteria_distant, catalog)
    # find the unique grid_idx
    grid_idx = distant_catalog_all['grid_idx'].unique()
    for j in tqdm(grid_idx):
        name = f'{dir_dic["latitude_residuals"]}/b_{(i+1)*10}_grid_{j}.npy'
        # check if the file exists
        if os.path.isfile(name):
            print(f'File {name} exists')
            continue
        residual_obs = process_single_grid_idx(j, catalog, distant_catalog_all)        
        np.save(f'{dir_dic["latitude_residuals"]}/b_{(i+1)*10}_grid_{j}.npy', residual_obs)
    # residual_final = np.vstack(residual)
    # composite_residual = st.composite(residual_final)
    # name = 'Residual in latitude range: ' + str(b)
    # name_composite = 'Composite residual in latitude range: ' + str(b)
    # np.save(f'latitude_sample/{name}.png', residual_final)
    # np.save(f'latitude_sample/{name_composite}.png', composite_residual)
    # sv(st_wave, residual_final, xlim=[5800, 6200], ylim=[0.8, 1.2], title='Residual in latitude range: ' + str(b))


# %%
