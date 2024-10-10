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
from rv_shift import resampling
from tqdm import tqdm
from server_shifting import server
from multiprocessing import Pool

serv, servers, dir_dic = server()
dir_dic['models'] = servers[serv] + 'Spectra/models_ver1'

st_wave = np.load(dir_dic['wave'] + 'segue_standard_wave.npy')
catalog = pd.read_csv(dir_dic['catalog'] + 'boan_standard_catalog_ver2.csv')

def model(local_catalog, grid, threshold = 30):
    # check if the model already exists
    if os.path.isfile(dir_dic['models'] + '/model_{}.npy'.format(grid)):
        return None
    local_catalog = local_catalog[local_catalog['grid_idx'] == grid].reset_index(drop = True)
    if len(local_catalog) <= threshold:
        return None
    readr_local = SpecReader(local_catalog, st_wave, path = dir_dic)
    c_local, f_local, i_local, m_local = readr_local.load_data(tqdm_disable=True)
    norm_local, _ = pn(st_wave, f_local, i_local, tqdm_disable=True)
    norm_local = resampling(st_wave, norm_local, c_local['RV'], reverse = False)
    model = np.median(norm_local, axis=0)
    np.save(dir_dic['models'] + '/model_{}.npy'.format(grid), model)
    return model

criteria = {'teff': [4000, 10000], 'dist_idx': [0, 1.5]}
local_catalog = catalog_by_criteria(criteria, catalog)
grid_idx = local_catalog['grid_idx'].unique()

def parallel_model(args):
    local_catalog, idx, threshold = args
    model(local_catalog, idx, threshold)

if __name__ == "__main__":
    grid_idx = range(len(grid_idx))
    num_processes = 8
    with Pool(num_processes) as pool:
        for _ in tqdm(pool.imap_unordered(parallel_model, [(local_catalog, idx, 30) for idx in grid_idx]), total=len(grid_idx)):
            pass
# %%
