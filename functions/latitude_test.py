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

distant_criteria = {'teff': [4000, 10000], 'dist_idx': [2, np.inf]}
distant_catalog = catalog_by_criteria(distant_criteria, catalog)
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
