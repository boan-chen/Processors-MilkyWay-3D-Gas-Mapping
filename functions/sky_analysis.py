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

def open_model(idx):
    file_name = dir_dic['models'] + '/model_{}.npy'.format(idx)
    model = np.load(file_name)
    return model

catalog = pd.read_csv(dir_dic['catalog'] + 'boan_standard_catalog_ver2.csv')
# %%
