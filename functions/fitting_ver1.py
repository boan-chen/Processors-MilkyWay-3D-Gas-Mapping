#%%
import numpy as np
import random
import json
import pandas as pd
import healpy as hp
import matplotlib.pyplot as plt
# go to parent directory and import continuum
import os
os.chdir('..')
from continuum import continuum
from SpecReader import SpecReader
from catalog_process import catalog_process, catalog_by_criteria
import SpecTools as st
from astropy.io import fits
from SpecNormalization import poly_normalization as pn
from Visualization import spectrum_visualization as sv
from Visualization import spectra_visualization as ssv
import scipy.stats as ss
import seaborn as sns
from rv_shift import resampling
from tqdm import tqdm
import sys

nside = 16

server = {'studio': '/Users/boanchen/Desktop/',
          'usb': '/Volumes/NTU_Astro/'}
serv = 'studio'
dir_dic = {'catalog': server[serv] + "Catalog/",
           'stellar_spectra' : server[serv] + 'Spectra/Stellar_Spectra/',
           'wave': server[serv] + 'Spectra/segue_wave/',
           'latitude_residuals': server[serv] + 'Spectra/latitude_residuals/',
           }
# if any of the dir_dic does not exist, switch to usb
for key in dir_dic:
    if os.path.isdir(dir_dic[key]) == False:
        serv = 'usb'
        dir_dic = {'catalog': server[serv] + "Catalog/",
                   'stellar_spectra' : server[serv] + 'Spectra/Stellar_Spectra/',
                   'wave': server[serv] + 'Spectra/segue_wave/',
                   'latitude_residuals': server[serv] + 'Spectra/latitude_residuals/',
                   }
        break

sys.path.append('..')

st_wave = np.load(dir_dic['wave'] + 'segue_standard_wave.npy')
catalog = pd.read_csv(dir_dic['catalog'] + 'boan_standard_catalog_ver2.csv')

grids = [4984,  4985,  5666,  4983,  6348,  5667,  4271,  4270,  4952,  5665,\
    2134,  1260,  2133,  1261,  1259,  2135,  2103,  1292,  1293,  2132]

grid_idx = 4984
criteria_local = {'dist_idx': [0, 1], 'grid_idx': grid_idx}
local_catalog = catalog_by_criteria(criteria_local, catalog)
local_readr = SpecReader(local_catalog, st_wave, path = dir_dic)
criteria_dist = {'dist_idx': [2, np.inf], 'grid_idx': grid_idx}
distant_catalog = catalog_by_criteria(criteria_dist, catalog)
distant_readr = SpecReader(distant_catalog, st_wave, path = dir_dic)
c_local, f_local, i_local, m_local = local_readr.load_data()
c_distant, f_distant, i_distant, m_distant = distant_readr.load_data()

#%%
import scipy.optimize as so
from scipy import ndimage

def black_body_radiation(wave, T):
    # wave in angstrom
    # T in K
    h = 6.626e-34
    c = 3e8
    k = 1.38e-23
    wave = wave * 1e-10
    return 2*h*c**2/wave**5/(np.exp(h*c/(wave*k*T))-1)

def polynomial(x, order, coeff):
    result = 0
    if len(coeff) != order + 1:
        print('Error: order and coeff not match')
        return None
    for i in range(0, order+1):
        result += coeff[i] * x**i
    return result

def continuum_model(wave, coeff):
    return polynomial(wave, len(coeff)-1, coeff)

def fit_func(x, *coeff):
    return continuum_model(x, coeff)


# fit sample with continuum model
def continuum(sample, st_wave):
    wave = st_wave.copy()
    flux_filtered = ndimage.percentile_filter(sample, percentile = 90, size = 51)
    finite = np.isfinite(flux_filtered)
    filtered, wave = flux_filtered[finite], st_wave[finite]
    popt, pcov = so.curve_fit(fit_func, wave, filtered, p0=[1]*8)
    conti = fit_func(st_wave, *popt)
    return conti

def normalization_all(st_wave, flux, catalog, ivar):
    normed = np.zeros_like(flux)
    std_normed = np.zeros_like(flux)
    print('Normalizing all fluxes...')
    for i in tqdm(range(len(flux))):
        # turn np.nan to 0
        sample, T =flux[i], catalog['teff'][i]
        sample[np.isnan(sample)] = 0
        conti = continuum(sample, st_wave)
        std = 1/np.sqrt(ivar[i])
        std[np.isnan(std)] = 0
        normed[i] = sample/conti
        std_normed[i] = std/conti
    return normed, std_normed

def percentile_filter(flux, percentile = 90, size = 51):
    norm = np.zeros_like(flux)
    for i in range(len(flux)):
        sample = flux[i]
        filtered = ndimage.percentile_filter(sample, percentile = percentile, size = size)
        # finite = np.isfinite(filtered)
        # filtered, wave = filtered[finite], st_wave[finite]
        norm[i] = sample/filtered
    return norm
# normed_local, std_normed_local = normalization_all(st_wave, f_local, c_local, i_local)
normed_local, std = normalization_all(st_wave, f_local, c_local, i_local)
# leave only finite values

#%%
corr_matrix = np.corrcoef(normed_local.T[:400])
# sns.heatmap(corr_matrix, cmap='coolwarm')

arr = np.array(find)
# k = empty set
k = []
set_arr = list(set(arr[0]))
n = len(set_arr)
for i in set_arr:
    target = arr[1][np.where(arr[0] == i)[0]]
    for j in range(0, n - i - 1):
        to_save = None
        step = 1
        if i - j >= 0 and i + j < n:
            element = np.arange(i - j, i + j + 1)
            neighbor_all_exist = np.all(np.isin(element, target))
            if j >= 1:
                to_save = element
                step = j
            if neighbor_all_exist == False:
                print(i, j)
                k.append(to_save)
                break

# find the arrays in list k that has len larger than 5, and turn all extract all the elements in it as a set
new_k = []
for i in k:
    if len(i) >= 5:
        new_k.append(i)




# based on the corr_matrix, find the highly correlated pairs in diagonal blocks


# %%
# use local as training set in pca, then apply pca to reconstruct distant
import sklearn.decomposition as skd
pca = skd.PCA(n_components=10)
pca.fit(normed_local)
pca_normed_local = pca.transform(normed_local)
normed_distant_shift = resampling(st_wave, normed_distant, c_distant['RV'], reverse = False)
normed_distant_shift[np.isnan(normed_distant_shift)] = 0
pca_normed_distant = pca.transform(normed_distant_shift)
recon_normed_distant = pca.inverse_transform(pca_normed_distant) 
recom_normed_distant = resampling(st_wave, recon_normed_distant, c_distant['RV'], reverse = True)
    
#%%
# for i in range(0, 20):
#     spectra = [normed_distant[i], recon_normed_distant[i]]
#     labels = ['original', 'reconstructed']
#     ssv(st_wave, spectra, label_list = labels, title = f'{i}')
    
residual = normed_distant/recon_normed_distant
plt.figure(figsize=(10, 5))
plt.plot(st_wave, np.median(residual, axis = 0))
plt.xlim(5800, 6000)
plt.ylim(0.95, 1.05)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Residual')
# plt.vlines([3933.66, 3968.47], 0.9, 1.1, color='red', alpha=0.5)
plt.show()
