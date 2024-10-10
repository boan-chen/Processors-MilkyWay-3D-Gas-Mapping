#%%
import numpy as np
import astropy.io.fits as pf
import pandas as pd
import os.path
import os
from tqdm import tqdm
import logging
import SpecDownload as sd
from scipy.interpolate import interp1d
from mpi4py import MPI

#%%
class SpecReader:
    def __init__(self, catalog, st_wave, path, directory = None):
        self.st_wave = st_wave
        self.path = path
        self.spectra_directory = path['stellar_spectra']
        self.directory = directory
        self.catalog = catalog
        self.spec_name = catalog['spec_name']
        self.plate_name, self.fiber_name, self.mjd_name = 'plate_1', 'fiber_1', 'mjd_1'
        data_formats = {
            pd.DataFrame: 'pandas',
            np.ndarray: 'numpy',
            pf.hdu.hdulist.HDUList: 'fits'
        }
        self.data_format = data_formats.get(type(catalog))
        if self.data_format is None:
            print("Unsupported catalog format.")
            return None
        
    def read(self, i):
        base_dir = self.spectra_directory
        if self.data_format == 'pandas':
            spectrum_name = self.spec_name[i]
            Path = base_dir
        else:
            plate = '{0:0>4}'.format(int(self.catalog[self.plate_name][i]))
            fiberid = '{0:0>4}'.format(int(self.catalog[self.fiber_name][i]))
            mjd = '{0:0>5}'.format(int(self.catalog[self.mjd_name][i]))
            spectrum_name = 'spec-' + plate + '-' \
                    + mjd + '-' + fiberid + '.fits'
            Path = base_dir + plate + '/' 
        if os.path.isfile(Path + spectrum_name) == False:
            sd.download_request(self.catalog, i, self.path, enforce = False)
        try:
            spectrum_read = pf.open(Path + spectrum_name)[1].data
        except (OSError, TypeError, IndexError) as e:
            # Log the reason for invalidity and continue the loop
            logging.error(f"File {i} is invalid: {str(e)}")
            sd.download_request(self.catalog, i, self.path, enforce = True)
            print(f"File{i} is downloaded")
            spectrum_read = pf.open(Path + spectrum_name)[1].data

        return spectrum_read

    def common_pix(self, log_wave, flux, ivar, model):
        wave = 10**log_wave
        st_wave = self.st_wave
        interp_flux_resampled = interp1d(wave, flux, kind='linear', bounds_error=False, fill_value=np.nan)
        flux_sampled = interp_flux_resampled(st_wave)
        interp_ivar_resampled = interp1d(wave, ivar, kind='linear', bounds_error=False, fill_value=np.nan)
        ivar_sampled = interp_ivar_resampled(st_wave)
        interp_model_resampled = interp1d(wave, model, kind='linear', bounds_error=False, fill_value=np.nan)
        model_sampled = interp_model_resampled(st_wave)
        return wave, flux_sampled, ivar_sampled, model_sampled

    
    def load_data(self, loop=None, tqdm_disable=False):
        catalog = self.catalog
        st_wave = self.st_wave
        
        if loop is None:
            if self.data_format == 'pandas':
                loop = len(catalog)
            elif self.data_format in ('numpy', 'fits'):
                loop = len(catalog[1].data)

        flux_array = np.empty((loop, len(st_wave)))
        ivar_array = np.empty((loop, len(st_wave)))
        model_array = np.empty((loop, len(st_wave)))
        self.invalid_file_list = []
        print("Loading data...")

        for i in tqdm(range(loop), disable = tqdm_disable):
            if self.data_format == 'pandas':
                if i not in catalog.index:
                    print("Catalog index is not reset.")
                    return None
            reader = self.read(i)
            log_wave, flux, ivar, model = reader['loglam'], reader['flux'], reader['ivar'], reader['model']
            output = self.common_pix(log_wave, flux, ivar, model)
            flux_array[i], ivar_array[i], model_array = output[1], output[2], output[3]

        if self.data_format == 'pandas':
            loaded_catalog = catalog.drop(catalog.index[self.invalid_file_list])
        elif self.data_format in ('numpy', 'fits'):
            loaded_catalog = catalog[1].data[np.delete(np.arange(loop), self.invalid_file_list)]

        flux_array = np.delete(flux_array, self.invalid_file_list, axis=0)
        ivar_array = np.delete(ivar_array, self.invalid_file_list, axis=0)
        model_array = np.delete(model_array, self.invalid_file_list, axis=0)

        return loaded_catalog, flux_array, ivar_array, model_array

    def save_data(self, catalog, flux_array, ivar_array, model_array):
        if os.path.isdir(self.directory + '/samples') == False:
            os.system("mkdir "+self.directory + '/samples')
        np.save(self.directory + '/samples/flux_array.npy', flux_array)
        np.save(self.directory + '/samples/ivar_array.npy', ivar_array)
        np.save(self.directory + '/samples/model_array.npy', model_array)
        if isinstance(catalog, pd.DataFrame):
            catalog.to_csv(self.directory + '/samples/catalog.csv', index=False)
        elif isinstance(catalog, np.ndarray):
            pf.writeto(self.directory + '/samples/catalog.fits', catalog)
        else:
            print("Unsupported catalog format.")
        return
    

#%% test
# import pandas as pd
# import random
# subcatalog = output[0]
# path.directory = "/Volumes/NTU_Astro/Workspace/sampling_v4"
# catalog = pd.read_csv(path.directory + '/catalog/catalog_sampled.csv')
# # choice = random.sample(range(0, len(catalog)), 3000)
# # subcatalog = catalog.iloc[choice].reset_index(drop = True)
# rdr = SpecReader(subcatalog, path)
# output = rdr.load_data()
# rdr.save_data(output[0], output[1], output[2], output[3])

# %%
