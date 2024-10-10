#%%
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
#%%

def resampling(st_wave, flux_array, radial_velocity_list, reverse = False, ivar_array=None):
    """
    Resample the flux and ivar arrays to rest-frame velocity using radial velocity information.

    Parameters:
        st_wave (numpy array): The original wavelength grid.
        flux_array (numpy array): The original flux array with shape (num_spectra, num_wavelengths).
        radial_velocity_list (list): List of radial velocities corresponding to each spectrum in km/s.
        ivar_array (numpy array, optional): The original inverse variance (ivar) array with shape (num_spectra, num_wavelengths).
                                           If None, only the flux will be resampled.

    Returns:
        resampled_flux (numpy array): The resampled flux array with shape (num_spectra, num_wavelengths).
        resampled_ivar (numpy array, optional): The resampled inverse variance (ivar) array with shape (num_spectra, num_wavelengths).
                                                Only returned if ivar_array is provided.
    """
    # Calculate rest-frame velocity
    c = 3e5  # Speed of light in km/s
    shift = 1 + np.array(radial_velocity_list) / c
    if reverse == True:
        shift = 1 - np.array(radial_velocity_list) / c

    # Calculate rest-frame wavelengths
    restframe_wavelengths = st_wave / shift[:, np.newaxis]

    # Initialize an array to store the interpolated flux spectra
    resampled_flux = np.empty_like(flux_array)
    for i in range(len(flux_array)):
        interp_function_resampled = interp1d(restframe_wavelengths[i], flux_array[i], kind='linear', bounds_error=False, fill_value=np.nan)
        resampled_flux[i] = interp_function_resampled(st_wave)
    
    if ivar_array is None:
        return resampled_flux
    else:
        # Initialize an array to store the interpolated var spectra

        # Initialize an array to store the resampled ivar spectra
        resampled_ivar = np.empty_like(ivar_array)

        # Interpolate each var spectrum back to the original st_wave
        for i in range(len(ivar_array)):
            interp_function_resampled = interp1d(restframe_wavelengths[i], ivar_array[i], kind='linear', bounds_error=False, fill_value=np.nan)
            resampled_ivar[i] = 1.0 / interp_function_resampled(st_wave)

        return resampled_flux, resampled_ivar
