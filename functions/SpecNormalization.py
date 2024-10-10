#%%
import numpy as np
import scipy.stats as ss
from scipy import ndimage
from scipy.signal import find_peaks
from tqdm import tqdm
#%%
def clipped_polyfit(st_wave, sample_spec):
    wave = st_wave
    spec = sample_spec
    wave = wave[~np.isnan(spec)]
    spec = spec[~np.isnan(spec)]
    for i in range(0, 3):
        fit = np.polyfit(wave, spec, 7)
        fitted = np.polyval(fit, wave)
        error = (spec - fitted)/fitted
        idx = np.where(abs(ss.zscore(error)) > 3)[0]
        if len(idx) == 0:
            break
        wave = np.delete(wave, idx)
        spec = np.delete(spec, idx)
    # fitted = np.polyval(fit, st_wave)
    filtered = ndimage.percentile_filter(spec, percentile = 50, size = 51)
    continuum = np.interp(st_wave, wave, filtered)
    return continuum

def normalization(st_wave, sample_spec, ivar):
    continuum = clipped_polyfit(st_wave, sample_spec)
    norm = sample_spec / continuum
    std = 1/np.sqrt(ivar)
    normed_std = std / continuum
    # turn np.nan to 0
    norm[np.isnan(norm)] = 0
    normed_std[np.isnan(normed_std)] = 0
    return norm, normed_std

def poly_normalization(st_wave, spec, ivar, tqdm_disable = False):
    normed = np.empty_like(spec)
    normed_std = np.empty_like(spec)
    print('normalizing spectra...')
    for i in tqdm(range(0, len(spec)), disable = tqdm_disable):
        normed[i], normed_std[i] = normalization(st_wave, spec[i], ivar[i])
    return normed, normed_std

def bias_filter(wave, flux, window_size=25, max_iterations=2, iteration=0, peak=True):
    if np.all(flux == 0):
        return np.repeat(np.nan, len(wave))
    flux_empty = np.zeros(len(wave))
    for i in range(0, len(flux) - window_size, 5):
        lFlux = flux[i: i + window_size]
        hst = np.histogram(lFlux, bins=50)
        find = max(np.where(hst[0] >= np.quantile(hst[0], 0.95))[0])
        st_bin = hst[1][find - 1 : find + 1]
        selected = np.where((lFlux >= st_bin[0]) & (lFlux <= st_bin[-1]))[0]
        flux_empty[i + selected] = lFlux[selected]
    filtered_flux = ndimage.percentile_filter(flux_empty[flux_empty != 0], 50, 25)
    if np.any(flux_empty != 0):
        filtered_flux = np.interp(wave, wave[flux_empty != 0], filtered_flux)
    if peak and iteration > 1:
        peak_flux = find_peaks(flux_empty[wave < 5500], height=0.95)[0]
        local_filter = flux[peak_flux]
        if len(local_filter) > 0:
            local_filter = np.interp(wave[wave < 5500], wave[peak_flux], local_filter)
            filtered_flux[wave < 5500] = (local_filter + filtered_flux[wave < 5500]) / 2
    else:
        pass
    new_norm = flux / filtered_flux
    if np.array_equal(new_norm, flux) or iteration >= max_iterations:
        return new_norm
    else:
        return bias_filter(wave, new_norm, window_size, max_iterations, iteration + 1, peak=peak)

def filtered_normalization(wave, flux, window_size=25, max_iterations=3, iteration=0, peak = True):
    continuum = clipped_polyfit(wave, flux)
    continuum[continuum <= 0] = 1e-4
    norm = flux / continuum
    flux_norm = bias_filter(wave, norm, window_size, max_iterations, iteration, peak = peak)
    return flux_norm
    
# %%
