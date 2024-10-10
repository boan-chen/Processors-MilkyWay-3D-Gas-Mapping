from scipy import ndimage
from scipy import optimize as opt
import numpy as np
from mpi4py import MPI
from tqdm import tqdm


class continuum:
    def __init__(self, st_wave, flux_array, ivar_array = None):
        self.continuum_array = np.zeros_like(flux_array)
        self.flux_array = flux_array
        self.ivar_array = ivar_array
        self.st_wave = st_wave
        if self.ivar_array is not None:
            self.ivar_array = ivar_array

    def polynomial_func(self, x, *coeffs):
        return np.polyval(coeffs, x)
    
    def curve_fit(self, wave, flux, ivar = None, degree = 7):
        initial_guess = np.ones(degree + 1)
        if ivar is None:
            fit_params, _ = opt.curve_fit(self.polynomial_func, wave, flux, p0=initial_guess)
        else:
            fit_params, _ = opt.curve_fit(self.polynomial_func, wave, flux, p0=initial_guess, sigma=1/np.sqrt(ivar))
        predicted_spectrum = self.polynomial_func(self.st_wave, *fit_params)
        return predicted_spectrum
    
    def normalize(self, flux, to_filter = True, filter_size = 25, percentile = 90):
        if to_filter == True:
            filtered = ndimage.percentile_filter(flux, percentile = percentile, size = filter_size)
            finite = np.isfinite(filtered)
            filtered, wave = filtered[finite], self.st_wave[finite]
            continuum = self.curve_fit(wave, filtered)
        else:
            continuum = self.curve_fit(self.st_wave, flux)
        return continuum
    
    def norm_all_flux(self, to_filter = True, filter_size = 25, percentile = 90):
        print("Normalizing all fluxes...")
        for i in tqdm(range(len(self.flux_array))):
            self.continuum_array[i] = self.normalize(self.flux_array[i], to_filter = to_filter, filter_size = filter_size, percentile = percentile)
        self.norm_array = self.flux_array / self.continuum_array
        return self.norm_array
    
    def norm_all_flux_parallel(self, to_filter = True, filter_size = 25, percentile = 90):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        for i in range(rank, len(self.flux_array), size):
            self.continuum_array[i] = self.normalize(self.flux_array[i], to_filter = to_filter, filter_size = filter_size, percentile = percentile)
        self.norm_array = self.flux_array / self.continuum_array
        return self.norm_array
    
    def ivar_propagate(self):
        self.ivar_array = self.ivar_array * self.continuum_array**2
        return self.ivar_array