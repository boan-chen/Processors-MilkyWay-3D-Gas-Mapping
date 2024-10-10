import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from astropy.io import ascii
from astropy.modeling import models, fitting
import scipy.stats as ss
import statistics 
import healpy as hp
from healpy.newvisufunc import projview, newprojplot
from astropy.modeling.models import Polynomial1D
from healpy.newvisufunc import projview, newprojplot
# import wget


from astropy import units as u
from astropy.coordinates import SkyCoord
import os
import time
from scipy.optimize import curve_fit
import csv
from sklearn.decomposition import PCA, IncrementalPCA
from scipy.signal import medfilt
from scipy import ndimage
from scipy.signal import find_peaks

CaII = [3934.78, 3969.59]
NaD = [5891.58, 5897.56]


def composite(A1):
    spectrum_median = np.zeros(len(A1[0]))
    for p in range(0, len(A1[0])):
        spectrum_median[p] = np.nanmedian(A1[:, p])
    return spectrum_median

def equivelent_width(A, sigma) :
    W0 = A * sigma * np.sqrt(2 * np.pi)
    return W0

def dbl_gauss(x, A1, mu1, sigma1, A2, mu2, sigma2, z):
    return A1*np.exp(-(x-mu1)**2/(2.*sigma1**2)) + A2*np.exp(-(x-mu2)**2/(2.*sigma2**2)) + z
def dbl_gauss_Ca(x, mu1, A1, sigma1, A2, sigma2, z):
    mu2 = CaII[1] * mu1/CaII[0] 
    return A1*np.exp(-(x-mu1)**2/(2.*sigma1**2)) + A2*np.exp(-(x-mu2)**2/(2.*sigma2**2)) + z
def dbl_gauss_Na(x, mu1, A1, sigma1, A2, sigma2, z):
    mu2 = NaD[1] * mu1/NaD[0] 
    return A1*np.exp(-(x-mu1)**2/(2.*sigma1**2)) + A2*np.exp(-(x-mu2)**2/(2.*sigma2**2)) + z

def dbl_gauss_fitting_Atoms(A1, obs_wave, error = None):
    Ca_X = np.where(abs(obs_wave - 3957) < 50)
    Ca_x = obs_wave[Ca_X[0]]
    Ca_error = error[Ca_X[0]]
    Ca_y = A1[Ca_X[0]]
    Na_X = np.where(abs(obs_wave - 5895) < 50)
    Na_x = obs_wave[Na_X[0]]
    Na_error = error[Na_X[0]]
    Na_y = A1[Na_X[0]]
    if np.max(A1) == 0:
        popt = np.zeros(6)
        popt1 = np.zeros(6)
        fitted = np.zeros(len(obs_wave))
    else:
        try:
            p0 = [3934, -0.2, 1, -0.2, 1, 1]
            popt = p0
            popt, pcov = curve_fit(dbl_gauss_Ca, Ca_x, Ca_y, bounds=([3932, -1, 0, -1, 0, 0.99], [3936, 0.2, 5, 0.2, 5, 1.05]), sigma = Ca_error, max_nfev=5100)
            p1 = [5891, -0.1, 1, -0.1, 1, 1]
            popt1 = p1
            popt1, pcov = curve_fit(dbl_gauss_Na, Na_x, Na_y, bounds=([5887, -1, 0, -1, 0, 0.995], [5895, 0.2, 4, 0.2, 4, 1.05]), sigma = Na_error, max_nfev=5100)
            fitted = np.ones(len(obs_wave))
            fitted[Ca_X[0]] = dbl_gauss_Ca(Ca_x, *popt)
            fitted[Na_X[0]] = dbl_gauss_Na(Na_x, *popt1)
        except RuntimeError:
            popt = np.zeros(6)
            popt1 = np.zeros(6)
            fitted = np.zeros(len(obs_wave))
            print('cannot fit')
    return popt, popt1, fitted

def R_EW_Atoms(v, r):
    W = np.zeros((1, 4))
    W[0, 0] = - (equivelent_width(v[1], v[2]))
    W[0, 1] = - (equivelent_width(v[3], v[4]))
    W[0, 2] = - (equivelent_width(r[1], r[2]))
    W[0, 3] = - (equivelent_width(r[3], r[4]))
    return W

def bootstrap(A2):
    L = len(A2)
    temp = np.zeros((L, len(A2[0, :])))
    if L > 0 :
        random = ss.randint.rvs(0, L, size = L)
        for r in range(0, L):
            temp[r, :] = A2[random[r], :]
    return temp