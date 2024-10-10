# Description: This file contains functions for visualization of the data
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
import healpy as hp

def visualize_stellar_count(catalog, nside, scale=[0, 100], title=None):
    # Counting the occurrences of each 'pix_idx' in the catalog
    pix_counts = np.zeros(hp.nside2npix(nside))
    for pix_idx in catalog['pix_idx']:
        pix_counts[pix_idx] += 1
    pix_series = np.maximum(pix_counts, 0)  # Ensuring no negative counts
    pix_series[pix_series == 0] = hp.UNSEEN  # Masking pixels with zero stars
    hp.projview(hp.pixelfunc.ma(pix_series, badval = hp.UNSEEN), \
            title = title,\
            graticule=True, graticule_labels=True, projection_type="mollweide",\
            cb_orientation="horizontal", alpha = 1, min = scale[0], max = scale[1])
    plt.show()
    return 

def spectrum_visualization(
    wave,
    flux,  
    flux_err = None,
    xlim = None,
    ylim = None,
    title = None,
    save = False,
    directory = None
):
    plt.figure(figsize=(10, 5))
    plt.plot(wave, flux, label='flux')
    if flux_err is not None:
        # plot the error as shaded region
        plt.fill_between(wave, flux - flux_err, flux + flux_err, alpha=0.2, label='flux_err')
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux')
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.title(title)
    if save == True:
        plt.savefig(f"{directory}{title}.png")
    plt.show()
    return
    
def spectra_visualization(
    wave,
    flux_list,
    flux_err_list = None,
    xlim = None,
    ylim = None,
    title = None,
    label_list = None,
    save = False, 
    directory = None
):
    plt.figure(figsize=(10, 5))
    for i in range(0, len(flux_list)):
        if flux_err_list is not None:
            flux_err = flux_err_list[i]
            plt.fill_between(wave, flux_list[i] - flux_err, flux_list[i] + flux_err, alpha=0.2, label=label_list[i])
        plt.plot(wave, flux_list[i], alpha = 0.5, label=label_list[i])
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux')
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.legend()
    plt.title(title)
    if save == True:
        plt.savefig(f"{directory}{title}.png")
    plt.show()
    return
    
    
    
