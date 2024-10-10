# PreProcessor Module Documentation

## Overview
The PreProcessor module contains various scripts and functions for preprocessing astronomical data. This includes downloading data, reading and saving data, normalizing spectra, and performing various analyses. This version is last updated in 2022 and is partially extracted from my private work that involves model details.

## Files and Functions

### **SpecNormalization.py**
This file contains functions for normalizing spectral data.

- **`clipped_polyfit(st_wave, sample_spec)`**:  
  Performs a polynomial fit with clipping.
  
- **`normalization(st_wave, sample_spec, ivar)`**:  
  Normalizes the spectrum.
  
- **`poly_normalization(st_wave, spec, ivar, tqdm_disable=False)`**:  
  Applies polynomial normalization to the spectrum.
  
- **`bias_filter(wave, flux, window_size=25, max_iterations=2, iteration=0, peak=True)`**:  
  Applies a bias filter to the spectrum.
  
- **`filtered_normalization(wave, flux, window_size=25, max_iterations=3, iteration=0, peak=True)`**:  
  Applies filtered normalization to the spectrum.

### **SpecDownload.py**
This file contains functions for downloading spectral data.

- **`plate_dir(catalog, path)`**:  
  Creates a directory for the plate data.
  
- **`download_request(catalog, i, path, enforce=False)`**:  
  Sends a download request for the data.
  
- **`download(catalog, i, path, enforce=False)`**:  
  Downloads the data.

### **SpecReader.py**
This file contains the `SpecReader` class for reading and saving spectral data.

- **`SpecReader.__init__(self, catalog, st_wave, path, directory=None)`**:  
  Initializes the `SpecReader` object.
  
- **`SpecReader.read(self, i)`**:  
  Reads the spectral data.
  
- **`SpecReader.common_pix(self, log_wave, flux, ivar, model)`**:  
  Finds common pixels in the data.
  
- **`SpecReader.load_data(self, loop=None, tqdm_disable=False)`**:  
  Loads the spectral data.
  
- **`SpecReader.save_data(self, catalog, flux_array, ivar_array, model_array)`**:  
  Saves the spectral data.

### **sky_analysis.py**
This file contains functions for analyzing the sky data.

- **`open_model(idx)`**:  
  Opens a model based on the index.
  
- **`catalog = pd.read_csv(dir_dic['catalog'] + 'boan_standard_catalog_ver2.csv')`**:  
  Loads the catalog data.

### **SpecTools.py**
This file contains various tools for spectral analysis.

- **`composite(A1)`**:  
  Creates a composite spectrum.
  
- **`equivelent_width(A, sigma)`**:  
  Calculates the equivalent width.
  
- **`dbl_gauss(x, A1, mu1, sigma1, A2, mu2, sigma2, z)`**:  
  Fits a double Gaussian model.
  
- **`dbl_gauss_Ca(x, mu1, A1, sigma1, A2, sigma2, z)`**:  
  Fits a double Gaussian model for Calcium.
  
- **`dbl_gauss_Na(x, mu1, A1, sigma1, A2, sigma2, z)`**:  
  Fits a double Gaussian model for Sodium.
  
- **`dbl_gauss_fitting_Atoms(A1, obs_wave, error=None)`**:  
  Fits a double Gaussian model to atomic lines.
  
- **`R_EW_Atoms(v, r)`**:  
  Calculates the equivalent width ratio.
  
- **`bootstrap(A2)`**:  
  Performs bootstrap resampling.

### **model_construction.py**
This file contains functions for constructing ISM models from the spectral data.

- **`model(local_catalog, grid, threshold=30)`**:  
  Constructs a model from the local catalog.
  
- **`parallel_model(args)`**:  
  Constructs models in parallel.

### **continuum.py**
This file contains the `continuum` class for continuum fitting.

- **`continuum`**:  
  Class for continuum fitting.

### **catalog_process.py**
This file contains functions for processing the catalog data.

- **`dust_map(nside, catalog, sky_map)`**:  
  Creates a dust map.
  
- **`catalog_by_criteria(criteria, catalog)`**:  
  Filters the catalog based on criteria.

### **server_shifting.py**
This file contains functions for server shifting.

- **`server()`**:  
  Initializes the server.

## Usage
To use the functions and classes in the PreProcessor module, import the necessary files and call the functions with the appropriate parameters. For example:

```python
from PreProcessor.SpecNormalization import poly_normalization
from PreProcessor.SpecReader import SpecReader

# Initialize SpecReader
spec_reader = SpecReader(catalog, st_wave, path)
# Load data
spec_reader.load_data()
# Normalize spectrum
normalized_spectrum = poly_normalization(st_wave, spec, ivar)
```