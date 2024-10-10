import numpy as np
import astropy.io.fits as pf
import pandas as pd
import os.path
import os
import requests

def plate_dir(catalog, path):
    plate_name, fiber_name, mjd_name = 'plate_1', 'fiber_1', 'mjd_1'
    base_dir = path['stellar_spectra']
    plates = list(set(catalog[plate_name]))
    for i in range(0,len(plate)):
        #print('tes')
        plate = '%04d' % (plates[i])
        plate_dir = base_dir + plate
        if os.path.isdir(plate_dir) == True:
            pass
        else:
            os.system("mkdir "+plate_dir)
    print('SpecReader is ready.')

def download_request(catalog, i, path, enforce=False):
    base_dir = path['stellar_spectra']
    if isinstance(catalog, pd.DataFrame):
        data_format = 'pandas'
    elif isinstance(catalog, np.ndarray):
        data_format = 'numpy'
    elif isinstance(catalog, pf.hdu.hdulist.HDUList):
        data_format = 'fits'
    else:
        print("Unsupported catalog format.")
        return None
    
    if data_format not in ['pandas', 'numpy', 'fits']:
        print("Unsupported catalog format.")
        return None

    if data_format == 'pandas':
        plate = catalog['spec_name'][i][:4]
        filename = catalog['spec_name'][i][5:]
    else:
        plate_name, fiber_name, mjd_name = 'plate_1', 'fiber_1', 'mjd_1'
        plate = '%04d' % (catalog[plate_name][i])
        filename = f'spec-{plate}-{catalog[mjd_name][i]}-{catalog[fiber_name][i]}.fits'

    if not enforce and os.path.exists(os.path.join(base_dir, plate, filename)):
        return

    plate_int = int(plate)
    plate_versions = [
        ('https://data.sdss.org/sas/dr17/sdss/spectro/redux/26/spectra/lite/', '26'),
        ('https://data.sdss.org/sas/dr17/sdss/spectro/redux/104/spectra/lite/', '104'),
        ('https://data.sdss.org/sas/dr17/sdss/spectro/redux/103/spectra/lite/', '103'),
    ]

    if plate_int < 3523:
        plate_versions.reverse()

    for version_url, version in plate_versions:
        dir_download = f'{version_url}{plate}/{filename}'
        response = requests.get(dir_download)

        if response.status_code == 200:
            os.makedirs(os.path.join(base_dir, plate), exist_ok=True)
            with open(os.path.join(base_dir, plate, filename), 'wb') as f:
                f.write(response.content)
            return

    print(f'{i}th file cannot be downloaded.')
    return


def download(catalog, i, path, enforce = False):
    base_dir = path['stellar_spectra']
    if isinstance(catalog, pd.DataFrame):
        data_format = 'pandas'
    elif isinstance(catalog, np.ndarray):
        data_format = 'numpy'
    elif isinstance(catalog, pf.hdu.hdulist.HDUList):
        data_format = 'fits'
    else:
        print("Unsupported catalog format.")
        return None
    
    if data_format == 'pandas':
        plate = catalog['spec_name'][i][:4]
        filename = catalog['spec_name'][i][5:]
    elif data_format == 'numpy' or data_format == 'fits':
        plate_name, fiber_name, mjd_name = 'plate_1', 'fiber_1', 'mjd_1'
        plate = '%04d' % (catalog[plate_name][i])
        if data_format == 'fits':
            filename = 'spec-%04d-%05d-%04d.fits' % (catalog[1].data[plate_name][i],catalog[1].data[mjd_name][i],catalog[1].data[fiber_name][i])
        else:
            filename = 'spec-%04d-%05d-%04d.fits' % (catalog[plate_name][i],catalog[mjd_name][i],catalog[fiber_name][i])
    else:
        print("Unsupported catalog format.")
        return None

    plate_int = int(plate)
    
    if os.path.isfile(base_dir + plate + '/' + filename)==True or enforce == True:
        os.system("rm "+base_dir + plate + '/' + filename)
    
    if os.path.isfile(base_dir + plate + '/' + filename)==False or enforce == True:
        if (plate_int < 3523):
            dir_short = 'https://data.sdss.org/sas/dr17/sdss/spectro/redux/26/spectra/lite/'		
            dir_download =dir_short + plate + '/' + filename
            os.system("wget "+dir_download)
            if os.path.isfile(filename)==True:		
                os.system("mv "+filename+' '+base_dir+plate)
                return 
            else:
                dir_short = 'https://data.sdss.org/sas/dr17/sdss/spectro/redux/104/spectra/lite/'
                dir_download =dir_short + plate + '/' + filename

                os.system("wget "+dir_download)
                if os.path.isfile(filename)==True:	
                    os.system("mv "+filename+' '+base_dir+plate)
                    return
                else:	
                    dir_short = 'https://data.sdss.org/sas/dr17/sdss/spectro/redux/103/spectra/lite/'
                    dir_download =dir_short + plate + '/' + filename

                    os.system("wget "+dir_download)
                    if os.path.isfile(filename)==True:	
                        os.system("mv "+filename+' '+base_dir+plate)
                        return
        else:
            dir_short = 'https://data.sdss.org/sas/dr17/sdss/spectro/redux/v5_13_2/spectra/lite/'
            dir_download =dir_short + plate + '/' + filename
            os.system("wget "+dir_download)
            if os.path.isfile(filename)==True:		
                os.system("mv "+filename+' '+base_dir+plate)
                return
    print(f'{i} th file cannot be downloaded.')
    return 0
    