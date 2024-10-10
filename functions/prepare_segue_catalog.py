#%%
import pandas as pd
import astropy.io.fits as pf

catalog = pf.open('/Volumes/NTU_Astro/Catalog/' + 'gaia_matched_catalog.fits')
catalog = catalog[1].data

teff = catalog['TEFF_ANNRR']
teff_err = catalog['TEFF_ANNRR_UNC']
logg = catalog['LOGG_ANNRR']
logg_err = catalog['LOGG_ANNRR_UNC']
feh = catalog['FEH_ANNRR']
feh_err = catalog['FEH_ANNRR_UNC']
parallax = catalog['parallax']
parallax_err = catalog['parallax_error']
rv = catalog['RV_ADOP']
RV_flag = catalog['RV_FLAG']
RV_error = catalog['RV_ADOP_UNC']
snr = catalog['SNR']
l = catalog['l']
b = catalog['b']
dist_sdss = catalog['DIST_ADOP']
segue1_target1 = catalog['SEGUE1_TARGET1']
segue1_target2 = catalog['SEGUE1_TARGET2']
segue2_target1 = catalog['SEGUE2_TARGET1']
segue2_target2 = catalog['SEGUE2_TARGET2']
program = catalog['PROGRAMNAME']
primary_target = catalog['PRIM_TARGET']
spec_type_hammer = catalog['SPECTYPE_HAMMER']
spec_type_subclass = catalog['SPECTYPE_SUBCLASS']
plate, fiber, mjd = catalog['plate_1'], catalog['fiber_1'], catalog['mjd_1']
dist_z = catalog['DIST_Z']
dist_AP = catalog['DIST_AP']
# check this website: https://data.sdss.org//datamodel/files/SSPP_REDUX/RERUN/PLATE4/output/param/ssppOut.html

# define spec_name
def output_specname(plate, fiber, mjd):
    spec_name = []
    for i in range(len(plate)):
        spec_name.append(f'{plate[i]:0>4}/spec-{plate[i]:0>4}-{mjd[i]:0>5}-{fiber[i]:0>4}.fits')
    return spec_name

spec_name = output_specname(plate, fiber, mjd)

catalog_light = pd.DataFrame(columns = [], index=range(len(catalog)))
catalog_light['teff'] = teff
catalog_light['logg'] = logg
catalog_light['feh'] = feh
catalog_light['teff_err'] = teff_err
catalog_light['logg_err'] = logg_err
catalog_light['feh_err'] = feh_err
catalog_light['parallax'] = parallax
catalog_light['parallax_err'] = parallax_err
catalog_light['RV'] = rv
catalog_light['RV_flag'] = RV_flag
catalog_light['RV_error'] = RV_error
catalog_light['SNR'] = snr
catalog_light['l'] = l
catalog_light['b'] = b
catalog_light['dist_adop'] = dist_sdss
catalog_light['dist_z'] = dist_z
catalog_light['segue1_target1'] = segue1_target1
catalog_light['segue1_target2'] = segue1_target2
catalog_light['segue2_target1'] = segue2_target1
catalog_light['segue2_target2'] = segue2_target2
catalog_light['program'] = program
catalog_light['primary_target'] = primary_target
catalog_light['spec_type_hammer'] = spec_type_hammer
catalog_light['spec_type_subclass'] = spec_type_subclass
catalog_light['ra'] = catalog['RA_1']
catalog_light['dec'] = catalog['DEC_1']
catalog_light['spec_name'] = spec_name

#save the catalog
catalog_light.to_csv('/Volumes/NTU_Astro/Catalog/' + 'gaia_matched_catalog_light.csv', index=False)
# %%
