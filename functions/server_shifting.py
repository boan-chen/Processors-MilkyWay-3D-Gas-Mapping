import os
def server():
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
    return serv, server, dir_dic