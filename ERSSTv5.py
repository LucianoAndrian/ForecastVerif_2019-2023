"""
Descarga ERRSSTv5 en .nc
"""
################################################################################
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'

# data from
# https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5
################################################################################
import os
################################################################################
update = True
# descarga #####################################################################
if update:
    try:
        url = 'https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/' \
              'sst.mnmean.nc'
        local_path = '/pikachu/datos/luciano.andrian/verif_2019_2023/' \
                     'sst.mnmean.nc'

        os.system(f'wget -O {local_path} {url}')
    except:
        print('Error al descargar ERSSTv5')

################################################################################
print('#######################################################################')
print('done')
################################################################################