"""
Descarga CMAP
"""
################################################################################
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
################################################################################
def update():
    url = 'https://downloads.psl.noaa.gov/Datasets/cmap/enh/precip.mon.mean.nc'
    local_path = '/pikachu/datos/luciano.andrian/verif_2019_2023/cmap/' \
                'pp_cmap.nc'
    os.system(f'wget -O {local_path} {url}')

if __name__ == "__main__":
    update()