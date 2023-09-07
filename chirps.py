"""
Descarga y procesado de CHIRPS
"""
################################################################################
import os
import glob
import xarray as xr
from datetime import datetime
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
# https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/netcdf/p25/
################################################################################
download_chirps_full = False
update = True
proc = True
proc_2019_2023 = True
anio = datetime.now().year
################################################################################
def download_full():
    print('DOWNLOAD = TRUE!!!')
    for y in range(1981, 2024):
        url = 'https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/' \
              'netcdf/p25/chirps-v2.0.' + str(y) + '.days_p25.nc'
        os.system('wget -O /pikachu/datos/luciano.andrian/verif_2019_2023/'
                  'chirps/ch_' + str(y) + '.nc ' + url)

    print('done download')

    print('proc = True')
    files = glob.glob('/pikachu/datos/luciano.andrian/verif_2019_2023/'
                      'chirps/ch_*.nc')
    files = sorted(files, key=lambda x: x.split()[0])

    data = xr.open_mfdataset(files)
    data = data.sel(time=slice('1990-01-01', '2020-12-01'))
    data = data.resample(time='1MS').mean('time')
    data.to_netcdf('/pikachu/datos/luciano.andrian/verif_2019_2023/'
                   'chirps/chirps_1990_2020_mmean.nc')

#------------------------------------------------------------------------------#
def update():
    print('UPDATE = TRUE!!!')
    y = anio
    url = 'https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/' \
          'netcdf/p25/chirps-v2.0.' + str(y) + '.days_p25.nc'
    os.system('wget -O /pikachu/datos/luciano.andrian/verif_2019_2023/'
              'chirps/ch_' + str(y) + '.nc ' + url)

    print('done download')

    print('proc_2019_2023 = True')
    files = glob.glob('/pikachu/datos/luciano.andrian/verif_2019_2023/'
                      'chirps/ch_*.nc')

    files = sorted(files, key=lambda x: x.split()[0])
    if anio != 2024:
        files = files[-6:]
        data = xr.open_mfdataset(files)

        data = data.resample(time='1MS').mean('time')
        data.to_netcdf('/pikachu/datos/luciano.andrian/verif_2019_2023/'
                       'chirps/chirps_2019_2023_mmean.nc')
        print('Done proc')
    else:
        print('CORREGIR files[-6:]!!!')

#------------------------------------------------------------------------------#
if __name__ == "__main__":
    update()
################################################################################