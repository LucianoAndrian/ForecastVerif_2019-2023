"""
Descarga y procesado de CHIRPS
"""
################################################################################
import os
import glob
import xarray as xr
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
# https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/netcdf/p25/
################################################################################
download = False
update = True
proc = True
proc_2019_2023 = True

################################################################################
if download:
    print('DOWNLOAD = TRUE!!!')
    for y in range(1981, 2024):
        url = 'https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/' \
              'netcdf/p25/chirps-v2.0.' + str(y) + '.days_p25.nc'
        os.system('wget -O /pikachu/datos/luciano.andrian/verif_2019_2023/'
                  'chirps/ch_' + str(y) + '.nc ' + url)

    print('done')

#------------------------------------------------------------------------------#
if update:
    print('UPDATE = TRUE!!!')
    for y in range(2023, 2024):
        url = 'https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/' \
              'netcdf/p25/chirps-v2.0.' + str(y) + '.days_p25.nc'
        os.system('wget -O /pikachu/datos/luciano.andrian/verif_2019_2023/'
                  'chirps/ch_' + str(y) + '.nc ' + url)

    print('done')

#------------------------------------------------------------------------------#
if proc:
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
if proc_2019_2023:
    print('proc_2019_2023 = True')
    files = glob.glob('/pikachu/datos/luciano.andrian/verif_2019_2023/'
                      'chirps/ch_*.nc')
    files = sorted(files, key=lambda x: x.split()[0])
    files = files[-5:]
    data = xr.open_mfdataset(files)

    data = data.resample(time='1MS').mean('time')
    data.to_netcdf('/pikachu/datos/luciano.andrian/verif_2019_2023/'
                      'chirps/chirps_2019_2023_mmean.nc')



# CMAP #########################################################################
url = 'https://downloads.psl.noaa.gov/Datasets/cmap/enh/precip.mon.mean.nc'
local_path = '/pikachu/datos/luciano.andrian/verif_2019_2023/cmap/pp_cmap.nc'

os.system(f'wget -O {local_path} {url}')
################################################################################
print('done')
################################################################################