"""
Descarga, ordena y seleccion el indice ONI
"""
################################################################################
oni_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/oni.data'
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'
# data from https://psl.noaa.gov/data/correlation/oni.data
################################################################################
import pandas as pd
import xarray as xr
import datetime as datet
from datetime import datetime
import sys
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
################################################################################
update = True
# descarga #####################################################################
if update:
    try:
        url = 'https://psl.noaa.gov/data/correlation/oni.data'
        local_path = '/pikachu/datos/luciano.andrian/verif_2019_2023/oni.data'

        os.system(f'wget -O {local_path} {url}')
    except:
        print('Error al descargar oni.data')
# apertura #####################################################################
cu_year = datet.datetime.now().year
n_filas = cu_year - 1950 + 1
try:
    aux_data = pd.read_csv(oni_dir, sep='\s+', header=None, skiprows=1,
                       nrows=n_filas, engine='python')
except:
    print('-----------------------')
    print('el directorio no existe')
    sys.exit(1)

#------------------------------------------------------------------------------#
aux_data = aux_data.loc[aux_data[aux_data.columns[0]]>=2019]

data = []
for index, row in aux_data.iterrows():
    year = int(row[0])
    for month in range(1, 13):
        data.append([datetime(year, month, 1), row[month]])

data = pd.DataFrame(data, columns=['time', 'oni'])
data = data.loc[data.index>-10]
data = data[data['oni'] !=-99.90]
#------------------------------------------------------------------------------#
data = data.to_xarray()
data.to_netcdf(out_dir + 'oni.nc')

################################################################################
print('#######################################################################')
print('done')
print('out_dir = ' + out_dir)
print('#######################################################################')
################################################################################