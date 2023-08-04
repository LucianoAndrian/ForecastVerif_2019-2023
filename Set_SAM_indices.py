sam_dir = '/home/luciano.andrian/verif_2019_2023/sam_monthly.csv'
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'

# data from
# https://www.cima.fcen.uba.ar/~elio.campitelli/asymsam/data/sam_monthly.csv
################################################################################
import pandas as pd
import xarray as xr
import sys
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
################################################################################
level = 700 #
sam_component_name = ['sam', 'ssam', 'asam']
test = False
################################################################################
try:
    data = pd.read_csv(sam_dir, sep=',', header=0)
except:
    print('-----------------------')
    print('el directorio no existe')
    sys.exit(1)
#------------------------------------------------------------------------------#

data['time'] = pd.to_datetime(data['time'])
data = data.loc[data.time.dt.year>=2019]

dates = xr.cftime_range(start='2019-01-01',
                        end=str(data.time[0]).split(' ')[0],
                        freq='MS')

for s in sam_component_name:
    sam_component = data.loc[(data['index'] == s) &
                             (data['lev'] == level)]

    sam_component = sam_component.drop(columns=['lev', 'index', 'time'])
    sam_component = sam_component.iloc[::-1]
    sam_component = sam_component.to_xarray()

    sam_component['index'] = dates
    sam_component = sam_component.rename({'index': 'time'})

    if test:
        print(s + '--------------------------')
        print(sam_component)
    else:
        sam_component.to_netcdf(out_dir + s + '.nc')

################################################################################