"""
Descarga, controla y separa en .nc las componentes del SAM en 700hPa
"""
################################################################################
sam_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/sam_monthly.csv'
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
update = True
test = False
sam_component_name = ['sam', 'ssam', 'asam']
# descarga #####################################################################
if update:
    try:
        url = 'https://www.cima.fcen.uba.ar/~elio.campitelli/asymsam/data/' \
              'sam_monthly.csv'
        os.system('wget -O /pikachu/datos/luciano.andrian/verif_2019_2023/' 
                  'sam_monthly.csv ' + url)
    except:
        print('Error al descargar sam_monthly.csv')
# apertura #####################################################################
try:
    data = pd.read_csv(sam_dir, sep=',', header=0)
except:
    print('-----------------------')
    print('el directorio no existe')
    sys.exit(1)
#------------------------------------------------------------------------------#
data['time'] = pd.to_datetime(data['time'])
data = data.loc[data.time.dt.year>=2019]

print('# Control #############################################################')
skip_anios = []
l_count = 0
j_count = 0
for s in sam_component_name:
    sam_component = data.loc[(data['index'] == s) &
                             (data['lev'] == level)]

    for i in range(0, data.time.dt.year[0]-2019):
        anioi = 2019 + i
        aniof = 2020 + i
        aux = sam_component.time.loc[(data.time.dt.year >= anioi) &
                                     (data.time.dt.year < aniof)]
        l = len(aux)

        # Todos los años tienen 12 meses?
        if l != 12:
            l_count = 1
            skip_anios.append(aniof)

        else:
            # Los que tienen los 12 meses
            # Estan en orden?
            for j in range(0, 11):
                m0 = pd.to_datetime(aux.values).month[j]
                m1 = pd.to_datetime(aux.values).month[j + 1]
                if m0 - m1 != 1:
                    j_count = 1
                    print(str(m0) + ' - ' + str(m1))
                    skip_anios.append(aniof)

if (l_count == 0) & (j_count == 0):
    print('Control OK')
    print('###################################################################')
    print('Creando NCs...')

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
else:
    print('Años incompletos: ' + str(skip_anios) )
    print('NCs desactualizados...')
################################################################################
print('#######################################################################')
print('done')
print('out_dir = ' + out_dir )
print('#######################################################################')
################################################################################