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
# descarga #####################################################################
def update():
    current_month = datetime.now().month
    data_last_trimester = xr.open_dataset(out_dir + 'oni.nc').time.values[-1].\
                              astype('datetime64[M]').astype(int) % 12 + 1

    # ONI actualizado hasta el el trimestre que conforma el mes anterior
    # al actual como 3er mes. Si es Sep, ONI hasta JJA
    if (current_month - 2) != data_last_trimester:
        print('ONI desactualizado')
        print('Actualizando y seteando ONI')

        try:
            url = 'https://psl.noaa.gov/data/correlation/oni.data'
            local_path = '/pikachu/datos/luciano.andrian/verif_2019_2023/' \
                         'oni.data'

            os.system(f'wget -O {local_path} {url}')
        except:
            print('Error al descargar oni.data')

         # apertura ############################################################

        cu_year = datet.datetime.now().year
        n_filas = cu_year - 1950 + 1
        try:
            aux_data = pd.read_csv(oni_dir, sep='\s+', header=None, skiprows=1,
                               nrows=n_filas, engine='python')
        except:
            print('-----------------------')
            print('el directorio no existe')
            sys.exit(1)

        #----------------------------------------------------------------------#
        aux_data = aux_data.loc[aux_data[aux_data.columns[0]]>=2019]

        data = []
        for index, row in aux_data.iterrows():
            year = int(row[0])
            for month in range(1, 13):
                data.append([datetime(year, month, 1), row[month]])

        data = pd.DataFrame(data, columns=['time', 'oni'])
        data = data.loc[data.index>-10]
        data = data[data['oni'] !=-99.90]
        #----------------------------------------------------------------------#
        data = data.to_xarray()
        data.to_netcdf(out_dir + 'oni.nc')

################################################################################
if __name__ == "__main__":
    update()
################################################################################