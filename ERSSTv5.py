"""
Descarga ERRSSTv5 en .nc
"""
################################################################################
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/'

# data from
# https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5
################################################################################
from datetime import datetime
import xarray as xr
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
################################################################################
def download(file):
    try:
        url = 'https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/' \
              'sst.mnmean.nc'
        os.system(f'wget -O {file} {url}')
    except:
        print('Error al descargar ERSSTv5')

# Update #######################################################################
def update():
    current_month = datetime.now().strftime('%m')
    file = out_dir + 'sst.mnmean.nc'
    file_month = datetime.utcfromtimestamp(
        os.path.getmtime(file)).strftime('%m')

    if current_month != file_month:
        print('ERSSTv5 desactualizado')
        print('Descargando versión más reciente')

        download(file)

    elif current_month == file_month:
        # normalmente con esta condición ya alcanzaria para verificar
        # que está actualizado
        data_last_month = xr.open_dataset(file).time.values[-1].\
                              astype('datetime64[M]').astype(int) % 12 + 1

        # en caso de no haberse actualizado pese al cambio de mes
        if data_last_month != int(current_month) - 1:
            # el último mes actualizado tiene que ser el anterior
            # al mes actual
            print('ERSSTv5 desactualizado')
            print('Descargando versión más reciente')
            download(file)

################################################################################
if __name__ == "__main__":
    update()
################################################################################