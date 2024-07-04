"""
Anomalias de PP observadas en cada trimestre solapado
"""
################################################################################
save = True
update = True
plot_full = True
################################################################################
nmme_pronos = '/pikachu/datos/luciano.andrian/verif_2019_2023/nmme_pronos/'
cmap_data = '/pikachu/datos/luciano.andrian/verif_2019_2023/cmap/'
chirps_data = '/pikachu/datos/luciano.andrian/verif_2019_2023/chirps/'
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'
dir_results = 'ppanoms_trim'

################################################################################
import numpy as np
import xarray as xr
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import colors
from Funciones import SelectFilesNMME, MakeMask, ChangeLons, ABNobs, \
    RPSO, RPSF, BSO, BSF, Plot, SameDateAs, DirAndFile, \
    CreateDirectory, OpenRegiones, Correlaciones, ColorBySeason
import set_indices, cmap, chirps, Scales_Cbars
from dateutil.relativedelta import relativedelta
from datetime import datetime
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
warnings.filterwarnings("ignore")

################################################################################
if save:
    dpi = 300
    CreateDirectory(out_dir, dir_results)
else:
    dpi = 100

seasons = [666,'DJF', 'JFM', 'FMA', 'MAM', 'AMJ', 'MJJ', 'JJA', 'JAS', 'ASO',
           'SON', 'OND', 'NDJ']
################################################################################
def Check_LastTrimPlot(endtime, data, dataset, out_dir, dir_results, full):

    seasons = [666, 'DJF', 'JFM', 'FMA', 'MAM', 'AMJ', 'MJJ',
               'JJA', 'JAS','ASO', 'SON', 'OND', 'NDJ']

    path = out_dir + dir_results + '/'
    try:
        file = 'pp_anom_trim_' + seasons[int(str(endtime)[5:7]) - 1] + '-' + \
               str(endtime)[:4] + '_' + dataset + '.jpg'
    except:
        file = 'pp_anom_trim_' + seasons[int(str(endtime)[5:7])] + '-' + \
               str(endtime)[:4] + '_' + dataset + '.jpg'


    file_path = os.path.join(path, file)

    if os.path.exists(file_path):
        dates_to_plot = []
    else:
        dates_to_plot = [data.time.values[-2]]
        if full:
            dates_to_plot = data.time.values[:-2]



    return dates_to_plot

################################################################################
date = datetime.now().replace(day=1, hour=0, minute=0, second=0, microsecond=0)
endtime = np.datetime64(date - relativedelta(months=1))
################################################################################
# Apertura de bases de datos
data = xr.open_dataset(
    chirps_data + 'chirps_1990_2020_mmean.nc').__mul__(365/12) #dia
data = ChangeLons(data, 'longitude')
data = data.rename({'latitude':'lat'})


data = data.sel(lon=slice(270,330), lat=slice(-60, 20))
chirps_clim = data.rolling(time=3, center=True).mean()
chirps_clim = chirps_clim.groupby('time.month').mean()

# 2019-2023
data = xr.open_dataset(
    chirps_data + 'chirps_2019_2023_mmean.nc').__mul__(365/12) #dia

# Ideem que con CMAP
if endtime > data.time.values[-2]:
    print('CHIRPS desactualizado')
    if update:
        print('Intentando actualizar CHIRPS')
        chirps.update()
        data = xr.open_dataset(
            chirps_data + 'chirps_2019_2023_mmean.nc').__mul__(365 / 12)  # dia
        if endtime == data.time.values[-2]:
            print('Actualizacion CHIRPS no disponible')
            chirps_updated = False
        else:
            chirps_updated = True
else:
    chirps_updated = True

data = ChangeLons(data, 'longitude')
data = data.rename({'latitude':'lat'})

# promedios trimestrales (VER)
data_verif_ch = data.rolling(time=3, center=True).mean()
data_verif_ch = data_verif_ch.sel(lon=slice(270,330), lat=slice(-60, 20))
data_verif_ch = data_verif_ch.sel(time=slice('2019-01-01', endtime))

aux = data.rolling(time=3, center=True).mean()
chirps_anom = aux.groupby('time.month') - chirps_clim
chirps_anom = chirps_anom*MakeMask(chirps_anom, 'precip')
chirps_anom = chirps_anom.sel(time=slice('2019-01-01', endtime))
################################################################################
if chirps_updated:
    dates_to_plot = Check_LastTrimPlot(endtime, chirps_anom, 'chirps',
                                       out_dir, dir_results, plot_full)

    scale = Scales_Cbars.get_scales('pp')
    cbar = Scales_Cbars.get_cbars('pp')

    for d in dates_to_plot:
        Plot(chirps_anom.sel(time=d), chirps_anom.sel(time=d).precip,
             scale, save, dpi,
             "PP' " + seasons[int(str(d)[5:7])] + '-' + str(d)[:4],
             out_dir + dir_results + '/pp_anom_trim_' +
             seasons[int(str(d)[5:7])] + '-' + str(d)[:4] + '_chirps',
             'k', cbar)
################################################################################