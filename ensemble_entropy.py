"""
Evolucion temporal de los índices y entropia del ensamble
Mapas de entropia del ensamble anuales
"""
################################################################################
nmme_pronos = '/pikachu/datos/luciano.andrian/verif_2019_2023/nmme_pronos/'
dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'
dir_results = 'Entropy'
################################################################################
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from Funciones import CreateDirectory, SelectFilesNMME, DirAndFile, \
    OpenRegiones, ColorBySeason, Entropy
from Tool_Poligonos import get_vertices, get_mask
from matplotlib.path import Path
import set_indices, nmme_update
from dateutil.relativedelta import relativedelta
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
warnings.filterwarnings("ignore")
################################################################################
save = True
test = False # solo va computar una region
lead = [0, 1, 2, 3]
CreateDirectory(out_dir, dir_results)
################################################################################
# some sets
if save:
    dpi = 300
else:
    dpi = 50

titulos, lon_regiones, lat_regiones = OpenRegiones('regiones_sa.csv')
try:
    if test:
        lon_regiones = lon_regiones[0]
        lat_regiones = lat_regiones[0]
        titulos = titulos[0]
        lead = [0, 1, 2]
        print('##########################################################')
        print('<<<<<<<<<<<<<<<<<<<<<<<<<< TEST >>>>>>>>>>>>>>>>>>>>>>>>>>')
        print('-----------------------Una sola region--------------------')
        print('--------------------------Lead 0--------------------------')
        print('##########################################################')
except:
    pass
################################################################################
def Proc(array):
    serie = np.reshape(array, array.size)
    serie_mean = np.nanmean(serie)
    return serie, serie_mean

def ComputeAndPlot(region, t, save):
    for l in lead:
        print('Lead: ' + str(l))

        entropy = Entropy(region, l)
        entropy = entropy.mean(['lon', 'lat'])

        # Plots ---------------------------------------------------------------#
        print('Plots by indices...')
        # Mean prob. ----------------------------------------------------------#
        dates = date_nmme[l::]
        plt.rcParams['date.converter'] = 'concise'
        fig = plt.figure(figsize=(10, 7), dpi=dpi)
        ax = fig.add_subplot(111)
        ax.xaxis.set_major_locator(
            mdates.AutoDateLocator(minticks=20, maxticks=26))
        ax2 = ax.twinx()

        # Colores para cada season
        for d in dates:
            color = ColorBySeason(d)
            aux = d.astype('M8[D]').astype('O')
            d2 = np.datetime64(aux + relativedelta(months=1))
            ax.axvspan(d, d2 , alpha=0.2, color=color)

        ax.plot(dates, entropy.values, color='k', marker='.',
                   label='Entropy', linewidth=1, zorder=2)

        max_entropy = xr.where(np.round(entropy.values, 3) > np.round(
            -1 * (0.33 * np.log(0.33) * 3), 3) - 0.01, entropy, np.nan).values

        ax.scatter(x=dates, y=max_entropy, color='red', marker='.',
                   label='Max. Entropy', zorder=3)

        # indices
        lndmi = ax2.plot(dates_index[l::], dmi[l::].values, label='DMI',
                         color='#289E64')

        lnn34 = ax2.plot(dates_index[l::], n34.oni[l::].values, label='N34',
                         color='#00C9ED')

        lnsam = ax2.plot(dates_index[l::], sam[l::].values, label='SAM',
                         color='#005EFF')

        lnasam = ax2.plot(dates_index[l::], asam[l::].values, label='A-SAM',
                          color='#960B00')

        lnssam = ax2.plot(dates_index[l::], ssam[l::].values, label='S-SAM',
                          color='#FF0088')

        ax.hlines(y=0, xmin=dates[0], xmax=dates[-1], color='gray')
        ax2.hlines(y=0, xmin=dates[0], xmax=dates[-1], color='gray')

        ax.grid(zorder=1)
        ymax = np.round(0.33*np.log(0.33)*3*-1,1) + 0.05
        ymin = np.round(min(entropy.values),1) - 0.1
        ax.set_ylim((ymin, ymax))
        ax2.set_ylim((-1.5, 7))
        ax.set_ylabel('Probabilidad', fontsize=10)
        ax2.set_ylabel('índices', fontsize=10)
        ax.set_title(
            'Prob.' + ' - ' + t + '\n' + 'Lead: ' + str(l),
            fontsize=15)

        lns = lndmi + lnn34 + lnsam + lnasam + lnssam
        ax.legend(lns, ['DMI', 'ONI', 'SAM', 'A-SAM', 'S-SAM'],
                  loc='upper left')

        if save:
            plt.savefig(DirAndFile(out_dir, dir_results, 'Entropy',
                                   [t, 'Lead', str(l)]), dpi=dpi)
            plt.close('all')
        else:
            plt.show()

################################################################################
print('Set Indices ###########################################################')
n34, dmi, sam, ssam, asam, endtime = set_indices.compute()
dates_index = n34.time.values

indices = [sam, asam, ssam, dmi, n34]
indices_name = ['SAM', 'A-SAM', 'S-SAM', 'DMI', 'Niño3.4']
################################################################################
# endtime determinado por el ONI, se actualiza al trimestre anterior
# e.g. al finalizar agosto actualiza ONI en JJA --> mm = 7
print('#######################################################################')
print('<<<<<<<<<<<<<<<<<<<<< indices hasta: ' + str(endtime).split('T')[0] +
      ' >>>>>>>>>>>>>>>>>>>>>>')
print('#######################################################################')
################################################################################
# NMME forecast
nmme_update.update()
files = SelectFilesNMME(nmme_pronos, 'prate', size_check=True)

# endtime de los pronos libre
endtime_nmme = files[-1].split('_')[-2]
endtime_nmme = f"{endtime_nmme[:4]}-{endtime_nmme[4:6]}-01"
endtime_nmme = np.datetime64(endtime_nmme)
print('#######################################################################')
print('<<<<<<<<<<<<<<<<<<<<<<< NMME hasta: ' + str(endtime_nmme).split('T')[0] +
      ' >>>>>>>>>>>>>>>>>>>>>>>>')
print('#######################################################################')

# y el ultimo target correspondiente
aux = endtime_nmme.astype('M8[D]').astype('O')
targetime = aux + relativedelta(months=6)
targetime = np.datetime64(targetime)

# pronos desde 201901
data_nmme = xr.open_mfdataset(files, decode_times=False, engine='netcdf4',
                              combine='nested', concat_dim='initial_time')
data_nmme = data_nmme.rename({'initial_time':'time'}) # para mas adelante
data_nmme['time'] = pd.date_range(start='2018-12-01', end=endtime_nmme,
                                  freq='M') + pd.DateOffset(days=1)
data_nmme['target'] = pd.date_range(start='2018-12-01', end=targetime,
                                  freq='M') + pd.DateOffset(days=1)

date_nmme = data_nmme.time.values
################################################################################
# Regiones SA
for ln, lt, t in zip(lon_regiones, lat_regiones, titulos):
    print(t)
    region = data_nmme.sel(lon=slice(ln[0], ln[1]), lat=slice(lt[1], lt[0]))

    ComputeAndPlot(region, t, save)

# -----------------------------------------------------------------------------#
# Regiones Arg
regiones_arg = 'hum_sur2', 'hum_norte2', 'centro2', 'patagonia_oeste', \
               'patagonia', 'noa'
for t in regiones_arg:
    print(t)
    path = Path(get_vertices(t))
    region = get_mask(data_nmme, path)

    ComputeAndPlot(region, t, save)


################################################################################
print('#######################################################################')
print('done')
print('out_dir = ' + out_dir + dir_results )
print('#######################################################################')
################################################################################

