"""
Evolucion temporal de los índices y "scatter" de las prob. de los puntos de
 grilla.
"""
"""
Scatter plot entre la magnitud de los índices y la probabilidad pronostica
"""
################################################################################
nmme_pronos = '/pikachu/datos/luciano.andrian/verif_2019_2023/nmme_pronos/'
dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/scatter/'
dir_results = 'EvolT'
################################################################################
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from Funciones import CreateDirectory, SelectFilesNMME, DMI, SameDateAs, \
    LeadMonth, DirAndFile, OpenRegiones
import set_indices, nmme_update
from dateutil.relativedelta import relativedelta
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
warnings.filterwarnings("ignore")
################################################################################
save = False
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

for ln, lt, t in zip(lon_regiones, lat_regiones, titulos):
    print(t)
    region = data_nmme.sel(lon=slice(ln[0], ln[1]), lat=slice(lt[1], lt[0]))

    # leads
    for l in lead:
        mean_below_probs = []
        mean_norm_probs = []
        mean_above_probs = []


        print('Lead: ' + str(l))
        # cada fecha
        for d in date_nmme:
            if lead != 0:
                d0 = LeadMonth(d, l)
            else:
                d0 = d

            try:
                aux = region.sel(time=d0, target=d)

                # below -------------------------------------------------------#
                below, below_mean = Proc(aux.prob_below.values)
                mean_below_probs.append(below_mean)

                # norm --------------------------------------------------------#
                norm, norm_mean = Proc(aux.prob_norm.values)
                mean_norm_probs.append(norm_mean)

                # above -------------------------------------------------------#
                above, above_mean = Proc(aux.prob_above.values)
                mean_above_probs.append(above_mean)

            except:
                print('Skip ' + np.datetime_as_string(d0, unit='M') + ' con' +
                      ' target ' + np.datetime_as_string(d, unit='M'))

        # Plots ---------------------------------------------------------------#
        print('Plots by indices...')
        # Mean prob. ----------------------------------------------------------#
        dates = date_nmme[l::]
        fig = plt.figure(figsize=(10, 7), dpi=dpi)
        ax = fig.add_subplot(111)
        ax2 = ax.twinx()

        ax.scatter(x=dates, y=mean_above_probs, color='dodgerblue', marker='^',
                   label='Above')
        ax.scatter(x=dates, y=mean_norm_probs, color='forestgreen', marker='s',
                   label='Normal')
        ax.scatter(x=dates, y=mean_below_probs, color='firebrick', marker='v',
                   label='Below')

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
        plt.rcParams['date.converter'] = 'concise'
        ax.xaxis.set_major_locator(
            mdates.AutoDateLocator(minticks=20, maxticks=26))

        ax.grid()
        ymax = np.round(np.max([mean_above_probs, mean_norm_probs,
                                mean_below_probs]) + .10, 1)
        ax.set_ylim((0.1, ymax))
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
            plt.savefig(DirAndFile(out_dir, dir_results, 'EvolT',
                                   [t, 'Lead', str(l)]), dpi=dpi)
            plt.close('all')
        else:
            plt.show()


        print('Done mean_prob plots')

################################################################################
print('#######################################################################')
print('done')
print('out_dir = ' + out_dir + dir_results )
print('#######################################################################')
################################################################################