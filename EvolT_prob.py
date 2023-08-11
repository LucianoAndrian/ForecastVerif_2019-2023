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
    LeadMonth
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

lon_regiones = [[296, 296 + 20], [296, 296 + 20], [296, 300 + 20],
                [295, 295 + 10], [290, 290 + 5]]
lat_regiones = [[-40, -40 + 20], [-40, -40 + 10], [-30, -30 + 17],
                [-40, -40 + 15], [-40, -40 + 20]]
titulos = ['SESA', 'S-SESA', 'N-SESA', 'NEA', 'NOA']
try:
    if test:
        lon_regiones = [[296, 296 + 20]]
        lat_regiones = [[-40, -40 + 20]]
        titulos = ['SESA']
        lead = [0, 2]
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
    serie_mean = serie.mean()
    return serie, serie_mean
################################################################################
# NMME forecast
files = SelectFilesNMME(nmme_pronos, 'prate')
files = files[:-8] # se descargaron tdo 2023 que no hay nada desde mayo
# pronos desde 201901 hasta 202304 (52)
data_nmme = xr.open_mfdataset(files, decode_times=False, engine='netcdf4',
                              combine='nested', concat_dim='initial_time')
data_nmme = data_nmme.rename({'initial_time':'time'}) # para mas adelante
data_nmme['time'] = pd.date_range(start='2018-12-01', end='2023-04-01',
                                  freq='M') + pd.DateOffset(days=1)
data_nmme['target'] = pd.date_range(start='2018-12-01', end='2023-10-01',
                                  freq='M') + pd.DateOffset(days=1)

################################################################################
# indices
print('Indices DMI, N34, SAM, S-SAM y A-SAM')
# SST actualizada
# https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/
dmi, aux, dmi_aux = DMI(filter_bwa=False, start_per=1920, end_per=2023)

#pendiente, arreglar Ninio3.4CPC
n34 = [0.7, 0.7, 0.7, 0.7, 0.5, 0.5, 0.3, 0.1, 0.2 ,0.3, 0.5 ,0.5, 0.5,	0.5,
0.4, 0.2, -0.1, -0.3, -0.4, -0.6, -0.9, -1.2, -1.3, -1.2, -1.0, -0.9, -0.8,
-0.7, -0.5, -0.4, -0.4, -0.5, -0.7, -0.8, -1.0, -1.0, -1.0, -0.9, -1.0, -1.1,
-1.0, -0.9, -0.8, -0.9, -1.0, -1.0, -0.9, -0.8, -0.7, -0.4, -0.1]

# SAM
sam = xr.open_dataset(dir + 'sam.nc')['mean_estimate']
asam = xr.open_dataset(dir + 'asam.nc')['mean_estimate']
ssam = xr.open_dataset(dir + 'ssam.nc')['mean_estimate']

sam = sam.sel(time=slice('2019-01-01', '2023-03-01'))
asam = SameDateAs(asam, sam)
ssam = SameDateAs(ssam, sam)
dmi = SameDateAs(dmi_aux, sam)

dates2 = sam.time.values

indices = [sam, asam, ssam, dmi, n34]
indices_name = ['SAM', 'A-SAM', 'S-SAM', 'DMI', 'Niño3.4']
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
        for d in dates2:
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
        dates = sam.time.values[l::]
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
        lndmi = ax2.plot(dates, dmi[l::].values, label='DMI',
                         color='forestgreen')

        lnn34 = ax2.plot(dates, n34[l::], label='N34', color='firebrick')

        lnsam = ax2.plot(dates, sam[l::].values, label='SAM', color='k')

        lnasam = ax2.plot(dates, asam[l::].values, label='A-SAM',
                          color='purple')

        lnssam = ax2.plot(dates, ssam[l::].values, label='S-SAM', color='lime')

        ax.hlines(y=0, xmin=dates[0], xmax=dates[-1], color='gray')
        ax2.hlines(y=0, xmin=dates[0], xmax=dates[-1], color='gray')
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y %b'))
        ax.grid()
        ax.set_ylim((0.1, .5))
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
            plt.savefig(
                out_dir + '/' + dir_results + '/' + 'EvolT' + '_' + t +
                '_lead_' + str(l) + '.jpg', dpi=dpi)
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