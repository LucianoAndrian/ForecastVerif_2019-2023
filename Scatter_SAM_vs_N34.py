"""
Scatter plot entre la magnitud de 2 índices y la probabilidad pronostica
"""
################################################################################
nmme_pronos = '/pikachu/datos/luciano.andrian/verif_2019_2023/nmme_pronos/'
dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/scatter/'
dir_results = 'index_vs_index'
################################################################################
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Funciones import SelectFilesNMME, DMI, SameDateAs, LeadMonth, \
    CreateDirectory, DirAndFile
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
    dpi = 100

lon_regiones = [[296, 296 + 20], [296, 296 + 20], [296, 300 + 20],
                [295, 295 + 10]]
lat_regiones = [[-40, -40 + 20], [-40, -40 + 10], [-30, -30 + 17],
                [-40, -40 + 15]]
titulos = ['SESA', 'S-SESA', 'N-SESA', 'NEA']
try:
    if test:
        lon_regiones = [[296, 296 + 20]]
        lat_regiones = [[-40, -40 + 20]]
        titulos = ['SESA']
        lead = [0]
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
    serie_mean = (np.array(np.nanmean(serie))*8)**5
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
n34 = np.array([0.7, 0.7, 0.7, 0.7, 0.5, 0.5, 0.3, 0.1, 0.2 ,0.3, 0.5 ,0.5, 0.5,
                0.5, 0.4, 0.2, -0.1, -0.3, -0.4, -0.6, -0.9, -1.2, -1.3, -1.2,
                -1.0, -0.9, -0.8, -0.7, -0.5, -0.4, -0.4, -0.5, -0.7, -0.8,
                -1.0, -1.0, -1.0, -0.9, -1.0, -1.1, -1.0, -0.9, -0.8, -0.9,
                -1.0, -1.0, -0.9, -0.8, -0.7, -0.4, -0.1])

# SAM
sam = xr.open_dataset(dir + 'sam.nc')['mean_estimate']
asam = xr.open_dataset(dir + 'asam.nc')['mean_estimate']
ssam = xr.open_dataset(dir + 'ssam.nc')['mean_estimate']

sam = sam.sel(time=slice('2019-01-01', '2023-03-01'))
asam = SameDateAs(asam, sam)
ssam = SameDateAs(ssam, sam)
dmi = SameDateAs(dmi_aux, sam)

dates = sam.time.values

indices = [sam, asam, ssam]
indices_name = ['SAM', 'A-SAM', 'S-SAM']
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
        for d in dates:
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

        mean_above_probs = np.array(mean_above_probs)
        mean_norm_probs = np.array(mean_norm_probs)
        mean_below_probs = np.array(mean_below_probs)

        # ---------------------------------------------------------------------#
        # Plots ---------------------------------------------------------------#
        print('Plots by indices...')
        # Mean prob. ----------------------------------------------------------#

        n34_aux = n34[l::]
        for i, ititle in zip(indices, indices_name):
            try:
                i = i.values
            except:
                pass
            i = i[l::]
            imax = np.round(max(abs(i)),1)
            fig = plt.figure(figsize=(6, 5), dpi=dpi)
            ax = fig.add_subplot(111)

            # plot sólo el mayor de los valores en cada punto
            max_prob = np.argmax([mean_above_probs, mean_norm_probs,
                                 mean_below_probs], axis=0)

            ax.scatter(n34_aux[max_prob == 0], i[max_prob == 0],
                       s=mean_above_probs[max_prob == 0], color='dodgerblue',
                       marker='^', label='Above', alpha=0.5)
            ax.scatter(n34_aux[max_prob == 1], i[max_prob == 1],
                       s=mean_norm_probs[max_prob == 1], color='forestgreen',
                       marker='s', label='Normal', alpha=0.5)
            ax.scatter(n34_aux[max_prob == 2], i[max_prob == 2],
                       s=mean_below_probs[max_prob == 2], color='firebrick',
                       marker='v',label='Below', alpha=0.5)

            ax.grid(True)
            plt.legend(markerscale=.5)
            ax.set_ylim(-imax-.2, imax+.2)
            plt.xlim((-1.75, 1.75))
            ax.set_ylabel(ititle, size=12)
            ax.set_xlabel('Niño3.4', size=12)

            plt.title(ititle + ' vs N34 - ' + t + '\n' + 'Lead: '
                      + str(l))

            plt.tight_layout()

            if save:
                plt.savefig(DirAndFile(out_dir, dir_results, 'Prob',
                                       [t, ititle, 'vs_N34', 'Lead', str(l)]))
            else:
                plt.show()


################################################################################
print('#######################################################################')
print('done')
print('out_dir = ' + out_dir + dir_results )
print('#######################################################################')
################################################################################