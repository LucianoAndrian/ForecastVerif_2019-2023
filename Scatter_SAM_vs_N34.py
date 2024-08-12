"""
Scatter plot entre la magnitud de 2 índices y la probabilidad pronostica
"""
################################################################################
save = True
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
    CreateDirectory, DirAndFile, OpenRegiones, RMean3
import set_indices, nmme_update
from dateutil.relativedelta import relativedelta
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
warnings.filterwarnings("ignore")
################################################################################
test = False # solo va computar una region
lead = [0, 1, 2, 3]
CreateDirectory(out_dir, dir_results)
################################################################################
# some sets
update = False
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
    pass
################################################################################
def Proc(array):
    serie = np.reshape(array, array.size)
    serie_mean = (np.array(np.nanmean(serie))*8)**5
    return serie, serie_mean
################################################################################
print('Set Indices ###########################################################')
n34, dmi, sam, ssam, asam, endtime = set_indices.compute()
dates = n34.time.values
# ---------------------------------------------------------------------------- #
# Hasta junio de 2024 debido al cambio de modelos en el ensamble canadiense
# CANSips en julio.
endtime = np.datetime64('2024-06-01T00:00:00.000000000')

indices = [sam, asam, ssam]
indices_name = ['SAM', 'A-SAM', 'S-SAM']
################################################################################
# endtime determinado por el ONI, se actualiza al trimestre anterior
# e.g. al finalizar agosto actualiza ONI en JJA --> mm = 7
endtime = n34.time.values[-1]
print('#######################################################################')
print('<<<<<<<<<<<<<<<<<<<<< Valores hasta: ' + str(endtime).split('T')[0] +
      ' >>>>>>>>>>>>>>>>>>>>>>')
print('#######################################################################')
################################################################################
# NMME forecast
if update:
    print('update = True !!!')
    nmme_update.update()
files = SelectFilesNMME(nmme_pronos, 'prate', size_check=True)

# para identificar el prono correspondiente
anio = endtime.astype('datetime64[Y]').astype(int) + 1970
mes = endtime.astype('datetime64[M]').astype(int) % 12 + 1
endtime_str = f"{anio}{mes:02d}"
# y el ultimo target correspondiente

aux = endtime.astype('M8[D]').astype('O')
targetime = aux + relativedelta(months=6)
targetime = np.datetime64(targetime)

# ultimo pronostico que que puede ser verificado
posf = [i for i, prono in enumerate(files) if endtime_str in prono][0]
files = files[0:posf+1]

# pronos desde 201901
data_nmme = xr.open_mfdataset(files, decode_times=False, engine='netcdf4',
                              combine='nested', concat_dim='initial_time')
data_nmme = data_nmme.rename({'initial_time':'time'}) # para mas adelante
data_nmme['time'] = pd.date_range(start='2018-12-01', end=endtime,
                                  freq='M') + pd.DateOffset(days=1)
data_nmme['target'] = pd.date_range(start='2018-12-01', end=targetime,
                                  freq='M') + pd.DateOffset(days=1)

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

        n34_aux = n34.oni[l::]
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