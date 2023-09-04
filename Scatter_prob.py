"""
Scatter plot entre la magnitud de los índices y la probabilidad pronostica
"""
################################################################################
nmme_pronos = '/pikachu/datos/luciano.andrian/verif_2019_2023/nmme_pronos/'
dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/scatter/'
dir_results = 'prob_vs_indices'
dir_results2 = 'prob_vs_indices_3d'
################################################################################
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Funciones import CreateDirectory, SelectFilesNMME, DMI, SameDateAs, \
    LeadMonth, DirAndFile, OpenRegiones

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
CreateDirectory(out_dir, dir_results2)
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
    pass
################################################################################
def Proc(array):
    serie = np.reshape(array, array.size)
    serie_mean = np.nanmean(serie)
    return serie, serie_mean
################################################################################
# endtime determinado por el SAM, último índice en actualizarse
endtime = xr.open_dataset(dir + 'sam.nc').time.values[-1]
################################################################################
# NMME forecast
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
# indices
print('Indices DMI, N34, SAM, S-SAM y A-SAM')
# ONI Descargado, para no cambiar tdo el codigo n34 = ONI
n34 = xr.open_dataset(dir + 'oni.nc')
dates = n34.time.values

# SAM
sam = xr.open_dataset(dir + 'sam.nc').mean_estimate
asam = xr.open_dataset(dir + 'asam.nc').mean_estimate
ssam = xr.open_dataset(dir + 'ssam.nc').mean_estimate

# DMI calculado a partir de ERSSTv5 actualizada
aux0, aux, dmi = DMI(filter_bwa=False, start_per=1920, end_per=anio)
dmi = SameDateAs(dmi, sam)

indices = [sam, asam, ssam, dmi, n34.oni]
indices_name = ['SAM', 'A-SAM', 'S-SAM', 'DMI', 'Niño3.4']
################################################################################

for ln, lt, t in zip(lon_regiones, lat_regiones, titulos):
    print(t)
    region = data_nmme.sel(lon=slice(ln[0], ln[1]), lat=slice(lt[1], lt[0]))

    # leads
    for l in lead:
        mean_below_probs = []
        below_probs = []
        mean_norm_probs = []
        norm_probs = []
        mean_above_probs = []
        above_probs = []

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
                below_probs.append(below)

                # norm --------------------------------------------------------#
                norm, norm_mean = Proc(aux.prob_norm.values)
                mean_norm_probs.append(norm_mean)
                norm_probs.append(norm)

                # above -------------------------------------------------------#
                above, above_mean = Proc(aux.prob_above.values)
                mean_above_probs.append(above_mean)
                above_probs.append(above)

            except:
                print('Skip ' + np.datetime_as_string(d0, unit='M') + ' con' +
                      ' target ' + np.datetime_as_string(d, unit='M'))

        # ---------------------------------------------------------------------#
        # Plots ---------------------------------------------------------------#
        print('Plots by indices...')
        # Mean prob. ----------------------------------------------------------#
        for i, ititle in zip(indices, indices_name):
            fig, ax = plt.subplots(dpi=dpi)

            i = i[l::]
            sd = np.std(i)
            max = np.max(i) + sd/2
            min = np.min(i) - sd/2

            ax.scatter(x=i, y=mean_above_probs, color='dodgerblue', marker='^',
                        label='Above')
            ax.scatter(x=i, y=mean_norm_probs, color='forestgreen', marker='s',
                        label='Normal')
            ax.scatter(x=i, y=mean_below_probs, color='firebrick', marker='v',
                        label='Below')

            plt.legend()
            plt.ylim((0.10,0.9))
            plt.xlim((min, max))
            ax.grid(True)

            fig.set_size_inches(6, 6)
            plt.xlabel(ititle, size=15)
            plt.ylabel('Prob.', size=15)


            plt.title(ititle + ' vs Prob. - ' + t + '\n' + 'Lead: ' + str(l))
            plt.tight_layout()
            if save:
                plt.savefig(DirAndFile(out_dir, dir_results, 'PROB',
                                       [t, ititle, 'Lead', str(l)]))
            else:
                plt.show()

        print('Done mean_prob plots')

################################################################################
print('#######################################################################')
print('done')
print('out_dir = ' + out_dir + dir_results )
print('#######################################################################')
################################################################################
# Descarte 3d
# # ---------------------------------------------------------------------#
# # All probs. ----------------------------------------------------------#
# for i, ititle in zip(indices, indices_name):
#     fig = plt.figure(dpi=dpi)
#     ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection='3d')
#     try:
#         i = i.values
#     except:
#         pass
#     i = i[l::]
#     sd = np.std(i)
#     max = np.max(i) + sd / 2
#     min = np.min(i) - sd / 2
#
#     for x, y_vals in zip(i, above_probs):
#         ax.scatter([x] * len(y_vals), [1] * len(y_vals), y_vals,
#                    color='dodgerblue', marker='^', s=10)
#
#     for x, y_vals in zip(i, norm_probs):
#         ax.scatter([x] * len(y_vals), [2] * len(y_vals), y_vals,
#                    color='forestgreen', marker='s', s=10)
#
#     for x, y_vals in zip(i, below_probs):
#         ax.scatter([x] * len(y_vals), [3] * len(y_vals), y_vals,
#                    color='firebrick', marker='v', s=10)
#
#     ax.scatter([], [], [], color='dodgerblue', marker='^',
#                label='Above')
#     ax.scatter([], [], [], color='forestgreen', marker='s',
#                label='Normal')
#     ax.scatter([], [], [], color='firebrick', marker='v', label='Below')
#     ax.view_init(elev=20, azim=-45)
#
#     ax.set_yticks([1, 2, 3])
#     ax.set_yticklabels(['Above', 'Normal', 'Below'])
#
#     plt.legend()
#     ax.set_zlim(0, .8)
#     plt.xlim((min, max))
#     ax.set_zlabel('Prob.', size=12)
#     ax.set_xlabel(ititle, size=12)
#
#     plt.title(ititle + ' vs Prob. - ' + t + '\n' + 'Lead: ' + str(l))
#     plt.tight_layout()
#     if save:
#         plt.savefig(DirAndFile(out_dir, dir_results2, 'Prob3D',
#                                [t, ititle, 'Lead', str(l)]))
#     else:
#         plt.show()