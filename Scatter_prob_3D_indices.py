"""
Scatter plot entre la magnitud de 2 índices y la probabilidad pronostica
"""
################################################################################
nmme_pronos = '/pikachu/datos/luciano.andrian/verif_2019_2023/nmme_pronos/'
dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/scatter/'
dir_results = '3d_index_vs_index'
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

lon_regiones = [[296, 296 + 20], [296, 296 + 20], [300, 300 + 20],
                [296, 296 + 8], [290, 290 + 5], [288, 288 + 8],
                [290, 290 + 5]]

lat_regiones = [[-40, -40 + 20], [-40, -40 + 10], [-30, -30 + 17],
                [-35, -35 + 13], [-35, -35 + 15], [-55, -55 + 15],
                [-40, -40 + 10]]
titulos = ['SESA', 'S-SESA', 'N-SESA', 'NEA', 'NOA', 'Patagonia', 'Cuyo']

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
    serie_mean = np.nanmean(serie)
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

dates = sam.time.values

indices = [sam, asam, ssam]
indices_name = ['SAM', 'A-SAM', 'S-SAM']

indices2 = [dmi, n34]
indices2_name = ['DMI', 'N34']
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

        # ---------------------------------------------------------------------#
        # Plots ---------------------------------------------------------------#
        print('Plots by indices...')
        # Mean prob. ----------------------------------------------------------#
        for i, ititle in zip(indices, indices_name):
            try:
                i = i.values
            except:
                pass
            i = i[l::]
            sd = np.std(i)
            max = np.max(i) + sd / 2
            min = np.min(i) - sd / 2

            for i2, i2title in zip(indices2, indices2_name):

                fig = plt.figure(dpi=dpi)
                ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection='3d')
                try:
                    i = i.values
                except:
                    pass

                # No entiendo xq ahora n34 no funciona igual q antes...
                try:
                    i2.mean()
                except:
                    i2 = np.array(i2)
                ###

                i2 = i2[l::]
                sd2 = np.std(i2)

                max2 = np.max(i2) + sd2 / 2
                min2 = np.min(i2) - sd2 / 2

                ax.scatter(i, i2, mean_above_probs, color='dodgerblue',
                           marker='^',
                           label='Above')
                ax.scatter(i, i2, mean_norm_probs, color='forestgreen',
                           marker='s',
                           label='Normal')
                ax.scatter(i, i2, mean_below_probs, color='firebrick',
                           marker='v',
                           label='Below')

                # ax.scatter([], [], [], color='dodgerblue', marker='^',
                #            label='Above')
                # ax.scatter([], [], [], color='forestgreen', marker='s',
                #            label='Normal')
                # ax.scatter([], [], [], color='firebrick', marker='v',
                #            label='Below')
                ax.view_init(elev=25, azim=-45)

                plt.legend()
                ax.set_ylim(min2, max2)
                plt.xlim((min, max))
                ax.set_zlim(0, .6)
                ax.set_ylabel(i2title, size=12)
                ax.set_zlabel('Prob.', size=12)
                ax.set_xlabel(ititle, size=12)

                plt.title(ititle + ' vs ' + i2title +'_' + t + '\n' + 'Lead: '
                          + str(l))

                plt.tight_layout()

                if save:
                    plt.savefig(DirAndFile(out_dir, dir_results, '3D',
                                           [t, ititle, 'vs', i2title, 'Lead',
                                            str(l)]))
                else:
                    plt.show()

################################################################################
print('#######################################################################')
print('done')
print('out_dir = ' + out_dir + dir_results )
print('#######################################################################')
################################################################################