nmme_pronos = '/pikachu/datos/luciano.andrian/verif_2019_2023/nmme_pronos/'
cmap_data = '/pikachu/datos/luciano.andrian/verif_2019_2023/cmap/'
chirps_data = '/pikachu/datos/luciano.andrian/verif_2019_2023/chirps/'
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'
dir_results = 'mapas_index'

#https://psl.noaa.gov/data/gridded/data.cmap.html
################################################################################
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import matplotlib.patches as mpatches
from matplotlib import colors
from Funciones import Nino34CPC, DMI # indices
from Funciones import SelectFilesNMME, MakeMask, ChangeLons, ABNobs, \
    RPSO, RPSF, BSO, BSF, CorrSP, Plot, SameDateAs, DirAndFile, CreateDirectory

import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
warnings.filterwarnings("ignore")
################################################################################
save = True
plot_mapas = True
mapa = True
correlaciones = True
test = False # = True solo va computar una region
lead = [0, 1, 2, 3]
if save:
    dpi = 300
    CreateDirectory(out_dir, dir_results)
else:
    dpi = 100
################################################################################
def ComputeAndPlot(index, correlaciones, dpi, save, lead=0, test=False):
    dates = pd.date_range(start='2018-12-01', end='2023-04-01',
                          freq='M') + pd.DateOffset(days=1)
    dates = dates[lead:-1]

    # lon_regiones = [[296, 296 + 20], [296, 296 + 20], [296, 300 + 20],
    #                 [295, 295 + 10], [290, 290 + 5]]
    # lat_regiones = [[-40, -40 + 20], [-40, -40 + 10], [-30, -30 + 17],
    #                 [-40, -40 + 15], [-40, -40 + 20]]
    #titulos = ['SESA', 'S-SESA', 'N-SESA', 'NEA', 'NOA']

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
            print('##########################################################')
            print('<<<<<<<<<<<<<<<<<<<<<<<<<< TEST >>>>>>>>>>>>>>>>>>>>>>>>>>')
            print('-----------------------Una sola region--------------------')
            print('##########################################################')

    except:
        pass

    if index.upper() == 'RPSS':
        index_cmap = RPSS_cmap
        index_chirps = RPSS_chirps
    elif index.upper() == 'BSS':
        index_cmap = BSS_cmap
        index_chirps = BSS_chirps

    print(index.upper() + ' por regiones... ')
    if correlaciones:
        print(' y correlaciones...')
    for ln, lt, t in zip(lon_regiones, lat_regiones, titulos):
        data_frames = True
        aux = index_cmap.sel(lon=slice(ln[0], ln[1]), lat=slice(lt[1], lt[0]))
        aux2 = index_chirps.sel(lon=slice(ln[0], ln[1]), lat=slice(lt[1],lt[0]))

        fig = plt.figure(figsize=(10, 7), dpi=dpi)
        ax = fig.add_subplot(111)
        ax2 = ax.twinx()

        # index score
        aux_lnscmap = aux.mean(['lon', 'lat']).mask[:-1]
        lnscmap = ax.plot(dates, aux_lnscmap,
                             color='#FF0003',
                             label=index.upper() + '_CMAP2.5', linewidth=2)

        aux_lnschirps = aux2.mean(['lon', 'lat']).mask[:-1]
        lnschirps = ax.plot(dates, aux_lnschirps,
                               color='#FFA500',
                               label=index.upper() + '_CHIRPS1', linewidth=2)

        # indices
        aux_lndmi = dmi_aux.sel(time=slice('2019-01-01', '2023-03-01'))[lead:]
        lndmi = ax2.plot(dates, aux_lndmi, label='DMI', color='#289E64')

        lnn34 = ax2.plot(dates, n34[lead:], label='N34', color='#00C9ED')

        aux_lnsam = sam.mean_estimate.sel(
            time=slice('2019-01-01', '2023-03-01'))[lead:]
        lnsam = ax2.plot(dates, aux_lnsam, label='SAM', color='#005EFF')

        aux_lnasam = asam.mean_estimate.sel(
            time=slice('2019-01-01', '2023-03-01'))[lead:]
        lnasam = ax2.plot(dates, aux_lnasam, label='A-SAM', color='#960B00')

        aux_lnssam = ssam.mean_estimate.sel(
            time=slice('2019-01-01', '2023-03-01'))[lead:]
        lnssam = ax2.plot(dates, aux_lnssam, label='S-SAM', color='#FF0088')

        if correlaciones:
            try:
                c = {'Region': [t],
                     index.upper() + '_CM-DMI': [CorrSP(aux_lnscmap.values,
                                                        aux_lndmi.values)],
                     index.upper() + '_CM-N34': [CorrSP(aux_lnscmap.values,
                                                        n34[lead:])],
                     index.upper() + '_CM-SAM': [CorrSP(aux_lnscmap.values,
                                                        aux_lnsam.values)],
                     index.upper() + '_CM-S_SAM': [CorrSP(aux_lnscmap.values,
                                                          aux_lnssam.values)],
                     index.upper() + '_CM-A_SAM': [CorrSP(aux_lnscmap.values,
                                                          aux_lnasam.values)],

                     index.upper() + '_CH-DMI': [
                         CorrSP(aux_lnschirps.values[1:],
                                aux_lndmi.values[1:])],
                     index.upper() + '-CH-N34': [
                         CorrSP(aux_lnschirps.values[1:], n34[lead + 1:])],
                     index.upper() + '_CH-SAM': [
                         CorrSP(aux_lnschirps.values[1:],
                                aux_lnsam.values[1:])],
                     index.upper() + '-CHIRPS-S_SAM': [
                         CorrSP(aux_lnschirps.values[1:],
                                aux_lnssam.values[1:])],
                     index.upper() + '_CH-A_SAM': [
                         CorrSP(aux_lnschirps.values[1:],
                                aux_lnasam.values[1:])]
                     }

                if t == 'SESA':
                    c_df = pd.DataFrame(c)
                else:
                    c_df = pd.concat([c_df, pd.DataFrame(c)], axis=0)

                c2 = {'Region': [t],
                      index.upper() + '_CM-DMI': [
                          CorrSP(aux_lnscmap.values, aux_lndmi.values, True)[
                              1]],
                      index.upper() + '_CM-N34': [CorrSP(aux_lnscmap.values,
                                                         n34[lead:], True)[1]],
                      index.upper() + '_CM-SAM': [
                          CorrSP(aux_lnscmap.values, aux_lnsam.values, True)[
                              1]],
                      index.upper() + '_CM-S_SAM': [CorrSP(aux_lnscmap.values,
                                                           aux_lnssam.values,
                                                           True)[1]],
                      index.upper() + '_CM-A_SAM': [CorrSP(aux_lnscmap.values,
                                                           aux_lnasam.values,
                                                           True)[1]],

                      index.upper() + '_CH-DMI': [
                          CorrSP(aux_lnschirps.values[1:], aux_lndmi.values[1:],
                                 True)[1]],
                      index.upper() + '-CH-N34': [
                          CorrSP(aux_lnschirps.values[1:], n34[1 + lead:],
                                 True)[1]],
                      index.upper() + '_CH-SAM': [
                          CorrSP(aux_lnschirps.values[1:], aux_lnsam.values[1:],
                                 True)[1]],
                      index.upper() + '-CHIRPS-S_SAM': [
                          CorrSP(aux_lnschirps.values[1:],
                                 aux_lnssam.values[1:], True)[1]],
                      index.upper() + '_CH-A_SAM': [
                          CorrSP(aux_lnschirps.values[1:],
                                 aux_lnasam.values[1:], True)[1]]
                      }

                if t == 'SESA':
                    c2_df = pd.DataFrame(c2)
                else:
                    c2_df = pd.concat([c2_df, pd.DataFrame(c2)], axis=0)
            except:
                print('NaN values in ' + t + ' with Lead' + str(lead))
                data_frames = False


        ax.hlines(y=0, xmin=dates[0], xmax=dates[-1], color='gray')
        ax2.hlines(y=0, xmin=dates[0], xmax=dates[-1], color='gray')
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y %b'))
        ax.grid()
        ax.set_ylim((-.4, .5))
        ax2.set_ylim((-1.5, 7))
        ax.set_ylabel(index.upper(), fontsize=10)
        ax2.set_ylabel('índices', fontsize=10)
        ax.set_title(index.upper() + ' - ' + t + '\n' + 'Lead: ' + str(lead),
                     fontsize=15)

        lns = lnscmap + lnschirps + lndmi + lnn34 + lnsam + lnasam + lnssam
        ax.legend(lns, [index.upper() + '_CMAP2.5',
                        index.upper() + '_CHIRPS1',
                        'DMI', 'ONI',
                        'SAM', 'A-SAM', 'S-SAM'],
                  loc='upper left')

        if save:
            plt.savefig(out_dir + index.upper() + '_' + t + '_lead_' +
                        str(lead) +  '.jpg', dpi=dpi)
            plt.close('all')
        else:
            plt.show()


        fig = plt.figure(figsize=(10, 7), dpi=dpi)
        ax = fig.add_subplot(111)
        ax2 = ax.twinx()

        aux_pp_CMAP = data_anom.sel(
            lon=slice(ln[0], ln[1]), lat=slice(lt[1], lt[0]))
        aux_pp_CMAP = SameDateAs(aux_pp_CMAP, index_cmap)

        aux_pp_CHIRPS = data_anom_CHIRPS.sel(
            lon=slice(ln[0], ln[1]), lat=slice(lt[0], lt[1]))
        aux_pp_CHIRPS = SameDateAs(aux_pp_CHIRPS, index_chirps)

        lnrpsscmap = ax.plot(dates, aux.mean(['lon', 'lat']).mask[:-1],
                             color='#FF0003',
                             label=index.upper() + '_CMAP2.5', linewidth=2)
        lnrpsschirps = ax.plot(dates, aux2.mean(['lon', 'lat']).mask[:-1],
                               color='#FFA500',
                               label=index.upper() + '_CHIRPS1', linewidth=2)

        ln_pp_CMAP = ax2.plot(dates,
                              aux_pp_CMAP.mean(['lon', 'lat']).precip[:-1],
                              color='gray',
                              label='CMAP2.5', linewidth=2)

        ln_pp_CHIRPS = ax2.plot(dates,
                                aux_pp_CHIRPS.mean(['lon', 'lat']).precip[:-1],
                                color='k',
                                label='CHIPRS', linewidth=2)

        ax.hlines(y=0, xmin=dates[0], xmax=dates[-1], color='gray')
        ax2.hlines(y=0, xmin=dates[0], xmax=dates[-1], color='gray')
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y %b'))
        ax.grid()
        ax.set_ylim((-.4, .5))
        ax2.set_ylim((-50, 150))
        ax.set_ylabel('RPSS', fontsize=10)
        ax2.set_ylabel('PP anom. [mm]', fontsize=10)
        ax.set_title(index.upper() + ' - ' + t + '\n' + 'Lead: ' + str(lead),
                     fontsize=15)

        lns = lnrpsscmap + lnrpsschirps + ln_pp_CMAP + ln_pp_CHIRPS
        ax.legend(lns, [index.upper() + '_CMAP2.5', index.upper() +
                        '_CHIRPS1', 'CMAP', 'CHIRPS'], loc='upper left')

        if save:
            plt.savefig(out_dir + index.upper() + '_' + t + '_lead_' +
                        str(lead) + '_wpp_anoms.jpg', dpi=dpi)
            plt.close('all')
        else:
            plt.show()

    if data_frames:
        c_df.to_csv(out_dir + index.upper() + '_correlaciones2_lead_' +
                    str(lead) + '.txt', sep='\t', index=False, header=True)

        c2_df.to_csv(out_dir + index.upper() + '_correlaciones2_pvalue_lead_' +
                     str(lead) + '.txt', sep='\t', index=False, header=True)



# Open and set #################################################################
# CMAP
data = xr.open_dataset(cmap_data + 'pp_cmap.nc').__mul__(365/12) #dia
lon_cmap = data.lon.values
lat_cmap = data.lat.values
# tiene resolucion 2.5...

data = data.sel(lon=slice(270,330), lat=slice(20, -60))
# promedios trimestrales (VER)
data = data.rolling(time=3, center=True).mean()
#data = data.interp(lon=data_nmme.lon.values, lat=data_nmme.lat.values)
# climatologia trimestral
data_clim_f = data.sel(time=slice('1990-01-01', '2020-12-01'))
data_clim = data_clim_f.groupby('time.month').quantile([.33, .66], dim='time')

data_verif = data.sel(time=slice('2019-01-01','2023-04-01'))

data_anom = data.groupby('time.month') - \
            data_clim_f.groupby('time.month').mean()
data_anom = data_anom*MakeMask(data_anom, 'precip')

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

data_nmme25 = data_nmme.interp(lon=lon_cmap, lat=lat_cmap)
data_nmme25 = data_nmme25.sel(lon=slice(270,330), lat=slice(20, -60))

data_nmme25.load()

################################################################################
# chirps
data = xr.open_dataset(
    chirps_data + 'chirps_1990_2020_mmean.nc').__mul__(365/12) #dia
data = ChangeLons(data, 'longitude')
data = data.rename({'latitude':'lat'})
data = data.interp(lon=data_nmme.lon,
                   lat=np.arange(-50,50,1))

data = data.sel(lon=slice(270,330), lat=slice(-60, 20))
# promedios trimestrales (VER)
data = data.rolling(time=3, center=True).mean()
# climatologia trimestral
data_clim_f_ch = data
data_clim_ch = data.groupby('time.month').quantile([.33, .66], dim='time')

# 2019-2023
data = xr.open_dataset(
    chirps_data + 'chirps_2019_2023_mmean.nc').__mul__(365/12) #dia
data = ChangeLons(data, 'longitude')
data = data.rename({'latitude':'lat'})
data = data.interp(lon=data_nmme.lon, lat=np.arange(-50,50,1))
# promedios trimestrales (VER)
data_verif_ch = data.rolling(time=3, center=True).mean()
data_verif_ch = data_verif_ch.sel(lon=slice(270,330), lat=slice(-60, 20))

data_anom_CHIRPS = data.groupby('time.month') - \
                   data_clim_f_ch.groupby('time.month').mean()
data_anom_CHIRPS = data_anom_CHIRPS*MakeMask(data_anom_CHIRPS, 'precip')

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
sam = xr.open_dataset(out_dir + 'sam.nc')
asam = xr.open_dataset(out_dir + 'asam.nc')
ssam = xr.open_dataset(out_dir + 'ssam.nc')
################################################################################

# En que categoría está lo observado según la climatología
print('Categorias observadas...')
abn_cmap = ABNobs(data_clim, data_verif)
abn_chirps = ABNobs(data_clim_ch, data_verif_ch)

# RPS ##########################################################################
# RPSo ------------------------------------------------------------------------#
print('RPS observado...')
RPSo_cmap = RPSO(abn_cmap, data_verif)
RPSo_chirps = RPSO(abn_chirps, data_verif_ch)

# BS ##########################################################################
# BSo--------------------------------------------------------------------------#
print('BS observado')
BSo_cmap = BSO(abn_cmap, data_verif)
BSo_chirps = BSO(abn_chirps, data_verif_ch)

for l in lead:
    print('Lead: ' + str(l))
    # RPS #######################################################################
    # RPSf  -------------------------------------------------------------------#
    print('RPS forecast...')
    print('CMAP')
    RPSf_cmap = RPSF(abn_cmap, data_nmme25, data_verif, l)
    print('CHIRPS')
    RPSf_chirps = RPSF(abn_chirps, data_nmme, data_verif_ch, l)

    ############################################################################
    # RPSS
    print('RPSS...')
    # SameDateAs es en el caso de usar lead ! = 0 RPSf_* pierde los primeros
    # lead tiempos y esta resta cuenta no puede hacerce
    RPSS_cmap = 1 - (RPSf_cmap/SameDateAs(RPSo_cmap, RPSf_cmap))
    RPSS_chirps = 1 - (RPSf_chirps/SameDateAs(RPSo_chirps, RPSf_chirps))

    # para plotear
    RPSS_cmap_nm = RPSS_cmap
    RPSS_chirps_nm = RPSS_chirps
    #
    RPSS_cmap = RPSS_cmap*MakeMask(RPSS_cmap)
    RPSS_chirps = RPSS_chirps*MakeMask(RPSS_chirps)

    # BS #######################################################################
    # BSf ---------------------------------------------------------------------#
    print('BS forecast')
    print('CMAP')
    BSf_cmap = BSF(abn_cmap, data_nmme25, data_verif, l)
    print('CHIRPS')
    BSf_chirps = BSF(abn_chirps, data_nmme, data_verif_ch, l)

    # BSS ---------------------------------------------------------------------#
    print('BSS...')
    aux_cmap = SameDateAs(BSo_cmap, BSf_cmap)
    BSS_cmap = - (BSf_cmap - aux_cmap) / aux_cmap
    aux_chirps = SameDateAs(BSo_chirps, BSf_chirps)
    BSS_chirps = - (BSf_chirps - aux_chirps) / aux_chirps

    # para plotear
    BSS_cmap_nm = BSS_cmap
    BSS_chirps_nm = BSS_chirps
    #
    BSS_cmap = BSS_cmap*MakeMask(BSS_cmap)
    BSS_chirps = BSS_chirps*MakeMask(BSS_chirps)

    ############################################################################
    ComputeAndPlot('rpss', True, dpi, save, l,test)
    ComputeAndPlot('bss', True, dpi, save, l, test)
    ############################################################################

    if plot_mapas:
        print('Mapas RPSS y BSS...')
        cbar = colors.ListedColormap(
            ['#5FFE9B', '#FEB77E', '#FE675C', '#CA3E72',
             '#782280', '#251255'])
        cbar.set_over('#251254')
        cbar.set_under('#5FFE9B')
        cbar.set_bad(color='white')

        rpss = [RPSS_cmap_nm, RPSS_chirps, BSS_cmap_nm, BSS_chirps]
        dataset = ['CMAP', 'CHIRPS', 'CMAP', 'CHIRPS']
        index = ['RPSS', 'RPSS', 'BSS', 'BSS']
        anio = [2019, 2020, 2021, 2022, 2023]
        # ---------------------------------------------------------------------#
        # Por periodos --------------------------------------------------------#
        for ds, ds_name, i in zip(rpss, dataset, index):
            for a in anio:
                try:
                    aux = ds.sel(time=ds.time.dt.year.isin(a)).mean('time')

                    titulo = i + ' - ' + str(a) + ' - ' + ds_name + ' Lead: ' +\
                             str(l)
                    name_fig = DirAndFile(out_dir, dir_results, 'MAPA',
                                   [str(a), i, 'Lead', str(l)])

                    try:
                        Plot(aux, aux['mask'], [-.1, 0, .1, .2, .3, .4, .5],
                             save,
                             dpi, titulo,
                             name_fig, 'gray', cbar)
                    except:
                        Plot(aux, aux, [-.1, 0, .1, .2, .3, .4, .5], save,
                             dpi, titulo, name_fig, 'k', cbar)

                except:
                    pass
        # All 7/2019-4/2023 ---------------------------------------------------#
        for ds, ds_name, i in zip(rpss, dataset, index):
            try:
                aux = ds.sel(time=slice('2019-07-01', '2023-12-01')).\
                    mean('time')
                titulo = i + ' - ' + '7/2019 - 3/2023' + ' - ' + ds_name  +\
                         ' Lead: ' + str(l)
                name_fig = DirAndFile(out_dir, dir_results, 'MAPA',
                                      ['2019_2023', i, 'Lead', str(l)])

                try:
                    Plot(aux, aux['mask'], [-.1, 0, .1, .2, .3, .4, .5], save,
                         dpi, titulo,
                         name_fig, 'gray', cbar )
                except:
                    Plot(aux, aux, [-.1, 0, .1, .2, .3, .4, .5], save,
                         dpi, titulo, name_fig, 'k', cbar)

            except:
                pass
################################################################################
################################################################################
print('done')
################################################################################
