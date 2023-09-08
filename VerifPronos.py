"""
Verifica los pronosticos desde 2019 hasta el presente con los índices RPSS
y BSS.

Se actualiza automaticamente cuando este disponible y/o sea necesario
La actualización de las bases de datos puede desactivarse con update = False

endtime_select = -1 | para verificar hasta el último trimestre disponible
en el índice ONI. Es posible que no se pueda verificar este trimestre
dependiendo de CMAP y CHIRPS. En ese caso, sólo se hará el grafico de
evolucion temporal del skill. En este caso el valor mas reciente de BSS no se
debe tener en cuenta.

endtime_select = -n verifica hasta el trimestre -n donde deberia estar tdo
actualizado.

los índices se actualizan en set_indices sólo si es necesario
los pronos del nmme se prueba descargar siempre que no esté ya descargo el
pronostico del mes actual o el archivo de ese mes sea menor a 1mb (vacio o
con algun error)
"""
################################################################################
save = False
update = True
plot_mapas = True
# para computar tdo hasta donde se pueda verificar, endtime_select = -2 -------#
endtime_select = -1
################################################################################
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
from matplotlib import colors
from Funciones import SelectFilesNMME, MakeMask, ChangeLons, ABNobs, \
    RPSO, RPSF, BSO, BSF, Plot, SameDateAs, DirAndFile, \
    CreateDirectory, OpenRegiones, Correlaciones
import set_indices, cmap, chirps, nmme_update
from dateutil.relativedelta import relativedelta
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
warnings.filterwarnings("ignore")

################################################################################
test = False # = True solo va computar una region
lead = [0, 1, 2, 3]
#------------------------------------------------------------------------------#
# Si explicitamente no se quiere plotear algun mapa, cambiar a False
plot_mapas_cmap = plot_mapas
correlaciones_cmap = plot_mapas
plot_mapas_chirps = plot_mapas
correlaciones_chirps = plot_mapas
#------------------------------------------------------------------------------#
if save:
    dpi = 300
    CreateDirectory(out_dir, dir_results)
else:
    dpi = 100
# Funciones ####################################################################
def ComputeAndPlot(index, correlaciones_cmap, correlaciones_chirps,
                   dpi, save, lead=0, test=False, region_name='regiones_sa'):
    dates = pd.date_range(start='2018-12-01', end=endtime,
                          freq='M') + pd.DateOffset(days=1)
    dates = dates[lead:]

    titulos, lon_regiones, lat_regiones = OpenRegiones(region_name + '.csv')

    try:
        if test:
            lon_regiones = lon_regiones[0]
            lat_regiones = lat_regiones[0]
            titulos = titulos[0]
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
    for ln, lt, t in zip(lon_regiones, lat_regiones, titulos):
        data_frames = True
        aux = index_cmap.sel(lon=slice(ln[0], ln[1]), lat=slice(lt[1], lt[0]))
        aux2 = index_chirps.sel(lon=slice(ln[0], ln[1]), lat=slice(lt[1],lt[0]))

        fig = plt.figure(figsize=(10, 7), dpi=dpi)
        ax = fig.add_subplot(111)
        ax2 = ax.twinx()

        # index score
        aux_lnscmap = aux.mean(['lon', 'lat']).mask
        lnscmap = ax.plot(dates, aux_lnscmap,
                             color='#FF0003',
                             label=index.upper() + '_CMAP2.5', linewidth=2)

        aux_lnschirps = aux2.mean(['lon', 'lat']).mask
        lnschirps = ax.plot(dates, aux_lnschirps,
                               color='#FFA500',
                               label=index.upper() + '_CHIRPS1', linewidth=2)

        # indices
        aux_lndmi = dmi[lead:]
        lndmi = ax2.plot(dates, aux_lndmi, label='DMI', color='#289E64')

        aux_lnn34 = n34[lead:]
        lnn34 = ax2.plot(dates, aux_lnn34, label='N34', color='#00C9ED')

        aux_lnsam = sam[lead:]
        lnsam = ax2.plot(dates, aux_lnsam, label='SAM', color='#005EFF')

        aux_lnasam = asam[lead:]
        lnasam = ax2.plot(dates, aux_lnasam, label='A-SAM', color='#960B00')

        aux_lnssam = ssam[lead:]
        lnssam = ax2.plot(dates, aux_lnssam, label='S-SAM', color='#FF0088')

        c_cmap, c2_cmap, c_chirps, c2_chirps = \
            Correlaciones(correlaciones_cmap, correlaciones_chirps, aux_lnscmap,
                          aux_lnschirps, lead, aux_lnn34, aux_lndmi, aux_lnsam,
                          aux_lnssam, aux_lnasam, index, t)

        if t == 'SESA':
            c_df_cmap = pd.DataFrame(c_cmap)
        else:
            c_df_cmap = pd.concat([c_df_cmap, pd.DataFrame(c_cmap)], axis=0)

        if t == 'SESA':
            c2_df_cmap = pd.DataFrame(c2_cmap)
        else:
            c2_df_cmap = pd.concat([c2_df_cmap, pd.DataFrame(c2_cmap)], axis=0)

        if t == 'SESA':
            c_df_chirps = pd.DataFrame(c_chirps)
        else:
            c_df_chirps = pd.concat([c_df_chirps, pd.DataFrame(c_chirps)],
                                    axis=0)

        if t == 'SESA':
            c2_df_chirps = pd.DataFrame(c2_chirps)
        else:
            c2_df_chirps = pd.concat([c2_df_chirps, pd.DataFrame(c2_chirps)],
                                     axis=0)


        ax.hlines(y=0, xmin=dates[0], xmax=dates[-1], color='gray')
        ax2.hlines(y=0, xmin=dates[0], xmax=dates[-1], color='gray')
        #ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y %b'))

        plt.rcParams['date.converter'] = 'concise'
        ax.xaxis.set_major_locator(
            mdates.AutoDateLocator(minticks=20, maxticks=26))
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

        ########################################################################
        fig = plt.figure(figsize=(10, 7), dpi=dpi)
        ax = fig.add_subplot(111)
        ax2 = ax.twinx()

        aux_pp_CMAP = data_anom.sel(
            lon=slice(ln[0], ln[1]), lat=slice(lt[1], lt[0]))
        aux_pp_CMAP = SameDateAs(aux_pp_CMAP, index_cmap)

        aux_pp_CHIRPS = data_anom_CHIRPS.sel(
            lon=slice(ln[0], ln[1]), lat=slice(lt[0], lt[1]))
        aux_pp_CHIRPS = SameDateAs(aux_pp_CHIRPS, index_chirps)

        lnrpsscmap = ax.plot(dates, aux.mean(['lon', 'lat']).mask,
                             color='#FF0003',
                             label=index.upper() + '_CMAP2.5', linewidth=2)
        lnrpsschirps = ax.plot(dates, aux2.mean(['lon', 'lat']).mask,
                               color='#FFA500',
                               label=index.upper() + '_CHIRPS1', linewidth=2)

        ln_pp_CMAP = ax2.plot(dates,
                              aux_pp_CMAP.mean(['lon', 'lat']).precip,
                              color='gray',
                              label='CMAP2.5', linewidth=2)

        ln_pp_CHIRPS = ax2.plot(dates,
                                aux_pp_CHIRPS.mean(['lon', 'lat']).precip,
                                color='k',
                                label='CHIPRS', linewidth=2)

        ax.hlines(y=0, xmin=dates[0], xmax=dates[-1], color='gray')
        ax2.hlines(y=0, xmin=dates[0], xmax=dates[-1], color='gray')

        plt.rcParams['date.converter'] = 'concise'
        ax.xaxis.set_major_locator(
            mdates.AutoDateLocator(minticks=20, maxticks=26))

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

    if correlaciones_cmap:
        c_df_cmap.to_csv(out_dir + index.upper() +
                         '_CMAP_correlaciones2_lead_' + str(lead) + '.txt',
                         sep='\t', index=False, header=True)

        c2_df_cmap.to_csv(out_dir + index.upper() +
                      '_CMAP_correlaciones2_pvalue_lead_' + str(lead) + '.txt',
                      sep='\t', index=False, header=True)

    if correlaciones_chirps:
        c_df_chirps.to_csv(out_dir + index.upper() +
                           '_chirps_correlaciones2_lead_' + str(lead) + '.txt',
                           sep='\t', index=False, header=True)

        c2_df_chirps.to_csv(out_dir + index.upper() +
                            '_chirps_correlaciones2_pvalue_lead_' + str(lead) +
                            '.txt', sep='\t', index=False, header=True)

################################################################################
print('Set Indices ###########################################################')
n34, dmi, sam, ssam, asam, endtime = set_indices.compute()

# ENDTIME ######################################################################
# endtime determinado por el ONI, se actualiza al trimestre anterior
# e.g. al finalizar agosto actualiza ONI en JJA --> mm = 7
endtime = n34.time.values[endtime_select]
n34 = n34.oni
if endtime_select != -1:
    n34 = n34[:endtime_select+1]
    time0 = dmi.time.values[0]
    dmi = dmi.sel(time=slice(time0, endtime))
    sam = SameDateAs(sam, dmi)
    asam = SameDateAs(asam, dmi)
    ssam = SameDateAs(ssam, dmi)
    update = False
print('#######################################################################')
print('<<<<<<<<<<<<<<<<<<< Verificación hasta: ' + str(endtime).split('T')[0] +
      ' >>>>>>>>>>>>>>>>>>>>')
print('#######################################################################')
# para identificar el prono correspondiente
anio = endtime.astype('datetime64[Y]').astype(int) + 1970
mes = endtime.astype('datetime64[M]').astype(int) % 12 + 1
endtime_str = f"{anio}{mes:02d}"
# y el ultimo target correspondiente
aux = endtime.astype('M8[D]').astype('O')
targetime = aux + relativedelta(months=6)
targetime = np.datetime64(targetime)

print('Open and set ##########################################################')
print('CMAP ------------------------------------------------------------------')
data = xr.open_dataset(cmap_data + 'pp_cmap.nc').__mul__(365/12) #dia
lon_cmap = data.lon.values
lat_cmap = data.lat.values
# tiene resolucion 2.5...

# El tiempo final de los datos debe ser un trimestre mas que el endtime del oni
# para poder calcular el promedio trimestral
# si el ONI tiene mm = 7 --> JJA,
# y el ultimo tiempo de CMAP  es del mes 7, el trimestre es MJJ
if endtime == data.time.values[-1]:
    print('CMAP desactualizado')
    if update:
        print('Intentando actualizar CMAP')
        cmap.update()
        data = xr.open_dataset(cmap_data + 'pp_cmap.nc').__mul__(
            365 / 12)  # dia
        if endtime == data.time.values[-1]:
            print('Actualizacion CMAP no disponible')
            print('Se computará la evolución temporal solamente')
            correlaciones_cmap = False
            plot_mapas_cmap = False
    else:
        correlaciones_cmap = False
        plot_mapas_cmap = False

data = data.sel(lon=slice(270,330), lat=slice(20, -60))
# promedios trimestrales (VER)
data = data.rolling(time=3, center=True).mean()
#data = data.interp(lon=data_nmme.lon.values, lat=data_nmme.lat.values)
# climatologia trimestral
data_clim_f = data.sel(time=slice('1990-01-01', '2020-12-01'))
data_clim = data_clim_f.groupby('time.month').quantile([.33, .66], dim='time')

data_verif = data.sel(time=slice('2019-01-01', endtime))

data_anom = data.groupby('time.month') - \
            data_clim_f.groupby('time.month').mean()
data_anom = data_anom*MakeMask(data_anom, 'precip')

################################################################################
print('NMME ------------------------------------------------------------------')
nmme_update.update()
files = SelectFilesNMME(nmme_pronos, 'prate', False)
# ultimo pronostico que que puede ser verificado
posf = [i for i, prono in enumerate(files) if endtime_str in prono][0]
files = files[0:posf+1]
data_nmme = xr.open_mfdataset(files, decode_times=False, engine='netcdf4',
                              combine='nested', concat_dim='initial_time')
data_nmme = data_nmme.rename({'initial_time':'time'}) # para mas adelante
data_nmme['time'] = pd.date_range(start='2018-12-01', end=endtime,
                                  freq='M') + pd.DateOffset(days=1)
data_nmme['target'] = pd.date_range(start='2018-12-01', end=targetime,
                                  freq='M') + pd.DateOffset(days=1)

data_nmme25 = data_nmme.interp(lon=lon_cmap, lat=lat_cmap)
data_nmme25 = data_nmme25.sel(lon=slice(270,330), lat=slice(20, -60))

data_nmme25.load()
################################################################################
print('CHIRPS ----------------------------------------------------------------')
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

# Ideem que con CMAP
if endtime == data.time.values[-1]:
    print('CHIRPS desactualizado')
    if update:
        print('Intentando actualizar CHIRPS')
        chirps.update()
        data = xr.open_dataset(
            chirps_data + 'chirps_2019_2023_mmean.nc').__mul__(365 / 12)  # dia
        if endtime == data.time.values[-1]:
            print('Actualizacion CHIRPS no disponible')
            print('Se computará la evolución temporal solamente')
            correlaciones_chirps = False
            plot_mapas_chirps = False
    else:
        correlaciones_cmap = False
        plot_mapas_cmap = False

data = ChangeLons(data, 'longitude')
data = data.rename({'latitude':'lat'})
data = data.interp(lon=data_nmme.lon, lat=np.arange(-50,50,1))
# promedios trimestrales (VER)
data_verif_ch = data.rolling(time=3, center=True).mean()
data_verif_ch = data_verif_ch.sel(lon=slice(270,330), lat=slice(-60, 20))
data_verif_ch = data_verif_ch.sel(time=slice('2019-01-01', endtime))

aux = data.rolling(time=3, center=True).mean()
data_anom_CHIRPS = aux.groupby('time.month') - \
                   data_clim_f_ch.groupby('time.month').mean()
data_anom_CHIRPS = data_anom_CHIRPS*MakeMask(data_anom_CHIRPS, 'precip')
data_anom_CHIRPS = data_anom_CHIRPS.sel(time=slice('2019-01-01', endtime))

print('#######################################################################')
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
    ComputeAndPlot('rpss', correlaciones_cmap, correlaciones_chirps, dpi, save,
                   l,test, 'regiones_sa')
    ComputeAndPlot('bss', correlaciones_cmap, correlaciones_chirps, dpi, save,
                   l, test, 'regiones_sa')
    ############################################################################

    if plot_mapas & (plot_mapas_cmap | plot_mapas_chirps):
        print('Mapas RPSS y BSS...')
        cbar = colors.ListedColormap(
            ['white', '#FEB77E', '#FE675C', '#EB0330', '#CA3E72',
             '#782280', '#251255'])
        cbar.set_over('#251254')
        cbar.set_under('white')
        cbar.set_bad(color='white')

        anio = [2019, 2020, 2021, 2022, 2023]

        if plot_mapas_chirps & plot_mapas_cmap:
            rpss = [RPSS_cmap_nm, RPSS_chirps, BSS_cmap_nm, BSS_chirps]
            dataset = ['CMAP', 'CHIRPS', 'CMAP', 'CHIRPS']
            index = ['RPSS', 'RPSS', 'BSS', 'BSS']
        elif plot_mapas_chirps & (plot_mapas_cmap==False):
            rpss = [RPSS_chirps, BSS_chirps]
            dataset = ['CHIRPS', 'CHIRPS']
            index = ['RPSS', 'BSS']
        elif (plot_mapas_chirps==False) & plot_mapas_cmap:
            rpss = [RPSS_cmap, BSS_cmap]
            dataset = ['CMAP', 'CMAP']
            index = ['RPSS', 'BSS']

        # ---------------------------------------------------------------------#
        # Por periodos --------------------------------------------------------#
        for ds, ds_name, i in zip(rpss, dataset, index):
            for a in anio:
                try:
                    aux = ds.sel(time=ds.time.dt.year.isin(a)).mean('time')

                    titulo = i + ' - ' + str(a) + ' - ' + ds_name + ' Lead: ' +\
                             str(l)
                    name_fig = DirAndFile(out_dir, dir_results, 'MAPA',
                                   [str(a), i, ds_name, 'Lead', str(l)])

                    try:
                        Plot(aux, aux['mask'],
                             [-.05, 0, .05, .1, .15, .25, .30, .35], save,
                             dpi, titulo, name_fig, 'gray', cbar)
                    except:
                        Plot(aux, aux, [-.05, 0, .05, .1, .15, .25, .30, .35],
                             save, dpi, titulo, name_fig, 'k', cbar)

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
                                      ['2019_2023', i, ds_name, 'Lead', str(l)])

                try:
                    Plot(aux, aux['mask'],
                         [-.05, 0, .05, .1, .15, .25, .30, .35], save,
                         dpi, titulo, name_fig, 'gray', cbar )
                except:
                    Plot(aux, aux,
                         [-.05, 0, .05, .1, .15, .25, .30, .35], save,
                         dpi, titulo, name_fig, 'k', cbar)

            except:
                pass
################################################################################
################################################################################
print('done')
################################################################################
