"""
Funciones generales para ENSO_IOD
"""
from itertools import groupby
import xarray as xr
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
import statsmodels.formula.api as sm
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import regionmask
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import glob
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

out_dir = '~/'
################################################################################
def SelectFilesNMME(dir, variable, size_check):
    files =  glob.glob(dir + variable + '_*_prob.adj.sea.nc')

    # Solo archivos con tamaño mayor a 0
    if size_check:
        files = [file for file in files if os.path.getsize(file) > 1]

    return sorted(files, key=lambda x: x.split()[0])
# -----------------------------------------------------------------------------#
def MakeMask(DataArray, dataname='mask'):
    """
    :param DataArray: campo que se quiere enmascarar (admite multiples tiempos)
    :param dataname: nombre que va recibir la variable
    :return: el dataarray enmascarada
    """
    import regionmask
    mask = regionmask.defined_regions.natural_earth_v5_0_0.countries_110.\
        mask(DataArray)
    mask = xr.where(np.isnan(mask), mask, 1)
    mask = mask.to_dataset(name=dataname)
    return mask

# -----------------------------------------------------------------------------#
def Detrend(xrda, dim):
    aux = xrda.polyfit(dim=dim, deg=1)
    try:
        trend = xr.polyval(xrda[dim], aux.var_polyfit_coefficients)
    except:
        trend = xr.polyval(xrda[dim], aux.polyfit_coefficients)
    dt = xrda - trend
    return dt
# -----------------------------------------------------------------------------#

def ChangeLons(data, lon_name='lon'):
    """
    Cambia el formato de las latitudes
    :param data:
    :param lon_name:
    :return:
    """
    data['_longitude_adjusted'] = xr.where(
        data[lon_name] < 0,
        data[lon_name] + 360,
        data[lon_name])

    data = (
        data
            .swap_dims({lon_name: '_longitude_adjusted'})
            .sel(**{'_longitude_adjusted': sorted(data._longitude_adjusted)})
            .drop(lon_name))

    data = data.rename({'_longitude_adjusted': 'lon'})

    return data
#------------------------------------------------------------------------------#
def Weights(data):
    weights = np.transpose(np.tile(np.cos(data.lat * np.pi / 180),
                                   (len(data.lon), 1)))
    data_w = data * weights
    return data_w

#------------------------------------------------------------------------------#
def LeadMonth(date, lead, suma=False):
    """
    suma o resta un mes (lead) a date
    :param date:
    :param lead:
    :return:
    """
    if suma:
        date2 = pd.Timestamp(date) + pd.DateOffset(months=lead)
    else:
        date2 = pd.Timestamp(date) - pd.DateOffset(months=lead)

    date2 = date2.to_numpy()
    return date2

#------------------------------------------------------------------------------#
def SameDateAs(data, datadate):
    """
    En data selecciona las mismas fechas que datadate
    :param data:
    :param datadate:
    :return:
    """
    return data.sel(time=datadate.time.values)

#------------------------------------------------------------------------------#
def ABNobs(data_clim, data_verif):
    """
    Categoriza las anomalias de una variable en 3 terciles
    :param data_clim:
    :param data_verif:
    :return:
    """
    data2 = xr.zeros_like(data_verif)
    first = True
    for d in data_verif.time.values:
        m = pd.to_datetime(d).month

        aux = xr.where(data_clim.sel(month=m, quantile=0.33)
                       < data_verif.sel(time=d), data_verif.sel(time=d), np.NAN)

        n = xr.where(aux < data_clim.sel(month=m, quantile=0.66), 1, 0)

        b = xr.where(data_verif.sel(time=d) <
                     data_clim.sel(month=m, quantile=0.33), 1, 0)

        a = xr.where(data_verif.sel(time=d) >
                     data_clim.sel(month=m, quantile=0.66), 1, 0)

        if first:
            first = False
            N = n
            B = b
            A = a
        else:
            N = xr.concat([N, n], dim='time')
            B = xr.concat([B, b], dim='time')
            A = xr.concat([A, a], dim='time')

    data2 = data2.assign(A=A.precip)
    data2 = data2.assign(B=B.precip)
    data2 = data2.assign(N=N.precip)
    data2 = data2.drop('quantile')
    return data2

#------------------------------------------------------------------------------#
def RPSO(data, data_verif):
    """
    Calculo del índice RPS observado
    :param data:
    :param data_verif:
    :return:
    """
    clim_prob = (.33333, .33334, .33333)
    first = True
    for d in data_verif.time.values:
        aux = data.sel(time=d)
        aux_bn = aux.B + aux.N
        aux_total = aux.B + aux.N + aux.A

        aux_rps = (clim_prob[0] - aux.B) ** 2 + \
                  (np.sum(clim_prob[0:2]) - aux_bn) ** 2 + \
                  (np.sum(clim_prob) - aux_total) ** 2

        if first:
            first = False
            RPSo = aux_rps
        else:
            RPSo = xr.concat([RPSo, aux_rps], dim='time')

    RPSo = xr.where(RPSo >= sum(np.cumsum(clim_prob) ** 2), np.NAN, RPSo)
    return RPSo

#------------------------------------------------------------------------------#
def RPSF(data, data_prono, data_verif, lead=0):
    """
    Calculo del RPS forecast
    :param data:
    :param data_prono:
    :param data_verif:
    :return:
    """
    first = True
    for d in data_verif.time.values:
        aux = data.sel(time=d)
        aux_bn = aux.B + aux.N
        aux_total = aux.B + aux.N + aux.A

        if lead != 0:
            d0 = LeadMonth(d, lead)
        else:
            d0 = d
        try:
            auxf = data_prono.sel(time=d0, target=d)
            auxf_bn = auxf.prob_below + auxf.prob_norm
            auxf_total = auxf.prob_below + auxf.prob_norm + auxf.prob_above

            aux_rps = (auxf.prob_below - aux.B) ** 2 + \
                      (auxf_bn - aux_bn) ** 2 + \
                      (auxf_total - aux_total) ** 2

            if first:
                first = False
                RPSf = aux_rps
            else:
                RPSf = xr.concat([RPSf, aux_rps], dim='time')

        except:
            print('Skip ' + np.datetime_as_string(d0, unit='M') + ' con' +
                  ' target ' + np.datetime_as_string(d, unit='M'))

    # para poder usar la dim time igual que las otras
    if lead != 0:
        RPSf = RPSf.rename({'target':'time'})

    return RPSf

#------------------------------------------------------------------------------#
def BSO(data, data_verif):
    """
    Calculo del índice BS observado
    :param data:
    :param data_verif:
    :return:
    """
    clim_prob = (.33333, .33334, .33333)
    first = True
    for d in data_verif.time.values:
        aux = data.sel(time=d)

        aux_bs = (clim_prob[0] - aux.B) ** 2 + \
                  (clim_prob[1] - aux.N) ** 2 + \
                  (clim_prob[2] - aux.A) ** 2

        if first:
            first = False
            BSo = aux_bs
        else:
            BSo = xr.concat([BSo, aux_bs], dim='time')

    BSo = xr.where(BSo >= sum(np.cumsum(clim_prob) ** 2), np.NAN, BSo)
    return BSo

#------------------------------------------------------------------------------#
def BSF(data, data_prono, data_verif, lead=0):
    """
    Calculo del índice BC forecast
    :param data:
    :param data_prono:
    :param data_verif:
    :return:
    """
    first = True
    for d in data_verif.time.values:
        aux = data.sel(time=d)

        if lead != 0:
            d0 = LeadMonth(d, lead)
        else:
            d0 = d

        try:
            auxf = data_prono.sel(time=d0, target=d)

            aux_bs = (auxf.prob_below - aux.B) ** 2 + \
                     (auxf.prob_norm - aux.N) ** 2 + \
                     (auxf.prob_above - aux.A) ** 2

            if first:
                first = False
                BSF = aux_bs
            else:
                BSF = xr.concat([BSF, aux_bs], dim='time')

        except:
            print('Skip ' + np.datetime_as_string(d0, unit='M') + ' con' +
                  ' target ' + np.datetime_as_string(d, unit='M'))

    # para poder usar la dim time igual que las otras
    if lead != 0:
        BSF = BSF.rename({'target':'time'})
    return BSF

#------------------------------------------------------------------------------#
def CorrSP(serie1, serie2, pvalue=False):
    """
    Correlación entre dos series temporales
    :param serie1:
    :param serie2:
    :param pvalue:
    :return:
    """
    from scipy.stats import pearsonr
    if pvalue:
        return np.round(pearsonr(serie1, serie2), 2)
    else:
        return np.round(pearsonr(serie1, serie2)[0],2)

#------------------------------------------------------------------------------#
def Plot(comp, comp_var, levels, save, dpi, title, name_fig,
         color_map, cmap):

    import matplotlib.pyplot as plt
    import cartopy.feature
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    import cartopy.crs as ccrs
    fig_size = (5, 6)
    extent = [270, 330, -60, 20]
    xticks = np.arange(270, 330, 10)
    yticks = np.arange(-60, 40, 20)
    crs_latlon = ccrs.PlateCarree()

    levels_contour = levels.copy()
    if isinstance(levels, np.ndarray):
        levels_contour = levels[levels != 0]
    else:
        levels_contour.remove(0)


    fig = plt.figure(figsize=fig_size, dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    ax.set_extent(extent, crs=crs_latlon)
    im = ax.contourf(comp.lon, comp.lat, comp_var, levels=levels,
                     transform=crs_latlon, cmap=cmap, extend='both')

    cb = plt.colorbar(im, fraction=0.042, pad=0.035, shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    ax.add_feature(cartopy.feature.LAND, facecolor='white', edgecolor=color_map)
    ax.add_feature(cartopy.feature.COASTLINE)
    #ax.add_feature(cartopy.feature.LAKES, color='skyblue')
    ax.add_feature(cartopy.feature.BORDERS, linestyle='-', color=color_map)
    #ax.add_feature(cartopy.feature.RIVERS, edgecolor='skyblue')
    ax.add_feature(cartopy.feature.STATES)
    ax.coastlines(color=color_map, linestyle='-', alpha=1)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(xticks, crs=crs_latlon)
    ax.set_yticks(yticks, crs=crs_latlon)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    ax.tick_params(labelsize=7)
    plt.title(title, fontsize=10)
    plt.tight_layout()

    if save:
        plt.savefig(name_fig + '.jpg')
        plt.close()
    else:
        plt.show()
#------------------------------------------------------------------------------#

# Niño3.4 & DMI ########################################################################################################
def MovingBasePeriodAnomaly(data, start=1920, end=2023):
    import xarray as xr
    # first five years
    start_num = start
    start = str(start)

    initial = data.sel(time=slice(start + '-01-01', str(start_num + 5) + '-12-31')).groupby('time.month') - \
              data.sel(time=slice(str(start_num - 14) + '-01-01', str(start_num + 5 + 10) + '-12-31')).groupby(
                  'time.month').mean('time')


    start_num = start_num + 6
    result = initial

    while (start_num != end-4) & (start_num < end-4):

        aux = data.sel(time=slice(str(start_num) + '-01-01', str(start_num + 4) + '-12-31')).groupby('time.month') - \
              data.sel(time=slice(str(start_num - 15) + '-01-01', str(start_num + 4 + 10) + '-12-31')).groupby(
                  'time.month').mean('time')

        start_num = start_num + 5

        result = xr.concat([result, aux], dim='time')

    if start_num > end - 4:
        start_num = start_num - 5

    aux = data.sel(time=slice(str(start_num) + '-01-01', str(start_num + 4) + '-12-31')).groupby('time.month') - \
          data.sel(time=slice(str(end-29) + '-01-01', str(end) + '-12-31')).groupby('time.month').mean('time')

    result = xr.concat([result, aux], dim='time')

    return (result)

def Nino34CPC(data, start=1920, end=2020):

    # Calculates the Niño3.4 index using the CPC criteria.
    # Use ERSSTv5 to obtain exactly the same values as those reported.

    #from Funciones import MovingBasePeriodAnomaly

    start_year = str(start-14)
    end_year = str(end)
    sst = data
    # N34
    ninio34 = sst.sel(lat=slice(4.0, -4.0), lon=slice(190, 240), time=slice(start_year+'-01-01', end_year + '-12-31'))
    ninio34 = ninio34.sst.mean(['lon', 'lat'], skipna=True)

    # compute monthly anomalies
    ninio34 = MovingBasePeriodAnomaly(data=ninio34, start=start, end=end)

    # compute 5-month running mean
    ninio34_filtered = np.convolve(ninio34, np.ones((3,)) / 3, mode='same')  #
    ninio34_f = xr.DataArray(ninio34_filtered, coords=[ninio34.time.values], dims=['time'])

    aux = abs(ninio34_f) > 0.5
    results = []
    for k, g in groupby(enumerate(aux.values), key=lambda x: x[1]):
        if k:
            g = list(g)
            results.append([g[0][0], len(g)])

    n34 = []
    n34_df = pd.DataFrame(columns=['N34', 'Años', 'Mes'], dtype=float)
    for m in range(0, len(results)):
        # True values
        len_true = results[m][1]

        # True values for at least 5 consecutive seasons
        if len_true >= 5:
            a = results[m][0]
            n34.append([np.arange(a, a + results[m][1]), ninio34_f[np.arange(a, a + results[m][1])].values])

            for l in range(0, len_true):
                if l < (len_true - 2):
                    main_month_num = results[m][0] + 1 + l
                    if main_month_num != 1210:
                        n34_df = n34_df.append({'N34': np.around(ninio34_f[main_month_num].values, 2),
                                            'Años': np.around(ninio34_f[main_month_num]['time.year'].values),
                                            'Mes': np.around(ninio34_f[main_month_num]['time.month'].values)},
                                           ignore_index=True)

    return ninio34_f, n34, n34_df

def DMIndex(iodw, iode, sst_anom_sd=True, xsd=0.5, opposite_signs_criteria=True):

    import numpy as np
    from itertools import groupby
    import pandas as pd

    limitsize = len(iodw) - 2

    # dipole mode index
    dmi = iodw - iode

    # criteria
    western_sign = np.sign(iodw)
    eastern_sign = np.sign(iode)
    opposite_signs = western_sign != eastern_sign



    sd = np.std(dmi) * xsd
    print(str(sd))
    sdw = np.std(iodw.values) * xsd
    sde = np.std(iode.values) * xsd

    valid_criteria = dmi.__abs__() > sd

    results = []
    if opposite_signs_criteria:
        for k, g in groupby(enumerate(opposite_signs.values), key=lambda x: x[1]):
            if k:
                g = list(g)
                results.append([g[0][0], len(g)])
    else:
        for k, g in groupby(enumerate(valid_criteria.values), key=lambda x: x[1]):
            if k:
                g = list(g)
                results.append([g[0][0], len(g)])


    iods = pd.DataFrame(columns=['DMI', 'Años', 'Mes'], dtype=float)
    dmi_raw = []
    for m in range(0, len(results)):
        # True values
        len_true = results[m][1]

        # True values for at least 3 consecutive seasons
        if len_true >= 3:

            for l in range(0, len_true):

                if l < (len_true - 2):

                    main_month_num = results[m][0] + 1 + l
                    if main_month_num != limitsize:
                        main_month_name = dmi[main_month_num]['time.month'].values  # "name" 1 2 3 4 5

                        main_season = dmi[main_month_num]
                        b_season = dmi[main_month_num - 1]
                        a_season = dmi[main_month_num + 1]

                        # abs(dmi) > sd....(0.5*sd)
                        aux = (abs(main_season.values) > sd) & \
                              (abs(b_season) > sd) & \
                              (abs(a_season) > sd)

                        if sst_anom_sd:
                            if aux:
                                sstw_main = iodw[main_month_num]
                                sstw_b = iodw[main_month_num - 1]
                                sstw_a = iodw[main_month_num + 1]
                                #
                                aux2 = (abs(sstw_main) > sdw) & \
                                       (abs(sstw_b) > sdw) & \
                                       (abs(sstw_a) > sdw)
                                #
                                sste_main = iode[main_month_num]
                                sste_b = iode[main_month_num - 1]
                                sste_a = iode[main_month_num + 1]

                                aux3 = (abs(sste_main) > sde) & \
                                       (abs(sste_b) > sde) & \
                                       (abs(sste_a) > sde)

                                if aux3 & aux2:
                                    iods = iods.append({'DMI': np.around(dmi[main_month_num].values, 2),
                                                        'Años': np.around(dmi[main_month_num]['time.year'].values),
                                                        'Mes': np.around(dmi[main_month_num]['time.month'].values)},
                                                       ignore_index=True)

                                    a = results[m][0]
                                    dmi_raw.append([np.arange(a, a + results[m][1]),
                                                    dmi[np.arange(a, a + results[m][1])].values])


                        else:
                            if aux:
                                iods = iods.append({'DMI': np.around(dmi[main_month_num].values, 2),
                                                    'Años': np.around(dmi[main_month_num]['time.year'].values),
                                                    'Mes': np.around(dmi[main_month_num]['time.month'].values)},
                                                   ignore_index=True)

    return iods, dmi_raw

def DMI(per = 0, filter_bwa = True, filter_harmonic = True,
        filter_all_harmonic=True, harmonics = [],
        start_per=1920, end_per=2020,
        sst_anom_sd=True, opposite_signs_criteria=True):

    western_io = slice(50, 70) # definicion tradicional

    start_per = str(start_per)
    end_per = str(end_per)

    if per == 2:
        movinganomaly = True
        start_year = '1906'
        end_year = '2020'
        change_baseline = False
        start_year2 = '1920'
        end_year2 = '2020_30r5'
        print('30r5')
    else:
        movinganomaly = False
        start_year = start_per
        end_year = end_per
        change_baseline = False
        start_year2 = '1920'
        end_year2 = end_per
        print('All')

    ##################################### DATA #####################################
    # ERSSTv5
    sst = xr.open_dataset("/pikachu/datos/luciano.andrian/verif_2019_2023/sst.mnmean.nc")
    dataname = 'ERSST'
    ##################################### Pre-processing #####################################
    iodw = sst.sel(lat=slice(10.0, -10.0), lon=western_io,
                       time=slice(start_year + '-01-01', end_year + '-12-31'))
    iodw = iodw.sst.mean(['lon', 'lat'], skipna=True)
    iodw2 = iodw
    if per == 2:
        iodw2 = iodw2[168:]
    # -----------------------------------------------------------------------------------#
    iode = sst.sel(lat=slice(0, -10.0), lon=slice(90, 110),
                   time=slice(start_year + '-01-01', end_year + '-12-31'))
    iode = iode.sst.mean(['lon', 'lat'], skipna=True)
    # -----------------------------------------------------------------------------------#
    bwa = sst.sel(lat=slice(20.0, -20.0), lon=slice(40, 110),
                  time=slice(start_year + '-01-01', end_year + '-12-31'))
    bwa = bwa.sst.mean(['lon', 'lat'], skipna=True)
    # ----------------------------------------------------------------------------------#

    if movinganomaly:
        iodw = MovingBasePeriodAnomaly(iodw)
        iode = MovingBasePeriodAnomaly(iode)
        bwa = MovingBasePeriodAnomaly(bwa)
    else:
        # change baseline
        if change_baseline:
            iodw = iodw.groupby('time.month') - \
                   iodw.sel(time=slice(start_year2 + '-01-01', end_year2 + '-12-31')).groupby('time.month').mean(
                       'time')

            iode = iode.groupby('time.month') - \
                   iode.sel(time=slice(start_year2 + '-01-01', end_year2 + '-12-31')).groupby('time.month').mean(
                       'time')

            bwa = bwa.groupby('time.month') - \
                  bwa.sel(time=slice(start_year2 + '-01-01', end_year2 + '-12-31')).groupby('time.month').mean(
                      'time')
            print('baseline: ' + str(start_year2) + ' - ' + str(end_year2))
        else:
            print('baseline: All period')
            iodw = iodw.groupby('time.month') - iodw.groupby('time.month').mean('time', skipna=True)
            iode = iode.groupby('time.month') - iode.groupby('time.month').mean('time', skipna=True)
            bwa = bwa.groupby('time.month') - bwa.groupby('time.month').mean('time', skipna=True)

    # ----------------------------------------------------------------------------------#
    # Detrend
    iodw_trend = np.polyfit(range(0, len(iodw)), iodw, deg=1)
    iodw = iodw - (iodw_trend[0] * range(0, len(iodw)) + iodw_trend[1])
    # ----------------------------------------------------------------------------------#
    iode_trend = np.polyfit(range(0, len(iode)), iode, deg=1)
    iode = iode - (iode_trend[0] * range(0, len(iode)) + iode_trend[1])
    # ----------------------------------------------------------------------------------#
    bwa_trend = np.polyfit(range(0, len(bwa)), bwa, deg=1)
    bwa = bwa - (bwa_trend[0] * range(0, len(bwa)) + bwa_trend[1])
    # ----------------------------------------------------------------------------------#

    # 3-Month running mean
    iodw_filtered = np.convolve(iodw, np.ones((3,)) / 3, mode='same')
    iode_filtered = np.convolve(iode, np.ones((3,)) / 3, mode='same')
    bwa_filtered = np.convolve(bwa, np.ones((3,)) / 3, mode='same')

    # Common preprocessing, for DMIs other than SY2003a
    iode_3rm = iode_filtered
    iodw_3rm = iodw_filtered

    #################################### follow SY2003a #######################################

    # power spectrum
    # aux = FFT2(iodw_3rm, maxVar=20, maxA=15).sort_values('Variance', ascending=False)
    # aux2 = FFT2(iode_3rm, maxVar=20, maxA=15).sort_values('Variance', ascending=False)

    # filtering harmonic
    if filter_harmonic:
        if filter_all_harmonic:
            for harmonic in range(15):
                iodw_filtered = WaveFilter(iodw_filtered, harmonic)
                iode_filtered = WaveFilter(iode_filtered, harmonic)
            else:
                for harmonic in harmonics:
                    iodw_filtered = WaveFilter(iodw_filtered, harmonic)
                    iode_filtered = WaveFilter(iode_filtered, harmonic)

    ## max corr. lag +3 in IODW
    ## max corr. lag +6 in IODE

    # ----------------------------------------------------------------------------------#
    # ENSO influence
    # pre processing same as before
    if filter_bwa:
        ninio3 = sst.sel(lat=slice(5.0, -5.0), lon=slice(210, 270),
                         time=slice(start_year + '-01-01', end_year + '-12-31'))
        ninio3 = ninio3.sst.mean(['lon', 'lat'], skipna=True)

        if movinganomaly:
            ninio3 = MovingBasePeriodAnomaly(ninio3)
        else:
            if change_baseline:
                ninio3 = ninio3.groupby('time.month') - \
                         ninio3.sel(time=slice(start_year2 + '-01-01', end_year2 + '-12-31')).groupby(
                             'time.month').mean(
                             'time')

            else:

                ninio3 = ninio3.groupby('time.month') - ninio3.groupby('time.month').mean('time', skipna=True)

            trend = np.polyfit(range(0, len(ninio3)), ninio3, deg=1)
            ninio3 = ninio3 - (trend[0] * range(0, len(ninio3)) +trend[1])

        # 3-month running mean
        ninio3_filtered = np.convolve(ninio3, np.ones((3,)) / 3, mode='same')

        # ----------------------------------------------------------------------------------#
        # removing BWA effect
        # lag de maxima corr coincide para las dos bases de datos.
        lag = 3
        x = pd.DataFrame({'bwa': bwa_filtered[lag:], 'ninio3': ninio3_filtered[:-lag]})
        result = sm.ols(formula='bwa~ninio3', data=x).fit()
        recta = result.params[1] * ninio3_filtered + result.params[0]
        iodw_f = iodw_filtered - recta

        lag = 6
        x = pd.DataFrame({'bwa': bwa_filtered[lag:], 'ninio3': ninio3_filtered[:-lag]})
        result = sm.ols(formula='bwa~ninio3', data=x).fit()
        recta = result.params[1] * ninio3_filtered + result.params[0]
        iode_f = iode_filtered - recta
        print('BWA filtrado')
    else:
        iodw_f = iodw_filtered
        iode_f = iode_filtered
        print('BWA no filtrado')
    # ----------------------------------------------------------------------------------#

    # END processing
    if movinganomaly:
        iodw_3rm = xr.DataArray(iodw_3rm, coords=[iodw.time.values], dims=['time'])
        iode_3rm = xr.DataArray(iode_3rm, coords=[iodw.time.values], dims=['time'])

        iodw_f = xr.DataArray(iodw_f, coords=[iodw.time.values], dims=['time'])
        iode_f = xr.DataArray(iode_f, coords=[iodw.time.values], dims=['time'])
        start_year = '1920'
    else:
        iodw_3rm = xr.DataArray(iodw_3rm, coords=[iodw2.time.values], dims=['time'])
        iode_3rm = xr.DataArray(iode_3rm, coords=[iodw2.time.values], dims=['time'])

        iodw_f = xr.DataArray(iodw_f, coords=[iodw2.time.values], dims=['time'])
        iode_f = xr.DataArray(iode_f, coords=[iodw2.time.values], dims=['time'])

    ####################################### compute DMI #######################################

    dmi_sy_full, dmi_raw = DMIndex(iodw_f, iode_f,
                                   sst_anom_sd=sst_anom_sd,
                                   opposite_signs_criteria=opposite_signs_criteria)

    return dmi_sy_full, dmi_raw, (iodw_f-iode_f)#, iodw_f - iode_f, iodw_f, iode_f

def DMI2(end_per=1920, start_per=2020, filter_harmonic=True, filter_bwa=False,
         sst_anom_sd=True, opposite_signs_criteria=True):

    # argumentos fijos ------------------------------------------------------------------------------------------------#
    movinganomaly = False
    change_baseline = False
    start_year2 = '6666'
    end_year2 = end_per
    #------------------------------------------------------------------------------------------------------------------#
    western_io = slice(50, 70)  # definicion tradicional
    start_per = str(start_per)
    end_per = str(end_per)

    start_year = start_per
    end_year = end_per
    ####################################################################################################################
    # DATA - ERSSTv5 --------------------------------------------------------------------------------------------------#
    sst = xr.open_dataset("/pikachu/datos/luciano.andrian/verif_2019_2023/sst.mnmean.nc")

    # Pre-processing --------------------------------------------------------------------------------------------------#
    iodw = sst.sel(lat=slice(10.0, -10.0), lon=western_io,
                       time=slice(start_year + '-01-01', end_year + '-12-31'))
    iodw = iodw.sst.mean(['lon', 'lat'], skipna=True)
    # -----------------------------------------------------------------------------------------------------------------#
    iode = sst.sel(lat=slice(0, -10.0), lon=slice(90, 110),
                   time=slice(start_year + '-01-01', end_year + '-12-31'))
    iode = iode.sst.mean(['lon', 'lat'], skipna=True)
    # -----------------------------------------------------------------------------------------------------------------#

    if movinganomaly:
        iodw = MovingBasePeriodAnomaly(iodw)
        iode = MovingBasePeriodAnomaly(iode)
    else:
        if change_baseline:
            iodw = iodw.groupby('time.month') - \
                   iodw.sel(time=slice(start_year2 + '-01-01', end_year2 + '-12-31')).groupby('time.month').mean(
                       'time')

            iode = iode.groupby('time.month') - \
                   iode.sel(time=slice(start_year2 + '-01-01', end_year2 + '-12-31')).groupby('time.month').mean(
                       'time')

            print('baseline: ' + str(start_year2) + ' - ' + str(end_year2))
        else:
            print('baseline: All period')
            iodw = iodw.groupby('time.month') - iodw.groupby('time.month').mean('time', skipna=True)
            iode = iode.groupby('time.month') - iode.groupby('time.month').mean('time', skipna=True)

    # Detrend ---------------------------------------------------------------------------------------------------------#
    iodw_trend = np.polyfit(range(0, len(iodw)), iodw, deg=1)
    iodw = iodw - (iodw_trend[0] * range(0, len(iodw)) + iodw_trend[1])
    #------------------------------------------------------------------------------------------------------------------#
    iode_trend = np.polyfit(range(0, len(iode)), iode, deg=1)
    iode = iode - (iode_trend[0] * range(0, len(iode)) + iode_trend[1])
    #------------------------------------------------------------------------------------------------------------------#
    # 3-Month running mean --------------------------------------------------------------------------------------------#
    iodw_filtered = np.convolve(iodw, np.ones((3,)) / 3, mode='same')
    iode_filtered = np.convolve(iode, np.ones((3,)) / 3, mode='same')

    # Filtering Harmonic ----------------------------------------------------------------------------------------------#
    if filter_harmonic:
        for harmonic in range(15):
            iodw_filtered = WaveFilter(iodw_filtered, harmonic)
            iode_filtered = WaveFilter(iode_filtered, harmonic)

    # Filter BWA #######################################################################################################
    if filter_bwa:
        bwa = sst.sel(lat=slice(20.0, -20.0), lon=slice(40, 110),
                      time=slice(start_year + '-01-01', end_year + '-12-31'))
        bwa = bwa.sst.mean(['lon', 'lat'], skipna=True)

        if movinganomaly:
            bwa = MovingBasePeriodAnomaly(bwa)
        else:
            bwa = bwa.groupby('time.month') - \
                  bwa.sel(time=slice(start_year2 + '-01-01', end_year2 + '-12-31')). \
                      groupby('time.month').mean('time')

        # Detrend -----------------------------------------̣̣------------------------------------------------------------#
        bwa_trend = np.polyfit(range(0, len(bwa)), bwa, deg=1)
        bwa = bwa - (bwa_trend[0] * range(0, len(bwa)) + bwa_trend[1])
        bwa_filtered = np.convolve(bwa, np.ones((3,)) / 3, mode='same')

        ninio3 = sst.sel(lat=slice(5.0, -5.0), lon=slice(210, 270),
                         time=slice(start_year + '-01-01', end_year + '-12-31'))
        ninio3 = ninio3.sst.mean(['lon', 'lat'], skipna=True)

        if movinganomaly:
            ninio3 = MovingBasePeriodAnomaly(ninio3)
        else:
            if change_baseline:
                ninio3 = ninio3.groupby('time.month') - \
                         ninio3.sel(time=slice(start_year2 + '-01-01', end_year2 + '-12-31'))\
                             .groupby('time.month').mean('time')
            else:
                ninio3 = ninio3.groupby('time.month') - ninio3.groupby('time.month').mean('time', skipna=True)

            trend = np.polyfit(range(0, len(ninio3)), ninio3, deg=1)
            ninio3 = ninio3 - (trend[0] * range(0, len(ninio3)) + trend[1])

        # 3-month running mean
        ninio3_filtered = np.convolve(ninio3, np.ones((3,)) / 3, mode='same')

        # -------------------------------------------------------------------------------------------------------------#
        # removing BWA effect
        # lag de maxima corr coincide para las dos bases de datos.
        lag = 3
        x = pd.DataFrame({'bwa': bwa_filtered[lag:], 'ninio3': ninio3_filtered[:-lag]})
        result = sm.ols(formula='bwa~ninio3', data=x).fit()
        recta = result.params[1] * ninio3_filtered + result.params[0]
        iodw_f = iodw_filtered - recta

        lag = 6
        x = pd.DataFrame({'bwa': bwa_filtered[lag:], 'ninio3': ninio3_filtered[:-lag]})
        result = sm.ols(formula='bwa~ninio3', data=x).fit()
        recta = result.params[1] * ninio3_filtered + result.params[0]
        iode_f = iode_filtered - recta
        print('BWA filtrado')
    else:
        iodw_f = iodw_filtered
        iode_f = iode_filtered

    ####################################################################################################################
    # END processing --------------------------------------------------------------------------------------------------#
    iodw_f = xr.DataArray(iodw_f, coords=[iodw.time.values], dims=['time'])
    iode_f = xr.DataArray(iode_f, coords=[iodw.time.values], dims=['time'])

    # Compute DMI ######################################################################################################
    dmi_sy_full, dmi_raw = DMIndex(iodw_f, iode_f,
                                   sst_anom_sd=sst_anom_sd,
                                   opposite_signs_criteria=opposite_signs_criteria)
    return dmi_sy_full, dmi_raw, (iodw_f - iode_f)
    ####################################################################################################################



def WaveFilter(serie, harmonic):

    import numpy as np

    sum = 0
    sam = 0
    N = np.size(serie)

    sum = 0
    sam = 0

    for j in range(N):
        sum = sum + serie[j] * np.sin(harmonic * 2 * np.pi * j / N)
        sam = sam + serie[j] * np.cos(harmonic * 2 * np.pi * j / N)

    A = 2*sum/N
    B = 2*sam/N

    xs = np.zeros(N)

    for j in range(N):
        xs[j] = A * np.sin(2 * np.pi * harmonic * j / N) + B * np.cos(2 * np.pi * harmonic * j / N)

    fil = serie - xs
    return(fil)
####################################################################################################################
def CreateDirectory(out_dir, *args):
    for arg in args:
        if arg is not None:
            if not os.path.exists(os.path.join(out_dir, str(arg))):
                os.mkdir(os.path.join(out_dir, str(arg)))

def DirAndFile(out_dir, dir_results, common_name, names):
    file_name = f"{'_'.join(names)}_{common_name}.jpg"
    path = os.path.join(out_dir, dir_results, file_name)
    return path

def OpenRegiones(name):
    out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'
    try:
        regiones = pd.read_csv(out_dir + name, sep=',', header=0)
    except:
        return print('Error al abrir el archivo')

    nrwos = regiones.shape[0]

    titulos = []
    lon_regiones = []
    lat_regiones = []
    for n in range(0, nrwos):
        titulos.append(regiones['region'][n])

        aux_lon = [regiones['loni'][n], regiones['lonf'][n]]
        aux_lat = [regiones['lati'][n], regiones['latf'][n]]

        lon_regiones.append(aux_lon)
        lat_regiones.append(aux_lat)

    return titulos, lon_regiones, lat_regiones


def RMean3(data):
    aux = data.rolling(time=3, center=True).mean()
    return aux[:-1]

def Correlaciones(correlaciones_cmap, correlaciones_chirps,
                  aux_lnscmap, aux_lnschirps, lead,
                  aux_lnn34, aux_lndmi, aux_lnsam, aux_lnssam, aux_lnasam,
                  index, t):
    c_cmap = []
    c2_cmap = []
    c_chirps = []
    c2_chirps = []
    if correlaciones_cmap:
        try:
            c_cmap = {'Region': [t],
                 index.upper() + '_CM-DMI': [CorrSP(aux_lnscmap.values,
                                                    aux_lndmi.values)],
                 index.upper() + '_CM-N34': [CorrSP(aux_lnscmap.values,
                                                    aux_lnn34)],
                 index.upper() + '_CM-SAM': [CorrSP(aux_lnscmap.values,
                                                    aux_lnsam.values)],
                 index.upper() + '_CM-S_SAM': [CorrSP(aux_lnscmap.values,
                                                      aux_lnssam.values)],
                 index.upper() + '_CM-A_SAM': [CorrSP(aux_lnscmap.values,
                                                      aux_lnasam.values)]
                 }

            c2_cmap = {'Region': [t],
                  index.upper() + '_CM-DMI': [
                      CorrSP(aux_lnscmap.values, aux_lndmi.values, True)[1]],
                  index.upper() + '_CM-N34': [CorrSP(aux_lnscmap.values,
                                                     aux_lnn34, True)[1]],
                  index.upper() + '_CM-SAM': [
                      CorrSP(aux_lnscmap.values, aux_lnsam.values, True)[
                          1]],
                  index.upper() + '_CM-S_SAM': [CorrSP(aux_lnscmap.values,
                                                       aux_lnssam.values,
                                                       True)[1]],
                  index.upper() + '_CM-A_SAM': [CorrSP(aux_lnscmap.values,
                                                       aux_lnasam.values,
                                                       True)[1]]
                  }

        except:
            print('CMAP')
            print('NaN values in ' + t + ' with Lead' + str(lead))

    if correlaciones_chirps:
        try:
            c_chirps = {'Region': [t],
                 index.upper() + '_CH-DMI': [
                     CorrSP(aux_lnschirps.values,
                            aux_lndmi.values)],
                 index.upper() + '-CH-N34': [
                     CorrSP(aux_lnschirps.values, aux_lnn34)],
                 index.upper() + '_CH-SAM': [
                     CorrSP(aux_lnschirps.values,
                            aux_lnsam.values)],
                 index.upper() + '-CHIRPS-S_SAM': [
                     CorrSP(aux_lnschirps.values,
                            aux_lnssam.values)],
                 index.upper() + '_CH-A_SAM': [
                     CorrSP(aux_lnschirps.values,
                            aux_lnasam.values)]
                 }

            c2_chirps = {'Region': [t],
                  index.upper() + '_CH-DMI': [
                      CorrSP(aux_lnschirps.values, aux_lndmi.values,
                             True)[1]],
                  index.upper() + '-CH-N34': [
                      CorrSP(aux_lnschirps.values, aux_lnn34.values,
                             True)[1]],
                  index.upper() + '_CH-SAM': [
                      CorrSP(aux_lnschirps.values, aux_lnsam.values,
                             True)[1]],
                  index.upper() + '-CHIRPS-S_SAM': [
                      CorrSP(aux_lnschirps.values,
                             aux_lnssam.values, True)[1]],
                  index.upper() + '_CH-A_SAM': [
                      CorrSP(aux_lnschirps.values,
                             aux_lnasam.values, True)[1]]
                  }

        except:
            print('CHIRPS')
            print('NaN values in ' + t + ' with Lead' + str(lead))

    return c_cmap, c2_cmap, c_chirps, c2_chirps

#------------------------------------------------------------------------------#
def Entropy(data_prono, lead=0):

    """
    Calculo del RPS forecast
    :param data:
    :return:
    """
    first = True
    for d in data_prono.time.values:

        if lead != 0:
            d0 = LeadMonth(d, lead)
        else:
            d0 = d

        try:
            auxf = data_prono.sel(time=d0, target=d)
            aux_entropy = -1* (auxf.prob_below*np.log(auxf.prob_below) +
                               auxf.prob_norm*np.log(auxf.prob_norm) +
                               auxf.prob_above * np.log(auxf.prob_above))

            if first:
                first = False
                entropy = aux_entropy
            else:
                entropy = xr.concat([entropy, aux_entropy], dim='time')

        except:
            print('Skip ' + np.datetime_as_string(d0, unit='M') + ' con' +
                  ' target ' + np.datetime_as_string(d, unit='M'))

    # para poder usar la dim time igual que las otras
    if lead != 0:
        entropy = entropy.drop('time')
        entropy = entropy.rename({'target':'time'})

    return entropy
#------------------------------------------------------------------------------#

def ColorBySeason(date):
    month = pd.Timestamp(date).month
    if 3 <= month <= 5:
        return '#E0E048'
    elif 6 <= month <= 8:
        return '#0049FF'
    elif 9 <= month <= 11:
        return '#3DFF80'
    else:
        return '#FF6003'
#------------------------------------------------------------------------------#