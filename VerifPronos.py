nmme_pronos = '/pikachu/datos/luciano.andrian/verif_2019_2023/nmme_pronos/'
cmap_data = '/pikachu/datos/luciano.andrian/observado/ncfiles/data_no_detrend/'
chirps_data = '/pikachu/datos/luciano.andrian/verif_2019_2023/chirps/'
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'

#https://psl.noaa.gov/data/gridded/data.cmap.html
################################################################################
import glob
import numpy as np
import xarray as xr
import pandas as pd
from Funciones import Nino34CPC, DMI
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import matplotlib.patches as mpatches
from matplotlib import colors
################################################################################
dpi=300
save=True
################################################################################
def SelectFilesNMME(dir, variable):
    files =  glob.glob(dir + variable + '_*_prob.adj.sea.nc')
    return sorted(files, key=lambda x: x.split()[0])

def MakeMask(DataArray, dataname='mask'):
    import regionmask
    mask = regionmask.defined_regions.natural_earth_v5_0_0.countries_110.mask(DataArray)
    mask = xr.where(np.isnan(mask), mask, 1)
    mask = mask.to_dataset(name=dataname)
    return mask

def ChangeLons(data, lon_name='lon'):
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

def ABNobs(data_clim, data_verif):
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

def RPSO(data, data_verif):
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

def RPSF(data, data_prono, data_verif):
    first = True
    for d in data_verif.time.values:
        aux = data.sel(time=d)
        aux_bn = aux.B + aux.N
        aux_total = aux.B + aux.N + aux.A

        # d2 = pd.Timestamp(d) + pd.DateOffset(months=1)
        # d2 = d2.to_numpy()
        auxf = data_prono.sel(time=d, target=d)
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

    return RPSf

def Plot(comp, comp_var, levels, save, dpi, title, name_fig, out_dir,
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
        plt.savefig(out_dir + name_fig + '.jpg')
        plt.close()
    else:
        plt.show()

# Open and set #################################################################
# CMAP
data = xr.open_dataset(cmap_data + 'pp_cmap.nc')
lon_cmap = data.lon.values
lat_cmap = data.lat.values
# tiene resolucion 2.5...
data = data.sel(lon=slice(270,330), lat=slice(20, -60))
# promedios trimestrales (VER)
data = data.rolling(time=3, center=True).mean()
#data = data.interp(lon=data_nmme.lon.values, lat=data_nmme.lat.values)
# climatologia trimestral
data_clim = data.sel(time=slice('1990-01-01', '2020-12-01'))
data_clim = data_clim.groupby('time.month').quantile([.33, .66], dim='time')

data_verif = data.sel(time=slice('2019-01-01','2023-04-01'))


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
data_clim_ch = data.groupby('time.month').quantile([.33, .66], dim='time')

data = xr.open_dataset(
    chirps_data + 'chirps_2019_2023_mmean.nc').__mul__(365/12) #dia
data = ChangeLons(data, 'longitude')
data = data.rename({'latitude':'lat'})
data = data.interp(lon=data_nmme.lon,
                   lat=np.arange(-50,50,1))

data = data.sel(lon=slice(270,330), lat=slice(-60, 20))
# promedios trimestrales (VER)
data_verif_ch = data.rolling(time=3, center=True).mean()
################################################################################
# En que categoría está lo observado según la climatología
abn_cmap = ABNobs(data_clim, data_verif)
abn_chirps = ABNobs(data_clim_ch, data_verif_ch)

# RPS ##########################################################################
# RPSo ------------------------------------------------------------------------#
RPSo_cmap = RPSO(abn_cmap, data_verif)
RPSo_chirps = RPSO(abn_chirps, data_verif_ch)

# RPSf ------------------------------------------------------------------------#
RPSf_cmap = RPSF(abn_cmap, data_nmme25, data_verif)
RPSf_chirps = RPSF(abn_chirps, data_nmme, data_verif_ch)

################################################################################
# RPSS
RPSS_cmap = 1 - (RPSf_cmap/RPSo_cmap)
RPSS_chirps = 1 - (RPSf_chirps/RPSo_chirps)

RPSS_cmap_nm = RPSS_cmap
RPSS_chirps_nm = RPSS_chirps

RPSS_cmap = RPSS_cmap*MakeMask(RPSS_cmap)
RPSS_chirps = RPSS_chirps*MakeMask(RPSS_chirps)
################################################################################
# indices
# SST actualizada
# https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/
dmi, aux, dmi_aux = DMI(filter_bwa=False, start_per=1920, end_per=2023)

#pendiente, arreglar Ninio3.4CPC
n34 = [0.7, 0.7, 0.7, 0.7, 0.5, 0.5, 0.3, 0.1, 0.2 ,0.3, 0.5 ,0.5, 0.5,	0.5,
0.4, 0.2, -0.1, -0.3, -0.4, -0.6, -0.9, -1.2, -1.3, -1.2, -1.0, -0.9, -0.8,
-0.7, -0.5, -0.4, -0.4, -0.5, -0.7, -0.8, -1.0, -1.0, -1.0, -0.9, -1.0, -1.1,
-1.0, -0.9, -0.8, -0.9, -1.0, -1.0, -0.9, -0.8, -0.7, -0.4, -0.2]

# incorporar SAM
#

# Plot ------------------------------------------------------------------------#
# Regiones
fig = plt.figure(figsize=(3, 4), dpi=100)
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
crs_latlon = ccrs.PlateCarree()
ax.set_extent([270, 330, -60, 20], crs_latlon)

ax.add_patch(mpatches.Rectangle(xy=[296, -40], width=20, height=20,
                                facecolor='gray',
                                alpha=0.5, edgecolor='black', linewidth=1,
                                transform=ccrs.PlateCarree())
)

ax.add_patch(mpatches.Rectangle(xy=[296, -40], width=20, height=10,
                                facecolor='red',
                                alpha=0.5, edgecolor='red', linewidth=1,
                                transform=ccrs.PlateCarree())
)

ax.add_patch(mpatches.Rectangle(xy=[300, -30], width=20, height=17,
                                facecolor='blue',
                                alpha=0.5, edgecolor='blue', linewidth=1,
                                transform=ccrs.PlateCarree())
)

ax.add_patch(mpatches.Rectangle(xy=[295, -40], width=10, height=15,
                                facecolor='magenta',
                                alpha=0.5, edgecolor='white', linewidth=1,
                                transform=ccrs.PlateCarree())
)

ax.add_patch(mpatches.Rectangle(xy=[290, -40], width=5, height=20,
                                facecolor='dodgerblue',
                                alpha=0.5, edgecolor='dodgerblue', linewidth=1,
                                transform=ccrs.PlateCarree())
)

ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey')
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS)
ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
ax.set_xticks(np.arange(270, 330, 10), crs=crs_latlon)
ax.set_yticks(np.arange(-60, 40, 20), crs=crs_latlon)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(labelsize=7)
plt.title('Regiones', fontsize=10)
plt.tight_layout()
if save:
    plt.savefig(out_dir + 'mapa.jpg', dpi=dpi)
else:
    plt.show()

dates = pd.date_range(start='2018-12-01', end='2023-04-01',
                                  freq='M') + pd.DateOffset(days=1)
dates = dates[:-1]

lon_regiones = [[296,296+20], [296,296+20], [296,300+20],
                [295,295+10], [290,290+5]]
lat_regiones = [[-40,-40+20],[-40,-40+10],[-30,-30+17],
                [-40,-40+15], [-40,-40+20]]
titulos = ['SESA', 'S-SESA', 'N-SESA', 'NEA', 'NOA']

for ln, lt, t in zip(lon_regiones, lat_regiones, titulos):
    aux = RPSS_cmap.sel(lon=slice(ln[0], ln[1]), lat=slice(lt[1], lt[0]))
    aux2 = RPSS_chirps.sel(lon=slice(ln[0], ln[1]), lat=slice(lt[1], lt[0]))
    fig = plt.figure(figsize=(10, 7), dpi=dpi)
    ax = fig.add_subplot(111)
    ax2 = ax.twinx()
    lnrpsscmap = ax.plot(dates, aux.mean(['lon', 'lat']).mask[:-1],
                         color='blue',
                         label='RPSS_CMAP2.5', linewidth=2)
    lnrpsschirps = ax.plot(dates, aux2.mean(['lon', 'lat']).mask[:-1],
                           color='orange',
                           label='RPSS_CHIRPS1', linewidth=2)
    lndmi = ax2.plot(dates, dmi_aux.sel(time=slice('2019-01-01', '2023-03-01')),
                     label='DMI', color='forestgreen')
    lnn34 = ax2.plot(dates, n34, label='ONI', color='firebrick')

    ax.hlines(y=0, xmin=dates[0], xmax=dates[-1], color='gray')
    ax2.hlines(y=0, xmin=dates[0], xmax=dates[-1], color='gray')
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y %b'))
    ax.grid()
    ax.set_ylim((-.4,.5))
    ax2.set_ylim((-1.5, 7))
    ax.set_ylabel('RPSS', fontsize=10)
    ax2.set_ylabel('índices', fontsize=10)
    ax.set_title(t, fontsize=15)

    lns = lnrpsscmap + lnrpsschirps + lndmi + lnn34
    ax.legend(lns, ['RPSS_CMAP2.5', 'RPSS_CHIRPS1', 'DMI', 'ONI'],
              loc='upper left')

    if save:
        plt.savefig(out_dir + t + '.jpg', dpi=dpi)
        plt.close('all')
    else:
        plt.show()
################################################################################
# Mapa RPSS
cbar = colors.ListedColormap(['#6FFE9B','#FEB77E', '#FE675C','#CA3E72',
                                  '#782281','#251255'])
cbar.set_over('#251255')
cbar.set_under('#6FFE9B')
cbar.set_bad(color='white')

rpss = [RPSS_cmap_nm, RPSS_chirps_nm]
dataset = ['CMAP', 'CHIRPS']
seasons = ['MAM', 'JJA', 'SON', 'DJF']
mm = [4, 7, 10, 1]
anio = [2019, 2020, 2021, 2022, 2023]
#------------------------------------------------------------------------------#
# por season y año ------------------------------------------------------------#
for ds, ds_name in zip(rpss, dataset):
    for m, s in zip(mm, seasons):
        for a in anio:
            try:
                if m==10:
                    aux = ds.sel(time=str(a) + '-0' + str(m) + '-01' )
                else:
                    aux = ds.sel(time=str(a) + '-' + str(m) + '-01')
            
                titulo = 'RPSS - ' + s + ' '+ str(a) + ' - ' + ds_name
                name_fig = 'rpss_' + s + '_'+ str(a) + '_' + ds_name

                Plot(aux, aux, [-.1,0,.1,.2,.3,.4,.5], save, dpi, titulo,
                     name_fig, out_dir, 'gray', cbar)
            except:
                pass

# Por periodos ----------------------------------------------------------------#
for ds, ds_name in zip(rpss, dataset):
    for a in anio:
        try:
            aux = ds.sel(time=ds.time.dt.year.isin(a)).mean('time')

            titulo = 'RPSS - ' + str(a) + ' - ' + ds_name
            name_fig = 'rpss_' + str(a) + '_' + ds_name

            Plot(aux, aux, [-.1, 0, .1, .2, .3, .4, .5], save, dpi, titulo,
                 name_fig, out_dir, 'gray', cbar)
        except:
            pass

# All 7/2019-4/2023 -----------------------------------------------------------#
for ds, ds_name in zip(rpss, dataset):
    try:
        aux = ds.mean('time')
        titulo = 'RPSS - ' + '7/2019 - 4/2023' + ' - ' + ds_name
        name_fig = 'rpss_' + ds_name

        Plot(aux, aux, [-.1, 0, .1, .2, .3, .4, .5], save, dpi, titulo,
             name_fig, out_dir, 'gray', cbar)
    except:
        pass
################################################################################