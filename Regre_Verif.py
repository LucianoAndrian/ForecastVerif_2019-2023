"""
Verificación de pronos con modelo lineal
"""
################################################################################
save = False
################################################################################
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/regreverif/'
chirps_data = '/pikachu/datos/luciano.andrian/verif_2019_2023/chirps/'
################################################################################
# import
import numpy as np
import xarray as xr
from Funciones import ChangeLons, MakeMask, Detrend, Weights, Nino34CPC
################################################################################
# Regre
def LinearReg(xrda, dim, deg=1):
    # liner reg along a single dimension
    aux = xrda.polyfit(dim=dim, deg=deg, skipna=True)
    return aux

def SpatialCorr(campo1, campo2):
    stacked_campo1 = campo1.stack(dim_0=('lon', 'lat'))
    stacked_campo2 = campo2.stack(dim_0=('lon', 'lat'))

    # mascara para identificar los NaN
    mask = np.isfinite(stacked_campo1) & np.isfinite(stacked_campo2)
    campo1_valid = stacked_campo1[mask]
    campo2_valid = stacked_campo2[mask]

    correlacion = np.corrcoef(campo1_valid, campo2_valid)[0, 1]

    return correlacion
######
sst_aux = xr.open_dataset("/pikachu/datos/luciano.andrian/verif_2019_2023/"
                      "sst.mnmean.nc")
sst_aux = sst_aux.sel(time=slice('1920-01-01', '2020-12-01'))
n34_or = Nino34CPC(sst_aux, start=1990, end=2020)[0]

print('CHIRPS ----------------------------------------------------------------')
data = xr.open_dataset(
    chirps_data + 'chirps_1990_2020_mmean.nc').__mul__(365/12) #dia
data = ChangeLons(data, 'longitude')
data = data.rename({'latitude':'lat'})
data = data.rename({'precip':'var'})
data = data.sel(lon=slice(270,330), lat=slice(-60, 20))
data = data.rolling(time=3, center=True).mean()

trimestres = ['JFM', 'FMA', 'MAM', 'AMJ', 'MJJ', 'JJA', 'JAS',
              'ASO', 'SON', 'OND', 'NDJ', 'DJF']
first = True
for mm, s_name in zip(np.arange(1,13), trimestres):
    aux = data.sel(time=data.time.dt.month.isin(mm))
    mmean = data.mean('time') # para sacar anoms de cada trimestre en la verif
    aux = Detrend(aux, 'time')

    # Cambio la dim "time" por los valores del N34 para luego usar la funcion
    # mas facil
    aux['time'] = n34_or.sel(time=n34_or.time.dt.month.isin(mm))
    # Regresion
    aux_regre = LinearReg(aux, 'time')

    if first:
        first = False
        m_means = mmean
        regre_patterns = aux_regre.var_polyfit_coefficients[0] + \
                         aux_regre.var_polyfit_coefficients[1]
    else:
        aux_m_eans = mmean
        m_means = xr.concat([m_means, aux_m_eans], dim='mm')

        aux_regre_pattern = aux_regre.var_polyfit_coefficients[0] + \
                            aux_regre.var_polyfit_coefficients[1]

        regre_patterns = xr.concat([regre_patterns, aux_regre_pattern],
                                   dim='pat')

del data

# Periodo de verificacion
n34_verif = xr.open_dataset('/pikachu/datos/luciano.andrian/verif_2019_2023/'
                            'salidas/oni.nc')
n34_verif = xr.Dataset(
    {
        "oni": ("time", n34_verif["oni"].values)
    },
    coords={
        "time": n34_verif["time"].values
    }
)
# 2019-2023-2024
data = xr.open_dataset(
    chirps_data + 'chirps_2019_2023_mmean.nc').__mul__(365/12) #dia
data = data.sel(time=slice('2019-01-01', None))
data = ChangeLons(data, 'longitude')
data = data.rename({'latitude':'lat'})
data = data.rename({'precip':'var'})
data = data.sel(lon=slice(270,330), lat=slice(-60, 20))

aux = data.rolling(time=3, center=True).mean()
first = True
for t in aux.time.values:
    mm = t.astype('datetime64[M]').astype(int) % 12
    aux_data_anom = data.sel(time=t) - m_means.sel(mm=mm)

    if first:
        first = False
        data_anom = aux_data_anom
    else:
        data_anom = xr.concat([data_anom, aux_data_anom], dim='time')

#########
rs=[]
for t in data_anom.time.values:
    mm = t.astype('datetime64[M]').astype(int) % 12

    # n34 del trimestre a verificar:
    n34_mm = n34_verif.sel(time=t)

    # patron de regresion de ese trimestre:
    regre_pat = regre_patterns.sel(pat=mm)*np.abs(n34_mm['oni'].values)
    regre_pat = regre_pat #- regre_pat.mean()

    # trimestre a verificar
    aux_data_anom = data_anom.sel(time=t)['var']
    aux_data_anom = aux_data_anom #- aux_data_anom.mean()

    aux_rs = SpatialCorr(regre_pat, aux_data_anom)
    rs.append(aux_rs)

    # fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    # contour1 = axes[0].contourf(aux_data_anom,
    #                             levels=np.linspace(-300, 300, 13),
    #                             cmap='RdBu', extend='both')
    # fig.colorbar(contour1, ax=axes[0])
    # axes[0].set_title('Data Anom')
    #
    # contour2 = axes[1].contourf(regre_pat,
    #                             levels=np.linspace(-150, 150, 13),
    #                             cmap='RdBu', extend='both')
    #
    # fig.colorbar(contour2, ax=axes[1])
    # axes[1].set_title('Regression Patterns')
    #
    # fig.suptitle(f"{t.astype('datetime64[M]')}: ONI: {n34_mm['oni'].values},"
    #              f" r = {round(aux_rs,3)}",
    #              fontsize=16)
    #
    # plt.tight_layout()
    # plt.show()

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
fig = plt.figure(figsize=(15, 7), dpi=100)
ax = fig.add_subplot(111)
ax2 = ax.twinx()
ax.xaxis.set_major_locator(
    mdates.AutoDateLocator(minticks=20, maxticks=26))
x = data_anom.time.values
ax.plot(x, rs,
        color='k',
        label='blue', linewidth=2)

ax2.plot(x, n34_verif['oni'], label='n34', color='red',
                 linewidth=2)
plt.show()

#
# plt.contourf(data_anom.sel(time=t)['var'],
#              levels=np.linspace(-300,300,13),
#              cmap='RdBu');plt.colorbar();plt.show()
#
# plt.contourf(regre_patterns.sel(pat=mm)*n34_mm['oni'].values,
#              levels=np.linspace(-300,300,13),
#              cmap='RdBu');plt.colorbar();plt.show()
#
# SpatialCorr(regre_patterns.sel(pat=mm)*n34_mm['oni'].values, data_anom.sel(time=t)['var'])