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
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from Funciones import ChangeLons, Detrend, Weights, Nino34CPC
import set_indices
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
# n34 ------------------------------------------------------------------------ #
sst_aux = xr.open_dataset("/pikachu/datos/luciano.andrian/verif_2019_2023/"
                      "sst.mnmean.nc")
sst_aux = sst_aux.sel(time=slice('1920-01-01', '2020-12-01'))
n34_or = Nino34CPC(sst_aux, start=1990, end=2020)[0]

# data ----------------------------------------------------------------------- #
data = xr.open_dataset(
    chirps_data + 'chirps_1990_2020_mmean.nc').__mul__(365/12) #dia
data = ChangeLons(data, 'longitude')
data = data.rename({'latitude':'lat'})
data = data.rename({'precip':'var'})
data = data.sel(lon=slice(270,330), lat=slice(-60, 15))
data = data.rolling(time=3, center=True).mean()
data = Weights(data)

# ---------------------------------------------------------------------------- #
trimestres = ['JFM', 'FMA', 'MAM', 'AMJ', 'MJJ', 'JJA', 'JAS',
              'ASO', 'SON', 'OND', 'NDJ', 'DJF']
first = True
for mm, s_name in zip(np.arange(1,13), trimestres):
    # seleccion de trimestre
    data_trim = data.sel(time=data.time.dt.month.isin(mm))

    # media para sacar anoms de cada trimestre en la verif
    mmean = data_trim.mean('time')
    data_trim = Detrend(data_trim, 'time')
    #data_trim = data_trim/data_trim.std('time')

    # Regresion. time = n34_or(en season) para facilitar el uso de LinearReg
    n34_aux = n34_or.sel(time=n34_or.time.dt.month.isin(mm))
    data_trim['time'] = n34_aux#/n34_aux.std('time')
    aux_regre = LinearReg(data_trim, 'time')

    aux_regre_pat = (aux_regre.var_polyfit_coefficients[0] +
                     aux_regre.var_polyfit_coefficients[1])

    if first:
        first = False
        m_means = mmean
        regre_patterns = aux_regre_pat
    else:
        m_means = xr.concat([m_means, mmean], dim='mm')
        regre_patterns = xr.concat([regre_patterns, aux_regre_pat],
                                   dim='pat')

# ---------------------------------------------------------------------------- #
# Periodo de verificacion ---------------------------------------------------- #
n34_verif = xr.open_dataset('/pikachu/datos/luciano.andrian/verif_2019_2023/'
                            'salidas/oni.nc')
n34_verif = xr.Dataset(
    {"oni": ("time", n34_verif["oni"].values)},
    coords={"time": n34_verif["time"].values})
#n34_verif = n34_verif/n34_verif.std('time')

n34, dmi, sam, ssam, asam, endtime = set_indices.compute()

# 2019-2023-2024 ------------------------------------------------------------- #
data = xr.open_dataset(
    chirps_data + 'chirps_2019_2023_mmean.nc').__mul__(365/12) #dia
data = data.sel(time=slice('2019-01-01', None))
data = ChangeLons(data, 'longitude')
data = data.rename({'latitude':'lat'})
data = data.rename({'precip':'var'})
data = data.sel(lon=slice(270,330), lat=slice(-60, 15))
data = data.rolling(time=3, center=True).mean()
#data = data/data.std('time')

# Anomalia de cada trimestre respecto a la diferencia
first = True
for t in data.time.values:
    mm = t.astype('datetime64[M]').astype(int) % 12
    aux_data_anom = data.sel(time=t) - m_means.sel(mm=mm)
    if first:
        first = False
        data_anom = aux_data_anom
    else:
        data_anom = xr.concat([data_anom, aux_data_anom], dim='time')

# ---------------------------------------------------------------------------- #
rs=[]
for t in data_anom.time.values: # para cada trimestre
    mm = t.astype('datetime64[M]').astype(int) % 12

    # n34 del trimestre a verificar
    n34_mm = n34_verif.sel(time=t)

    # patron de regresion de ese trimestre * n34
    regre_pat = regre_patterns.sel(pat=mm) * np.abs(n34_mm['oni'].values)
    regre_pat = regre_pat - regre_pat.mean()

    # trimestre a verificar
    aux_data_anom = data_anom.sel(time=t)['var']
    aux_data_anom = aux_data_anom - aux_data_anom.mean()

    # Correlacion espacial
    aux_rs = SpatialCorr(regre_pat, aux_data_anom)
    rs.append(aux_rs)

# ---------------------------------------------------------------------------- #
# plots
# ---------------------------------------------------------------------------- #
# Evolucion temporal
plt.rcParams['date.converter'] = 'concise'
fig = plt.figure(figsize=(12, 7), dpi=100)
ax = fig.add_subplot(111)
ax2 = ax.twinx()
ax.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=20, maxticks=26))
x = data_anom.time.values
rline = ax.plot(x, rs, color='k', label='r', linewidth=2.5)
oniline = ax2.plot(x, n34_verif['oni'], label='ONI', color='red', linewidth=2.5)
dmiline = ax2.plot(x, dmi, label='DMI', color='#007F22', linewidth=1)
samline = ax2.plot(x, sam, label='SAM', color='#FF9B00', linewidth=1)
ax.hlines(y=0, xmin=x[0], xmax=x[-1], color='k', linestyle='-')
ax.grid()
lns = rline + oniline + dmiline + samline
ax.legend(lns, ['r', 'ONI', 'DMI', 'SSAM'], loc='lower right')
ax.set_ylabel('PP pattern correlation', fontsize=12)
ax2.set_ylabel('ONI', fontsize=12)
ax.set_ylim(-0.7,0.7)
ax2.set_ylim(-2.5,2.5)
plt.tight_layout()
if save:
    plt.savefig(f"{out_dir}evolR_indices.jpg")
    plt.close()
else:
    plt.show()

# ---------------------------------------------------------------------------- #
# Scatter
fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(111)

formatted_dates = np.array([f"{str(date)[2:4]}-{str(date)[5:7]}" for date in x])
formatted_dates2 = np.array([f"{str(date)[5:7]}" for date in x])

mm_pos = []
mm_year = []
for y in ['19', '20', '21', '22', '23', '24']:
    aux_mm_pos = []
    aux_mm_year = []
    for t_pos, t in enumerate(formatted_dates):
        if y in t:
            aux_mm_pos.append(t_pos)
            aux_mm_year.append(t[3:5])

    mm_pos.append(aux_mm_pos)
    mm_year.append(aux_mm_year)

ax.set_ylim(-0.6,0.7)
ax.set_xlim(-1.5,2.5)
ax.grid()
ax.hlines(y=0, xmin=-2.5, xmax=2.5, color='k', linestyle='--')
ax.vlines(x=0, ymin=-0.8, ymax=0.8, color='k', linestyle='--')

for pos, mes, c, y in zip(mm_pos, mm_year, ['green', 'blue', 'red', 'orange',
                                            '#006D7B', '#FF00AC'],
                       ['2019', '2020', '2021', '2022', '2023', '2024']):

    ax.scatter(n34_verif['oni'].values[pos], rs[pos[0]:pos[-1]+1], s=9,
               color=c, marker='o', alpha=1, label=y)

    for i, m in zip(pos, mes):
        ax.text(n34_verif['oni'].values[i], rs[i], m, fontsize=9,
            alpha=1, color=c,
            verticalalignment='bottom', horizontalalignment='center')

plt.legend()
plt.xlabel('ONI')
plt.ylabel('PP pattern correlation')
plt.tight_layout()
if save:
    plt.savefig(f"{out_dir}scatter_full.jpg")
    plt.close()
else:
    plt.show()

################################################################################
# ---------------------------------------------------------------------------- #
# a modo de ejemplo:
# ---------------------------------------------------------------------------- #
# a veeeer
# max positivo
mmax = np.where(rs==max(rs[1:-1]))[0][0]
#mmax = 15+12
t = data_anom.time.values[mmax]
mm = t.astype('datetime64[M]').astype(int) % 12
n34_mm = n34_verif.sel(time=t)

aux_rs = SpatialCorr(regre_patterns.sel(pat=mm)*np.abs(n34_mm['oni'].values),
                     data_anom.sel(time=t)['var'])
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
contour1 = axes[0].contourf(data_anom.sel(time=t)['var'],
                            levels=np.linspace(-75, 75, 13),
                            cmap='BrBG', extend='both')
fig.colorbar(contour1, ax=axes[0])
axes[0].set_title('Data Anom')

contour2 = axes[1].contourf(regre_patterns.sel(pat=mm)*np.abs(n34_mm['oni'].values),
                            levels=np.linspace(-75, 75, 13),
                            cmap='BrBG', extend='both')

fig.colorbar(contour2, ax=axes[1])
axes[1].set_title('Regression Patterns')

fig.suptitle(f"{t.astype('datetime64[M]')}: ONI: {n34_mm['oni'].values},"
             f" r = {round(aux_rs,3)}",
             fontsize=16)
plt.tight_layout()
if save:
    plt.savefig(f"{out_dir}ej_max.jpg")
    plt.close()
else:
    plt.show()

# min negativo
mmax = np.where(rs==min(rs[1:-1]))[0][0]
#mmax = 15+12
t = data_anom.time.values[mmax]
mm = t.astype('datetime64[M]').astype(int) % 12
n34_mm = n34_verif.sel(time=t)

aux_rs = SpatialCorr(regre_patterns.sel(pat=mm)*np.abs(n34_mm['oni'].values),
                     data_anom.sel(time=t)['var'])
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
contour1 = axes[0].contourf(data_anom.sel(time=t)['var'],
                            levels=np.linspace(-75, 75, 13),
                            cmap='BrBG', extend='both')
fig.colorbar(contour1, ax=axes[0])
axes[0].set_title('Data Anom')

contour2 = axes[1].contourf(regre_patterns.sel(pat=mm)*np.abs(n34_mm['oni'].values),
                            levels=np.linspace(-75, 75, 13),
                            cmap='BrBG', extend='both')

fig.colorbar(contour2, ax=axes[1])
axes[1].set_title('Regression Patterns')

fig.suptitle(f"{t.astype('datetime64[M]')}: ONI: {n34_mm['oni'].values},"
             f" r = {round(aux_rs,3)}",
             fontsize=16)
plt.tight_layout()
if save:
    plt.savefig(f"{out_dir}ej_min.jpg")
    plt.close()
else:
    plt.show()
################################################################################