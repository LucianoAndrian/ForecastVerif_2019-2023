"""
Mapa de regiones
"""
################################################################################
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'
dir_results = 'mapas_index'
################################################################################
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import matplotlib.patches as mpatches
from Funciones import CreateDirectory
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
warnings.filterwarnings("ignore")
################################################################################
save = True
if save:
    dpi = 300
    CreateDirectory(out_dir, dir_results)
else:
    dpi = 100
################################################################################
lon_regiones_sa = [[296, 296 + 10], [285, 285 + 8]]
lat_regiones_sa = [[-39, -39 + 14], [-40, -40 + 10]]
names_regiones_sa = ['SESA', 'Cuyo-Chile']

################################################################################
# Plot regiones SA ------------------------------------------------------------#
print('plot regiones')
fig = plt.figure(figsize=(3, 4), dpi=dpi)
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
crs_latlon = ccrs.PlateCarree()
ax.set_extent([270, 330, -60, 20], crs_latlon)

for r, rname in enumerate(names_regiones_sa):
    w = lon_regiones_sa[r][1] - lon_regiones_sa[r][0]
    h = np.abs(lat_regiones_sa[r][0]) - np.abs(lat_regiones_sa[r][1])
    ax.add_patch(mpatches.Rectangle(xy=[lon_regiones_sa[r][0],
                                        lat_regiones_sa[r][0]],
                                    width = w, height=h,
                                facecolor='None', alpha=1, edgecolor='k',
                                linewidth=2, transform=ccrs.PlateCarree()))

    d = {'region': [rname],
         'loni': [lon_regiones_sa[r][0]],
         'lonf': [lon_regiones_sa[r][1]],
         'lati': [lat_regiones_sa[r][0]],
         'latf': [lat_regiones_sa[r][1]]}

    if r == 0:  # la primera
        df = pd.DataFrame(d)
    else:
        df = df.append(pd.DataFrame(d), ignore_index=True)

ax.add_feature(cartopy.feature.LAND, facecolor='white')
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
    plt.savefig(out_dir + 'mapa_regiones.jpg', dpi=dpi)
else:
    plt.show()

df.to_csv(out_dir + 'regiones_sa.csv', index=False)
# Plot regiones ARG------------------------------------------------------------#
################################################################################
print('#######################################################################')
print('done')
print('out_dir = ' + out_dir + dir_results )
print('#######################################################################')
################################################################################