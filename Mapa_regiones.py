"""
Mapa de regiones
"""
################################################################################
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'
dir_results = 'mapas_index'
################################################################################
import numpy as np
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
save = False
if save:
    dpi = 300
    CreateDirectory(out_dir, dir_results)
else:
    dpi = 100
################################################################################
# Plot regiones SA ------------------------------------------------------------#
print('plot regiones')
fig = plt.figure(figsize=(3, 4), dpi=100)
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
crs_latlon = ccrs.PlateCarree()
ax.set_extent([270, 330, -60, 20], crs_latlon)

# SESA idem que en paper NMME
ax.add_patch(mpatches.Rectangle(xy=[296, -39], width=10, height=14,
                                facecolor='None', alpha=1, edgecolor='k',
                                linewidth=2, transform=ccrs.PlateCarree()))

# ~ Cujo - Chile
ax.add_patch(mpatches.Rectangle(xy=[285, -40], width=8, height=10,
                                facecolor='None', alpha=1, edgecolor='k',
                                linewidth=2, transform=ccrs.PlateCarree()))

# # Este Brasil ?
# ax.add_patch(mpatches.Rectangle(xy=[310, -24], width=10, height=12,
#                                 facecolor='None', alpha=1, edgecolor='k',
#                                 linewidth=2, transform=ccrs.PlateCarree()))

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

# Plot regiones ARG------------------------------------------------------------#
################################################################################