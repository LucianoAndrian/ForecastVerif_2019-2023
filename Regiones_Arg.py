"""
Mapa de regiones Arg
"""
################################################################################
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'
dir_results = 'mapas_index'
data_arg = '/pikachu/datos/luciano.andrian/verif_2019_2023/' \
           'gadm41_ARG_shp/gadm41_ARG_1.shp'
################################################################################
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.patches as mpatches
from Funciones import CreateDirectory
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
warnings.filterwarnings("ignore")
from Tool_Poligonos import get_vertices
################################################################################
save = True
if save:
    dpi = 300
    CreateDirectory(out_dir, dir_results)
else:
    dpi = 100
################################################################################
regiones_arg = 'hum_sur2', 'hum_norte2', 'centro2', 'patagonia_oeste', \
               'patagonia', 'noa'#, 'centro_norte'

colores = ['firebrick', 'forestgreen', 'navy', 'dodgerblue', 'yellow',
           'orange', 'pink', 'purple', 'lime']
################################################################################
# Plot regiones SA ------------------------------------------------------------#
print('plot regiones')
fig = plt.figure(figsize=(3, 4), dpi=dpi)
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
crs_latlon = ccrs.PlateCarree()
ax.set_extent([280, 310, -60, -20], crs_latlon)

for r, c in zip(regiones_arg, colores):
    vertices = get_vertices(r)
    x, y = zip(*vertices)
    ax.fill(y, x, color=c, edgecolor='k', linewidth=0.5, alpha=0.5,
            transform=crs_latlon)

adm1_shapes = list(shpreader.Reader(data_arg).geometries())
ax.add_geometries(adm1_shapes, ccrs.PlateCarree(),
                  edgecolor='black', facecolor='none', alpha=1)

ax.add_feature(cartopy.feature.LAND, facecolor='white')
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS)
ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
ax.set_xticks(np.arange(280, 310, 5), crs=crs_latlon)
ax.set_yticks(np.arange(-60, -20, 5), crs=crs_latlon)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(labelsize=7)

plt.tight_layout()
if save:
    plt.savefig(out_dir + 'mapa_regiones_arg.jpg', dpi=dpi)
else:
    plt.show()
################################################################################
print('#######################################################################')
print('done')
print('out_dir = ' + out_dir + dir_results )
print('#######################################################################')
################################################################################

