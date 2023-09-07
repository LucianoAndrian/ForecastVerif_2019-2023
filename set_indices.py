"""
Abre y setea los trimestre de los indices y calcula el endtime usado para nmme
y los datos observados
"""
################################################################################
nmme_pronos = '/pikachu/datos/luciano.andrian/verif_2019_2023/nmme_pronos/'
dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/scatter/'
dir_results = 'index_vs_index'
################################################################################
import xarray as xr
from Funciones import DMI, SameDateAs, RMean3
################################################################################
def compute():
    # ONI Descargado, para no cambiar tdo el codigo n34 = ONI -----------------#
    n34 = xr.open_dataset(dir + 'oni.nc')

    endtime = n34.time.values[-1]
    anio = endtime.astype('datetime64[Y]').astype(int) + 1970

    # SAM ---------------------------------------------------------------------#
    sam = SameDateAs(RMean3(xr.open_dataset(dir + 'sam.nc').mean_estimate), n34)
    asam = SameDateAs(RMean3(xr.open_dataset(dir + 'asam.nc').mean_estimate),
                      n34)
    ssam = SameDateAs(RMean3(xr.open_dataset(dir + 'ssam.nc').mean_estimate),
                      n34)

    # DMI calculado a partir de ERSSTv5 actualizada ---------------------------#
    aux0, aux, dmi = DMI(filter_bwa=False, start_per=1920, end_per=anio)
    dmi = SameDateAs(dmi, n34)

    return n34, dmi, sam, ssam, asam, endtime


if __name__ == "__main__":
    n34, dmi, sam, ssam, asam, endtime = compute()