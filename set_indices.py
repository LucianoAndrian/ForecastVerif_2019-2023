"""
Abre y setea los trimestre de los indices y calcula el endtime usado para nmme
y los datos observados
Comprueba si estan actualizados en cada caso
"""
################################################################################
nmme_pronos = '/pikachu/datos/luciano.andrian/verif_2019_2023/nmme_pronos/'
dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/'
out_dir = '/pikachu/datos/luciano.andrian/verif_2019_2023/salidas/scatter/'
dir_results = 'index_vs_index'
################################################################################
import xarray as xr
import numpy as np
from dateutil.relativedelta import relativedelta
from Funciones import DMI, SameDateAs, RMean3
import Set_ONI as ONI
import Set_SAM_indices as SAM
import ERSSTv5
################################################################################
def compute():
    # ONI Descargado, para no cambiar tdo el codigo n34 = ONI -----------------#
    ONI.update()
    n34 = xr.open_dataset(dir + 'oni.nc')

    endtime = n34.time.values[-1]
    anio = endtime.astype('datetime64[Y]').astype(int) + 1970

    # DMI calculado a partir de ERSSTv5 actualizada ---------------------------#
    ERSSTv5.update()
    aux0, aux, dmi = DMI(filter_bwa=False, start_per=1920, end_per=anio)
    dmi = SameDateAs(dmi, n34)


    # SAM ---------------------------------------------------------------------#
    SAM.update()
    try:
        sam = SameDateAs(RMean3(xr.open_dataset(dir + 'sam.nc').mean_estimate), n34)
        asam = SameDateAs(RMean3(xr.open_dataset(dir + 'asam.nc').mean_estimate),
                          n34)
        ssam = SameDateAs(RMean3(xr.open_dataset(dir + 'ssam.nc').mean_estimate),
                          n34)
    except KeyError as e:
        if "not all values found in index 'time'" in str(e):
            print('SAM desactualizado un trimestre')
            print('Seteando todos los indices al mes anterior')

            n34 = n34.isel(index=slice(None, -1))

            try:
                sam = SameDateAs(
                    RMean3(xr.open_dataset(dir + 'sam.nc').mean_estimate), n34)
                asam = SameDateAs(
                    RMean3(xr.open_dataset(dir + 'asam.nc').mean_estimate), n34)
                ssam = SameDateAs(
                    RMean3(xr.open_dataset(dir + 'ssam.nc').mean_estimate), n34)
                dmi = SameDateAs(dmi, n34)
            except:
                print('SAM 2 trimestre desactualizado')
                return
    return n34, dmi, sam, ssam, asam, endtime

################################################################################
if __name__ == "__main__":
    n34, dmi, sam, ssam, asam, endtime = compute()
################################################################################