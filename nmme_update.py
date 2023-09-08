"""
Chequea si el prono del mes actual ya está descargado sino lo descarga
"""
################################################################################
ruta_nmme = '/pikachu/datos/luciano.andrian/verif_2019_2023/nmme_pronos/'
################################################################################
import subprocess
from Funciones import SelectFilesNMME
from datetime import datetime
################################################################################
def update():
    current_month = datetime.now().strftime('%m')
    # con size_check = True, abre todos los archivos que su tamaño sea > 1mb
    # los pronos pesan 5.3mb
    files = SelectFilesNMME(ruta_nmme, 'prate', True)
    last_f = files[-1]
    last_f_month = last_f.split('_')[-2][4:7]

    if current_month != last_f_month:
        subprocess.run(
            ['bash', '/pikachu/datos/luciano.andrian/verif_2019_2023/'
                     'nmme_download/update_nmme.sh'])

################################################################################
if __name__ == "__main__":
    update()
################################################################################