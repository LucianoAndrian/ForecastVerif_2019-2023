"""
Barras de colores y escalas varias
"""
################################################################################
from matplotlib import colors
import numpy as np
################################################################################
################################################################################
# Colorbars ####################################################################
def get_cbars(VarName):
    cbar_pp = colors.ListedColormap(['#003C30', '#004C42', '#0C7169',
                                     '#79C8BC', '#B4E2DB',
                                     'white',
                                     '#F1DFB3', '#DCBC75', '#995D13',
                                     '#6A3D07', '#543005', ][::-1])
    cbar_pp.set_under('#3F2404')
    cbar_pp.set_over('#00221A')
    cbar_pp.set_bad(color='white')

    if VarName.lower() == 'hgt200':
        return
    elif VarName.lower() == 'pp':
        return cbar_pp
    elif VarName.lower() == 't':
        return

# scales #######################################################################
def get_scales(VarName):
    scale_reg_hgt = [-150, -100, -75, -50, -25, -15, 0, 15, 25, 50, 75, 100, 150]

    scale_reg_pp = [-45,-30, -15, -5,0,5, 15, 30, 45]

    scale_reg_t =  [-.6,-.4,-.2,-.1,-.05,0,0.05,0.1,0.2,0.4,0.6]

    if VarName.lower() == 'hgt200':
        return
    elif VarName.lower() == 'pp':
        return scale_reg_pp
    elif VarName.lower() == 't':
        return

################################################################################
if __name__ == "__main__":
    cbar = get_cbars(VarName)
    scale = get_scales(VarName)
################################################################################