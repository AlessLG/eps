'''
pandas==2.3.3
numpy==2.2.6

Este código se divide en tres:
1. Extracción de datos
2. La física
3. Generación de gráficos
'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import os

control_calidad = 0.05 # Mientras más bajo, mejor. Si es igual a cero, se ignora
total_galaxias = 0
galaxias_consideradas = 0

# EXTRACCIÓN DE DATOS
ruta = os.getcwd()
archivos_dat = glob.glob(os.path.join(ruta, '../Rotmod_LTG/*.dat'))

def calculo_discrepancia_masa(archivo):
    # Obtiene los datos de un archivo .dat determinado
    datos = pd.read_table(archivo, skiprows=3, sep=r'\s+', names=['Rad', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk',	'SBbul'])

    # Guarda los datos en variables 
    radio, vobs, errv, vgas, vdisk, vbul = datos['Rad'], datos['Vobs'], datos['errV'], datos['Vgas'], datos['Vdisk'], datos['Vbul']
    # print(vobs, vgas, vdisk, vbul)

    # FÍSICA
    # kpc -> m
    radio = radio * 3.086e19
    # v: km/s -> m/s
    vobs, vgas, vdisk, vbul, errv = 1e3*vobs, 1e3*vgas, 1e3*vdisk, 1e3*vbul, 1e3*errv

    # Calcula la discrepancia de masa (D)
    #vbar_2 = vbar^2
    vbar_2 = vgas**2 + vdisk**2 + vbul**2
    D = vobs**2/vbar_2
    errD = 2*D*errv/vobs # Asumiendo que los errores en vbar son cero
    a = vobs**2/radio
    g_n = vbar_2/radio

    if control_calidad > 0:
        error = np.average(errv/vobs)
    return D, radio, a, g_n, error, errD


# GENERACIÓN DE GRÁFICOS
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10), dpi=300)

for archivo in archivos_dat:
    D, radio, a, g_n, error, errD = calculo_discrepancia_masa(archivo)
    total_galaxias += 1
    if error <= control_calidad:
        ax1.errorbar(radio, D, errD, fmt='k.', alpha=0.2) 
        ax2.errorbar(a, D, errD, fmt='k.', alpha=0.2)
        ax3.errorbar(g_n, D, errD, fmt='k.', alpha=0.2)



        galaxias_consideradas += 1

for ax in (ax1, ax2, ax3):
    ax.set_ylim(-1.5, 20)
    ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.set_ylabel(r'$\mathcal{D}=(V/V_b)^2$')
    ax.axhline(y=1, xmin=0, xmax=1, color='black', alpha=0.5)

ax1.set_xlabel(r'$r$ [m]')
ax2.set_xlabel(r'$a = V^2/r$ [m s$^{-2}$]')
ax3.set_xlabel(r'$g_N = V_b^2/r$ [m s$^{-2}$]')

fig.suptitle(f'Error' r'$\leq$' f'{control_calidad*100}% - {galaxias_consideradas} galaxias consideradas ({(galaxias_consideradas/total_galaxias)*100:.1f}%)')
fig.tight_layout()
fig.savefig('figuras/D_vs_radio-a-gn.jpg')

print(f'De {total_galaxias} galaxias, {galaxias_consideradas} fueron consideradas ({(galaxias_consideradas/total_galaxias)*100:.1f}%). ')
plt.close()