import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import os

a_0 = 1.2e-10

ruta = os.getcwd()
archivos_dat = glob.glob(os.path.join(ruta, '../Rotmod_LTG/*.dat'))

def calculo_discrepancia_masa(archivo, controlCalidad):
    # Obtiene los datos de un archivo .dat determinado
    datos = pd.read_table(archivo, skiprows=3, sep=r'\s+', names=['Rad', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk',	'SBbul'])

    # Guarda los datos en variables 
    vobs, vgas, vdisk, vbul, errv, radio = datos['Vobs'], datos['Vgas'], datos['Vdisk'], datos['Vbul'], datos['errV'], datos['Rad']
    # print(vobs, vgas, vdisk, vbul)

    # kpc -> m
    radio = radio * 3.086e19
    # v: km/s -> m/s
    vobs, vgas, vdisk, vbul, errv = 1e3*vobs, 1e3*vgas, 1e3*vdisk, 1e3*vbul, 1e3*errv

    # Calcula la discrepancia de masa (D)
    #vbar_2 = vbar^2
    vbar_2 = vgas**2 + vdisk**2 + vbul**2
    D = vobs**2/vbar_2
    a = vobs**2/radio
    x = a/a_0
    g_n = a * (x/np.sqrt(1 + x**2))

    if controlCalidad > 0:
        error = np.average(errv/vobs)
    return D, radio, a, g_n, error

controlCalidad = 0.05

plt.figure(figsize=[8,4], dpi=300)
for archivo in archivos_dat:
    D, radio, a, g_n, error = calculo_discrepancia_masa(archivo, controlCalidad)
    if error <= controlCalidad:
    # print(D)
    # plt.scatter(np.log(datos['Rad']), D)
        plt.scatter(g_n, D)
plt.ylim(-0.1, 20)
plt.xscale('log')
plt.xlabel(r'$g_N = V_b^2/r$ (m s$^{-2}$)')
plt.ylabel(r'Discrepancia de masa $= (V/V_b)^2$')
plt.axhline(y=1, xmin=0, xmax=1, color='black', alpha=0.5)
plt.title(f'Error <= {controlCalidad*100}%')
plt.tight_layout()
plt.savefig('figuras/D_vs_gn_v2.jpg')
plt.close()