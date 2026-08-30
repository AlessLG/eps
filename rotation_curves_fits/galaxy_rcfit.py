# pandas==2.3.3
# numpy==2.2.6

# print(datos)
# print(datos.dtypes)
# datos.info()
# print(datos['Vobs'].shape)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import os

ruta = os.getcwd()
archivos_dat = glob.glob(os.path.join(ruta, '../Rotmod_LTG/*.dat'))

def ajuste(archivo):
    # Obtiene los datos de un archivo .dat determinado
    datos = pd.read_table(archivo, skiprows=3, sep=r'\s+', names=['Rad', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk',	'SBbul'])

    # Guarda los datos en variables 
    vobs, vgas, vdisk, vbul, errv, radio = datos['Vobs'], datos['Vgas'], datos['Vdisk'], datos['Vbul'], datos['errV'], datos['Rad']
    # print(vobs, vgas, vdisk, vbul)

    error = 0

    # Calcula la velocidad debido a las componentes bariónicas
    vbar_2 = vgas**2 + vdisk**2 + vbul**2

    # kpc -> m
    # radio = radio * 3.086e19
    # v: km/s -> m/s
    # vobs, vgas, vdisk, vbul, errv = 1e3*vobs, 1e3*vgas, 1e3*vdisk, 1e3*vbul, 1e3*errv

    error = np.average(errv/vobs)

    return vobs, np.sqrt(vbar_2), radio, errv, error

controlCalidad = 0.05
# nombre_galaxia = 'UGC02953'
# nombre_galaxia = 'F563-V1'
# # LSB
# nombre_galaxia = 'F563-1'
# nombre_galaxia = 'F563-V2'

for archivo in archivos_dat:
    nombre_galaxia = os.path.basename(archivo).split('_')[0]
    vobs, vbar, radio, errv, error = ajuste(f'../Rotmod_LTG/{nombre_galaxia}_rotmod.dat')

    plt.figure(figsize=[6,3], dpi=300)
    plt.errorbar(radio, vobs, errv, fmt='o-k', alpha=0.3, label=r'$v_{obs}$')
    plt.scatter(radio, vbar, label=r'$v_{bar}$')
    plt.ylim(bottom=0)
    plt.xlabel(r'Radio [kpc]')
    plt.ylabel(r'$V$ [km/s]')
    # plt.title(f'{nombre_galaxia}', 'Error:' r'$\leq$' f'{error*100}%')
    plt.title(f'{nombre_galaxia} | Error: ' f'{error*100:.1f}%')
    plt.legend()
    plt.tight_layout()
    if error <= 0.05:
        plt.savefig(f'figuras/0err5/{nombre_galaxia}.jpg')
    elif error > 0.05 and error <= 0.10:
        plt.savefig(f'figuras/5err10/{nombre_galaxia}.jpg')
    else:
        plt.savefig(f'figuras/10err100/{nombre_galaxia}.jpg')
    plt.close()

print('Terminado')

# plt.figure(figsize=[8,4], dpi=300)
# for archivo in archivos_dat:
#     D, radio, a, g_n, error = ajuste(archivo, controlCalidad)
#     if error <= controlCalidad:
#     # print(D)
#     # plt.scatter(np.log(datos['Rad']), D)
#     # radio = radio/3.086e16
#         plt.scatter(radio, D)
# plt.ylim(-0.1, 20)
# plt.xscale('log')
# plt.xlabel(r'Radio (m)')
# plt.ylabel(r'Discrepancia de masa $= (V/V_b)^2$')
# plt.axhline(y=1, xmin=0, xmax=1, color='black', alpha=0.5)
# plt.title(f'Error <= {controlCalidad*100}%')
# plt.tight_layout()
# plt.savefig('figuras/D_vs_radio.jpg')
# plt.close()

# plt.figure(figsize=[8,4], dpi=300)
# for archivo in archivos_dat:
#     D, radio, a, g_n, error = calculo_discrepancia_masa(archivo, controlCalidad)
#     if error <= controlCalidad:
#     # print(D)
#     # plt.scatter(np.log(datos['Rad']), D)
#         plt.scatter(a, D)
# plt.ylim(-0.1, 20)
# plt.xscale('log')
# plt.xlabel(r'Aceleración $= V^2/r$ (m s$^{-2}$)')
# plt.ylabel(r'Discrepancia de masa $= (V/V_b)^2$')
# plt.axhline(y=1, xmin=0, xmax=1, color='black', alpha=0.5)
# plt.title(f'Error <= {controlCalidad*100}%')
# plt.tight_layout()
# plt.savefig('figuras/D_vs_a.jpg')
# plt.close()

# plt.figure(figsize=[8,4], dpi=300)
# for archivo in archivos_dat:
#     D, radio, a, g_n, error = calculo_discrepancia_masa(archivo, controlCalidad)
#     if error <= controlCalidad:
#     # print(D)
#     # plt.scatter(np.log(datos['Rad']), D)
#         plt.scatter(g_n, D)
# plt.ylim(-0.1, 20)
# plt.xscale('log')
# plt.xlabel(r'$g_N = V_b^2/r$ (m s$^{-2}$)')
# plt.ylabel(r'Discrepancia de masa $= (V/V_b)^2$')
# plt.axhline(y=1, xmin=0, xmax=1, color='black', alpha=0.5)
# plt.title(f'Error <= {controlCalidad*100}%')
# plt.tight_layout()
# plt.savefig('figuras/D_vs_gn.jpg')
# plt.close()

