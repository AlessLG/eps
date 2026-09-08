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

ignorar_valores_negativos = True

def ajuste(archivo):
    # Obtiene los datos de un archivo .dat determinado
    datos = pd.read_table(archivo, skiprows=3, sep=r'\s+', names=['Rad', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk',	'SBbul'])

    # Guarda los datos en variables 
    vobs, vgas, vdisk, vbul, errv, radio = datos['Vobs'], datos['Vgas'], datos['Vdisk'], datos['Vbul'], datos['errV'], datos['Rad']
    # print(vobs, vgas, vdisk, vbul)

    # FÍSICA
    # kpc -> m
    # radio = radio * 3.086e19
    # v: km/s -> m/s
    # vgas, vdisk, vbul, errv = 1e3*vgas, 1e3*vdisk, 1e3*vbul, 1e3*errv

    error = 0
    error = np.average(errv/vobs)

    if ignorar_valores_negativos:
        mascara_negativos = vgas < 0
        if mascara_negativos.any():
            ultimo_indice_vgas_negativo = np.where(mascara_negativos)[0][-1]
            # Los valores de la velocidad bariónica se vuelven nan para vgas negativos
            # .iloc modifica las filas por su posición (desde 0 hasta el índice)
            vgas.iloc[:ultimo_indice_vgas_negativo + 1] = np.nan
            vdisk.iloc[:ultimo_indice_vgas_negativo + 1] = np.nan
            vbul.iloc[:ultimo_indice_vgas_negativo + 1] = np.nan

    Y_disk = 0.5 # Según Lelli et al. (2016), valor más realista para el Near-InfraRed (NIR)
    Y_bul = 1.4 * Y_disk # Lelli et al. (2016)

    vbar_2 = abs(vgas)*vgas + Y_disk * abs(vdisk) * vdisk + Y_bul * abs(vbul) * vbul
    D = vobs**2/vbar_2

    mascara_menor_uno = D < 1
    '''Si en algún punto de la curva de rotación, D < 1 (es decir, vbar > vobs),
    se recalcula vbar con un valor de Y_disk menor.'''
    while mascara_menor_uno.any():
        Y_disk -= 0.01
        Y_bul = 1.4 * Y_disk
        vbar_2 = abs(vgas)*vgas + Y_disk * abs(vdisk) * vdisk + Y_bul * abs(vbul) * vbul
        D = vobs**2/vbar_2
        mascara_menor_uno = D < 1

    # m -> kpc
    # radio = radio * 3.241e-20
    # v: m/s -> km/s
    # vbar, vgas, vdisk, vbul, errv = 1e-3*vbar_2, 1e-3*vgas, 1e-3*vdisk, 1e-3*vbul, 1e-3*errv
   
    return vobs, np.sqrt(vbar_2), radio, errv, error, Y_disk

controlCalidad = 0.05
# nombre_galaxia = 'UGC02953'
# nombre_galaxia = 'F563-V1'
# # LSB
# nombre_galaxia = 'F563-1'
# nombre_galaxia = 'F563-V2'
# nombre_galaxia = 'UGC03205'

i = 0
for archivo in archivos_dat:
    nombre_galaxia = os.path.basename(archivo).split('_')[0]
    vobs, vbar, radio, errv, error, Y_disk = ajuste(f'../Rotmod_LTG/{nombre_galaxia}_rotmod.dat')
    plt.figure(figsize=[6,3], dpi=300)
    plt.errorbar(radio, vobs, errv, fmt='o-k', alpha=0.3, label=r'$v_{obs}$')
    plt.scatter(radio, vbar, label=r'$v_{bar}$')
    plt.ylim(bottom=0)
    plt.xlabel(r'Radio [kpc]')
    plt.ylabel(r'$V$ [km/s]')
    # plt.title(f'{nombre_galaxia}', 'Error:' r'$\leq$' f'{error*100}%')
    plt.title(f'{nombre_galaxia} | Error: ' f'{error*100:.1f}% | ' r'$\Upsilon_{\star}$: ' f'{Y_disk:.2f}')
    plt.legend()
    plt.tight_layout()
    if error <= 0.05:
        plt.savefig(f'figuras/0err5/{nombre_galaxia}.jpg')
    elif error > 0.05 and error <= 0.10:
        plt.savefig(f'figuras/5err10/{nombre_galaxia}.jpg')
    else:
        plt.savefig(f'figuras/10err100/{nombre_galaxia}.jpg')
    plt.close()

    i += 1
    if i % 10 == 0: print(f'{(i/len(archivos_dat)*100):.1f}%')

print('Terminado')
