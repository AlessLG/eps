'''
pandas==2.3.3
numpy==2.2.6

Este código se divide en tres:
1. Extracción de datos
2. Control de calidad
3. La física
4. Generación de gráficos

Este código:
Dibuja las curvas de rotación de v_obs y v_bar para toda la muestra del catálogo SPARC (175 galaxias). Se basa fuertemente en el paper de Lelli et al. (2016) para:
- Justificar Y_disk = 0.5 a lo largo de toda la muestra.
- Justificar Y_bulb = 1.4 * Y_disk.
- Calcular v_bar. v_bar = sqrt(abs(vgas)*vgas + Y_disk*abs(vdisk)*vdisk + Y_bul*abs(vbul)*vbul)

El gráfico de cada galaxia se guarda en una de tres carpetas dependiendo del error de sus datos: error = np.average(errv/vobs), donde vobs (velocidad observada) y errv (error en vobs) son vectores. Los rangos son: [0, 5], (5, 10] y (10, infinito). 
'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import os

ruta = os.getcwd()
archivos_dat = glob.glob(os.path.join(ruta, '../Rotmod_LTG/*.dat'))

ignorar_valores_negativos = True

controlCalidad = 0.05

i = 0
for archivo in archivos_dat:
    # Obtiene los datos de un archivo .dat determinado
    datos = pd.read_table(archivo, skiprows=3, sep=r'\s+', names=['Rad', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk',	'SBbul'])
    # Guarda los datos en variables 
    radio, vobs, errv, vgas, vdisk, vbul = datos['Rad'], datos['Vobs'], datos['errV'], datos['Vgas'], datos['Vdisk'], datos['Vbul']
    nombre_galaxia = os.path.basename(archivo).split('_')[0]
    
    error = 0
    error = np.average(errv/vobs)

    # Los valores de la velocidad bariónica se vuelven nan si corresponden a valores de vgas negativos
    if ignorar_valores_negativos:
        mascara_negativos = vgas < 0
        if mascara_negativos.any():
            ultimo_indice_vgas_negativo = np.where(mascara_negativos)[0][-1]
            # .iloc modifica las filas por su posición (desde 0 hasta el índice)
            vgas.iloc[:ultimo_indice_vgas_negativo + 1] = np.nan

    Y_disk = 0.5 # Según Lelli et al. (2016), valor más realista para el Near-InfraRed (NIR)
    Y_bul = 1.4 * Y_disk # Lelli et al. (2016)

    #vbar_2 = vbar^2
    vbar_2 = abs(vgas)*vgas + Y_disk*abs(vdisk)*vdisk + Y_bul*abs(vbul)*vbul
    vbar = np.sqrt(vbar_2)

    '''Se comentó el siguiente bloque de código, al no ser este un método
    científico para estimar correctamente el valor de Y_disk.'''
    # D = vobs**2/vbar_2
    # mascara_menor_uno = D < 1
    # Si en algún punto de la curva de rotación, D < 1 (es decir, vbar > vobs),
    # se recalcula vbar con un valor de Y_disk menor.
    # while mascara_menor_uno.any():
    #     Y_disk -= 0.01
    #     Y_bul = 1.4 * Y_disk
    #     vbar_2 = abs(vgas)*vgas + Y_disk * abs(vdisk) * vdisk + Y_bul * abs(vbul) * vbul
    #     D = vobs**2/vbar_2
    #     mascara_menor_uno = D < 1

    plt.figure(figsize=[6,3], dpi=300)
    plt.errorbar(radio, vobs, errv, fmt='o-k', alpha=0.3, label=r'$v_{obs}$')
    plt.scatter(radio, vbar, label=r'$v_{bar}$')
    plt.ylim(bottom=0)
    plt.xlabel(r'Radio [kpc]')
    plt.ylabel(r'$V$ [km/s]')
    # plt.title(f'{nombre_galaxia} | Error: ' f'{error*100:.1f}% | ' r'$\Upsilon_{\star}$: ' f'{Y_disk:.2f}')
    plt.title(f'{nombre_galaxia} | Error: ' f'{error*100:.1f}%')
    plt.legend()
    plt.tight_layout()
    # if error <= 0.05:
    #     plt.savefig(f'figuras/0err5/{nombre_galaxia}.jpg')
    # elif error > 0.05 and error <= 0.10:
    #     plt.savefig(f'figuras/5err10/{nombre_galaxia}.jpg')
    # else:
    #     plt.savefig(f'figuras/10err100/{nombre_galaxia}.jpg')
    if ((vbar > vobs + errv).any() == False): # Si ningún valor de vbar es mayor que vobs + errv, se dispara la condición.
        plt.savefig(f'figuras/bp/{nombre_galaxia}_{error*100:.1f}.jpg')
    else:
        plt.savefig(f'figuras/mp/{nombre_galaxia}_{error*100:.1f}.jpg')
    plt.close()

    i += 1
    if i % 10 == 0: print(f'{(i/len(archivos_dat)*100):.1f}%')

print('Terminado')
