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

def calculo_discrepancia_masa(archivo):
    # Obtiene los datos de un archivo .dat determinado
    datos = pd.read_table(archivo, skiprows=3, sep=r'\s+', names=['Rad', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk',	'SBbul'])

    # Guarda los datos en variables 
    vobs, vgas, vdisk, vbul, radio = datos['Vobs'], datos['Vgas'], datos['Vdisk'], datos['Vbul'], datos['Rad']
    # print(vobs, vgas, vdisk, vbul)

    # Calcula la discrepancia de masa
    vbar_2 = vgas**2 + vdisk**2 + vbul**2 #vbar_2 = vbar^2
    D = vobs**2/vbar_2

    return D, radio

plt.figure(figsize=[8,4], dpi=300)
for archivo in archivos_dat:
    D, radio = calculo_discrepancia_masa(archivo)
    # print(D)
    # plt.scatter(np.log(datos['Rad']), D)
    plt.scatter(radio, D)
plt.ylim(-0.1, 20)
plt.xscale('log')
plt.xlabel('Radio (kpc)')
plt.ylabel(r'Discrepancia de masa $(V/V_b)^2$')
plt.tight_layout()
plt.savefig('figuras/D_vs_radio.jpg')
plt.close()