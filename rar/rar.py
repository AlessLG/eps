import matplotlib.pyplot as plt
import numpy as np
import glob
import os

# Vectores para almacenado de datos
A_OBS = np.array([], dtype='f')
A_BAR = np.array([], dtype='f')

a_0 = 1.2e-10
Y_disk = 0.5 # Según Lelli et al. (2016), valor más realista para el Near-InfraRed (NIR)
Y_bul = 1.4 * Y_disk # Lelli et al. (2016)

archivos = glob.glob(os.path.join(os.getcwd(), '../Rotmod_LTG/*.dat'))

# for archivo in archivos[1:2]:
for archivo in archivos[:]:
    gal = np.loadtxt(archivo)
    radio, vobs, errv, vgas, vdisk, vbul = gal[:,0], gal[:,1], gal[:,2], gal[:,3], gal[:,4], gal[:,5]

    aobs = vobs**2/radio
    abar = abs(vgas)*vgas/radio + Y_disk*abs(vdisk)*vdisk/radio + Y_bul*abs(vbul)*vbul/radio

    mascara_negativos = vgas < 0
    if mascara_negativos.any():
        ultimo_indice_vgas_negativo = np.where(mascara_negativos)[0][-1]
        aobs = aobs[ultimo_indice_vgas_negativo + 1:]
        abar = abar[ultimo_indice_vgas_negativo + 1:]
    
    A_OBS = np.append(A_OBS, aobs)
    A_BAR = np.append(A_BAR, abar)

    # print(aobs)
    # print(abar)
    # print(len(aobs))
    # print(len(abar))

    # nombre_galaxia = os.path.basename(archivo).split('_')[0]
    # if (np.isnan(abar)).any():
    #     print(nombre_galaxia)
# print(len(A_OBS))
# print(len(A_BAR))

plt.figure(figsize=(8,8), dpi=300)
plt.hist2d(A_BAR, A_OBS)
plt.xscale('log')
plt.yscale('log')
plt.savefig('hist2d.pdf')
