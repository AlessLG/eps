'''
numpy==2.2.6

Este código:
Dibuja las curvas de rotación de v_obs y v_bar para toda la muestra del catálogo SPARC (175 galaxias). Se basa fuertemente en el paper de Lelli et al. (2016) para:
- Justificar Y_disk = 0.5 a lo largo de toda la muestra.
- Justificar Y_bulb = 1.4 * Y_disk.
- Calcular v_bar. v_bar = sqrt(abs(vgas)*vgas + Y_disk*abs(vdisk)*vdisk + Y_bulb*abs(vbul)*vbul)

El gráfico de cada galaxia se guarda en:
 - Una de tres carpetas dependiendo del error de sus datos: error = np.average(errv/vobs), donde vobs (velocidad observada) y errv (error en vobs) son vectores. Los rangos son: [0, 5], (5, 10] y (10, infinito). 
 - Una de dos carpetas ningún valor de vbar es mayor que vobs + errv (para el mismo radio): bp (galaxias bien portadas) si esto no ocurre y mp (mal portadas) si sí.
'''

import matplotlib.pyplot as plt
import numpy as np
import glob
import os

archivos = glob.glob(os.path.join(os.getcwd(), '../Rotmod_LTG/*.dat'))

a_0 = 1.2e-10
Y_disk = 0.5 # Según Lelli et al. (2016), valor más realista para el Near-InfraRed (NIR)
Y_bul = 1.4 * Y_disk # Lelli et al. (2016)

error = 0
i = 0
ppi = 100

def nu_simple_bar(r, v_b, a0):
    y = v_b**2 / (r * a0) # y = a_bar/a_0
    return (1 + np.sqrt(1 + 4/y))/2

def nu_estandar_bar(r, v_b, a0):
    y = v_b**2 / (r * a0) # y = a_bar/a_0
    return np.sqrt((1 + np.sqrt(1 + 4/y**2))/2)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=[6,5], dpi=ppi, sharex=True, gridspec_kw={'height_ratios': [1, 0.3], 'hspace': 0})
# fig.subplots_adjust(hspace=.0)

for archivo in archivos[:]:
    gal = np.loadtxt(archivo)
    radio, vobs, errv, vgas, vdisk, vbul = gal[:,0], gal[:,1], gal[:,2], gal[:,3], gal[:,4], gal[:,5]

    nombre_galaxia = os.path.basename(archivo).split('_')[0]
    
    error = np.average(errv/vobs)

    # identifica el índice del último valor de vgas < 0, convierte a nan todos los valores de vgas desde 0 hasta dicho índice
    mascara_negativos = vgas < 0
    if mascara_negativos.any():
        ultimo_indice_vgas_negativo = np.where(mascara_negativos)[0][-1]
        vgas[:ultimo_indice_vgas_negativo + 1] = np.nan

    # FÍSICA
    # kpc -> m
    radio = radio * 3.086e19
    # v: km/s -> m/s
    vobs, vgas, vdisk, vbul, errv = 1e3*vobs, 1e3*vgas, 1e3*vdisk, 1e3*vbul, 1e3*errv

    #vbar_2 = vbar^2
    vbar_2 = abs(vgas)*vgas + Y_disk*abs(vdisk)*vdisk + Y_bul*abs(vbul)*vbul
    vbar = np.sqrt(vbar_2)
    vmond = vbar * np.sqrt(nu_simple_bar(radio, vbar, a_0))
    vmond2 = vbar * np.sqrt(nu_estandar_bar(radio, vbar, a_0))

    radio = radio * (1/3.086e19)
    vobs, vbar, vmond, vmond2, errv = 1e-3*vobs, 1e-3*vbar, 1e-3*vmond, 1e-3*vmond2, 1e-3*errv

    '''Se comentó el siguiente bloque de código, al no ser este un método científico para estimar correctamente el valor de Y_disk.'''
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
    # plt.title(f'{nombre_galaxia} | Error: ' f'{error*100:.1f}% | ' r'$\Upsilon_{\star}$: ' f'{Y_disk:.2f}')

    ax1.clear()
    ax2.clear()
    ax1.errorbar(radio, vobs, errv, fmt='o-k', label=r'$v_{obs}$', capsize=3)
    ax1.plot(radio, vbar, label=r'$v_{bar}$', marker='o')
    ax1.plot(radio, vmond, label=r'$v_{mond, sim}$', marker='o', color='red', linestyle='dashed')
    ax1.plot(radio, vmond2, label=r'$v_{mond, est}$', marker='o', color='blue', linestyle='dashed')
    ax1.set_ylim(bottom=0)
    # ax1.get_xaxis().set_visible(False)
    ax1.set_ylabel(r'$V$ [km/s]')
    ax1.set_title(f'{nombre_galaxia} | Error: ' f'{error*100:.1f}%')
    ax1.legend()

    ax2.scatter(radio, (vmond - vobs), color='red')
    ax2.scatter(radio, (vmond2 - vobs), color='blue')
    ax2.set_xlabel(r'Radio [kpc]')
    ax2.set_ylabel('Residuales')
    ax2.axhline(y=0, xmin=0, xmax=1, color='black', alpha=0.5, linestyle='dotted')

    fig.tight_layout()
    # if error <= 0.05:
    #     fig.savefig(f'figuras/0err5/{nombre_galaxia}.jpg')
    # elif error > 0.05 and error <= 0.10:
    #     fig.savefig(f'figuras/5err10/{nombre_galaxia}.jpg')
    # else:
    #     fig.savefig(f'figuras/10err100/{nombre_galaxia}.jpg')
    if ((vbar > vobs + errv).any() == False): # Si ningún valor de vbar es mayor que vobs + errv, se dispara la condición.
        fig.savefig(f'figuras/bp/{nombre_galaxia}_{error*100:.1f}.jpg')
    else:
        fig.savefig(f'figuras/mp/{nombre_galaxia}_{error*100:.1f}.jpg')

    i += 1
    if i % 10 == 0: print(f'{(i/len(archivos)*100):.1f}%')

print('Terminado')