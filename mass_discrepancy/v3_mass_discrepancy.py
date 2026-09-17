'''
pandas==2.3.3
numpy==2.2.6

Este código se divide en tres:
1. Extracción de datos
2. Control de calidad
3. La física
4. Generación de gráficos

Este código:
Grafica la discrepancia de masa para las galaxias que cumplan cuyas incertezas sean menores al criterio establecido por 'control_calidad'; respecto al radio, aceleración observacional y aceleración bariónica.

Convierte a NAN los datos de vbar (velocidad bariónica) que correspondan a valores negativos de vgas (velocidad del gas), puesto que estos no son representativos de la dinámica galáctica general. Posteriormente, realiza diversos cálculos y guarda los resultados en vectores. Con el propósito de luego ajustarles una curva datos dada por $\mu (x)^{-1}$.
'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import os
from scipy.optimize import curve_fit

# Controles del análisis
control_calidad = 0.05 # Menor valor, criterio más estricto. Si es igual a cero, se ignora
ignorar_valores_negativos = True

total_galaxias = 175
galaxias_consideradas = 0

# Arrays para guardado de datos
Discrepancia = np.array([], dtype='f')
Error_discrepancia = np.array([], dtype='f')
Radio = np.array([], dtype='f')
Aceleracion_observable = np.array([], dtype='f')
Aceleracion_barionica = np.array([], dtype='f')

# mu_simple(x) = x/1+x
# mu_estandar(x) = x/sqrt(1+x^2); x = a/a_0
# a es la variable independiente
def mu_inv_simple(a, a_0): return (a/a_0)/(1 - (a/a_0))
def mu_inv_estandar(a, a_0): return (a/a_0)/np.sqrt(1 - (a/a_0)**2)


# EXTRACCIÓN DE DATOS
ruta = os.getcwd()
archivos_dat = glob.glob(os.path.join(ruta, '../Rotmod_LTG/*.dat'))

for archivo in archivos_dat:
    # Obtiene los datos de un archivo .dat determinado
    datos = pd.read_table(archivo, skiprows=3, sep=r'\s+', names=['Rad', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk',	'SBbul'])

    radio, vobs, errv, vgas, vdisk, vbul = datos['Rad'], datos['Vobs'], datos['errV'], datos['Vgas'], datos['Vdisk'], datos['Vbul']

    # CONTROL DE CALIDAD
    error = np.average(errv/vobs)
    if error > control_calidad: continue
    galaxias_consideradas += 1

    # Los valores de la velocidad bariónica se vuelven nan si corresponden a valores de vgas negativos
    if ignorar_valores_negativos:
        mascara_negativos = vgas < 0
        if mascara_negativos.any():
            ultimo_indice_vgas_negativo = np.where(mascara_negativos)[0][-1]
            # .iloc modifica las filas por su posición (desde 0 hasta el índice)
            vgas.iloc[:ultimo_indice_vgas_negativo + 1] = np.nan
            # vdisk.iloc[:ultimo_indice_vgas_negativo + 1] = np.nan
            # vbul.iloc[:ultimo_indice_vgas_negativo + 1] = np.nan
    
    # FÍSICA
    # kpc -> m
    radio = radio * 3.086e19
    # v: km/s -> m/s
    vobs, vgas, vdisk, vbul, errv = 1e3*vobs, 1e3*vgas, 1e3*vdisk, 1e3*vbul, 1e3*errv

    Y_disk = 0.5 # Según Lelli et al. (2016), valor más realista para el Near-InfraRed (NIR)
    Y_bul = 1.4 * Y_disk # Lelli et al. (2016)

    # vbar_2 = vbar^2
    vbar_2 = abs(vgas)*vgas + Y_disk*abs(vdisk)*vdisk + Y_bul*abs(vbul)*vbul
    D = vobs**2/vbar_2
    errD = 2*D*errv/vobs # Asumiendo que los errores en vbar son cero
    a_obs = vobs**2/radio
    a_bar = vbar_2/radio
    # return D, radio, a_obs, a_bar, error, errD

    Discrepancia            = np.append(Discrepancia, D)
    Error_discrepancia      = np.append(Error_discrepancia, errD)
    Radio                   = np.append(Radio, radio)
    Aceleracion_observable  = np.append(Aceleracion_observable, a_obs)
    Aceleracion_barionica   = np.append(Aceleracion_barionica, a_bar)

# Discrepancia            = Discrepancia[~np.isnan(Discrepancia)]
# Error_discrepancia      = Error_discrepancia[~np.isnan(Error_discrepancia)]
# Aceleracion_barionica   = Aceleracion_barionica[~np.isnan(Aceleracion_barionica)]


# print(len(Discrepancia          ))
# print(len(Error_discrepancia    ))
# print(len(Radio                 ))
# print(len(Aceleracion_observable))
# print(len(Aceleracion_barionica ))

popt_aobs, pcov_aobs = curve_fit(mu_inv_simple, Aceleracion_observable, Discrepancia, bounds=(10**(-11), 10**(-9)), nan_policy='omit')
popt_abar, pcov_abar = curve_fit(mu_inv_simple, Aceleracion_barionica, Discrepancia, bounds=(10**(-11), 10**(-9)), nan_policy='omit')

print(popt_aobs, pcov_aobs)
print(popt_abar, pcov_abar)

# GENERACIÓN DE GRÁFICOS
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10), dpi=300)
x = np.linspace(10**(-8), 10**(-11))

ax1.errorbar(Radio, Discrepancia, Error_discrepancia, fmt='k.', alpha=0.2) 

ax2.errorbar(Aceleracion_observable, Discrepancia, Error_discrepancia, fmt='k.', alpha=0.2)
ax2.plot(x, mu_inv_simple(x, popt_aobs), color='red')

ax3.errorbar(Aceleracion_barionica, Discrepancia, Error_discrepancia, fmt='k.', alpha=0.2)
ax3.plot(x, mu_inv_simple(x, popt_abar), color='red')

for ax in (ax1, ax2, ax3):
    ax.set_ylim(-1.5, 20)
    ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.set_ylabel(r'$\mathcal{D}=(V/V_b)^2$')
    ax.axhline(y=1, xmin=0, xmax=1, color='black', alpha=0.5)

ax1.set_xlabel(r'$r$ [m]')
ax2.set_xlabel(r'$a_{obs} = V_{obs}^2/r$ [m s$^{-2}$]')
ax3.set_xlabel(r'$a_{bar} = V_{bar}^2/r$ [m s$^{-2}$]')

fig.suptitle(f'Error' r'$\leq$' f'{control_calidad*100}% - {galaxias_consideradas} galaxias consideradas ({(galaxias_consideradas/total_galaxias)*100:.1f}%)')
fig.tight_layout()
fig.savefig('figuras/D_vs_radio-aobs-abar_nf.jpg') #nf = nueva fórmula

print(f'De {total_galaxias} galaxias, {galaxias_consideradas} fueron consideradas ({(galaxias_consideradas/total_galaxias)*100:.1f}%). ')
plt.close()