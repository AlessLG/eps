import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from scipy.optimize import curve_fit

# Vectores para almacenado de datos
A_OBS = np.array([], dtype='f')
A_BAR = np.array([], dtype='f')

a_0 = 1.2e-10
Y_disk = 0.5 # Según Lelli et al. (2016), valor más realista para el Near-InfraRed (NIR)
Y_bul = 1.4 * Y_disk # Lelli et al. (2016)

def nu_simple_bar(a_bar, a_0):
    y = a_bar / a_0
    return 0.5 + 0.5 * np.sqrt(1 + 4/y)

# Para el ajuste a los residuales
def gaussiana(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2)) / np.sqrt(2 * np.pi * sigma**2)

archivos = glob.glob(os.path.join(os.getcwd(), '../Rotmod_LTG/*.dat'))

# for archivo in archivos[1:2]:
for archivo in archivos[:]:
    gal = np.loadtxt(archivo)
    radio, vobs, errv, vgas, vdisk, vbul = gal[:,0], gal[:,1], gal[:,2], gal[:,3], gal[:,4], gal[:,5]

    # kpc -> m
    radio = radio * 3.086e19
    # v: km/s -> m/s
    vobs, vgas, vdisk, vbul, errv = 1e3*vobs, 1e3*vgas, 1e3*vdisk, 1e3*vbul, 1e3*errv

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

ejes = np.logspace(start=-12.5, stop=-8, num=75, base=10.0)
amond = ejes * nu_simple_bar(ejes, a_0)
# discrepancia = np.array([])
# for j in range(len(A_OBS)):
#     d = A_OBS[j] - A_BAR[j] * nu_simple_bar(A_BAR[j], a_0)
#     discrepancia = np.append(discrepancia, d)
plt.figure(figsize=(5, 5), dpi=200)
plt.plot(ejes, amond)

plt.hist2d(A_BAR, A_OBS, bins=[ejes, ejes], cmap='Blues')
plt.axline((1, 1), slope=1)
plt.xscale('log')
plt.yscale('log')
# plt.xlabel(r'$\log{(a_{bar})}$')
# plt.ylabel(r'$\log{(a_{obs})}$')
plt.xlabel(r'$a_{bar}$')
plt.ylabel(r'$a_{obs}$')
plt.tight_layout()
plt.gca().set_box_aspect(1)
plt.savefig('hist2d.pdf')

mask = (A_OBS > 0) & (A_BAR > 0)
res = np.log10(A_OBS[mask]) - np.log10(A_BAR[mask] * nu_simple_bar(A_BAR[mask], a_0))

filtro = np.abs(res) <= 0.7
res = res[filtro]

bins = np.linspace(-0.7, 0.7, 40)

plt.figure(figsize=(5, 5), dpi=200)

conteos, bordes, _ = plt.hist(res, bins=bins, color='cornflowerblue', edgecolor='k', linewidth=0.3)

centros = (bordes[:-1] + bordes[1:]) / 2

p0_inicial = [np.max(conteos), 0.0, 0.13]
popt, pcov = curve_fit(gaussiana, centros, conteos, p0=p0_inicial)

A_fit, mu_fit, sigma_fit = popt
sigma_fit = abs(sigma_fit)

x_suave = np.linspace(-0.7, 0.7, 200)
plt.plot(x_suave, gaussiana(x_suave, *popt), 'k', linewidth=1.5)
plt.text(0.05, 0.90, rf'$\sigma = {sigma_fit:.2f}$ dex', transform=plt.gca().transAxes)
plt.xlabel(r'Residuales [dex]')
plt.ylabel('Mediciones')
plt.tight_layout()
plt.savefig('residuales.pdf')