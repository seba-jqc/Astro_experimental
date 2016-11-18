#!/usr/bin/env python
# -*- coding: utf-8 -*-


import scipy as sp
import scipy.ndimage as spn
from astropy.io import fits
import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from zscale import zscale

wasp_list = glob.glob("wasp43/wasp43*fits")
master_bias = fits.open("MasterBias.fits")
master_bias = master_bias[0].data
master_flat = fits.open("MasterFlat.fits")
master_flat = master_flat[0].data
wasp = fits.open(wasp_list[0])
wasp = wasp[0].data
science = np.divide(wasp - master_bias, master_flat)*np.mean(wasp)


def stamp(data,x,y,r):
    """
    Crea una estampilla centrada en x,y de radio SR
    data: (N,)array_like
    Arreglo en 2D que representa los datos de una imagen fits
    x: int
    posicion en el eje x sobre el cual centrar la estampilla
    y: int
    posicion en el eje y sobre el cual centrar la estampilla
    r: int
    Radio de la estampilla
    """
    stamp=data[y-r:y+r, x-r:x+r]
    return stamp


def centroid(stamp):
    cx, cy = spn.measurements.center_of_mass(stamp)
    return int(cx), int(cy)


def centroid2(stamp):
    largo_x = np.shape(stamp)[1]
    largo_y = np.shape(stamp)[0]
    vect_x = np.linspace(0,largo_x-1,largo_x)
    vect_y = np.linspace(0,largo_y-1,largo_y)
    centro_y = np.median(vect_y)
    flujo_tot_x = np.sum(vect_x*stamp[centro_y, :])
    cx = flujo_tot_x/np.sum(stamp[centro_y, :])
    flujo_tot_y  = np.sum(vect_y*stamp[:, cx])
    cy = flujo_tot_y/np.sum(stamp[:, cx])
    return cx, cy


def perfil_radial(stamp, cx, cy):
    largo_x = len(stamp[1])
    pixel_x = np.arange(0,largo_x-cx)
    flujo_x = stamp[cy,cx:largo_x]
    plt.plot(pixel_x,flujo_x)
    plt.xlabel('Distancia al CM [pixeles]')
    plt.ylabel('Cuentas[ADU]')
    plt.title('Perfil Radial estrella de referencia')
    plt.savefig('perfilradialref.png')
    plt.show()


def distancia(x1, x2, y1, y2):
    dist = np.sqrt((x2-x1)**2+(y2-y1)**2)
    return dist


def ap_phot(stamp, cx, cy, apert, sk1, sk2):
    largo = np.shape(stamp)[1] #imagen es cuadrada
    vect = np.linspace(0,largo-1,largo)
    estrella = np.array([])
    cielo = np.array([])
    for i in range(len(vect)):
        for j in range(len(vect)):
            if distancia(cx,i,cy,j) <= apert:
                estrella = np.append(estrella, stamp[j,i])
            if distancia(cx,i,cy,j) >= sk1 and distancia(cx,i,cy,j) <= sk2:
                cielo = np.append(cielo,stamp[j,i])
    flujo_cielo = np.mean(cielo)
    flujo_est = np.sum(estrella-flujo_cielo)
    gain = 1
    N_A = len(estrella)
    N_B = len(cielo)
    error_b = np.std(cielo)
    error = np.sqrt(flujo_est/gain + N_A*error_b**2 +(N_A**2/N_B)*error_b**2)
    return flujo_est, error

#wasp43 (1064,991)
#ref_star (1810,1160)


stam=stamp(science,1064,991,50)
stam2=stamp(science,786,1384,50)
mn,mx = zscale(stam2)
plt.figure()
plt.imshow(stam2,vmin=mn,vmax=mx)
plt.savefig('estrella_ref2.png')
plt.show()
#ganancia = 1
#cx, cy = centroid(stam2)
#perfil_radial(stam2, cx, cy)
#flujo_est,err  = ap_phot(stam, cx, cy, 10, 30, 40)
#print flujo_est
#print err
#apert = 15
#sk1 = 30
#sk2 = 45
