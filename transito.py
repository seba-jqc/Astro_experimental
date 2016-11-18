#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import glob
from codigo import *

wasp_list = glob.glob("wasp43/wasp43*fits")
master_bias = fits.open("MasterBias.fits")
master_bias = master_bias[0].data
master_flat = fits.open("MasterFlat.fits")
master_flat = master_flat[0].data

flujos_wasp43 = np.array([])
flujos_ref = np.array([])
err_wasp43 = np.array([])
err_ref = np.array([])
flujos_final = np.array([])
err_final = np.array([])
Julian_day = np.array([])

for i in range(len(wasp_list)):
    HDUlist = fits.open(wasp_list[i])
    scidata = HDUlist[0].data
    JD = HDUlist[0].header['JD']
    Julian_day = np.append(Julian_day, JD)
    science = np.divide(scidata - master_bias, master_flat)*np.mean(master_flat)
    wasp43 = stamp(science, 1064, 991, 50)
    ref = stamp(science, 786, 1384, 32)
    cxw, cyw = centroid(wasp43)
    cxr, cyr = centroid(ref)
    f_wasp, e_wasp = ap_phot(wasp43, cxw, cyw, 10, 30, 40)
    f_ref, e_ref = ap_phot(ref, cxr, cyr, 8, 12, 22)
    flujos_wasp43 = np.append(flujos_wasp43, f_wasp)
    flujos_ref = np.append(flujos_ref, f_ref)
    err_wasp43 = np.append(err_wasp43, e_wasp)
    err_ref = np.append(err_ref, e_ref)
    f_final = f_wasp/f_ref
    e_final = (f_wasp/f_ref)*np.sqrt((e_wasp/f_wasp)**2 + (e_ref/f_ref)**2)
    flujos_final = np.append(flujos_final, f_final)
    err_final = np.append(err_final, e_final)

print err_ref

plt.figure(1)
plt.clf()
plt.errorbar(Julian_day, flujos_wasp43, yerr=err_wasp43, ls='None', marker='.',
             capsize=0, label='wasp43')
plt.errorbar(Julian_day, flujos_ref, yerr=err_ref, ls='None', marker='.',
             capsize=0, label='ref_star')
plt.xlabel('Tiempo [dias julianos]')
plt.ylabel('Flujo total [ADU]')
plt.legend(loc = 'upper left')
plt.title('comparacion de flujos')
plt.savefig('comp_flujos.png')
plt.show()


plt.figure(2)
plt.clf()
plt.errorbar(Julian_day, flujos_ref, yerr=err_ref, ls='None', marker='.',
             capsize=0)
plt.xlabel('Tiempo [dias julianos]')
plt.ylabel('Flujo total [ADU]')
plt.title('Registro de flujo Estrella de referencia')
plt.savefig('flujo_referencia.png')
plt.show()


plt.figure(3)
plt.clf()
plt.errorbar(Julian_day, flujos_wasp43, yerr=err_wasp43, ls='None', marker='.',
             capsize=0)
plt.xlabel('Tiempo [dias julianos]')
plt.ylabel('Flujo total [ADU]')
plt.title('Registro de flujo wasp43')
plt.savefig('flujo_wasp43.png')
plt.show()

plt.figure(4)
plt.clf()
plt.errorbar(Julian_day, flujos_final, yerr=err_final, ls='None', marker='.',
             capsize=0)
plt.xlabel('Tiempo [dias julianos]')
plt.ylabel('Flujo total [ADU]')
plt.title('Cociente de flujo wasp43/estrella de referencia')
plt.savefig('cociente_flujo.png')
plt.show()
