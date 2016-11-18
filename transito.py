#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

flujos_ref = np.array([])
err_flujo = np.ara
Julian_day = np.array([])

for i in range(len(wasp_list)):
    wasp = fits.open(wasp_list[i])
    wasp = wasp[0].data
    science = np.divide(wasp - master_bias, master_flat)*np.mean(master_flat)
    wasp43 = stamp(science,786,1384,32)
    ref = stamp(786,1384,32)
    cxw, cyw = centroid(wasp43)
    cxr, cyr = centroid(ref)
    f_wasp, e_wasp = ap_phot(wasp43, cxw, cyw, 10, 30, 40)
    f_ref, e_ref = ap_phot(ref, cxr, cyr, 8, 12, 22)
    f_plot = f_wasp/f_ref
    e_plot = #agregar formula propagacion de errore
    flujo_ref = np.append(flujo_ref,f_plot)
    err_flujo = np.append(err_flujo,e_plot)
