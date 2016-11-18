import matplotlib.pyplot as plt
import numpy as np
import pyfits as pf
import glob as glob
import scipy as sp

Bias = glob.glob('wasp43/Bias_*.fits')
all_bias = sp.array ([pf.open(f)[0].data for f in Bias])
mbias= sp.mean(all_bias,0)
all_bias=0 #ahorra memoria
pf.writeto("MasterBias.fits",mbias)
print mbias
