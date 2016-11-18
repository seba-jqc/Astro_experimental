
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np



Flat = glob.glob('wasp43/SkyFlatI_*.fits')
master_bias = fits.open("MasterBias.fits")
master_bias = master_bias[0].data
all_flats = sp.array ([fits.open(f)[0].data for f in Flat])
mflats= sp.mean(all_flats-master_bias,0)
all_flats=0 #ahorra memoria
fits.writeto("MasterFlat2.fits",mflats)
print mflats
