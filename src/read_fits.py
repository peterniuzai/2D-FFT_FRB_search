import readGBT
import pyfits
import numpy as np

f_dir = '/home/nch/work/burst_data/'
f_name   = 'crane_1711858_1.00_0067.fits'
hdulist = pyfits.open(f_dir + f_name)
data ,tx, ra, dec, az, el, freq = readGBT.read_fits(hdulist)
print data.shape
