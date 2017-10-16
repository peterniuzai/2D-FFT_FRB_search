import numpy as np
import matplotlib.pyplot as plt
import sys,os,time
import scipy.signal as signal
from DM_calculate import DM_calculate
def Signal_finding(data, ang_min = 5,ang_max = 85, pixel=2, T=0, nbin=0, fy=0):
         deg   = np.linspace(ang_min,ang_max,data.shape[1])
         data  = abs(data)
         lo    = np.where(data == np.max(data))
         d_max = data.max()
#         for i in np.arange(-pixel,pixel):
#             for j in np.arange(-pixel,pixel):
#                d_max += data[lo[0][0]+i,lo[1][0]+j]
         #print data.max(),data.min(),data.std(),'*/*/*/*/*/'
         snr   = (d_max - data.mean())/data.std()
	 angle_2ndFFT = deg[lo[1][0]]
	 DM = DM_calculate(fy,abs(angle_2ndFFT),nbin,T)
	 if snr > 6 and 1500 > DM > 200:
		snr	= snr
		DM	= DM
	 else:
		snr	= 0
		DM	= 0
	  
    	 return snr , DM

if __name__ == '__main__':
        ang_min = -90
        ang_max = 0
        pixel = 3
        data  = np.load('/home/ycli/data/burst_data/filtered_short.npy')
        data  = np.load('/home/nch/plot_raw/wigglez_found/data_pick/0_0.npy')
        if   data.shape[1] == 4:
             data = [data[:,0,:]]
        else:
             data = [data]
        SNR , location = Signal_finding(data, pixel)
        print SNR,location
