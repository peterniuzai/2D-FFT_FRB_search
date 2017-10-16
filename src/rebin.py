import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
#import gc

def rebin(data,freq, nbin,tx=0):
#Change the freqency axis to wave squre axis
#The signal line will become a straight line after this process.
    fy           = freq**-2
    f_axis       = np.linspace(fy.min(),fy.max(),nbin)
    b_aray,b_edg = np.histogram(fy,bins=nbin)
#    data	 = np.nan_to_num(data)
    data1	 = np.zeros((nbin,data.shape[1]),data.dtype)
    for i in np.arange(nbin):
            length  =  b_aray[i]
            if   i == 0:
                 index = 0
            elif i == 1:
                 index = b_aray[0]
            else:
                 index = b_aray[0:i].sum()#when i = 1 ,the bin_array[0:1] doesn't include the bin_array[1]
            for ii in np.arange(length):
                 data1[i,:]+=data[ii+index,:]
    tem = np.ones(data.shape[1])
    bin_weight = b_aray[:,None] * tem[None,:]
    data1 = data1 / 1.0 / bin_weight
    data1 = np.nan_to_num(data1)
    #for i in np.arange(5):
    #           y_s     =  data1.sum( axis = 1 )
    #           y_max   =  np.argmax(y_s)
    #           data1[y_max-10:y_max+11,:]=0
#    del tem,data
#    gc.collect()

    return data1, f_axis, fy

def rebin_inter(data, freq, nbin, t_axis):
    fy           = freq**-2
    fyy 	 = fy.reshape(-1,1)
    f_axis       = np.linspace(fy.min(),fy.max(),nbin)
    tx		 = t_axis
    data    	 = data.reshape(-1)
#begin to interpolate the data
    c    = np.broadcast_arrays(fyy,tx) #create the coordinates of the grid, the c[0].shape = c[1].shape = [4096,2048]
    f    = c[0].reshape(-1)
    t    = c[1].reshape(-1)
    d    = np.nan_to_num(data)
    points         = (f,t)
    grid_f, grid_t = np.mgrid[f.min():f.max():nbin*1j, t.min():t.max():nbin*1j]
    data1          = griddata(points,d,(grid_f,grid_t),method='nearest')
    data1          = np.nan_to_num(data1)
    #for i in np.arange(5):
    #           y_s     =  data1.sum( axis = 1 )
    #           y_max   =  np.argmax(y_s)
    #           data1[y_max-10:y_max+11,:]=0

    return data1, f_axis, fy


if __name__ == '__main__':

	freq  = np.load('../data/freq.npy')
	data  = np.load('../data/filtered_short.npy')
	if   data.shape[1] == 4:
	     data = [data[:,0,:]]
	else:
	     data = [data]
        t_sets    = [np.arange(data[0].shape[1])] 
	re_sets , f_axis ,nbin = rebin(data ,t_sets, freq)
	print re_sets[0].shape ,len(f_axis)
        plt.pcolormesh(re_sets[0])
        plt.colorbar()
        plt.show()
        np.save('../data/rebin.npy',re_sets[0])
