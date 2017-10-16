import numpy as np


def FFT(data, Dim = 2 , msk_cycle = 0 ):
    if Dim == 2:   #1st 2-D FFT
           data   = np.fft.fft2(data)
           data   = np.nan_to_num(data)
           data   = np.fft.fftshift(data)/np.sqrt(data.size)

           #data   =  shift[1:shift.shape[0]/2+1,shift.shape[1]/2:]#oringinal
	   data	  = data[1:data.shape[0]/2+1,data.shape[1]/2:data.shape[1]/2+data.shape[0]/2]
           data[-1, :] = 0
           data[ :, 0] = 0
#	   if msk_cycle >0:
#           	for i in np.arange(msk_cycle):
#              	 	x_sum   =  np.abs(data[-100:,:]).sum(axis = 0)
#               		y_sum   =  np.abs(data[    :,:]).sum(axis = 1)
#              	 	x_max   =  np.argmax(x_sum)
#              		y_max   =  np.argmax(y_sum)
#               		data[:,	  x_max] = 0
#               		data[y_max,  :  ] = 0
#               y_max      =  np.argmax(data_short)/data_short.shape[1]
#               x_max      =  np.argmax(data_short)%data_short.shape[1]
#               if data.shape[0]-20  < y_max:
#                       data[-40:,x_max-1:x_max+2]=0
    elif Dim == 1: #2nd 1D-FFT along radius
           data   = np.fft.fft(data,axis=0)
           data   = data/np.sqrt(data.shape[0])
           data   = np.fft.fftshift(data,axes=0)
	
    return data

if __name__ == '__main__':
      
        msk_cycle = 2
        data  = np.load('/home/ycli/data/burst_data/filtered_short.npy')
#        data  = np.load('/home/nch/plot_raw/wigglez_found/data_pick/0_0.npy')
        if   data.shape[1] == 4:
             data = [data[:,0,:]]
        else:
             data = [data]
        re_sets  = FFT(data, 2 , msk_cycle)
        print re_sets[0].shape 

