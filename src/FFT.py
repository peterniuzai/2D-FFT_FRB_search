import numpy as np
import matplotlib.pyplot as plt

def FFT(L_fft,data, Dim = 2 , msk_cycle = 0):
    if Dim == 2:   #1st 2-D FFT
#Window Function
	   blackman_x    = np.blackman(data.shape[1])
	   blackman_y    = np.blackman(data.shape[0])

	   window  = blackman_y[:,None] * blackman_x[None,:]
#	   data   = data*window
	   data   = np.fft.rfft2(data,norm = 'ortho')
#	   data	  = data[int(-L_fft):, 1:int(L_fft)+1]
           data   = data[data.shape[0]/2-1:,1:data.shape[0]/2+1]
#	   data   = data[data.shape[0]/2:,:]
#           data   = data[data.shape[0]/4:, 1:data.shape[0]/4]
#	   exit(1)


           data[-1, :] = data.mean()
           data[ :, 0] = data.mean()
	   data[-3:,:4] = data.mean()


	   #data[-1,0] = 0#data.max()
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
           data   = np.fft.fft(data,axis=0,norm = 'ortho')
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

