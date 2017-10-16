import numpy as np
from sigpyproc.Readers import FilReader
#import readGBT
#import pyfits
#import h5py
def read_data(f_dir, f_name, t_len, nbin, comm_size=10):

#load data from GBT, chunked them to data sets with shape :t_length * t_length (Default :2048*2048)
#Then distribute them to multiprocess, 
#Finally,return the data sets

	 f	= FilReader(f_dir + f_name)
	 hdr	= f.header
	 ftop	= hdr['ftop']
	 fbot	= hdr['fbottom']
	 nch	= hdr['nchans']
	 t_rsl  = hdr['tsamp']*1000 # unit (ms)
	 nsamp	= hdr['nsamples']
         if t_len != 0:
		 t_len	= t_len
	 else:
		 t_len	= nch
	 if nbin != 0:
		 nbin = nbin
	 else:
		 nbin = nch
	 freq	= np.linspace(ftop,fbot,nch)
         num    = nsamp/t_len  #Total number of chunks to process 
         p_n    = num/comm_size #Number of chunks for each thread
	# rad_grid  =  1         #radius resolutiongrid size for interpolate in polar-coordin transform.
     	# ang_grid  = (1./2**0.5/data.shape[1])*180/np.pi   #angle resolution for interpolate in polar-coordin transform.
         return f, num, p_n, freq, t_rsl, t_len, nbin, nch

if   __name__ == '__main__':
     f_dir  = '../data/'
     f_name = 'data_2017-08-30_17-35-36.fil'
     f, num,p_n,freq, t_rsl, t_len = read_data(f_dir ,f_name ,0,0,10)
     print freq.shape
     print p_n
     b = f.readBlock(0,100)
     print b.shape
     print 'load over!'
     exit()
