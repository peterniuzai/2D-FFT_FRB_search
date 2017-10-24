import numpy as np
from sigpyproc.Readers import FilReader
from DM_calculate import angle_range
from DM_calculate import time_delay
from DM_calculate import length_calculate
#import readGBT
#import pyfits
#import h5py
def read_data(f_dir, f_name, t_len, nbin, comm_size=10,DM_range=[100,1500],Wp=3,ang=[0,0]):

#load data,chunked them to data sets with shape :nbin * t_length (Default got from DM range). Wp means pulse width in (ms).

	 fil	= FilReader(f_dir + f_name)
	 hdr	= fil.header
	 ftop	= hdr['ftop']
	 fbot	= hdr['fbottom']
	 nch	= hdr['nchans']
	 t_rsl  = hdr['tsamp']*1000 # unit (ms)
	 nsamp	= hdr['nsamples']
	 C	= 4.148908e6 # (ms)
	 t_de	= time_delay(DM_range,fbot,ftop)
         #print t_de,'ms'
	 t_gulp	= round(t_de / t_rsl)

	 if nbin != 0:
		 nbin = nbin
	 else:
		 nbin = 2**(int(np.log2(nch)))

         if t_len != 0:
                 t_len  = t_len
         else:
                 t_len  = nbin

		

	 freq	= np.linspace(ftop,fbot,nch)
	 fy	= freq**-2

#Multi processors parameter
         num    = int(nsamp/t_len)  #Total number of chunks to process 
         p_n    = int(num/comm_size) #Number of chunks for each thread to process
	 T	= t_len * t_rsl
#Image parameter
	 FFT_rs = ((1./nbin)**2 + (1./t_len)**2)**0.5
	 Ang_rs	= ((nbin**2 + t_len**2)**-0.5)/FFT_rs #Angle resolution constraind by signal line length
	 Rad_rs	= 1
         angle	= angle_range(fy,DM_range,nbin,T)
	 if (ang[0] != 0) or (ang[1] !=0):
		angle = ang

	 if int(Wp/t_rsl) == 0:
		Wp = 1
	 else:
		Wp = int(Wp/t_rsl) 

	 L_fft	= length_calculate(fy,t_rsl ,DM_range,nbin,Wp, FFT_rs)
	 part	= (angle[1]-angle[0])/90.
	 t_gulp = round(t_gulp/t_len)+1
	 ang_rs = (2**0.5/L_fft)*180/np.pi
 
         print 'Angle:', angle
         print 'L_fft:', L_fft
         print 'Rad_rs', Rad_rs
         print 'Width resolution:Ang_rs', Ang_rs
         print 'Matrix Ang resoulution',ang_rs

	 if Ang_rs > ang_rs :
		Ang_rs = ang_rs
		print 'Using matrix angle resolution..'
	 
	# Radius resolutiongrid size for interpolate in polar-coordin transform.
     	# Angle resolution for interpolate in polar-coordin transform.
	
         return fil, num, p_n, freq, t_rsl, t_len, nbin, nch , T , fy,angle, Ang_rs,Rad_rs, L_fft,t_gulp

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
