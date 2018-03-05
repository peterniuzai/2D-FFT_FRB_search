import numpy as np
import matplotlib.pyplot as plt
import sys,os,time
import scipy.signal as signal
from DM_calculate import DM_calculate
def Signal_finding(data,angle , pixel=2, T=0, nbin=0, fy=0, DM_range=[0,1500],seq=0):
	 ang_min	= angle[0]
	 ang_max	= angle[1]
         tan_min	= np.tan(ang_min*np.pi/180)
         tan_max	= np.tan(ang_max*np.pi/180)
         tan_grid	= np.linspace(tan_min,tan_max,data.shape[1])
         deg		= np.arctan(tan_grid)/np.pi*180
	 #deg		= np.linspace(ang_min,ang_max,data.shape[1])
         data  = abs(data)
	 d_max = 0
         n_point = 4
	 edge_t	 = 0
	 lo_s = np.where(data == np.max(data))
	 if lo_s[1][0] < data.shape[1]/2:
		 std_clean = data[:,data.shape[1]/2:].std()
		 mean_clean = data[:,data.shape[1]/2:].mean()
	 else :
		 std_clean = data[:,:data.shape[1]/2].std()
                 mean_clean = data[:,:data.shape[1]/2].mean()
         for ii in range(n_point):
   	     lo    = np.where(data == np.max(data))
	     if abs(lo[1][0] -lo_s[1][0]) < pixel :
		 if (pixel  <= abs(lo[0][0] -lo_s[0][0])) or (ii ==0):
	   	 
	 #           if  0 < ii < 4:	 
	 #		 print "\n/////////"
	#	    	 print "ii:",ii,"\n"
	#	    	 print "TRUE",lo_s[1][0],lo[1][0]," seq:",seq
        #	         print "TRUE",lo_s,lo," seq:",seq
	#	    	 print "TRUE, location error:",lo[1][0] -lo_s[1][0]
	#		 print "#########\n"
		    #    exit(0)
        	    for i in np.arange(-pixel,pixel):
	             	  for j in np.arange(-pixel,pixel):
			     if 0 < lo[0][0]+i < data.shape[0] and 0 < lo[1][0]+j < data.shape[1]:
	        	         d_max += data[lo[0][0]+i,lo[1][0]+j]
			 	 data[lo[0][0]+i,lo[1][0]+j]=0
			
			     else:
				 edge_t = 1
				 break
	
	 if edge_t ==1:
		d_max =0
	#	print "\n####\nseq,d_max at edge:",seq, '\n#####'
	 
#	 d_max = d_max/1.0/((pixel*2)**2)/ii
#         snr   = (d_max - data.mean())/data.std()
	 snr	= (d_max-mean_clean)/std_clean
	 angle_2ndFFT = deg[lo_s[1][0]]
	 DM = DM_calculate(fy,abs(angle_2ndFFT),nbin,T)
	
	 if  snr > 5 and DM_range[1] > DM > DM_range[0]:
		snr	= snr
		DM	= DM
	 else:
		snr	= 0
		DM	= 0
	  
    	 return snr , DM , angle_2ndFFT

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
