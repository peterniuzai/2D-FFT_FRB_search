import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys,os,time
#import readGBT
#import pyfits
#import h5py
import mpi4py.MPI as MPI
#import scipy.signal as signal
#from scipy.interpolate import griddata
from readfile import read_data
from dir_create import dir_create
#from calibrated import calibration
from rebin import rebin, rebin_inter
from FFT import FFT
from polar_transform import polar_coordinates_convert, polar_coordinates_convert_inter
from Signal_finding import Signal_finding
from DM_calculate import DM_calculate
from plot_all import plot
import warnings
warnings.filterwarnings("ignore")


if __name__ == '__main__':
     # instance for invoking MPI relatedfunctions
     comm = MPI.COMM_WORLD
     # the node rank in the whole community
     comm_rank = comm.Get_rank()
     # the size of the whole community, i.e.,the total number of working nodes in the MPI cluster
     comm_size = comm.Get_size()


     if ('-h' in sys.argv) and comm_rank ==0:
                print "USAGE :mpirun -n <number of processor> python __main__.py "
                quit = 1
     else:
                quit = 0
     quit = comm.bcast(quit,root=0)
     comm.barrier()
     if quit == 1:
                exit()


############
#Arguments #
############


#     f_name	= '2011-02-20-01:52:19.fil'
#     f_dir	= '/data0/FRB_Parkes_data/FRB110220/'    
#     f_name    = 'data_2017-08-30_17-35-36.fil'
     f_name    = 'fake_test.fil'
#     f_name    = 'BJ0009_02551.fil'
#     f_name	= 'FRB110626.fil'
#     f_name	= 'FRB010621.fil'
#     f_name	= 'FRB110220.fil'
#     f_name	= 'PM0141_017A1.fil'
#     f_name    = '1.fil'
     f_dir     = '../data/'
     plot_dir  = '../graph/' + f_name[:-4] + '/'
     plot_proc = '2ndFFT_3D,polar_sets_3D,1stFFT,raw,rebin,polar_sets_2D,2ndFFT_2D'
     # Plot_proc: list  processes we  want to plot.

     if comm_rank == 0:
                dir_create(plot_dir)
                print 'directory build complete!'
     comm.barrier()

     t_len     = 0	#time length for each smallest unit to process.
     DM_range  = [100,5000]	#Min and Max DM
     Wp	       = 4		#Wp means pulse width in (ms)
     nbin      = 0
     ang     = [0,0] #Angle range for search
     msk_cycle = 5	#the number of channels to be zeros in 2D-FFT(Noise remove).
     pixel     = 2	#the number of pixel to sum in 2ndFFT3D SNR compute.
     SNR_l     = []
     DM_l      = []   
 
######################
#Read Data from file #
######################

     time_1    = time.time()

     if comm_rank == 0:	 print 'Begin to load data from ' + f_name 
     fil, num, p_n, freq, t_rsl, t_len, nbin ,nch, T,fy,angle,Ang_rs,Rad_rs,L_fft = read_data(f_dir, f_name ,t_len, nbin, comm_size,DM_range,Wp,ang)
     
###################
#Begin to search  #
###################
     time_s = time.time()
     for  i_ch in range(p_n):  #i_chunk 
     #for  i_ch in range(1):  #i_chunk
	     t_p    = comm_rank*p_n   #the thread position in total time in unit(chunk)
	     data   = fil.readBlock(t_len*(i_ch+t_p),t_len)
#	     data[:220,:]=0
	     data   = np.nan_to_num(data)
	     t_ch_s = t_len*(i_ch+t_p)*t_rsl   #time of chunck start.
	     t_ch_e = t_len*(i_ch+t_p+1)*t_rsl 
	     t_axis = np.linspace(t_ch_s,t_ch_e,t_len) 
	     
	     if comm_rank == 0:    print 'Load Data Over, Datasize for each chunk:',data.shape,'\nBegin to rebin... '

	     re_data, f_axis  =  rebin_inter(data, fy, nbin,t_axis)
	
	     if comm_rank == 0:    print 'Rebin over. \nBegin to do 1st 2-D FFT on rebin data...'
	
	     FFT1st_data = FFT(L_fft, re_data, 2, msk_cycle=0)
	
	     if comm_rank == 0:    print '1st FFT over.\nBegin to transform rectangular coordinates into polar coordinates...'
	
	     polar_data  = polar_coordinates_convert_inter( FFT1st_data, angle, Ang_rs)
	
	     if comm_rank == 0:    print 'Polar transform over.\nBegin to do the 2nd 1-D FFT along radius direction...'
	
	     FFT2nd_data = FFT(L_fft, polar_data, 1 )# 1 means 1 Dimension FFT
	
	     if comm_rank == 0:    print '2nd FFT over.\nBegin to locate the signal and calculate SNR...'
	
	     SNR , DM, A_f = Signal_finding(FFT2nd_data,  angle , pixel, T, nbin, fy, DM_range)
	     SNR_l.append(SNR)
	     DM_l.append(DM)
	     if comm_rank == 0:    print 'Searching Over. '
	     time_2  =  time.time()
	     consume =  time_2 - time_1
	
	
	     if comm_rank == 0:
                         print '\n#############'
                         print 'Process matrix size:',nch,' * '+str(t_len)
                         print 'Time cost:', consume ,'seconds,  ','equal',consume/60.,'minutes.'
                         print'Process:',i_ch,' of ',p_n,' for total:',p_n*t_len*comm_size,'samples'
			 print 'Angle range:',angle
                         print 'SNR:',SNR,';DM: ',DM
                         print '###############\n\nBegin to plot...'
#	     plot(comm_rank,t_axis,data,re_data,polar_data,FFT1st_data,FFT2nd_data,plot_proc,freq,f_axis,Rad_rs,Ang_rs,plot_dir,pixel,angle,i_ch,p_n,SNR,DM,A_f)
             if comm_rank == 0: 	print 'Plot Over...\n\n'	
#########################################
# gather the results from all processes #
#########################################
     SNR_l = np.array(SNR_l)
     DM_l  = np.array(DM_l)
     comm.barrier()
     time_e= time.time()
     if comm_rank ==0:   print 'Total consume:',(time_e-time_s)/60.,'mins'
     combine_SNR      = comm.gather(SNR_l,root=0)
     combine_SNR      = np.array(combine_SNR).reshape(-1)
     combine_DM = comm.gather(DM_l,root=0)
     combine_DM = np.array(combine_DM).reshape(-1)
     lo	= np.where(combine_SNR == combine_SNR.max() )
     DM = combine_DM[lo[0]]
     SNR= combine_SNR.max()
    
     if comm_rank == 0:   
		print '*********'
		print 'Found FRB with DM:',DM[0],'wit SNR:',SNR,' in ',lo[0],'file.'
		print '********'
                print 'multiprocess plot over....'
	#	exit()
                N_cut = len(combine_SNR)
                plt.ylabel('SNR')
                plt.xlabel(str(N_cut)+'files')
                plt.title('SNR of each '+str(N_cut)+' files')
                plt.plot(np.arange(N_cut),combine_SNR,'ro',label='SNR of 2nd FFT in 2-D map')
                plt.grid()
                plt.savefig(plot_dir+'SNR')
  #              plt.show()
                plt.close()

                N_cut = len(combine_DM)
                plt.ylabel('Result DM')
                plt.xlabel(str(N_cut)+'files')
                plt.title('Found DM:'+str(DM)+'with SNR:'+str(SNR))
                plt.plot(np.arange(N_cut),combine_DM,label='Location of 2nd FFT in 2-D map')
		plt.plot(lo[0],DM,'ro')
                plt.grid()
                plt.savefig(plot_dir + 'DM')
 #               plt.show()
                plt.close()

       	        print 'plot over.... rank:',comm_rank
