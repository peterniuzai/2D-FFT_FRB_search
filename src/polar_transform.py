import scipy.signal as signal
from scipy import ndimage
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")


def polar_coordinates_convert(data, angle ):
    #Transform the rectanular coordinates into polar coordinates
    #Take histogram2D method
	  print angle
	  ang_min = angle[0]
	  ang_max = angle[1]
          rang    = data.shape
          ang_rsl = (1./rang[1])*180/np.pi
#	  ang_rsl =  
	  rad_rsl =  1
	  
	  #take the left top of the matrix as center
          row     = np.arange(-rang[0]+1,1,dtype=np.float32) 
          rank    = np.arange(rang[1],dtype=np.float32)

	  # calculate the angle of each pixe
          angle   = np.nan_to_num(np.arctan(row[:,None]/rank[None,:])) / np.pi * 180 

	  # calculate radius of each pixel
          radius  = np.sqrt(row[:,None]**2 + rank[None,:]**2)
          ang     = -angle.reshape(-1)
          rad     = radius.reshape(-1)
          data    = data.reshape(-1)
	  rad_bin = np.arange(1,rang[1],rad_rsl)
	  ang_bin = np.arange(ang_min,ang_max,ang_rsl)
	  data_r,bin_eage_x,bin_eage_y   = np.histogram2d(rad, ang, weights=np.real(data), bins= (rad_bin,ang_bin))
	  data_i, bin_eage_x,bin_eage_y  = np.histogram2d(rad, ang, weights=np.imag(data), bins= (rad_bin,ang_bin))
	  polar_data   = data_r + data_i * 1j

          return  polar_data,ang_rsl,rad_rsl

def polar_coordinates_convert_inter(data, angle,Ang_rs):
	#Transform the rectanular coordinates into polar coordinates
	#Take interpolation method

	  rang     = data.shape
	  ang_min  = angle[0]
	  ang_max  = angle[1]
	  tan_min  = np.tan(ang_min*np.pi/180)
	  tan_max  = np.tan(ang_max*np.pi/180)
	  
	  # Make grid for interpolation
	#  ang_rsl  = Ang_rs
	  #ang_rsl  = (1.0/rang[1])/np.pi*180.
#	  rad_rsl  = 1# 2**0.5
	  
#	  rad_grid = np.arange(4,rang[0],rad_rsl)
#	  n_rad	   = 512*2#round((rang[1]-4.0)/rad_rsl)
#	  rad_grid = np.linspace(0,rang[0],n_rad)
          #n_rad    = rang[0]#int(np.sqrt(data.size))
	  #n_rad	   = 2**(round(np.log2(n_rad)))
#	  n_deg    = round(np.sqrt(data.size))
#	  n_deg	   = 2**(round(np.log2(n_deg)))
#	  n_deg    = 1024
#	  n_deg	   = int(np.sqrt(data.size))#(ang_max-ang_min)*20
#	  n_deg	   = round(data.size/1.0/len(rad_grid))
#	  n_deg    = #len(rad_grid)
	  n_deg	  = (ang_max-ang_min)*20
	  n_rad   = int(data.size/n_deg)
	  n_rad   = 2**(round(np.log2(n_rad)))
	  n_rad	  = 512*2
	  rad_grid = np.linspace(1,rang[0]-1,n_rad)
#	  print "n_rad",n_rad
	  #exit(1)
#	  ang_grid = np.linspace(ang_min,ang_max,n_deg)
          tan_min  = np.tan(ang_min*np.pi/180)
          tan_max  = np.tan(ang_max*np.pi/180)
	  tan_grid  = np.linspace(tan_min,tan_max,n_deg)
	  ang_grid  = np.arctan(tan_grid)/np.pi*180
	  #print "ang_gird:",ang_grid.min(),ang_grid.max()
	  #plt.plot(ang_grid)
	  #plt.grid()
	  #plt.show()
	  #exit(1)
	  grid_a,grid_r = np.meshgrid(ang_grid,rad_grid)
	  x_p     = grid_r * np.cos(grid_a*np.pi/180.) -1
	  y_p     = rang[0] - grid_r * np.sin(grid_a*np.pi/180.) 
	  x_p     = x_p.reshape(-1)
	  y_p     = y_p.reshape(-1)
	  cord    = [y_p,x_p]
	  polar_matrix_r  = ndimage.map_coordinates(np.real(data),cord,order=0)
	  polar_matrix_i  = ndimage.map_coordinates(np.imag(data),cord,order=0)
	  polar_data    = polar_matrix_r+polar_matrix_i*1j
 	  polar_data    = polar_data.reshape(grid_r.shape)
#	  polar_data[:4,:]=0
	  print polar_data.shape


	  return polar_data


if __name__ == '__main__':

        rad_grid = 2400
        ang_grid = 600
        data  = np.load('/home/ycli/data/burst_data/filtered_short.npy')
        data  = np.load('/home/nch/plot_raw/wigglez_found/data_pick/0_0.npy')
        if   data.shape[1] == 4:
             data = [data[:,0,:]]
        else:
             data = [data]
        polar_sets  = polar_coordinates_convert( data,rad_grid ,ang_grid )
    
        print polar_sets[0].shape 


