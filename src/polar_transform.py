import scipy.signal as signal
from scipy import ndimage
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")


def polar_coordinates_convert(data, ang_min = 5 ,ang_max = 85 ):
    #Transform the rectanular coordinates into polar coordinates
    #Take histogram2D method
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

def polar_coordinates_convert_inter(data, ang_min = 5 ,ang_max = 85 ):
	#Transform the rectanular coordinates into polar coordinates
	#Take interpolation method

	  rang     = data.shape
	  
	  # Make grid for interpolation
	  ang_rsl  = (1./data.shape[1]/2**0.5)*180/np.pi
	  rad_rsl  = 1
	  rad_grid = np.arange(5,rang[1],rad_rsl)
	  ang_grid = np.arange(ang_min,ang_max,ang_rsl)

	  grid_a,grid_r = np.meshgrid(ang_grid,rad_grid)
	  x_p      = grid_r * np.cos(grid_a*np.pi/180)
	  y_p      = rang[0]-grid_r * np.sin(grid_a*np.pi/180)-1
	  x_p     = x_p.reshape(-1)
	  y_p     = y_p.reshape(-1)
	  cord    = [y_p,x_p]
	  polar_matrix_r  = ndimage.map_coordinates(np.real(data),cord,order=1)
	  polar_matrix_i  = ndimage.map_coordinates(np.imag(data),cord,order=1)
	  polar_data    = polar_matrix_r+polar_matrix_i*1j
 	  polar_data    = polar_data.reshape(grid_r.shape)
	  return polar_data,ang_rsl,rad_rsl


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


