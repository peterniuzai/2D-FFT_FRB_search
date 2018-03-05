import numpy as np
from sigpyproc.Readers import FilReader
import sys
def DM_calculate(fy,location,nbin,T):
    '''Calculate DM from angle get'''
    C     = 4.148908e6 # (ms)
    k2    = np.tan(location/180.*np.pi)
    f_rsl = (fy[-1] - fy[0])
    unit  = f_rsl / T
    DM    = k2 / C / unit
    return DM

def length_calculate(fy,t_rsl ,DM_range,nbin,Wp,FFT_rs):
    '''Calculate signal Length in re-bin map'''

    C		= 4.148908e6 # (ms)
    f_rsl	= (fy[-1] - fy[0])/nbin
    unit	= f_rsl / t_rsl
    k1_l	= 1. / DM_range[1] / C / unit #K1 with max DM in rebin map
    deg_min	= np.degrees( np.arctan(k1_l) ) 
    L_fft	= 1./(Wp*np.sin(deg_min*np.pi/180.)) 
    #L_fft	= 1.0/Wp
    L_fft	= L_fft / FFT_rs
    return L_fft

def angle_range(fy,DM_range,nbin,T):
    '''Calculate angle range according to given DM range'''
    C    	= 4.148908e6 # (ms)
    f_rsl	= (fy[-1] - fy[0])
    unit	= f_rsl / T
    k2_b    	= DM_range[0] * C * unit #K2 with min DM, close to horizontal
    k2_t    	= DM_range[1] * C * unit #K2 with max DM  far from horizontal
    theta2_b	= np.arctan(k2_b) *180 / np.pi
    theta2_t	= np.arctan(k2_t) *180 / np.pi
    angle = [theta2_b,theta2_t]
    return angle

def time_delay(DM_range,fbot,ftop):
    '''Calculate time delay according to top frequency and bottom frequency with max DM'''
    C     = 4.148908e6 # (ms)
    t_delay   = C * DM_range[1] * (fbot**-2  -  ftop**-2) #Time delay betweent top and bottom frequency with max DM value. (ms)
    return t_delay


if __name__ == '__main__':
    ang     =  0.3357#np.float(sys.argv[1])
    DM_range=  [100,1000]#np.float(sys.argv[1])
    f_dir   = '../data/'
#    f_name  = 'fake_test.fil'
    f_name = '1024.fil'
    f = FilReader(f_dir + f_name)
    hdr    = f.header
    ftop   = hdr['ftop']
    fbot   = hdr['fbottom']
    nch    = hdr['nchans']
    t_rsl  = hdr['tsamp']*1000. # unit (ms)
    N_s_chunck = nch * 4
#    t_len  = 250/1000.
    freq   = np.linspace(ftop,fbot,nch)
#    freq   = np.linspace(1450,1050,4096)
    fy     = freq**-2
    nbin   = nch
    T	   = N_s_chunck*t_rsl
    Nsamp  = time_delay(DM_range,fbot,ftop)/t_rsl
    Nchunck  = time_delay(DM_range,fbot,ftop)/t_rsl/N_s_chunck
    dm = DM_calculate(fy,ang,nbin,T)
    degree = angle_range(fy,DM_range,nbin,t_rsl * N_s_chunck)
    print 'DM is : ',dm ,'pc*cm^-3 at',ang,'deg'
    print 'Chunks samples:',N_s_chunck
    print 'The degree is :',degree, ' degree'
    print 'Load file from:',f_name
    print 'Delay within samples :',Nsamp
    print "Delay within chuncks:",Nchunck
