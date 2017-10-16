import numpy as np
from sigpyproc.Readers import FilReader
import sys
def DM_calculate(fy,location,nbin,T):
    d     = 4.148908e6 # (ms)
#    loc   = np.float32(90 - location)
    k2    = np.tan(location/180.*np.pi)
#    k1    = np.tan(loc/180*np.pi)
    f_rsl = (fy[-1] - fy[0])
    unit  = f_rsl / T
    DM    = k2 /d / unit
    return DM

def degree_calculate(fy,t_rsl ,DM,nbin):
    d      = 4.148908e6 # (ms)
#    d      = d / 1000

    f_rsl  = (fy[-1] - fy[0])/nbin
    unit   = f_rsl / t_rsl
    k      = 1. / DM / d
    ang    = k / unit
    degree = np.degrees( np.arctan(ang) )     
    degree = 90 - degree
    return degree

if __name__ == '__main__':
    ang     =  52#np.float(sys.argv[1])
    DM      =  800#np.float(sys.argv[1])
    f_dir   = '../data/'
    f_name  = 'fake_test.fil'
    f = FilReader(f_dir + f_name)
    hdr    = f.header
    ftop   = hdr['ftop']
    fbot   = hdr['fbottom']
    nch    = hdr['nchans']
    t_rsl  = hdr['tsamp']*1000. # unit (ms)
    freq   = np.linspace(ftop,fbot,nch)
    fy     = freq**-2
    nbin   = nch
    T	   = nch*t_rsl
    dm = DM_calculate(fy,ang,nbin,T)
    degree = degree_calculate(fy,t_rsl,DM,nbin)
    print 'the DM is : ',dm ,'pc*cm^-3'
    print 'The degree is :',degree, ' degree'
