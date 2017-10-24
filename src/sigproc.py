import numpy as np
import matplotlib.pyplot as plt
from sigpyproc.Readers import FilReader
import time
import sys

f_dir	= '../data/'
f_name	= 'fake_test.fil'
f	= FilReader(f_dir+f_name)
DM	= float(sys.argv[1])
t1	= time.time()
d	= f.dedisperse(DM)
t2	= time.time()
t_c	= (t2-t1)/60
d0	= f.dedisperse(0)
snr	= (d.max()-d.mean())/d.std()
lo	= np.where(d == d.max())


print 'SNR:',snr
print 'DM:',DM
print 'Location:',lo[0]
print 'Dedispersion cost:',t_c,'(min)'
print 'Data length(at DM:'+str(DM)+'): ',d.shape[0]
print 'Data length(at DM0): ',d0.shape[0]
plt.figure(figsize=[14,8])
plt.plot(d,label='DM at '+str(DM)+', SNR:'+str(int(snr)))
plt.plot(d0,label='DM at 0')
plt.plot(lo[0],d[lo[0]],'ro',label='Max loc:'+str(lo[0]))
plt.grid()
plt.legend(loc='best')
plt.show()
