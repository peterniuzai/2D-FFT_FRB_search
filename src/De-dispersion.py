import numpy as np
import matplotlib.pyplot as plt
from sigpyproc.Readers import FilReader
import time

t1	= time.time()
f_dir	= '../data/'
f_name	= 'out.fil'
f	= FilReader(f_dir + f_name )
#dm = np.linspace(1,2000,1000)
DM	= 100
snr	= []
dd	= f.dedisperse(DM)
snr	= (dd.max()-dd.mean())/dd.std()
t2	= time.time()
t_tol	= t2 -t1 
print "Time Total:",t_tol," s = ",t_tol/60. ," m"

plt.title("SNR:"+str(snr))
plt.xlabel("SRC:"+f.header['source_name'])
plt.plot(dd)
plt.show()

