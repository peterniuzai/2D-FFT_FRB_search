import numpy as np
import matplotlib.pyplot as plt
from sigpyproc.Readers import FilReader
import time

t1	= time.time()
f_dir	= '../data/'
f_name	= 'out.fil'
f = FilReader(f_dir + f_name )
dm = np.linspace(1,2000,1000)
snr= []
for i in range(len(dm)): 
	DM	= dm[i]
	dd	= f.dedisperse(DM)
	snr.append((dd.max()-dd.mean())/dd.std())

snr	= np.array(snr)
t2	= time.time()
t_tol	= t2 -t1 
print "Time Total:",t_tol," s = ",t_tol/60. ," m"

plt.title("SNR:"+str(snr.max()))
plt.xlabel("SRC:"+f.header['source_name'])
plt.plot(dm,snr)
plt.show()

