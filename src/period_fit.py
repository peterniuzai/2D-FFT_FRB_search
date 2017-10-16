import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.optimize import leastsq

len  = 32
time = 1000*len
raw  = np.load('/home/nch/raw_data/raw.npy')
freq = np.load('/home/nch/pulsar_data/freq.npy')
raw  = np.nan_to_num(raw)
data = raw[:,0,0:time]
data[2078:2150,:] = data[1654:1725,:] = data[800:921,:] =data[542:644,:] = data[410:424,:] =  0
d_on  = []
d_off = []
for i in range(data.shape[1]/len/2):
      for j in range(len):
          d_on.append(  data[:,j + 2*i*len ])
          d_off.append( data[:,j + (2*i+1)*len])
d_on    = np.array(d_on)
d_off   = np.array(d_off)
print d_on.dtype,'d_on'
s_on    = d_on.sum(axis=0)/time*2
print s_on.dtype,'s_on'
s_off   = d_off.sum(axis=0)/time*2
fit_on  = signal.medfilt(s_on,15)
fit_on  = np.float32(fit_on)
print fit_on.dtype,'fit_on'
fit_off = signal.medfilt(s_off,15)
fit_off = np.float32(fit_off)



def SquareWave(size = 100 , amplitude=1, length=10,counter = 5):
    
    num_steps   =  size
    s = []
    state = 1

    for n in range(num_steps):

        value = state * amplitude
        s.append( value )

        counter += 1

        if counter == length:
            counter = 0
            state *= -1

    return np.array(s)


def func(p,x):
    
    a,l,c= p
    
    y = SquareWave( len(x), a , l, c)

    return y

def residuals(p,x,obs):
    return obs - func(p,x)








a = 4000
l = 25
c = 5
#p = a ,l, c
#x = np.arange(fit.size)
#fit_coe  = leastsq(residuals,p,args=(x, fit))
#fitmap   = func(fit_coe[0],x)
#y_initial= SquareWave(len(x), a , l, c)

plt.figure(figsize = (12,12))
plt.subplot(2,2,1)
plt.title('Noise on')
plt.plot(freq,s_on,'r',label='sum')
plt.plot(freq,fit_on,'black',label = 'filtered by medium')
plt.ylabel('Intensity')
plt.xlabel('Frequency channel')
plt.legend( loc = 'upper left')
plt.grid()

plt.subplot(2,2,2)
plt.title('Noise off')
plt.plot(freq,s_off,'r',label='sum')
plt.plot(freq,fit_off,'black',label = 'filtered by medium')
plt.ylabel('Intensity')
plt.xlabel('Frequency channel')
plt.legend(loc = 'upper left')
plt.grid()

plt.subplot(2,2,3)
plt.title('Substraction')
plt.plot(freq,s_on,'r',label='sum_on')
plt.plot(freq,s_off,'black',label = 'sum_off')
sub = s_on - s_off
s_off_c = s_off / sub
s_on_c  = s_on  / sub
s_off_c = np.nan_to_num(s_off_c)
s_on_c = np.nan_to_num(s_on_c)
plt.plot(freq,s_off_c,'r',label='cali_off')
plt.plot(freq,s_on_c,label='cali_on')
plt.plot(freq,sub,'green',label = 'substraction')
plt.ylabel('Intensity')
plt.xlabel('Frequency channel')
plt.plot(freq[np.argmax(sub)],sub.max())#'ro',label = str(sub.max()))
plt.plot(freq[np.argmin(sub)],sub.min())#,'bo',label = str(sub.min()))
plt.legend(loc = 'upper left')
plt.grid()

plt.subplot(2,2,4)
plt.title('Calibration')
plt.plot(freq,fit_on,'r',label='fit_on')
plt.plot(freq,fit_off,'black',label = 'fit_off')
sub = fit_on - fit_off
s_off_c = s_off / sub
print fit_on.dtype,sub.dtype,s_off_c.dtype,'fit,sub,s_off_c'
s_on_c  = s_on  / sub 
s_off_c = np.nan_to_num(s_off_c)
s_on_c = np.nan_to_num(s_on_c)
plt.plot(freq,sub,'green',label = 'substraction')
plt.plot(freq,s_off_c,'blue',label = 'cal_off')
plt.plot(freq,s_on_c,'r',label='cali_on')
plt.ylabel('Intensity')
plt.xlabel('Frequency channel')
plt.plot(freq[np.argmax(s_off)],s_off.max(),'ro')
plt.plot(freq[np.argmax(s_off_c)],s_off_c.max(),'bo')
plt.legend(loc = 'upper left')
plt.grid()


np.save('/home/nch/FFT_search/src/chan_equaliz',sub)
print sub.shape,sub.dtype,'sub'
plt.savefig('/home/nch/FFT_search/graph/Noise_cal_raw')
plt.show()
plt.close()

plt.plot(sub)
plt.show()

