import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j0,j1,itj0y0
from numpy.fft import rfft,rfftfreq,irfft

# load data
data = np.loadtxt("seis-2")
t = data[:,0]
nt = len(t)

# get frequency range
dt = t[1] - t[0]
freqs = rfftfreq(nt,dt)
f1 = 0.01
f2 = 2
nf1 = np.where(freqs>=f1)[0][0]
nf2 = np.where(freqs >=f2)[0][0]

nsta = data.shape[1] - 1
dist = 0.5 + np.arange(nsta) * 0.5
#for i in range(nsta):
#    maxd = np.max(data[:,i+1])
#    data[:,i+1] /= maxd

# change to freq domain
cross = rfft(data[:,1:],axis=0)

# compute intensity
nc = 100
phaseV = np.linspace(3,5,nc)
spectra = np.zeros((nc,nf2-nf1+1),dtype=np.complex)
for i,c in enumerate(phaseV):
    print(i)
    for j in range(nf1,nf2+1):
        w = freqs[j] * 2. * np.pi
        k = w / c
        s = 0.0 + 0.0j
        for k in range(1,nsta):
            rk = dist[k]
            rk1 = dist[k-1]
            bk = (cross[j,k] - cross[j,k-1]) / (rk - rk1)
            #ak = cross[j,k-1] - rk1 * bk

            #temp1 = ak /k* ( rk*j1(k*rk) - rk1*j1(k*rk1))
            #temp2 = bk /k**2* ( k*rk**2*j1(k*rk) - k*rk1**2*j1(k*rk1))
            #temp3 = bk / k**2 *(rk * j0(k * rk) - rk1 * j0(k*rk1))
            #temp4 = -bk / k**2 * (itj0y0(k * rk)[0] - itj0y0(k * rk1)[0])

            s = s + cross[j,k] / k * rk * j1(k * rk) + \
                bk/k**3  *(k * rk * j0(k * rk) - itj0y0(k * rk)[0])
            
            s= s - (cross[j,k-1] / k * rk1 * j1(k * rk1) + \
                bk/k**3  *(k * rk1 * j0(k * rk1) - itj0y0(k * rk1)[0]) )
        #s = temp1 + temp2 + temp3 + temp4
        spectra[i,j-nf1] = s

plt.contourf(freqs[nf1:nf2+1],phaseV,abs(spectra))
plt.show()