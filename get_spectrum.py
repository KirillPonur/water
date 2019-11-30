import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

from water import *


df = pd.read_csv('surface1.txt', sep="\t", header = 0 ).T



def fft(t, x, ftype = None):
    if ftype == 'real':
        N = x.size
        m = 1
        X = np.fft.rfft(x, n = m*N)
        freq = np.fft.rfftfreq(n= m*N, d = (t[-1]-t[0])/t.size)
    else:
        N = x.size
        freq = np.linspace(-0.1,0.1,10**3)
        # freq = np.logspace(np.log10(0.06), np.log10(150), 100)
        X = np.zeros(freq.size, dtype='complex')
        progress_bar = tqdm(total = freq.size)
        for j in range(freq.size):
            for i in range(t.size):
                X[j] += x[i]*np.exp(2j*np.pi*freq[j]*t[i])  
            progress_bar.update(1)
        progress_bar.close()

    return freq,X


fig,ax = plt.subplots(nrows = 1, ncols = 2)

t = np.linspace(0,2000,10**5)
psi = np.random.uniform(0,2*pi)
a = 0.20373118
k= 0.0683977725
heights = a*np.sin(k*t+psi) 

ax[0].plot(t,heights)
freq,S = fft(t,heights,'real')
ax[1].plot(freq,np.abs(S)/t.size)


t = np.array(df.iloc[0])[::100]
heights = np.array(df.iloc[1])[::100]

ax[0].plot(t,heights)
freq,S = fft(t,heights,'real')
ax[1].plot(freq,np.abs(S)/t.size)


plt.show()