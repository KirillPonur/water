import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tqdm

x = np.linspace(0,10,1000)
k = 2*np.pi
a = 0.1
z = a * np.cos(k*x)

N = z.size
hat_z = np.fft.fft(z, n = N)
freq = np.fft.fftfreq(n = N, d = 1/10)
D = np.fft.ifft(1j*np.sign(freq)*hat_z)
# D = - a * np.sin(k*x)
# plt.plot(freq,np.abs(hat_z)/z.size)
plt.plot(x+D,z)
# plt.plot(x,D)

plt.show()
