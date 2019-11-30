import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm
import water
from matplotlib import rc
plt.rcParams['axes.labelsize'] = 20
rc('text',usetex=True)
rc('text.latex',preamble=[r'\usepackage[russian]{babel}'])
rc('font',family = 'serif')

N = 1000
M = 1
surface = water.Surface(N=N,M=M)
print(surface.k_m)
print(surface.sigma_sqr)
x0 = np.linspace(0,50,10**5)
y0 = 0
T = [0]

fig,ax = plt.subplots(nrows = 1, ncols = 1)
x, y = np.meshgrid(x0, y0)
for t in T:
    z_real = surface.model([x,y],t)[0]
    Dx,Dy = surface.D([x,y],t)
    ax.plot(x0+Dx[0],z_real,label='CWM',color='darkblue')
    ax.plot(x0,z_real,'--r',label='Стандартный метод')
    ax.set_xlabel(r'X, м')
    ax.set_ylabel(r'Z, м')
# # ax.plot(ans)
plt.legend()
plt.show()

