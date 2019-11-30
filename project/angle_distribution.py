import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm
import water
from matplotlib import rc


surface = water.Surface(wind=np.pi/2)


# rcParams['figure.figzise'] = [8,8]
plt.rcParams['axes.labelsize'] = 20
rc('text',usetex=False)
rc('text.latex',preamble=[r'\usepackage[russian]{babel}'])
rc('font',family = 'serif')

k_m = surface.k_m
x = np.linspace(-np.pi,np.pi,1000)
temp = [k_m/2,k_m,2*k_m]

plt.figure()
for i in temp:
    y = surface.Phi(i,x)
    plt.polar(x,y/k_m,label=r'$k/k_m$='+str(round(i/k_m,3)))
plt.legend()
temp = [10*k_m,50*k_m,100*k_m]
plt.figure()
for i in temp:
    y = surface.Phi(i,x)
    plt.polar(x,y/k_m,label=r'$k/k_m$='+str(round(i/k_m,3)))
plt.legend()
plt.show()
