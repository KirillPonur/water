import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm
import water
from matplotlib import rc
plt.rcParams['axes.labelsize'] = 20
rc('text',usetex=True)
rc('text.latex',preamble=[r'\usepackage[russian]{babel}',r'\usepackage{amsmath}'])
rc('font',family = 'serif')

N = 256
M = 100

t = 0


def plot_surface(x , y, t):
    fig,ax = plt.subplots(nrows = 1, ncols = 1)
    surface = water.Surface(N=N,M=M,U10=5,wind= np.pi/6)
    x, y = np.meshgrid(x, y)
    z = surface.model([x,y],t)
    print(z.shape)
    from matplotlib.cm import winter
    plt.contourf(z,100,cmap=winter)
    plt.colorbar()
    plt.ylabel('Y, м',fontsize=16)
    plt.xlabel('X, м',fontsize=16)
    plt.savefig('/home/kannab/documents/water/poster/fig/water5.png',  pdi=10**6,transparent=True)

x0 = np.linspace(0,200,200)
y0 = x0
plot_surface(x0,y0,t)
plt.show()
