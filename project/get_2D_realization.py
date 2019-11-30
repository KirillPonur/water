import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm
import water

N = 512
M = 256

surface = water.Surface(N=N,M=M)
x0 = np.linspace(0,400,200)
y0 = np.linspace(0,400,200)
t = 0
x, y = np.meshgrid(x0, y0)
z = surface.model([x,y],t)
fig,ax = plt.subplots(nrows = 1, ncols = 1)





with open('surface2D.txt','w') as f:
    f.write('##### Two dimensional realization ##### \n')
    f.write('U10 \t N \t M \t k_begin \t k_end \t x,y_begin \t x,y_end \n')
    f.write(str(surface.U10)+'\t'+ str(N) + '\t' + str(M) + '\t' 
    + str(round(surface.KT[0],3)) +'\t'
    + str(round(surface.KT[1],1)) + '\t'+ 
    '0 \t 200 \t' + '\n')
    f.write('####################################### \n')
    for i in range(z.shape[0]):
        f.write('\n')
        for j in range(z.shape[1]):
            f.write(str(z[i][j]) +'\t')
f.close()

from matplotlib.cm import winter
plt.contourf(z,100,cmap=winter)
plt.colorbar()
plt.ylabel('Y, м',fontsize=16)
plt.xlabel('X, м',fontsize=16)
plt.show()
