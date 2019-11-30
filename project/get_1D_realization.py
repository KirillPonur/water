import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm
import water

N = 2049
M = 1

surface = water.Surface(N=N,M=M)
x0 = np.linspace(0,200,10**5)
y0 = 0
t = 0
x, y = np.meshgrid(x0, y0)
z = surface.model([x,y],t)[0]
fig,ax = plt.subplots(nrows = 1, ncols = 1)
ax.plot(x0,z)




with open('surface1D.txt','w') as f:
    f.write('##### One dimensional realization ##### \n')
    f.write('U10 \t N \t M \t k_begin \t k_end \n')
    f.write(str(surface.U10)+'\t'+ str(N-1) + '\t' + str(M) + '\t' 
    + str(round(surface.KT[0],3)) +'\t'
    + str(round(surface.KT[1],1)) + '\n')
    f.write('####################################### \n')
    for i in range(z.size):
        f.write(str(x0[i]) + '\t' + str(z[i]) +'\n')
f.close()

plt.show()