from water import *
import os
N=1024
water = water(N=N, M=2,)
surface = water.model(water.k, water.phi, water.t)
x, y  = np.meshgrid(water.x0,water.y0)
z = surface([x,y],water.t)
plt.plot(water.x0, z[0])

s = water.spectrum(water.k)
a = water.amplitude(water.k)

print(a)
with open('surface'+str(N-1)+'.txt','w') as f:
    f.write(str(water.k[0]) + '\t' + str(a[0])+ '\t')
    f.write('\n')
    for i in range(z[0].size):
        f.write(str(water.x0[i]) + '\t' + str(z[0][i]) +'\n')
f.close()

plt.show()