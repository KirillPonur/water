import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import root,minimize
from scipy.misc import derivative
a = 0.1
k = 2*np.pi


x0_ = np.linspace(0,2,200)

def x0(x):
    f = lambda x0_, x_: x0_ - a*np.sin(k*x0_) - x_
    x0 = root(f, x0 = x, args=(x,) ).x 
    return x0

def z(x_):
    return a*np.cos(k*x0(x_))


fig, ax = plt.subplots(1,3)
y1 = [z(x0_[i]) for i in range(x0_.size)]
print('Среднее \t Реализация '+ str(np.mean(y1)) + '\t Теория '+str(-a**2/2*k))
print('STD', np.std(y1))
y1 = np.array(y1).T[0]
y2 = a*np.cos(k*x0_)
yy1  =derivative( z, x0 = x0_ , dx = 0.01)  
yy2  =derivative( lambda x: a*k*np.cos(k*x), x0 = x0_, dx = 0.01)  

ax[0].plot(x0_,y1)
ax[0].plot(x0_,y2)

f = -a*k*np.sin(k*x0_)/(1 - a*k*np.cos(k*x0_))
ax[1].plot(x0_,yy1)
ax[1].plot(x0_,yy2)
ax[1].plot(x0_,f)

ax[2].plot(y1,yy1)
ax[2].plot(y2,yy2)
ax[2].plot(y2,f)

# def covariation1(x0_):
#     covariation = np.zeros(x0_.size)
#     # x = x0(np.linspace(0,1,200))
#     x = np.linspace(0,1,200)
#     m = np.mean(z(x))
#     for i in range(x0_.size):
#         covariation[i] = ( 
#                     np.trapz( z(x)*z(x+x0_[i]) \
#                     # * (1 + a*k*np.cos(k*x)) 
#                     )     
#                          )
#     return covariation

# def covariation2(x0):
#     covariation = np.zeros(x0.size)
#     z = lambda x: a*np.cos(k*x)
#     x = np.linspace(0,1,200)
#     for i in range(x0.size):
#         covariation[i] = np.trapz( z(x)*z(x+x0[i]) )  
#     return covariation
# # N = len(y1)
# # z1 = np.fft.rfft(y1, n = N)
# # z2 = np.fft.rfft(y2, n = N)
# # freq = np.fft.rfftfreq(n = N)
# fig, ax = plt.subplots(1,1)

# C1 = covariation1(x0_)
# ax.plot(x0_,C1)

# C1 = covariation2(x0_)
# ax.plot(x0_,C1)
# # print(C.shape)

plt.show()