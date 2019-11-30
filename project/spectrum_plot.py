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

# for U in range(5,20,5):
    # surface = water.Surface(U10=U)
    # k = surface.k0
    # spectrum = surface.spectrum(k)
    # plt.loglog(k,spectrum,label=r'$U_{10}=$'+str(U)+r' м/с')
    # plt.xlim(0.007,50)
# plt.legend()
# plt.show()

for X in range(5000,25000,5000):
    surface = water.Surface(x=X)
    k = surface.k0
    spectrum = surface.spectrum(k)
    plt.loglog(k,spectrum,label=r'$\tilde x=$'+str(X) )
plt.legend()
plt.show()
