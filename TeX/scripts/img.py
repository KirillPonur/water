from pylab import *
from functions import parsing
import os.path as path
from scipy.signal import medfilt
from scipy import interpolate 
img = path.abspath('..'+'\\img\\DSC_0021 (1).txt')
x,y=parsing(img,0,1)
y=medfilt(y,3)
plot(x,y)
show()