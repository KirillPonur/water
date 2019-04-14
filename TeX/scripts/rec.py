from pylab import *
from matplotlib import rc
import os.path as path
import sys
from scipy import interpolate
from pandas import read_excel as read
rc('text', usetex=True)
rc('text.latex', preamble=[r'\usepackage[russian]{babel}',
                         r'\usepackage{amsmath}',
                        r'\usepackage{amssymb}'])


rc('font', family='serif')

rec = path.abspath('..'+'\\rec\\rec.xlsx')

df=read(rec, sheet_name='1')

figure('Задание 1')
z=df['Цена деления, В']
x=array(df['f'])
y=array(df['U'])/max(df['U'])
g = interpolate.interp1d(x,y, 'quadratic' )
x=linspace(x[0],x[-1],1000)
y=g(x)
plot(x,y,label='интерполяция')
plot(df['f'],df['U']/max(df['U']),'ro',label='эксперимент')
xlabel(r'$f$, \text{кГц}',fontsize=16)
ylabel(r'$K(f)=\frac{U_{\text{вых}}}{U^{m}_{\text{вых}}}$, ',fontsize=16)
grid(which='major', linestyle='-')
grid(which='minor', linestyle=':')
minorticks_on()
legend()
savefig(path.abspath('..'+'\\fig\\task1.pdf'))


show()