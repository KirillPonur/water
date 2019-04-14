from pylab import *
from numpy import random
from scipy import interpolate
from scipy  import integrate
import os
 
rcParams['figure.figsize'] = [7, 6]
rcParams['axes.labelsize'] = 20
import pandas as pd
import time

def find_decision(omega):
    P = 9.8*1000.0/0.074;
    Q = -1000.0*omega**2/0.074;
    x1= -Q/2.0 + sqrt( (Q/2)**2 + (P/3)**3 )
    x2= -Q/2.0 - sqrt( (Q/2)**2 + (P/3)**3 )
    k=x1**(1/3)-(-x2)**(1/3)
    return k

#    Функция возвращает детерминант при переходе от частоты к 
#волновым числам по полному дисперсионному
def det(k): 
    det=(9.8+3*k**2*0.074/1000)/(2*sqrt(9.8*k+k**3*0.074/1000) )
    return det

global g
g=9.81
def k_max(omega_max):
    k_max=omega_max**2/g
    return k_max
    
# Пересчет волнового числа в частоту по полному дисперсионному
def omega_k(k): 
    omega_k=(g*k+0.074*k**3/1000)**(1/2) 
    return omega_k  

def full_spectrum(k,x=20170, long_calculate=True):
    def JONSWAP(k):
        if k<=k_m:
            sigma=0.07
        else:
            sigma=0.09
        Sw=(
            Alpha(x)/2*k**(-3)*exp(-1.25*(k_m/k)**2 )*
            Gamma(x)**(exp(- ( sqrt(k/k_m)-1)**2 / (2*sigma**2) ))
           )
        return Sw

    def Gamma(x):
        if x>=20170:
            return 1
        gamma=(
               +5.253660929
               +0.000107622*x
               -0.03778776*sqrt(x)
               -162.9834653/sqrt(x)
               +253251.456472*x**(-3/2)
              )
        return gamma

    def Alpha(x):
        if x>=20170:
            return 0.0081
        alpha=array([],dtype='float64')
        alpha=[( 
               +0.0311937
               -0.00232774*log(x)
               -8367.8678786/x**2
               +4.5114599e+300*exp(-x)*1e+300*1e+17
    #            +4.5114599e+17*exp(-x)
              )]
        return alpha[0]

    def Omega(x): #Вычисление безразмерной частоты по безразмерному разгону
        if x>=20170:
            return 0.835
        omega_tilde=(0.61826357843576103 
                     + 3.52883010586243843e-06*x
                     - 0.00197508032233982112*sqrt(x)
                     + 62.5540113059129759/sqrt(x)
                     - 290.214120684236224/x
        )
        return omega_tilde

    def spectrum1(k):


        omega0=omega_k(limit_k[0])
        beta0= JONSWAP(limit_k[0])*omega0**4/det(limit_k[0])

        omega0=omega_k(k)

        return beta0/omega0**4*det(k)

    def spectrum2(k):


        omega0=omega_k(limit_k[1])
        beta0= spectrum1(limit_k[1])*omega0**5/det(limit_k[1])

        omega0=omega_k(k)

        return beta0/omega0**5*det(k)

    def spectrum3(k):


        omega0=omega_k(limit_k[2])
        beta0= spectrum2(limit_k[2])*omega0**2.7/det(limit_k[2])

        omega0=omega_k(k)

        return beta0/omega0**2.7*det(k)

    def spectrum4(k):


        omega0=omega_k(limit_k[3])
        beta0= spectrum3(limit_k[3])*omega0**5/det(limit_k[3])

        omega0=omega_k(k)

        return beta0*det(k)/omega0**5
  
    
    gamma,alpha,omega_m=Gamma(x),Alpha(x),Omega(x)
    
    try:
        full_spectrum=zeros(len(k))
    except:
        full_spectrum=[0]
        k=[k]
    if long_calculate==False:
        k=logspace(log10(k[0]),log10(k[-1]),1000)
        
    omega=omega_m*g/U10
    global k_m
    k_m=k_max(omega)
    limit_1= 1.2 
    limit_2=(
             +0.371347584096022408 
             + 0.290241610467870486*U10
             + 0.290178032985796564/U10
            )
    limit_3= 270.0 
    limit_4= 1020.0 

    limit_k=np.zeros(4)
    limit_k[0]=find_decision(limit_1*omega)
    limit_k[1]=find_decision(limit_2*omega)
    limit_k[2]=limit_3
    limit_k[3]=limit_4

    for i in range(len(k)):
        if k[i] <= limit_k[0]:
            full_spectrum[i] =  JONSWAP(k[i])
        elif k[i] <= limit_k[1]:
            full_spectrum[i] = spectrum1(k[i])
        elif k[i] <= limit_k[2]:
            full_spectrum[i] = spectrum2(k[i])
        elif k[i] <= limit_k[3]:
            full_spectrum[i] = spectrum3(k[i])
        else:
            full_spectrum[i] = spectrum4(k[i])

    return full_spectrum

def B(k):
    def b(k):
        b=(
            -0.28+0.65*exp(-0.75*log(k/k_m))
            +0.01*exp(-0.2+0.7*log10(k/k_m))  
          )          
        return b
    B=10**b(k)
    return B

def Normalization(B):
    Normalization=B/arctan(sinh(2*pi*B))
    return Normalization

def Phi(k,phi):
    try:
        Phi=zeros((len(k),len(phi)))
        for i in range(len(k)):
            B0=B(k[i])
            A0=Normalization(B0)
            Phi[i]=A0/cosh(2*B0*phi)
    except:
        B0=B(k)
        A0=Normalization(B0)
        Phi=A0/cosh(2*B0*phi)
    return Phi

def angle(k,phi):
    angle=sqrt( 2*pi/100 * Phi(k,phi) )
    return angle

def amplitude(k):
    # k-- выбранный диапазон волновых чисел
    # N-- количество моделируемых гармоник
    N=len(k)
    S=full_spectrum(k)
    dk=zeros(N)
#     dk[0]=k[0]
    dS=zeros(N)
    if len(k)==1:
        dk[0]=k[0]
        dS[0]=S[0]
    for i in range(1,N):
        dk[i]=(k[i]-k[i-1])
        dS[i]=S[i]
    amplitude=sqrt(2*dk * dS)
    return amplitude



def model(k,phi,N):
    U10=3
    rho=linspace(0,400,400)
    A=amplitude(k)
    F=angle(k,phi)
    psi=array([
    [ random.uniform(0,2*pi) for j in range(100)] 
      for i in range(N)              ]) 
    
    def water(r,phi,t=0):
        model=0
        for n in range(N):
            for m in range(100):
                model+=A[n]*F[n][m]* cos( k[n]*(r[0]*cos(phi[m])+r[1]*sin(phi[m]))+psi[n][m] )
        return model
    
    return water


def interspace(k,N):
    y=lambda k: full_spectrum(k)
    sigma=trapz(y(k),k)
    b0=sigma/(N)
    k_new=k[0]
    a=k[-1]
    k=zeros(N+1)
    k[0]=k_new
    I=zeros(N)
    err=zeros(N)
    epsabs=1.49e-12
    epsrel=1.49e-12
    for i in range(N-1):
        integral=0
        n=1
        m=2
        while integral<b0:
            k_new*=1+10**(-m)
            if k_new>a:
                k_new=a
                break

            integral,error=integrate.quad(y,k[i],k_new,
            limit=50,epsabs=epsabs, epsrel=epsrel)
            if integral>b0 and abs(integral-b0)>epsabs:
                k_new*=1-10**(-m)
                m+=1
                integral,error=integrate.quad(y,k[i],k_new,
                limit=50,epsabs=epsabs, epsrel=epsrel)
            n+=1
        k[i+1]=k_new
    k=k[0:N+1]
    k[-1]=a
    return k,b0,err

def nodes(ki,b0):
    y=lambda k: k**2*full_spectrum(k)
    nodes=zeros(len(ki))
    epsabs=1.49e-12
    epsrel=1.49e-12
    for i in range(1,len(ki)):
        integral,error=integrate.quad(y,ki[i-1],ki[i],
        limit=50,epsabs=epsabs, epsrel=epsrel)
        B=(sqrt(integral/b0))
        nodes[i-1]=B
#         print(B)
    return nodes[:-1]

if __name__ == "__main__":

	U10=6
	x=20170
	KT=[0.02,1999.99]


	k=arange(KT[0]/10,KT[-1]+1,0.01)
	S=full_spectrum(k)
	full_spectrum=interpolate.interp1d(k,S)

	rho=linspace(0,400,400)
	# N=find_garmonics(k,rho)
	N=200
	k=logspace(log10(KT[0]),log10(KT[-1]),10**2)
	k,b0,err=interspace(k,N)
	k=nodes(k,b0)
	phi=linspace(-pi,pi,100)

	sigma=model(k,phi,N)

	x=linspace(0,400,400)
	y=linspace(0,400,400)
	rcParams['figure.figsize'] = [7, 6]
	rcParams['axes.labelsize'] = 20
	x, y = np.meshgrid(x, y)


	z=sigma([x,y],phi)
	contourf(z,100,cmap=cm.winter)
	colorbar()
	ylabel(r'Y',fontsize=16)
	xlabel(r'X',fontsize=16)
	savefig('water'+str(U10)+'.png',pdi=2500)
	show()
