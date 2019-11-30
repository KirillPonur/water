
import numpy as np
from scipy import interpolate, integrate
from numpy import pi
class spectrum:
    def __init__(self, U10 = 5, x = 20170,):
        # ускорение свободного падения.
        self.g = 9.81
        # скорость ветра на высоте 10 м над уровнем моря.
        self.U10 = U10
        # коэффициент gamma (см. спектр JONSWAP)
        self.gamma = self.Gamma(x)
        # коэффициент alpha (см. спектр JONSWAP)
        self.alpha = self.Alpha(x)
        # координата пика спектра по частоте
        self.omega_m = self.Omega(x) * self.g/self.U10
        # координата пика спектра по волновому числу
        self.k_m = self.k_max( self.omega_m )
        # массив с границами моделируемого спектра.
        self.KT = np.array([self.k_m/4,self.k_m*500])
        # k0 -- густая сетка, нужна для интегрирования и интерполирования
        self.k0= np.logspace(np.log10(self.KT[0]), np.log10(self.KT[-1]), 10**4)

    def get_spectrum(self):
        # интерполируем смоделированный спектр
        self.spectrum = self.interpolate()
        self.sigma_sqr = integrate.quad(self.spectrum, self.KT[0],self.KT[-1])[0]
        return self.spectrum

    def find_decision(self,omega):
        P = 9.8 * 1000.0/0.074
        Q = -1000.0*omega**2/0.074
        x1= -Q/2.0 + np.sqrt( (Q/2)**2 + (P/3)**3)
        x2= -Q/2.0 - np.sqrt( (Q/2)**2 + (P/3)**3)
        k=x1**(1/3)-(-x2)**(1/3)
        return k

    def det(self,k):
    #        Функция возвращает Якобиан при переходе от частоты к
    #    волновым числам по полному дисперсионному уравнению
        det=(self.g+3*k**2*0.074/1000)/(2*np.sqrt(self.g*k+k**3*0.074/1000) )
        return det

    def k_max(self,omega_max):
        # k_max -- координата пика спектра
        k_max=omega_max**2/self.g
        return k_max

    def omega_k(self,k):
        #    Пересчет волнового числа в частоту по полному дисперсионному
        # уравнению
        omega_k=(self.g*k+0.074*k**3/1000)**(1/2)
        return omega_k

    def JONSWAP(self,k):
        if k<=self.k_m:
            sigma=0.074
        else:
            sigma=0.09
        Sw=(
            self.alpha/2*k**(-3)*np.exp(-1.25*(self.k_m/k)**2 )*
            self.gamma**(np.exp(- ( np.sqrt(k/self.k_m)-1)**2 / (2*sigma**2) ))
           )
        return Sw

    # Безразмерный коэффициент Gamma
    def Gamma(self,x):
        if x>=20170:
            return 1
        gamma = (
               +5.253660929
               +0.000107622*x
               -0.03778776*np.sqrt(x)
               -162.9834653/np.sqrt(x)
               +253251.456472*x**(-3/2)
                )
        return gamma

    # Безразмерный коэффициент Alpha
    def Alpha(self,x):
        if x >= 20170:
            return 0.0081
        alpha = np.array( [],dtype = 'float64')
        alpha = [(
               +0.0311937
               -0.00232774 * np.log(x)
               -8367.8678786/x**2
               +4.5114599e+300*np.exp(-x)*1e+300*1e+17
    #            +4.5114599e+17*exp(-x)
              )]
        return alpha[0]

    #Вычисление безразмерной частоты Omega по безразмерному разгону x
    def Omega(self,x):
        if x>=20170:
            return 0.835
        omega_tilde=(0.61826357843576103
                     + 3.52883010586243843e-06*x
                     - 0.00197508032233982112*np.sqrt(x)
                     + 62.5540113059129759/np.sqrt(x)
                     - 290.214120684236224/x
        )
        return omega_tilde

    def spectrum0(self,n,k,spectrum_type = 'Karaev'):
        if spectrum_type == 'Karaev':
            power = [0,4,5,2.7,5]
        if n==0:
            return self.JONSWAP(k)
        else:
            omega0 = self.omega_k(self.limit_k[n-1])
            beta0  = self.spectrum0(n-1,self.limit_k[n-1]) * \
                        omega0**power[n]/self.det(self.limit_k[n-1])
            omega0 = self.omega_k(k)
            return beta0/omega0**power[n]*self.det(k)


    def full_spectrum(self,k,x=20170):
        #    Спектр JONSWAP.
        #    По совместительству, граница моделируюмого спектра #0
        # 0< omega < 1.2*omega_max
        # См. функции spectrum_{1-4}(k).
        # limit_{1-4} -- это соответствующие границы
        self.limit_1 = 1.2
        self.limit_2 =(
                 + 0.371347584096022408
                 + 0.290241610467870486 * self.U10
                 + 0.290178032985796564 / self.U10
                     )
        self.limit_3 = self.omega_k(270.0)
        self.limit_4 = self.omega_k(1020.0)
        self.limit_k = np.zeros(4)
        self.limit_k[0] = self.find_decision(self.limit_1 * self.omega_m)
        self.limit_k[1] = self.find_decision(self.limit_2 * self.omega_m)
        self.limit_k[2] = 270.0
        self.limit_k[3] = 1020.0
        try:
            full_spectrum = np.zeros(len(k))
        except:
            full_spectrum = [0]
            k = [k]

        for i in range(len(k)):
            if k[i] <= self.limit_k[0]:
                full_spectrum[i] =  self.spectrum0(0,k[i])
            elif k[i] <= self.limit_k[1]:
                full_spectrum[i] = self.spectrum0(1,k[i])
            elif k[i] <= self.limit_k[2]:
                full_spectrum[i] = self.spectrum0(2,k[i])
            elif k[i] <= self.limit_k[3]:
                full_spectrum[i] = self.spectrum0(3,k[i])
            else:
                full_spectrum[i] = self.spectrum0(4,k[i])
        return full_spectrum


    def interpolate(self):
        # Интерполируем наш спектр.
        # Позволяет не думать над различными размерами массивов при счете
        full_spectrum = interpolate.interp1d(self.k0,
                                             self.full_spectrum(self.k0))
        return full_spectrum



class surface(spectrum):
    def __init__(self, N=256, M=100, space='log', phi_optimize=0, random_phases = 1):
        spectrum = spectrum()
        self.M = M
        self.N = N
        KT = self.KT
        if space =='log':
            self.k_log = np.logspace(np.log10(KT[0]), np.log10(KT[-1]),N)
            self.k = self.k_log
        elif space =='lin':
            self.k_lin = np.linspace(KT[0], KT[-1], N)
            self.k = self.k_lin

        if phi_optimize == 0:
            self.M = M
            self.phi=np.linspace(-pi,pi,self.M)

        # случайные фазы
        if random_phases == 0:
            self.psi = np.array([
                    [0 for m in range(self.M) ] for n in range(self.N) ])
        elif random_phases == 1:
            self.psi = np.array([
                [ np.random.uniform(0,2*pi) for m in range(self.M)]
                            for n in range(self.N) ])
        # массив с амплитудами i-ой гармоники
        self.A = self.amplitude(self.k)
        # угловое распределение
        self.F = self.angle(self.k,self.phi)

    def B(self,k):
          def b(k):
              b=(
                  -0.28+0.65*np.exp(-0.75*np.log(k/self.k_m))
                  +0.01*np.exp(-0.2+0.7*np.log10(k/self.k_m))
                )
              return b
          B=10**b(k)
          return B
    # Перепиши ради христа эту фигню. Не читабельно
    def Phi(self,k,phi):
        # Функция углового распределения
        N = self.N
        M = self.M
        normalization = lambda B: B/np.arctan(np.sinh(2*pi*B))
        B0 = self.B(k)
        A0 = normalization(B0)
        try:
            Phi = [[ A0[i]/np.cosh(2*B0[i]*phi[j]) for j in range(M) ] for i in range(N)]
            Phi=np.array(Phi)
        except:
            Phi = A0/np.cosh(2*B0*phi)
        return Phi


    def angle(self,k,phi):
        N, M = self.N, self.M
        Phi = self.Phi(k,phi)
        dphi = np.array([phi[j] - phi[j-1]  for j in range(1,M) ] ) 
        dF = np.delete(Phi,0,1)
        F = ([ [np.sqrt(dF[n][m]*dphi[m]) for m in range(M-1)]
                                            for n in range(N) ])
        for i in range(N):
                F[i].append(0)
        F=np.array(F)
        return F
    
    def amplitude(self, k):
        N = len(k)
        S = self.spectrum
        integral = np.array([ integrate.quad(S,k[i-1],k[i])[0] for i in range(1,N) ])
        amplitude = np.sqrt(2 *integral )
        amplitude = np.concatenate((amplitude,[0]))
        return amplitude

    def model(self, k, phi,time=[0]):
        def water(r,t):
            N = self.N
            M = self.M
  
            A = self.A
            F = self.F
            psi = self.psi
            self.surface = 0
            progress_bar = tqdm( total = N*int(np.mean(M)),  leave = False )
#            progress_bar.set_description("Processing %s" % t)
            for n in range(N):
                for m in range(M):
                  self.surface += A[n] * \
                  np.cos(
                      +k[n]*(r[0]*np.cos(phi[m])+r[1]*np.sin(phi[m]))
                      +psi[n][m]
                      +self.omega_k(k[n])*t) \
                      * F[n][m]
                  progress_bar.update(1)
            progress_bar.close()
            progress_bar.clear()

            return self.surface
        return water