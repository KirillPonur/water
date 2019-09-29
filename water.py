import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm



class spectrum:
    def __init__(self, KT = [0.05, 2000], \
                 U10 = 5, x = 20170):
        # ускорение свободного падения.
        self.g = 9.81
        # скорость ветра на высоте 10 м над уровнем моря.
        self.U10 = U10
        # массив с границами моделируемого спектра.
        self.KT = np.array(KT)
        # k0 -- густая сетка, нужна для интегрирования и интерполирования
        self.k0= np.logspace(np.log10(KT[0]), np.log10(KT[-1]), 10**4)
        # коэффициент gamma (см. спектр JONSWAP)
        self.gamma = self.Gamma(x)
        # коэффициент alpha (см. спектр JONSWAP)
        self.alpha = self.Alpha(x)
        # координата пика спектра по частоте
        self.omega_m = self.Omega(x) * self.g/self.U10
        # координата пика спектра по волновому числу
        self.k_m = self.k_max( self.omega_m )
        
        # интерполируем смоделированный спектр
        self.spectrum = self.interpolate()
        self.sigma_sqr = integrate.quad(self.spectrum, KT[0],KT[-1])[0]

        
    def find_decision(self,omega):
        P = 9.8 * 1000.0/0.074;
        Q = -1000.0*omega**2/0.074;
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
            sigma=0.07
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

    def spectrum1(self,k):
        # Граница моделируюмого спектра #1
        # 1.2*omega_max < omega < a_m omega_m
        omega0 = self.omega_k(self.limit_k[0])
        beta0  = self.JONSWAP(self.limit_k[0]) * \
                    omega0**4/self.det(self.limit_k[0])
        omega0 = self.omega_k(k)
        return beta0/omega0**4*self.det(k)

    def spectrum2(self,k):
        # Граница моделируюмого спектра #2
        # a_m(U_10) omega_m < omega < omega_gk simeq 60  рад/с
        omega0 = self.omega_k(self.limit_k[1])
        beta0  = self.spectrum1(self.limit_k[1]) * \
                    omega0**5/self.det(self.limit_k[1])
        omega0 = self.omega_k(k)
        return beta0/omega0**5*self.det(k)

    def spectrum3(self,k):
        # Граница моделируюмого спектра #3
        #  omega_gk < omega < omega_h simeq 290 рад/с
        omega0 = self.omega_k(self.limit_k[2])
        beta0  = self.spectrum2(self.limit_k[2]) * \
                    omega0**2.7/self.det(self.limit_k[2])
        omega0 = self.omega_k(k)
        return beta0/omega0**2.7*self.det(k)

    def spectrum4(self,k):
        # Граница моделируюмого спектра #4
        #  omega_h < omega
        omega0 = self.omega_k(self.limit_k[3])
        beta0  = self.spectrum3(self.limit_k[3]) * \
                    omega0**5/self.det(self.limit_k[3])
        omega0 = self.omega_k(k)
        return beta0*self.det(k)/omega0**5

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
        self.limit_3 = 270.0
        self.limit_4 = 1020.0
        self.limit_k = np.zeros(4)
        self.limit_k[0] = self.find_decision(self.limit_1 * self.omega_m)
        self.limit_k[1] = self.find_decision(self.limit_2 * self.omega_m)
        self.limit_k[2] = self.limit_3
        self.limit_k[3] = self.limit_4
        try:
            full_spectrum = np.zeros(len(k))
        except:
            full_spectrum = [0]
            k = [k]

        # Цикл ниже сшивает спектр на границах.
        for i in range(len(k)):
            if k[i] <= self.limit_k[0]:
                full_spectrum[i] =  self.JONSWAP(k[i])
            elif k[i] <= self.limit_k[1]:
                full_spectrum[i] = self.spectrum1(k[i])
            elif k[i] <= self.limit_k[2]:
                full_spectrum[i] = self.spectrum2(k[i])
            elif k[i] <= self.limit_k[3]:
                full_spectrum[i] = self.spectrum3(k[i])
            else:
                full_spectrum[i] = self.spectrum4(k[i])
        return full_spectrum

    def interpolate(self):
        # Интерполируем наш спектр.
        # Позволяет не думать над различными размерами массивов при счете
        full_spectrum = interpolate.interp1d(self.k0,
                                             self.full_spectrum(self.k0))
        return full_spectrum
    

    
    
class whitening():
# Экспериментальный способ. Метод <<отбеливания>> спектра
    def interspace(self, k, N, power=2):
        #     Функция разбивает заданный k-интервал на N участков, в каждом
        # из которых интеграл по спектру принимает одно и тоже значение.
        full_spectrum=self.get_full_spectrum()
        y = lambda k: full_spectrum(k)*k**power
        sigma = self.sigma_sqr
        b0 = sigma/(N)
        k_new = k[0]
        k_end = k[-1]
        k = np.zeros(N+1)
        k[0] = k_new
#        err = np.zeros(N)
        epsabs = 1.49e-8
#        epsrel = 1.49e-8
        progress_bar = tqdm(total = N-1 )
        for i in range(N-1):
            integral = 0
            n = 1
            m = 1
            while integral < b0:
                k_new*= 1 + 10**(-m)
                integral,error = integrate.quad(y, k[i], k_new)
                if integral> b0 and abs(integral-b0) > epsabs:
                    k_new*= 1 - 10**(-m)
                    m+= 1
                    integral, error = integrate.quad(y, k[i], k_new)
                n+= 1
            progress_bar.update(1)
            k[i+1] = k_new
        k = k[0:N+1]
        k[-1] = k_end
        # Возвращает k размерности N+1, т.о. имеем как раз N интервалов.
        # err -- ошибка вычисения интегралов
        progress_bar.close()
        return k

    def nodes(self,ki, power=2):
        # Функции interspace и nodes используются последовательно.
        # Первая возвращает интервалы, а вторая вычисляет в возвращенных
        # интервалах координаты узлов
        sigma = self.sigma_sqr
        b0 = sigma/(N)
        progress_bar = tqdm(total = len(ki)-1)
        full_spectrum=self.get_full_spectrum()

        y=lambda k: k**(2+power)*full_spectrum(k)

        nodes=np.zeros(len(ki))
        A=np.sqrt(2*b0)

        for i in range(1,len(ki)):
            integral,error=integrate.quad(y,ki[i-1],ki[i])
            B=(np.sqrt(integral/A**2))
            nodes[i-1]=B
            progress_bar.update(1)

        return nodes[:-1]

class surface:
    def __init__(self,N=256, M=80, space='log',phi_optimize=0,
                 random_phases = 1):
#        water()
        self.M = M
        self.N = N    
        KT = self.KT
        if space =='log':
            self.k_log = np.logspace(np.log10(KT[0]), np.log10(KT[-1]),N)
            self.k = self.k_log
        elif space =='lin':
            self.k_lin = np.linspace(KT[0], KT[-1], N)
            self.k = self.k_lin
#        elif space=='white':
#            ki,b0=self.interspace(self.k0,N)
#            self.k=nodes(ki,b0)

        if phi_optimize == 0:
            self.M = M*np.ones(self.N, dtype='int')
            self.phi=([np.linspace(-pi,pi,self.M[i]) for i in range(N)])

        # случайные фазы
        if random_phases == 0:
            self.psi = np.array([
                    [0 for m in range(self.M[n]) ] for n in range(N) ])
        elif random_phases == 1:
            self.psi = np.array([
                [ np.random.uniform(0,2*pi) for m in range(self.M[n])]
                            for n in range(N) ])
        # массив с амплитудами i-ой гармоники
        self.A = self.amplitude(self.k)
        # угловое распределение
        self.F = self.angle(self.k,self.phi)
        
        water.k = self.k
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
            Phi = ([
                    [ A0[i]/np.cosh(2*B0[i]*phi[i][j]) for j in range(M[i]) ]
                                                    for i in range(N)])
            for i in range(len(Phi)):
                for j in range(max(M)-len(Phi[i])):
                    Phi[i].append(0)
            Phi=np.array(Phi)
        except:
            try:
                Phi = [[ A0[i]/np.cosh(2*B0[i]*phi[j]) for j in range(M) ]
                                                    for i in range(N)]
                Phi=np.array(Phi)
            except:
                Phi = A0/np.cosh(2*B0*phi)
        return Phi


    def angle(self,k,phi):
        N, M = self.N, self.M
        Phi = self.Phi(k,phi)
        dphi = np.array([[phi[i][j] - phi[i][j-1] for j in range(1,M[i])]
                                                    for i in range(N)])
        dF = np.delete(Phi,0,1)
        F = ([ [np.sqrt(dF[n][m]*dphi[n][m]) for m in range(M[n]-1)]
                                            for n in range(N) ])
        for i in range(len(F)):
            for j in range(max(M)-len(F[i])):
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
                for m in range(M[n]):
                  self.surface += A[n] * \
                  np.cos(
                      +k[n]*(r[0]*np.cos(phi[n][m])+r[1]*np.sin(phi[n][m]))
                      +psi[n][m]
                      +self.omega_k(k[n])*t) \
                      * F[n][m]
                  progress_bar.update(1)
            progress_bar.close()
            progress_bar.clear()

            return self.surface
        return water
      
class correlation:
    def angles(self,rho):
        try:
            return self.correlation_slopes(rho)
        except:
            progress_bar=tqdm(total = len(rho) + 1)
            progress_bar.set_description("Вычисление функции корреляции")
            k= np.logspace(np.log10(self.KT[0]), np.log10(self.KT[-1]), 5*10**5)
            integral=np.zeros(len(rho))
            y=lambda k: k**2*self.spectrum(k)
            for i in range(len(rho)):
                integral[i]=np.trapz(y(k)*np.cos(k*rho[i]),x=k)
                progress_bar.update(1)
            progress_bar.set_description("Сохранение функции корреляции")
            self.correlation_slopes = interpolate.interp1d(rho,integral)
            progress_bar.update(1)
            progress_bar.close()
            return integral
    

    def angles_sum(self,k,rho):
        f=0
        A=self.amplitude(k)
        f=np.zeros(len(rho))
        for j in range(len(rho)):
                f[j]=sum( k**2*A**2/2*np.cos(k*rho[j]) )
        return f
    
    def height_restore(self, rho, k, fourier = 'real'):
        if fourier == 'real':
            self.K=self.height(rho,k)
            integral=np.zeros(len(k))
            for i in range(len(k)):
                integral[i]=integrate.simps(self.K*np.cos(k[i]*rho),x=rho)
            return integral/np.sqrt(np.pi)
        
        elif fourier == 'complex':
            self.K=self.height(rho,k,'complex')
            integral=np.fft.ifft(self.K,n=k.size)
            return integral.real
            
    def height(self,rho,k, fourier = 'real'):
        if fourier == 'real':
            S=self.spectrum(k)
            integral=np.zeros(len(rho))
            for i in range(len(rho)):
                integral[i]=np.trapz(S*np.cos(k*rho[i]),x=k)
            self.correlation_height = interpolate.interp1d(rho,integral)
            return integral/max(integral)
        
        if fourier == 'complex':
            S=self.spectrum(k)
            integral=np.fft.fft(S,n=k.size)
            return integral
    def height_sum(self,k,rho):
        f=0
        A=self.amplitude(k)
        f=np.zeros(len(rho))
        for j in range(len(rho)):
                f[j]=sum( A**2/2*np.cos(k*rho[j]) )
        return f
    
class water(spectrum,correlation,surface):
    # Этот класс содержит некоторые инструменты для работы с моделью
    def __init__(self,N=256,KT=[0.05,2000]):
        self.x=np.linspace(0,200,200)
        self.y=np.linspace(0,200,200)
        self.t=np.array([0])
        self.N=N
        self.KT=np.array(KT)
        spectrum.__init__(self,KT=self.KT)
        surface.__init__(self,N=self.N)
        self.rho=np.linspace(0,100,1000)
        self.rho0=np.linspace(0,100,self.k0.size)
    def plot(self,fig):
        plt.figure()
        if fig=='slopes' or fig=='s':
            self.plot_slopes(self.k0)
        elif fig=='heights' or fig=='h':
            self.plot_heights(self.k0)
        elif fig=='correlation slopes' or fig=='cs':
            self.plot_correlation_slopes(self.rho)
        elif fig=='correlation heights' or fig=='ch':
            self.plot_correlation_heights(self.rho)
        elif fig=='surface':
            self.plot_surface(self.x,self.y,self.t)
            
    def plot_correlation_slopes(self,rho):
        plt.plot(self.rho,self.angles(self.rho))
        plt.ylabel(r'$K_q(\rho)$',fontsize=16)
        plt.xlabel(r'$\rho$',fontsize=16)
        
    def plot_correlation_heights(self,rho):
        plt.plot(self.rho,self.height(self.rho,self.k0))
        plt.ylabel(r'$K(\rho)$',fontsize=16)
        plt.xlabel(r'$\rho$',fontsize=16)
        
    def plot_heights(self,k):
        plt.loglog(self.k0,self.spectrum(self.k0))
        plt.ylabel(r'$S_q(k)=S(k)k^2$',fontsize=16)
        plt.xlabel(r'$k$',fontsize=16)
        
    def plot_slopes(self,k):
        plt.loglog(self.k0,self.k0**2*self.spectrum(self.k0))    
        plt.ylabel(r'$S(k)$',fontsize=16)
        plt.xlabel(r'$k$',fontsize=16)
        
    def plot_surface(self, x , y, t):
        surface = self.model(self.k, self.phi, t)
        x, y = np.meshgrid(x, y)
        self.z=surface([x,y],t)
        from matplotlib.cm import winter
        plt.contourf(self.z,100,cmap=winter)
        plt.colorbar()
        plt.ylabel(r'Y, \text{м}',fontsize=16)
        plt.xlabel(r'X, \text{м}',fontsize=16)
    def plot_restore(self,fourier='real'):
        rho=self.rho
        k=self.k0
        plt.loglog(k,self.spectrum(k),label='Исходный спектр')
        plt.loglog(self.k,self.height_restore(rho,self.k,fourier),label='Восстановленный спектр')
        plt.legend()
    def plot_restore1(self):
        rho=self.rho0
        k=self.k0
        plt.plot(rho,self.height(rho,k,'real'),label='Вещественный')
        plt.plot(rho,self.height(rho,k,'complex'),label='Комплексный')
        plt.legend()       
        #   savefig(path.abspath('..'+'\\water\\anim\\'+'water'+str(i)+'.png'),
        #             pdi=10**6,bbox_inches='tight')
#        show()
water=water()
