import numpy as np
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm
from water.spectrum import Spectrum

class Surface(Spectrum):
    def __init__(self,  N=256, M=100, space='log', 
        random_phases = 1, whitening = 0, wind = 0,**kwargs):
        Spectrum.__init__(self,**kwargs)
        self.get_spectrum()
        self.N = N
        self.M = M
        KT = self.KT
        self.wind = wind # Направление ветра
        self.k = np.logspace(np.log10(KT[0]), np.log10(KT[-1]),self.N + 1)
        if whitening == 1:
            interspace = self.interspace(self.k, N//2, power=0)
            # print(interspace)
            self.k_heights = self.nodes(interspace,power=0)
            interspace = self.interspace(self.k, N, power=2)
            self.k_slopes = self.nodes(interspace,power=2)  
            # self.k_slopes = np.logspace(np.log10(KT[0]), np.log10(KT[-1]),self.N)
            self.k = np.hstack((self.k_heights,self.k_slopes))
            self.k = np.sort(self.k)

        self.phi=np.linspace(0,2*np.pi,self.M + 1) 
        self.phi_c = self.phi + self.wind

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

    def Phi(self,k,phi):
        # Функция углового распределения
        phi = phi -self.wind
        normalization = lambda B: B/np.arctan(np.sinh(2* (pi)*B))
        B0 = self.B(k)
        A0 = normalization(B0)
        Phi = A0/np.cosh(2*B0*(phi) )
        return Phi


    def angle(self,k,phi):
        M = self.M
        N = self.N
        # print(k.size)
        Phi = lambda phi,k: self.Phi(k,phi)
        integral = np.zeros((N,M))
        # print(integral.shape)
        for i in range(N):
            for j in range(M):
                # integral[i][j] = integrate.quad( Phi, phi[j], phi[j+1], args=(k[i],) )[0]
                integral[i][j] = np.trapz( Phi( phi[j:j+2],k[i] ), phi[j:j+2])
        amplitude = np.sqrt(2 *integral )
        return amplitude
    
    def amplitude(self, k):
        N = len(k)
        S = self.spectrum
        integral = np.array([ integrate.quad(S,k[i-1],k[i])[0] for i in range(1,N) ])
        amplitude = np.sqrt(2 *integral )
        return np.array(amplitude)

    def model(self,r,t):
        N = self.N 
        M = self.M
        self.k = self.k[:N]
        k = self.k
        phi = self.phi
        A = self.A
        F = self.F
        psi = self.psi
        self.surface = 0
        self.amplitudes = np.array([ A[i]*sum(F[i])  for i in range(N)])
        progress_bar = tqdm( total = N*M,  leave = False )
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
        print()
        return self.surface

    def D(self,r,t):
        N = self.N
        M = self.M
        k = self.k
        phi = self.phi
        A = self.A
        F = self.F
        psi = self.psi
        Dx = 0
        Dy = 0
        for n in range(N):
            for m in range(M):
                Dx += A[n]  * np.cos(phi[m]) *\
                np.cos(
                        +k[n]*(r[0]*np.cos(phi[m])+r[1]*np.sin(phi[m]))
                        +psi[n][m] + np.pi/2+self.omega_k(k[n])*t)  \
                        * F[n][m]

                Dy += A[n]  * np.sin(phi[m]) *\
                np.cos(
                        +k[n]*(r[0]*np.cos(phi[m])+r[1]*np.sin(phi[m]))
                        +psi[n][m] + np.pi/2+self.omega_k(k[n])*t)  \
                        * F[n][m]


        return np.array(Dx),np.array(Dy)

# Экспериментальный способ. Метод <<отбеливания>> спектра
    def interspace(self, k, N, power=0):
        #     Функция разбивает заданный k-интервал на N участков, в каждом
        # из которых интеграл по спектру принимает одно и тоже значение.
        full_spectrum=self.spectrum
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

    def nodes(self,ki, power=0):
        # Функции interspace и nodes используются последовательно.
        # Первая возвращает интервалы, а вторая вычисляет в возвращенных
        # интервалах координаты узлов
        sigma = self.sigma_sqr
        b0 = sigma/(len(ki)+1)
        progress_bar = tqdm(total = len(ki)-1)
        full_spectrum=self.spectrum

        y=lambda k: k**(2+power)*full_spectrum(k)

        nodes=np.zeros(len(ki))
        A=np.sqrt(2*b0)

        for i in range(1,len(ki)):
            integral,error=integrate.quad(y,ki[i-1],ki[i])
            B=(np.sqrt(integral/A**2))
            nodes[i-1]=B
            progress_bar.update(1)
        nodes = nodes[:-1]
        if nodes[-1]>self.KT[-1]:
            nodes[-1]=self.KT[-1]

        return nodes