% Testing wave simulations of cosinusoidal wave superposition method
%              date: 2009-4-1
%--------------------------------------------------------------------------
clear all
U10=5;
gravity=9.81;
x1=20170;
LL=0.1;   % radar wavelength (L = 0.008, 0.021, 0.03, 0.055, 0.1, 0.23m)
KM=calkm(U10,x1); 
LM=2*pi/KM,
Kmin=KM/4, 
Kmax= KL_LIMIT(U10,LL), 
Wmin=k2w(Kmin);
Wmax=k2w(Kmax);

Nx=2^9;Ny=2^9;Nt=2^0;
dx=pi/Kmax,dy=dx;dt=pi/Wmax,

x=(0:Nx-1)*dx;y=(0:Ny-1)*dy;t=(0:Nt-1)*dt;
%--------------------------------------------------------------------------
%        simulaion
%--------------------------------------------------------------------------
Nw = 100;                    % 把频域分割成500份
Ntheta = 100;                % 把方向分割成500份
d_w=(Wmax-Wmin)/Nw;          % Δω
d_phi=2*pi/Ntheta;           % Δφ
w1=linspace(Wmin,Wmax,Nw+1);  % generates Nw+1 linearly equally spaced points between 0 and w_max;
w_m=zeros(1,Nw);
for i=1:Nw
     w_m(i)=(w1(i)+w1(i+1))/2;
end
RAN=rand(Ntheta,Nw); 
w=zeros(Ntheta,Nw);phi=zeros(Ntheta,Nw);a=zeros(Ntheta,Nw);
for  m=1:Ntheta
     for  n=1:Nw
         w(m,n)=w_m(n)-0.5*d_w+(m-1+RAN(m,n))*d_w/Ntheta;        %ω(m)         
         phi(m,n)=-pi+(m-0.5)*d_phi;                             %φ(n)
         S=ModWavSpc(w(m,n),U10,x1)*spread(w(m,n),phi(m),U10,x1);% S(ω,φ)
         a(m,n)=sqrt(2*S*d_w*d_phi);                   % amplitude a(ω,φ)
      end
end
 
eta=zeros(Ny,Nx);
randphase=2*pi*rand(Ntheta,Nw); % random initial phase,uniform distribution between (0,2π)
for i=136:Ny
    for j=1:Nx
%         for k=1:Nt
          h=a.*cos(0-w2k(w).*(x(j)*cos(phi)+y(i)*sin(phi))+randphase);
          eta(i,j)=sum(sum(h));   %ξ(y,x,t)
%         end
     end
end
eta=squeeze(eta);
%--------------------------------------------------------------------------
%        show sea surface 
%--------------------------------------------------------------------------
figure(1);clf
colormap('winter')
surf(x,y,eta)
shading interp
axis([x(1) x(end) y(1) y(end) min(min(eta)) max(max(eta))])
axis('square')
view([0 90])
colorbar
%--------------------------------------------------------------------------
%     compare spectrum with input spectrum 
%--------------------------------------------------------------------------       
% Spec=(abs(fft2(eta-mean(mean(eta))))).^2*dx*dy/(Nx*Ny*2*pi*2*pi); 
%                                         % calculate spectrum from ξ(y,x,t)
% 
% Sd=windspec(U10,x1,Wmin,Wmax,Nw,Ntheta); 
% Snew = spec2spec(Sd,'k2d');  % S(ω,φ)--> F(Kx,Ky)
% figure,subplot(1,2,1),
% contour(Snew.k,Snew.k2,Snew.S,5),% wspecplot(Snew),
% xlabel('kx'),ylabel('ky'),title('Input wave spectrum')
% caxis([min(min(Snew.S)) max(max(Snew.S))]),h=colorbar('location','southoutside');
% xlabel('kx'),ylabel('ky'),title('Input wave spectrum for simulation')
% axis([-2*KM 2*KM -2*KM 2*KM]),axis square
% grid on
% 
% subplot(1,2,2),
% kx=2*pi*(-Nx/2:Nx/2-1)/Nx/dx;
% ky=2*pi*(-Ny/2:Ny/2-1)'/Ny/dy;
% figure,contour(kx,ky,fftshift(squeeze(Spec)),5)
% caxis([min(min(Snew.S)) max(max(Snew.S))]),h=colorbar('location','southoutside');
% xlabel('kx'),ylabel('ky'),title('Wave spectrum from simulation')
% axis([-2*KM 2*KM -2*KM 2*KM]),axis square
% grid on
% 
% HeightVariance=var(reshape(eta,Nx*Ny,1))
% [fx,fy]=gradient(eta);fx=fx/dx;fy=fy/dy;
% SlopVar_x=var(reshape(fx,Nx*Ny,1))
% SlopVar_y=var(reshape(fy,Nx*Ny,1))
    
%--------------------------------------------------------------------------
% Results：
%--------------------------------------------------------------------------

% LM =
%   22.989040764296604
% Kmin =
%    0.068328050000000
% Kmax =
%   19.838488169839607
% dx =
%    0.158358471003146
% dt =
%    0.224977031698681
% HeightVariance =
%    0.016227130855940
% SlopVar_x =
%    0.007065869028233
% SlopVar_y =
%    0.003773346104886