% Testing wave simulations of FFT method
%        
% Figure 1: sea surface 
% Figure 2: compare data spectrum with input spectrum
% Figure 3: comparison of correlation functions
% Figure 4: spectrum of height,slop and orbital velocity
%              date: 2009-4-7
%--------------------------------------------------------------------------

clear all
U10=5;    % wind speed in 10m height is 5m/s;
x1=20170; % fully dveloped wind waves
LL=0.1;   % radar wavelength (L = 0.008, 0.021, 0.03, 0.055, 0.1, 0.23m)
KM=calkm(U10,x1); 
LM=2*pi/KM;
Kmin=KM/4; 
Kmax= KL_LIMIT(U10,LL); 
%--undirectional wave spectrum--
S=createspec('freq','w');  
Nw=2^10;  % N¦Ø 
S.w=linspace(k2w(Kmin),k2w(Kmax),Nw);
S.S=ModWavSpc(S.w,U10,x1);
%--------------------------------------------------------------------------
%        simulaion
%--------------------------------------------------------------------------
dt=0.25;Nt=200*60*4;t=(0:Nt-1)*dt;
[xs, xsder]=spec2sdat(S,[Nt 1],dt); % 1D wave simulation
                                    % xs --> X(t), xsder --> X'(t)
%--------------------------------------------------------------------------
%        plot sea surface 
%--------------------------------------------------------------------------
% SWH=4*sqrt(spec2mom(S,0)) % answer>> SWH = 0.576295156780311
figure,plot(xs(:,1),xs(:,2))
xlabel(['times [s]','   TM=',num2str(2*pi/k2w(KM)),'s']);
ylabel( 'wave height[m]');
saveas(gca,'TestSimu_New2_surface','fig')  % save figure
%--------------------------------------------------------------------------
%     compare spectrum with input spectrum 
%--------------------------------------------------------------------------
figure,wspecplot(S,1,'b')           % plot input wave spectrum
Sest = dat2spec2(xs,300);           % calculate wave spectrum from X(t)
hold on, wspecplot(Sest,1,'r')
xlabel('Angular frequency [rad/s]')
saveas(gca,'TestSimu_New2_ComSpec','fig')  % save figure
%--------------------------------------------------------------------------
%     comparison of covariance functions
%--------------------------------------------------------------------------
R_S= spec2cov2(S,2,400,0.25);            % covariance from input spectrum
figure,plot((0:400)*0.25,R_S(:,1),'b')
R_xs=dat2cov(xs(:,[1 2]),400,'unbiased');% covariance from X(t)
hold on,plot(R_xs.t,R_xs.R,'r')
axis([0 100 -0.02 0.025])
xlabel('Time Lag [sec]')
ylabel('ACF')
title('Auto Covariance Function (ACF)')
legend('From spectrum','From simulation data')
saveas(gca,'TestSimu_New2_ComCov','fig')  % save figure
%--------------------------------------------------------------------------
%     spectrums 
%--------------------------------------------------------------------------
figure,semilogy(Sest.w/2/pi,Sest.S*2*pi,'k') % plot spectrum of height
OrbVel = dat2spec2(xsder,300);               % spectrum of orbital velocity
hold on,semilogy(OrbVel.w/2/pi,OrbVel.S*2*pi,'r')
xlabel('Frequency [1/s]')
ylabel('S(f) [m^2 s]')
title('Spectral density')

Sd=createspec('dir','w');           % 2D directional wave spectrum
Sd.w=linspace(k2w(Kmin),k2w(Kmax),Nw);
Sd.theta=linspace(-pi,pi,20)';
Sd.S=ones(size(Sd.theta))*ModWavSpc(Sd.w,U10,x1)/(2*pi);
     % omnidirectional wave spectrum
     % wspecplot(Sd) % plot 2D spectrum in polar coordinates
Nx=2^10;Ny=2^0;Nt=2^0;dx=pi/Kmax;dy=dx;dt=0.25;
H=seasim(Sd,Nx,Ny,Nt,dx,dy,dt,2,0); % wave simulation using FFT method
fx=gradient(H.Z,dx);                % slop
slopx=[(0:Nx-1)'*dx fx'];
SlopxS = dat2spec2(slopx,300);
SlopxS.w=k2w(SlopxS.w);
hold on,semilogy(SlopxS.w/2/pi,SlopxS.S*2*pi,'b') % plot spectrum of orbital velocity
legend('height','orbital velocity','slop')
saveas(gca,'TestSimu_New2_Specs','fig')           % save figure