% Testing wave simulations of cosinusoidal wave superposition method
%        
% Figure 1: sea surface 
% Figure 2: compare data spectrum with input spectrum
% Figure 3: comparison of correlation functions
% Figure 4: spectrum of height,slop and orbital velocity
%              date: 2009-4-7
%--------------------------------------------------------------------------
clear all
gravity=9.8
g=gravity;

U10=5;    % wind speed in 10m height is 5m/s;
x1=20170; % fully dveloped wind waves
LL=0.1;   % radar wavelength (L = 0.008, 0.021, 0.03, 0.055, 0.1, 0.23m)
KM=calkm(U10,x1); 
LM=2*pi/KM,
Kmin=KM/4, 
Kmax= KL_LIMIT(U10,LL), 
Wmin=k2w(Kmin);
Wmax=k2w(Kmax);

dr=pi/Kmax,dt=pi/Wmax
Nr=2^10;Nt=2^10; 
r=(0:Nr-1)*dr;t=(0:Nt-1)*dt;
%--------------------------------------------------------------------------
%        simulaion
%--------------------------------------------------------------------------
Nw = 500;                        % N¦Ø
d_w=(Wmax-Wmin)/Nw;              % ¦¤¦Ø
w1=linspace(Wmin,Wmax,Nw+1);  % generates Nw+1 linearly equally spaced points between 0 and w_max;
w_m=zeros(1,Nw);
for i=1:Nw
     w_m(i)=(w1(i)+w1(i+1))/2;
end
RAN=rand(Nw);
w=zeros(1,Nw);a=zeros(1,Nw);
for  n=1:Nw
     w(n)=w_m(n)-0.5*d_w+RAN(n)*d_w;   %¦Ø(m)         
     S=ModWavSpc(w(n),U10,x1);         % S(¦Ø)
     % S=ModWavSpc(w(n),U10,x1)*spread(w(m,n),phi(m),U10,x1)  %S(¦Ø,¦Õ)
     a(n)=sqrt(2*S*d_w);               % amplitude a(¦Ø)
end
 
eta=zeros(Nr,Nt);          % initialization wave height
randphase=2*pi*rand(1,Nw); % random phase,uniform distribution between (0,2¦Ð)
for i=1:Nr
   i
   for k=1:Nt     
     h=a.*cos(w*t(k)-w2k(w).*r(i)+randphase);
     eta(i,k)=sum(sum(h));   %¦Î(r,t)
    end
end
%--------------------------------------------------------------------------
%        plot sea surface 
%--------------------------------------------------------------------------
figure,plot(r,eta(:,1),'k');
h1=gca;
set(h1,'XAxisLocation','top','box','off','XColor','k','YColor','k');
p=get(h1,'Position');
xlabel(['space [m]','   LM=',num2str(2*pi/KM),'m']);
ylabel( 'wave height[m]');
xlim([min(r) max(r)])
ylim([-0.5 0.5])

h2=axes('Position',p);
set(h2,'XAxisLocation','bottom',...
    'YAxisLocation','right','Color','none','XColor','r');
hold on,plot(t,eta(1,:),'r');
xlabel(['times [s]','   TM=',num2str(2*pi/k2w(KM)),'s']);
xlim([min(t) max(t)])
ylim([-0.5 0.5])
saveas(gca,'TestSimu_New1_surface','fig')  % save figure
    % OR:
    % xlabels{1} = 'space [m]';
    % xlabels{2} = 'times [s]';
    % ylabels{1} = 'wave height[m]';
    % ylabels{2} = '';
    % figure,[ax,h1,h2]=plotxx(r,eta(:,1),t,eta(1,:),xlabels,ylabels); % plot sea surface
    % yLim=[-0.5 0.5];
    % set(ax(1),'YLim',yLim),set(ax(1),'xlim',[min(r) max(r)])
    % set(ax(2),'YLim',yLim),set(ax(2),'xlim',[min(t) max(t)])
    % title(['LM=',num2str(2*pi/KM) ,'m  TM=',num2str(2*pi/k2w(KM)),'s'])
%--------------------------------------------------------------------------
%     compare spectrum with input spectrum 
%--------------------------------------------------------------------------
xs=[t' eta'];
L=200;                      %  maximum lag size of the parzen window
Sest = dat2spec2(xs,L);     %  calculate wave spectrum from X(t)
figure,wspecplot(Sest,1,'r') 

S=createspec('freq','w');   % undirectional wave spectrum
S.w=linspace(k2w(Kmin),k2w(Kmax),Nw);
S.S=ModWavSpc(S.w,U10,x1);
                             % SWH=4*sqrt(spec2mom(S,0))= 0.576295156780311
hold on, wspecplot(S,1,'b')  % plot input wave spectrum
xlabel('Angular frequency [rad/s]')
saveas(gca,'TestSimu_New1_ComSpec','fig')  % save figure
%--------------------------------------------------------------------------
%     comparison of correlation functions
%--------------------------------------------------------------------------
R_S= spec2cov2(S,2,400,0.25);             % covariance from input spectrum
figure,plot((0:400)*0.25,R_S(:,1),'b')
R_xs=dat2cov(xs(:,[1 2]),400,'unbiased'); % covariance from X(t)
hold on,plot(R_xs.t,R_xs.R,'r')
axis([0 100 -0.02 0.025])
xlabel('Time Lag [sec]')
ylabel('ACF')
title('Auto Covariance Function (ACF)')
legend('From spectrum','From simulation data')
saveas(gca,'TestSimu_New1_ComCov','fig')  % save figure
%--------------------------------------------------------------------------
%     spectrums 
%--------------------------------------------------------------------------
%-- spectrum of height--
figure,semilogy(Sest.w/2/pi,Sest.S*2*pi,'k') 
xlabel('Frequency [1/s]')
ylabel('S(f) [m^2 s]')
title('Spectral density')
%--  spectrum of slop --
[orbvel slop]=gradient(eta,dt,dr);
slop=[r' slop];
SlopS=dat2spec2(slop,L);
SlopS.w=k2w(SlopS.w);
hold on,semilogy(SlopS.w/2/pi,SlopS.S*2*pi,'b')
%--  spectrum of orbital velocity --
orbvel=[t' orbvel'];
OrbVel=dat2spec2(orbvel,L);
hold on,semilogy(OrbVel.w/2/pi,OrbVel.S*2*pi,'r');
legend('height','slop','orbital velocity')
saveas(gca,'TestSimu_New1_Specs','fig')  % save figure
%--------------------------------------------------------------------------
%  function of dat2spec2(xs,L)  
%--------------------------------------------------------------------------
% L=200;
% dt=xs(2,1)-xs(1,1);Nt=length(xs);
% xs(:,2)=xs(:,2)-mean(mean(xs(:,2)));
% R=dat2cov(xs(:,[1 2]));
% r=R.R(:);
% win=parzen(2*L-1);  % Does not give negative estimates of the spectral density
%   
% nf=2^nextpow2(2*L-2);  %  Interpolate the spectrum with rate = 2
% nfft=2*nf;
% 
% S=createspec('freq');
% S.tr=9.8;
% S.L=L;
% S.S=zeros(nf+1,1);
% 
% add a nugget effect to ensure that round off errors
% do not result in negative spectral estimates
% 
% r=r+eps;
% rwin=r(1:L).*win(L:(2*L-1)); 
% Rper=real(fft([rwin; zeros(nfft-(2*L-1),1); rwin(L:-1:2)],nfft));
% 
% S.w=[0:(nf)]'/nf*pi/dt;              % (rad/s)
% S.S=abs(Rper(1:(nf+1),1))*dt/pi;     % (m^2*s/rad) one-sided spectrum