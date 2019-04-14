% Testing wave simulations of FFT method
%              date: 2009-4-1
%--------------------------------------------------------------------------
clear all
U10=5;
x1=20170;
LL=0.1;   % radar wavelength (L = 0.008, 0.021, 0.03, 0.055, 0.1, 0.23m)
KM=calkm(U10,x1); 
LM=2*pi/KM,
Kmin=KM/4, 
Kmax= KL_LIMIT(U10,LL)/2, 

Nw=2^4;
Ntheta=2^4; 
Sd=windspec(U10,x1,k2w(Kmin),k2w(Kmax),Nw,Ntheta); % S(¦Ø,¦Õ)
    % figure,wspecplot(Sd), % plot Sd
Nx=2^10;Ny=2^10;Nt=2^0;dt=1;
dx=pi/Kmax
dy=dx; 
Spec1=zeros(10,Ny,Nx);SWH=4*sqrt(spec2mom(Sd,0))
for i=1:10                              % simulation 10 times 
    H=seasim(Sd,Nx,Ny,Nt,dx,dy,dt,2,0);
        % HeightVariance=var(reshape(H.Z,Nx*Ny,1))
        % [fx,fy]=gradient(H.Z);fx=fx/dx;fy=fy/dy;
        % SlopVar_x=var(reshape(fx,Nx*Ny,1))
        % SlopVar_y=var(reshape(fy,Nx*Ny,1))
        % randn('state',1);
        % Hs=spectrum.cov;
        % psd(Hs,H.Z,'NFFT',512)
    Spec1(i,:,:)=(abs(fft2(H.Z-mean(mean(H.Z))))).^2;
end
%--------------------------------------------------------------------------
%        show sea surface 
%--------------------------------------------------------------------------
figure
colormap('winter')
surf(H.x,H.y,H.Z)
shading interp
axis([H.x(1) H.x(end) H.y(1) H.y(end) min(min(H.Z)) max(max(H.Z))])
axis('square')
view([0 90])
colorbar
saveas(gca,'TestSimu_seasim_surface','fig')  % save figure
%--------------------------------------------------------------------------
%     compare spectrum with input spectrum 
%--------------------------------------------------------------------------

figure,Snew = spec2spec(Sd,'k2d'); % S(¦Ø,¦Õ)--> F(Kx,Ky)
subplot(1,2,1),
contour(Snew.k,Snew.k2,Snew.S,5),
caxis([min(min(Snew.S)) max(max(Snew.S))]),
colorbar('location','southoutside');
xlabel('kx'),ylabel('ky'),title('Input wave spectrum')
axis([-2*KM 2*KM -2*KM 2*KM]),axis square
grid on

subplot(1,2,2),
Spec=fftshift(squeeze(mean(Spec1,1)*dx*dy/(Nx*Ny*2*pi*2*pi))); 
kx=2*pi*(-Nx/2:Nx/2-1)/Nx/dx;
ky=2*pi*(-Ny/2:Ny/2-1)'/Ny/dy;        

h = fspecial('gaussian',[3,3],1.5); 
Spec=filter2(h,Spec);

%         [in1,in2]=find(Spec==max(max(Spec)));
%         Km=sqrt((kx(in2)).^2+(ky(in1)).^2)
%         Lm=2*pi/km

contour(kx,ky,Spec,5)
caxis([min(min(Snew.S))/2 max(max(Snew.S))/2]),
colorbar('location','southoutside');
xlabel('kx'),ylabel('ky'),title('Wave spectrum from simulation')
axis([-2*KM 2*KM -2*KM 2*KM]),axis square
grid on
saveas(gca,'TestSimu_seasim_ComSpec','fig')  % save figure

simpson(Sd.w,simpson(Sd.theta,Sd.S,1))
simpson(Snew.k,simpson(Snew.k2,Snew.S,1))
%--------------------------------------------------------------------------
% Results£º
%--------------------------------------------------------------------------

% LM =
%   22.989040764296604
% Kmin =
%    0.068328050000000
% calling KL_LIMIT function
% Kmax =
%    9.919244084919804
% dx =
%    0.316716942006291
% SWH =
%    0.576109682690251
% ans =
%    0.020743897905591
% ans =
%    0.020749879584320
