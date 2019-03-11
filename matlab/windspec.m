function S=windspec(U10,x1,w_min,w_max,Nw,Ntheta)
% -------------------------------------------------------------------------
% calculate directional wave spectrum
% U10   is  wind speed
% x1    is  dimensionless wind fetch
% w_min is  minimal angular frequency
% w_max is  maximal angular frequency
% Nw    is  number of angular frequency
% Ntheta is number of angle
% -------------------------------------------------------------------------
km=calkm(U10,x1);
S=createspec('dir');  % undirectional wave spectrum
S.w=linspace(w_min,w_max,Nw);
S.theta=linspace(-pi,pi,Ntheta)';
S1=ModWavSpc(S.w,U10,x1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=-0.28+0.65*exp(-0.75*log(w2k(S.w)./km))+0.01*exp(-0.2+0.7*log10(w2k(S.w)./km));
B=10.^b;
A=B./atan(sinh(2*pi*B));
S2=zeros(Ntheta,Nw);
for ii= 1:Ntheta % the length of ¦Õ
    for jj=1:Nw  % Nw=length(KT)
        S2(ii,jj)=A(jj)*2./(exp(2*B(jj)*S.theta(ii))+exp(-2*B(jj)*S.theta(ii))); % spreading(¦Õ)
        S.S(ii,jj)=S1(jj)*S2(ii,jj);
    end
end
index=isnan(S.S);S.S(index)=0;
% figure,wspecplot(S)
if ~isempty(find(abs(simpson(S.theta,S2,1)-1)>0.01, 1))
    error ('Angular distribution is not normalized.')
end
return
