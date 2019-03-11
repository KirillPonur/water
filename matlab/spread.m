function Sp=spread(w,phi,U10,x1)  
% angular distribution ¦Õ£¨w£¬¦Õ£©
% w is angular frequency
% phi-->¦Õ is angular
% U10 is wind speed
% x1 is dimensionless wind fetch
%---------------------------------------------------------
km=calkm(U10,x1);  
b=-0.28+0.65*exp(-0.75*log(w2k(w)./km))+0.01*exp(-0.2+0.7*log10(w2k(w)./km));
B=10.^b;
A=B./atan(sinh(2*pi*B));
Sp=A*2./(exp(2*B*phi)+exp(-2*B*phi)); % spreading(¦Õ)
index=isnan(Sp);Sp(index)=0;
return 