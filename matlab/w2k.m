function k=w2k(w)
% W2K Translates from frequency to wave number
%     using the dispersion relation w = sqrt(g*k+0.074*k^3/1000)   (  0 <   w   < inf)
%
% CALL:  k=w2k(w)%
%     k   = wave numbers
%     w   = angular frequency   

p= 9.8*1000.0/0.074;
q= -1000.0*w.^2/0.074;
k=(-q/2+(q.^2/4+p.^3/27).^(1/2)).^(1/3)+... 
-(q/2+(q.^2/4+p.^3/27).^(1/2)).^(1/3); 


