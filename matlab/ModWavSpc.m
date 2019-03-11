function Spc=ModWavSpc(w,U10,x1)
%--------------------------
% this function used for calculating 1D wave spectrum S(w)
% w is angular frequency, 1D or 2D matrix
% x1 is dimensionless wind fetch
% U10 is a number----wind speed in 10m
%------------------------------
index=(w==0);w(index)=eps;
gravity=9.81
if x1<1430
    error('dimensionless wind fetch must bigger than 1430')
elseif  x1>20170 || x1==20170
    w_dev= 0.835;alpha=0.0081;gamma=1;
else 
    w_dev=0.61826357843576103 + 3.52883010586243843e-06*x1...
                          - 0.00197508032233982112*sqrt(x1)...
                          + 62.5540113059129759/sqrt(x1)...
                           - 290.214120684236224/x1;
    gamma=5.25366092888284494 + 0.000107621806575821253*x1...
                         - 0.0377677617998704552*sqrt(x1)...
                         - 162.983465281012658/sqrt(x1)...
                         + 253251.456471558638/x1^1.5;
    alpha=0.0311936666714071662 - 0.0023277357121814146*log(x1)...
                          - 8367.86787861323644/x1^2 ...
                         + 4.51145991018772117*10^(617-x1/log(10));                    
end
w_m=w_dev*gravity/U10;

a_m=0.3713+0.29024*U10+0.2902/U10;
Spc=zeros(size(w));

if 1.2*w_m==1
   spec_w12=alpha*gravity^2*(1.2*w_m)^(-5)*exp(-1.25/(1.2)^4)*gamma^(exp(-(1.2*w_m-w_m)^2/(2*0.09^2*w_m^2)));
else
   spec_w12=alpha*gravity^2*(1.2*w_m)^(-5-1.25*(w_m/(1.2*w_m))^4/log(1.2*w_m))*gamma^(exp(-(1.2*w_m-w_m)^2/(2*0.09^2*w_m^2)));
end

for t1=1:size(w,1)
    for t2=1:size(w,2)
        if w(t1,t2)<1.2*w_m || w(t1,t2)==1.2*w_m
            if w(t1,t2)==1;
               Spc(t1,t2)=alpha*gravity^2*w(t1,t2)^(-5)*exp(-1.25*(w_m/w(t1,t2))^4)*gamma^(exp(-(w(t1,t2)-w_m)^2/(2*(0.07*(w(t1,t2)<w_m | w(t1,t2)==w_m)+0.09*(w(t1,t2)>w_m))^2*w_m^2)));
            else
               Spc(t1,t2)=alpha*gravity^2*w(t1,t2)^(-5-1.25*(w_m/w(t1,t2))^4/log(w(t1,t2)))*gamma^(exp(-(w(t1,t2)-w_m)^2/(2*(0.07*(w(t1,t2)<w_m | w(t1,t2)==w_m)+0.09*(w(t1,t2)>w_m))^2*w_m^2)));
            end
        elseif w(t1,t2)>1.2*w_m && (w(t1,t2)<a_m*w_m || w(t1,t2)==a_m*w_m)
            alpha2=spec_w12*(1.2*w_m)^4;
            Spc(t1,t2)=alpha2/w(t1,t2)^4;
        elseif w(t1,t2)>a_m*w_m && (w(t1,t2)<64 || w(t1,t2)==64)
            alpha3=spec_w12*(1.2*w_m)^4*a_m*w_m;
            Spc(t1,t2)=alpha3/w(t1,t2)^5;
        elseif w(t1,t2)>64 && (w(t1,t2)<298 || w(t1,t2)==298)
            alpha4=spec_w12*(1.2*w_m)^4*a_m*w_m/(64^2.3);
            Spc(t1,t2)=alpha4/w(t1,t2)^2.7;
        else
            alpha5=spec_w12*(1.2*w_m)^4*a_m*w_m/(64^2.3)*298^2.3;
            Spc(t1,t2)=alpha5/w(t1,t2)^5;
        end        
    end
end
return