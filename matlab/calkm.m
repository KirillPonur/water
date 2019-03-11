function KM=calkm(U10,x1)
% CALKM for calulate km from U10 and x1
if x1<1430
    error('dimensionless wind fetch must bigger than 1430')
elseif  x1>20170 || x1==20170
    w_dev= 0.835;
else 
    w_dev=0.61826357843576103 + 3.52883010586243843e-06*x1...
                          - 0.00197508032233982112*sqrt(x1)...
                          + 62.5540113059129759/sqrt(x1)...
                           - 290.214120684236224/x1;                    
end
w_m=w_dev*9.8/U10;
KM=(w_m)^2/9.8;