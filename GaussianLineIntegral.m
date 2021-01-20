function sum = GaussianLineIntegral(fx,xa,ya,za,xb,yb,zb)

%Positive zeros of 12th order Legendre polynomial
A = [0.1252334085,  0.3678314989,  0.5873179542, ...
     0.7699026741,  0.9041172563,  0.9815606342];
%Gaussian Integration coefficients for a 12th order polynomial
WT =[0.2491470458,  0.2334925365,  0.2031674267, ...
     0.1600783285,  0.1069393259,  0.0471753363]; 
 
dx = (xb-xa)/2.0;
dy = (yb-ya)/2.0;
dz = (zb-za)/2.0;
xbar = (xa+xb)/2.0;
ybar = (ya+yb)/2.0;
zbar = (za+zb)/2.0;
% 
sum = 0;
% sumk = sum(fx(xbar+A.*dx).*WT + fx(xbar-A.*dx).*WT)/2;

for k=1:length(A)
    sum = sum +  WT(k).* ( fx(xbar+A(k).*dx, ybar+A(k).*dy, zbar+A(k).*dz) + ...
                           fx(xbar-A(k).*dx, ybar-A(k).*dy, zbar-A(k).*dz) );
end

sum = sum/2;
end