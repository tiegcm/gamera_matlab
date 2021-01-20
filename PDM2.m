function [f_left,f_right]=PDM2(f0,f1,f2,f3,f,PDMB)

% This subroutine splits one interface value into L/R states using the PDM
% limiter, which moniters the maximum amount of flux that can be pushed
% across the interface. The numerical diffusion is controlled by the
% parameter "PDMB", which is a function of the CFL condition. The detailed
% derivation of the PDM limiter algorithm is in Kain [1987], JCP:
%            https://doi.org/10.1016/0021-9991(87)90110-0.
%
% INPUT: Four values (f0,f1,f2,f3) alone a stencil (interface = 1/2)
%        High-order interpolated value (f)
%        PDMB: parameter controls numerical diffusion (PDMB ->0 then donor cell)

% first clipping the interpolated value to make sure it's between f1 and f2
maxf = max(f1,f2);
minf = min(f1,f2);
f = max(minf,min(f,maxf));

% the amount that can be pushed across interfaces
df0 = PDMB.*(f1-f0);
df1 = PDMB.*(f2-f1);
df2 = PDMB.*(f3-f2);

% based on the changes in gradient, determine whether the flux should be
% limited
s0 = sign(df0);
s1 = sign(df1);
s2 = sign(df2);

df0 = abs(df0);
df1 = abs(df1);
df2 = abs(df2);

q0 = abs(s0+s1);
q1 = abs(s1+s2);

df_left = f - f1;
df_righ = f2 - f;

% calculate the left- and right-state
f_left = f - s1.*max(0,abs(df_left) - q0.*df0);
f_right= f + s1.*max(0,abs(df_righ) - q1.*df2); 

end
