function [f_left]=PDMU7(f1,f2,f3,f4,f5,f_itp,PDMB,clipping)

% This subroutine is basically the same as the PDM2() function, but only
% computes the left-state based on the shifted, upwind interpolation
% The limiter step is the same as the default PDM algorithm, whichmoniters 
% the maximum amount of flux that can be pushed across the interface. 
% The numerical diffusion is controlled by the parameter "PDMB", which is 
% a function of the CFL condition. See Kain [1987], JCP for details.
%            https://doi.org/10.1016/0021-9991(87)90110-0.
%
% INPUT: (f0,f1,f2,f3): Four values alone a stencil (interface = 1/2)
%        (f)          : High-order interpolated value 
%        PDMB:        : parameter controls numerical diffusion (PDMB ->0 then donor cell)
%        clipping     : (optional) non-clipping choice (default off)

if (nargin<8)
    clipping=0; % default non-clipping off
end

% clipping
maxf = max(f3,f4);
minf = min(f3,f4);
f = max(minf,min(f_itp,maxf));

% find the PDM value for left-state
df0 = PDMB.*(f3-f2);
df1 = PDMB.*(f4-f3);

s0 = sign(df0);
s1 = sign(df1);
df0 = abs(df0);
q0 = abs(s0+s1);

df_left = f - f3;
f_pdm = f - s1.*max(0,abs(df_left) - q0.*df0);

if(clipping==0)
    isSmooth = 0; % smoothness indicator, default = 0
else
    % compute the deltas
    D2 = f2 - f1;
    D3 = f3 - f2;
    D4 = f4 - f3;
    D5 = f5 - f4;
    % compute the local extremum indicators (LEI)
    SmC=(D2>0)&(D3>0)&(D4<0)&(D5<0)&(abs(D3)<abs(D2))&(abs(D4)<abs(D5)); % locap peak
    SpC=(D2<0)&(D3<0)&(D4>0)&(D5>0)&(abs(D3)<abs(D2))&(abs(D4)<abs(D5)); % local dip
    % reset the smoothness indicator    
    isSmooth = sign(SmC+SpC);
end

% adjust the solution based on the local smoothness indicator
f_left = f_pdm.*(1-isSmooth) + isSmooth.*f_itp;

end
