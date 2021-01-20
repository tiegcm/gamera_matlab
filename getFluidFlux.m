function [Frho_p, FrhoVx_p,FrhoVy_p,FrhoVz_p,Feng_p, ...
          Frho_n, FrhoVx_n,FrhoVy_n,FrhoVz_n,Feng_n] = getFluidFlux(rho,Vx,Vy,Vz,p,gamma)

% This subroutine computes the interface flux based on Gas-Kinetic flux
% functions. The detailed derivations can be found in Xu et al. [1998], JCP
% The interface fluxes are basically calculated from the moment integrals of 
% a distribution function described by macroscopic fluid variables:
%                      1st-moment : mass flux
%                      2nd-moment : momentum flux
%                      3rd-moment : energy flux
% The normal velocity here is assumed to be the x-component, y- and z-components 
% are tengential components. 
%
% INPUT  : rho, vx (interface-normal), vy, vz, p on either side of the face
% OUTPUT : the positive (v>0) and the negative (v<0) moment integral
 
    eng= 0.5.*rho.*(Vx.^2+Vy.^2+Vz.^2)+p./(gamma-1); % plasma energy
    lamda = rho./(2*p);  % temperature - the width of a Maxwellian

    Vx0_p = 0.5*erfc(-sqrt(lamda).*Vx); % zeroth-moment integral for v>0
    Vx0_n = 0.5*erfc(+sqrt(lamda).*Vx); % zeroth-moment integral for v<0
    Vx1_p = Vx.*Vx0_p + 0.5.*exp(-lamda.*Vx.^2)./sqrt(pi.*lamda); % 1st-moment integral for v>0
    Vx1_n = Vx.*Vx0_n - 0.5.*exp(-lamda.*Vx.^2)./sqrt(pi.*lamda); % 1st-moment integral for v<0
   
    % Flux in the positive x direction (use positive moment integrals)
    Frho_p   = rho  .* Vx1_p;
    FrhoVx_p = rho.*Vx .* Vx1_p + p .* Vx0_p;
    FrhoVy_p = rho.*Vy .* Vx1_p;
    FrhoVz_p = rho.*Vz .* Vx1_p;
    Feng_p   = (eng + 0.5*p) .* Vx1_p + (0.5*p.*Vx).*Vx0_p;
    
    % Flux in the negative x direction (use negative moment integrals)
    Frho_n   = rho  .* Vx1_n;
    FrhoVx_n = rho.*Vx .* Vx1_n + p .* Vx0_n;
    FrhoVy_n = rho.*Vy .* Vx1_n;
    FrhoVz_n = rho.*Vz .* Vx1_n;
    Feng_n   = (eng + 0.5*p) .* Vx1_n + (0.5*p.*Vx).*Vx0_n;
    
end
