function [rhovx,rhovy,rhovz,eng] = getConservedVariables(rho,vx,vy,vz,p,gamma)

% from primitive to conserved variables
rhovx = rho.*vx; % x-momentum
rhovy = rho.*vy; % y-momentum
rhovz = rho.*vz; % z-momentum
eng= 0.5.*rho.*(vx.^2+vy.^2+vz.^2)+p./(gamma-1); % plasma energy

end
