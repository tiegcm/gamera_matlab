function dt = getDT(rho,vx,vy,vz,bx_total,by_total,bz_total,p,gamma,di,dj,dk,CFL)

global CA ic_act jc_act kc_act

% calculate wave speeds for dt
Valfvn = sqrt(bx_total.^2+by_total.^2+bz_total.^2)./sqrt(rho); % Alfven speed VA
Va_eff = Valfvn.*CA./sqrt(Valfvn.^2+CA.^2); % Boris-corrected (limited) VA
Vsound = sqrt(gamma.*p./rho);               % sound speed
Vfluid = sqrt(vx.^2+vy.^2+vz.^2);           % fluid speed

% calculate dt
VCFL = Vfluid+sqrt(Va_eff.^2+Vsound.^2);    % total speed within cells
dtCFL = CFL ./(VCFL./di+VCFL./dj+VCFL./dk); % dt calculation based on the CFL condition

dt = min(min(min(dtCFL(ic_act,jc_act,kc_act)))); % only use the active region

end