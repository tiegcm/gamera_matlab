function [drho,dpxF,dpyF,dpzF,deng,dpxB,dpyB,dpzB]=Hydro(rho_h,vx_h,vy_h,vz_h,p_h,bi_h,bj_h,bk_h,bx_h,by_h,bz_h,TRANS_i,TRANS_j,TRANS_k,...
                                                         fxi,fyi,fzi,fnormi,fxj,fyj,fzj,fnormj,fxk,fyk,fzk,fnormk,...
                                                         i_face_area,j_face_area,k_face_area,...
                                                         CA,Va_eff,volume,PDMB,limiter,gamma,n,dt)

% This subroutine computes the volume integrated fluid stresses: 
%                     drho: volume-integrated mass flux
%                     dpxF: volume-integrated x-momentum flux 
%                     dpyF: volume-integrated y-momentum flux 
%                     dpzF: volume-integrated z-momentum flux 
%                     deng: volume-integrated energy flux 
% and magnetic stresses:
%                     dpxB: volume-integrated x-component of Lorentz force 
%                     dpyB: volume-integrated y-component of Lorentz force 
%                     dpzB: volume-integrated z-component of Lorentz force 
%
% MAIN ALGORITHMS 
%   Step 1. Start from the i-direction
%   Step 2. Reconstruct along stencils to get left- and right-states
%   Step 3. Transform left- and right-velocity into face-normal coordinates using the TRANS_ matrix
%   Step 4. Evaluate interface fluid fluxes using Gas-Kinetic flux functions (including diffusive correction - HOGS flux)
%   Step 5. Inverse transform momentum stresses back to the (x,y,z) system
%   Step 6. Evaluate interface magnetic stresses using magneto-Gas-Kinetic flux functions (diffusive correction - HOGS flux - only need in Boris
%   Step 7. Loop over j-direction and k-direction
%   Step 8. Face-integrals for the flux terms for finite volume update
%
% NOTE: 1. The reconstruction algorithm applies to the conserved variables,
%          e.g., reconstruct mass and then divided by interface volume to
%          get high-order estimations for interface density;
%       2. The reconstruction of interface magnetic fields follow similar
%          ideas, e.g., reconstruct bx*volume and then divide by interface volume
%       3. The HOGS flux are only needed when using the semi-conservative
%          MHD equations; not needed for the fully-conservative equations.

global ic_act jc_act kc_act if_act jf_act kf_act

% Step 0. preliminary step
    % estimate the interface volume for reconstruction in non-uniform grids
        [vol_left, vol_right] = reconstruct_3DV(volume,if_act,jf_act,kf_act,PDMB,1,'8th');  
        vol_iface = 0.5.*(vol_left + vol_right);
        [vol_left, vol_right] = reconstruct_3DV(volume,if_act,jf_act,kf_act,PDMB,2,'8th');  
        vol_jface = 0.5.*(vol_left + vol_right);
        [vol_left, vol_right] = reconstruct_3DV(volume,if_act,jf_act,kf_act,PDMB,3,'8th');  
        vol_kface = 0.5.*(vol_left + vol_right);
        
% Step 1. Starting from the i-direction   
% Step 2. reconstruct fluid variables in the i-direction
    [rho_left, rho_right] = reconstruct_3DV(rho_h,if_act,jf_act,kf_act,PDMB,1,limiter,volume,vol_iface);  
    [vx_left, vx_right] = reconstruct_3DV(vx_h,if_act,jf_act,kf_act,PDMB,1,limiter,volume,vol_iface);   
    [vy_left, vy_right] = reconstruct_3DV(vy_h,if_act,jf_act,kf_act,PDMB,1,limiter,volume,vol_iface);       
    [vz_left, vz_right] = reconstruct_3DV(vz_h,if_act,jf_act,kf_act,PDMB,1,limiter,volume,vol_iface);        
    [p_left, p_right] = reconstruct_3DV(p_h,if_act,jf_act,kf_act,PDMB,1,limiter,volume,vol_iface);   

% Step 3. Transform left- and right-velocity into face-normal coordinates     
    % transform the interface velocities (vx_, vy_, vj_) to the face-normal 
    % coordinate systems (vi_ vj_ vk_). NOTE vi_ is the interface normal velocity
    % left-state velocities
    vi_left(if_act,jc_act,kc_act) = vx_left(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,1) + ...
                                    vy_left(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,2) + ...
                                    vz_left(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,3);
    vj_left(if_act,jc_act,kc_act) = vx_left(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,4) + ...
                                    vy_left(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,5) + ...
                                    vz_left(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,6);
    vk_left(if_act,jc_act,kc_act) = vx_left(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,7) + ...
                                    vy_left(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,8) + ...
                                    vz_left(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,9);
    % right-state velocities
    vi_right(if_act,jc_act,kc_act) = vx_right(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,1) + ...
                                     vy_right(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,2) + ...
                                     vz_right(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,3);
    vj_right(if_act,jc_act,kc_act) = vx_right(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,4) + ...
                                     vy_right(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,5) + ...
                                     vz_right(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,6);
    vk_right(if_act,jc_act,kc_act) = vx_right(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,7) + ...
                                     vy_right(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,8) + ...
                                     vz_right(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,9);

% Step 4. Evaluate the fluid fluxes at i-interfaces, using the Gas-Kinetic flux functions  
    % right-integral of the left distribution function
    [Frho_p(if_act,jc_act,kc_act),FrhoVx_p(if_act,jc_act,kc_act),FrhoVy_p(if_act,jc_act,kc_act),FrhoVz_p(if_act,jc_act,kc_act),Feng_p(if_act,jc_act,kc_act), ...
     ~,~,~,~,~] = getFluidFlux(rho_left(if_act,jc_act,kc_act),vi_left(if_act,jc_act,kc_act),vj_left(if_act,jc_act,kc_act),vk_left(if_act,jc_act,kc_act),p_left(if_act,jc_act,kc_act),gamma);
    % left-integral of the right distribution function
    [~,~,~,~,~, ...
     Frho_n(if_act,jc_act,kc_act),FrhoVx_n(if_act,jc_act,kc_act),FrhoVy_n(if_act,jc_act,kc_act),FrhoVz_n(if_act,jc_act,kc_act),Feng_n(if_act,jc_act,kc_act)] = ...
                  getFluidFlux(rho_right(if_act,jc_act,kc_act),vi_right(if_act,jc_act,kc_act),vj_right(if_act,jc_act,kc_act),vk_right(if_act,jc_act,kc_act),p_right(if_act,jc_act,kc_act),gamma); 
    % sum the contributions from the two interface distributions
    rho_flux_i = Frho_p + Frho_n;
    vx_flux_x = FrhoVx_p + FrhoVx_n;
    vy_flux_x = FrhoVy_p + FrhoVy_n;
    vz_flux_x = FrhoVz_p + FrhoVz_n;
    eng_flux_i = Feng_p + Feng_n;              
    
    % add the HOGS fix to the interface fluid flux in i-direction. The HOGS
    % flux is needed because the diffusion speed from the pure Gas-Kinetic
    % flux function is only related to the entropy and thermal speeds,
    % which is not enough for magnetically-dominated regions with beta<<1.
    % The HOGS flux is basically a local LxF flux with a diffusive speed VA
    dif_factor = 0.25; 
    diff_va = 0*Va_eff;
    diff_va(if_act,jc_act,kc_act) = dif_factor*(Va_eff(if_act-1,jc_act,kc_act) + Va_eff(if_act,jc_act,kc_act));
    rho_flux_i(if_act,jc_act,kc_act) = rho_flux_i(if_act,jc_act,kc_act) - diff_va(if_act,jc_act,kc_act).*(rho_right(if_act,jc_act,kc_act) - rho_left(if_act,jc_act,kc_act));
    vx_flux_x(if_act,jc_act,kc_act)  = vx_flux_x(if_act,jc_act,kc_act)  - diff_va(if_act,jc_act,kc_act).*(rho_right(if_act,jc_act,kc_act).*vi_right(if_act,jc_act,kc_act) - rho_left(if_act,jc_act,kc_act).*vi_left(if_act,jc_act,kc_act));
    vy_flux_x(if_act,jc_act,kc_act)  = vy_flux_x(if_act,jc_act,kc_act)  - diff_va(if_act,jc_act,kc_act).*(rho_right(if_act,jc_act,kc_act).*vj_right(if_act,jc_act,kc_act) - rho_left(if_act,jc_act,kc_act).*vj_left(if_act,jc_act,kc_act));
    vz_flux_x(if_act,jc_act,kc_act)  = vz_flux_x(if_act,jc_act,kc_act)  - diff_va(if_act,jc_act,kc_act).*(rho_right(if_act,jc_act,kc_act).*vk_right(if_act,jc_act,kc_act) - rho_left(if_act,jc_act,kc_act).*vk_left(if_act,jc_act,kc_act));
    eng_flux_i(if_act,jc_act,kc_act) = eng_flux_i(if_act,jc_act,kc_act) - diff_va(if_act,jc_act,kc_act).*(0.5.*rho_right(if_act,jc_act,kc_act).*(vi_right(if_act,jc_act,kc_act).^2+vj_right(if_act,jc_act,kc_act).^2+vk_right(if_act,jc_act,kc_act).^2) + p_right(if_act,jc_act,kc_act)./(gamma-1) ...
                                                                                                        - 0.5.*rho_left(if_act,jc_act,kc_act).*(vi_left(if_act,jc_act,kc_act).^2+vj_left(if_act,jc_act,kc_act).^2+vk_left(if_act,jc_act,kc_act).^2) - p_left(if_act,jc_act,kc_act)./(gamma-1));
    
% Step 5. inverse transform the momentum fluxes from face-normal (i,j,k) back to base (x,y,z)
    vx_flux_i(if_act,jc_act,kc_act) = vx_flux_x(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,1) + ...
                                      vy_flux_x(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,4) + ...
                                      vz_flux_x(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,7);
    vy_flux_i(if_act,jc_act,kc_act) = vx_flux_x(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,2) + ...
                                      vy_flux_x(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,5) + ...
                                      vz_flux_x(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,8);
    vz_flux_i(if_act,jc_act,kc_act) = vx_flux_x(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,3) + ...
                                      vy_flux_x(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,6) + ...
                                      vz_flux_x(if_act,jc_act,kc_act).*TRANS_i(if_act,jc_act,kc_act,9);

% Step 6. Evaluate interface magnetic stresses using magneto-Gas-Kinetic flux 
% functions (diffusive correction - HOGS flux - only need when the Boris
% correction is significant)
   % (this is a part of step 2) reconstruct left- and right- interface
   % magnetic fields, similar as fluid variables considering volume effect
    [bx_left, bx_right] = reconstruct_3DV(bx_h,if_act,jf_act,kf_act,PDMB,1,limiter,volume,vol_iface);   
    [by_left, by_right] = reconstruct_3DV(by_h,if_act,jf_act,kf_act,PDMB,1,limiter,volume,vol_iface);   
    [bz_left, bz_right] = reconstruct_3DV(bz_h,if_act,jf_act,kf_act,PDMB,1,limiter,volume,vol_iface);   
    
    % magneto-Gas-Kinetic flux evaluation, similar to the fluid flux function
    bnormi = bi_h./i_face_area;% normal field strength needed in the flux function
    % integrating v>0 for the left distribution function
    [Bstress_x_p(if_act,jc_act,kc_act), Bstress_y_p(if_act,jc_act,kc_act), Bstress_z_p(if_act,jc_act,kc_act), ...
     ~, ~, ~] = getMagneticStress(rho_left(if_act,jc_act,kc_act),vi_left(if_act,jc_act,kc_act),p_left(if_act,jc_act,kc_act),...
                                           bx_left(if_act,jc_act,kc_act),by_left(if_act,jc_act,kc_act),bz_left(if_act,jc_act,kc_act),bnormi(if_act,jc_act,kc_act),...
                                           fxi(if_act,jc_act,kc_act),fyi(if_act,jc_act,kc_act),fzi(if_act,jc_act,kc_act),fnormi(if_act,jc_act,kc_act),CA,...
                                           TRANS_i(if_act,jc_act,kc_act,1),TRANS_i(if_act,jc_act,kc_act,2),TRANS_i(if_act,jc_act,kc_act,3));
    % integrating v<0 for the right distribution function
     [~, ~, ~, ...
     Bstress_x_n(if_act,jc_act,kc_act), Bstress_y_n(if_act,jc_act,kc_act), Bstress_z_n(if_act,jc_act,kc_act)] ...
              = getMagneticStress(rho_right(if_act,jc_act,kc_act),vi_right(if_act,jc_act,kc_act),p_right(if_act,jc_act,kc_act),...
                                           bx_right(if_act,jc_act,kc_act),by_right(if_act,jc_act,kc_act),bz_right(if_act,jc_act,kc_act),bnormi(if_act,jc_act,kc_act),...
                                           fxi(if_act,jc_act,kc_act),fyi(if_act,jc_act,kc_act),fzi(if_act,jc_act,kc_act),fnormi(if_act,jc_act,kc_act),CA,...
                                           TRANS_i(if_act,jc_act,kc_act,1),TRANS_i(if_act,jc_act,kc_act,2),TRANS_i(if_act,jc_act,kc_act,3));
   
    % sum the contribution from the two interface distributions
    BstressX_i = Bstress_x_p + Bstress_x_n;
    BstressY_i = Bstress_y_p + Bstress_y_n;
    BstressZ_i = Bstress_z_p + Bstress_z_n; 
    
    % now add the magnetic HOGS flux (if CA>>1, then HOGS flux == 0)
    bsqi(if_act,jc_act,kc_act) = 0.5*((bx_left(if_act,jc_act,kc_act)+fxi(if_act,jc_act,kc_act)).^2+(by_left(if_act,jc_act,kc_act)+fyi(if_act,jc_act,kc_act)).^2+(bz_left(if_act,jc_act,kc_act)+fzi(if_act,jc_act,kc_act)).^2 +...
                                      (bx_right(if_act,jc_act,kc_act)+fxi(if_act,jc_act,kc_act)).^2+(by_right(if_act,jc_act,kc_act)+fyi(if_act,jc_act,kc_act)).^2+(bz_right(if_act,jc_act,kc_act)+fzi(if_act,jc_act,kc_act)).^2)/2/CA;
    BstressX_i(if_act,jc_act,kc_act) = BstressX_i(if_act,jc_act,kc_act) - bsqi(if_act,jc_act,kc_act).*(vx_right(if_act,jc_act,kc_act) - vx_left(if_act,jc_act,kc_act));
    BstressY_i(if_act,jc_act,kc_act) = BstressY_i(if_act,jc_act,kc_act) - bsqi(if_act,jc_act,kc_act).*(vy_right(if_act,jc_act,kc_act) - vy_left(if_act,jc_act,kc_act));
    BstressZ_i(if_act,jc_act,kc_act) = BstressZ_i(if_act,jc_act,kc_act) - bsqi(if_act,jc_act,kc_act).*(vz_right(if_act,jc_act,kc_act) - vz_left(if_act,jc_act,kc_act));

    %~~~~~~~~~END OF I-DIRECTION FLUX CALCULATIONS~~~~~~~~~
    
%   Step 7. Now do the j-direction flux calculations
    % reconstruct in the j-direction 
    [rho_left, rho_right] = reconstruct_3DV(rho_h,if_act,jf_act,kf_act,PDMB,2,limiter,volume,vol_jface);  
    [vx_left, vx_right] = reconstruct_3DV(vx_h,if_act,jf_act,kf_act,PDMB,2,limiter,volume,vol_jface);    
    [vy_left, vy_right] = reconstruct_3DV(vy_h,if_act,jf_act,kf_act,PDMB,2,limiter,volume,vol_jface);         
    [vz_left, vz_right] = reconstruct_3DV(vz_h,if_act,jf_act,kf_act,PDMB,2,limiter,volume,vol_jface);         
    [p_left, p_right] = reconstruct_3DV(p_h,if_act,jf_act,kf_act,PDMB,2,limiter,volume,vol_jface);    
          
    % transform the velocities to i,j,k coordinate
    vi_left(ic_act,jf_act,kc_act) = vx_left(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,1) + ...
                                    vy_left(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,2) + ...
                                    vz_left(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,3);
    vj_left(ic_act,jf_act,kc_act) = vx_left(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,4) + ...
                                    vy_left(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,5) + ...
                                    vz_left(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,6);
    vk_left(ic_act,jf_act,kc_act) = vx_left(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,7) + ...
                                    vy_left(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,8) + ...
                                    vz_left(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,9);

    vi_right(ic_act,jf_act,kc_act) = vx_right(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,1) + ...
                                     vy_right(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,2) + ...
                                     vz_right(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,3);
    vj_right(ic_act,jf_act,kc_act) = vx_right(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,4) + ...
                                     vy_right(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,5) + ...
                                     vz_right(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,6);
    vk_right(ic_act,jf_act,kc_act) = vx_right(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,7) + ...
                                     vy_right(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,8) + ...
                                     vz_right(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,9);

    % Evaluate the i-interface fluxes
    [Frho_p(ic_act,jf_act,kc_act),FrhoVx_p(ic_act,jf_act,kc_act),FrhoVy_p(ic_act,jf_act,kc_act),FrhoVz_p(ic_act,jf_act,kc_act),Feng_p(ic_act,jf_act,kc_act), ...
     ~,~,~,~,~] = getFluidFlux(rho_left(ic_act,jf_act,kc_act),vi_left(ic_act,jf_act,kc_act),vj_left(ic_act,jf_act,kc_act),vk_left(ic_act,jf_act,kc_act),p_left(ic_act,jf_act,kc_act),gamma);
    [~,~,~,~,~, ...
     Frho_n(ic_act,jf_act,kc_act),FrhoVx_n(ic_act,jf_act,kc_act),FrhoVy_n(ic_act,jf_act,kc_act),FrhoVz_n(ic_act,jf_act,kc_act),Feng_n(ic_act,jf_act,kc_act)] ...
     = getFluidFlux(rho_right(ic_act,jf_act,kc_act),vi_right(ic_act,jf_act,kc_act),vj_right(ic_act,jf_act,kc_act),vk_right(ic_act,jf_act,kc_act),p_right(ic_act,jf_act,kc_act),gamma); 
 
    rho_flux_j = Frho_p + Frho_n;
    vx_flux_x = FrhoVx_p + FrhoVx_n;
    vy_flux_x = FrhoVy_p + FrhoVy_n;
    vz_flux_x = FrhoVz_p + FrhoVz_n;
    eng_flux_j = Feng_p + Feng_n;              

    % ~~~~~~~~~~~~~add the HOGS fix, j-direction
    diff_va = 0*Va_eff;
    diff_va(ic_act,jf_act,kc_act) = dif_factor*(Va_eff(ic_act,jf_act-1,kc_act) + Va_eff(ic_act,jf_act,kc_act));
    rho_flux_j(ic_act,jf_act,kc_act) = rho_flux_j(ic_act,jf_act,kc_act) - diff_va(ic_act,jf_act,kc_act).*(rho_right(ic_act,jf_act,kc_act) - rho_left(ic_act,jf_act,kc_act));
    vx_flux_x(ic_act,jf_act,kc_act)  = vx_flux_x(ic_act,jf_act,kc_act)  - diff_va(ic_act,jf_act,kc_act).*(rho_right(ic_act,jf_act,kc_act).*vi_right(ic_act,jf_act,kc_act) - rho_left(ic_act,jf_act,kc_act).*vi_left(ic_act,jf_act,kc_act));
    vy_flux_x(ic_act,jf_act,kc_act)  = vy_flux_x(ic_act,jf_act,kc_act)  - diff_va(ic_act,jf_act,kc_act).*(rho_right(ic_act,jf_act,kc_act).*vj_right(ic_act,jf_act,kc_act) - rho_left(ic_act,jf_act,kc_act).*vj_left(ic_act,jf_act,kc_act));
    vz_flux_x(ic_act,jf_act,kc_act)  = vz_flux_x(ic_act,jf_act,kc_act)  - diff_va(ic_act,jf_act,kc_act).*(rho_right(ic_act,jf_act,kc_act).*vk_right(ic_act,jf_act,kc_act) - rho_left(ic_act,jf_act,kc_act).*vk_left(ic_act,jf_act,kc_act));
    eng_flux_j(ic_act,jf_act,kc_act) = eng_flux_j(ic_act,jf_act,kc_act) - diff_va(ic_act,jf_act,kc_act).*(0.5.*rho_right(ic_act,jf_act,kc_act).*(vi_right(ic_act,jf_act,kc_act).^2+vj_right(ic_act,jf_act,kc_act).^2+vk_right(ic_act,jf_act,kc_act).^2) + p_right(ic_act,jf_act,kc_act)./(gamma-1) ...
                                                                                                        - 0.5.*rho_left(ic_act,jf_act,kc_act).*(vi_left(ic_act,jf_act,kc_act).^2+vj_left(ic_act,jf_act,kc_act).^2+vk_left(ic_act,jf_act,kc_act).^2) - p_left(ic_act,jf_act,kc_act)./(gamma-1));
    
    % transform the momentum fluxes back to x,y,z
    vx_flux_j(ic_act,jf_act,kc_act) = vx_flux_x(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,1) + ...
                                      vy_flux_x(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,4) + ...
                                      vz_flux_x(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,7);
    vy_flux_j(ic_act,jf_act,kc_act) = vx_flux_x(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,2) + ...
                                      vy_flux_x(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,5) + ...
                                      vz_flux_x(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,8);
    vz_flux_j(ic_act,jf_act,kc_act) = vx_flux_x(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,3) + ...
                                      vy_flux_x(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,6) + ...
                                      vz_flux_x(ic_act,jf_act,kc_act).*TRANS_j(ic_act,jf_act,kc_act,9);
                  
   % calculate magnetic stress in the j direction
    [bx_left, bx_right] = reconstruct_3DV(bx_h,if_act,jf_act,kf_act,PDMB,2,limiter,volume,vol_jface);  
    [by_left, by_right] = reconstruct_3DV(by_h,if_act,jf_act,kf_act,PDMB,2,limiter,volume,vol_jface);   
    [bz_left, bz_right] = reconstruct_3DV(bz_h,if_act,jf_act,kf_act,PDMB,2,limiter,volume,vol_jface);   
 
    bnormj = bj_h./j_face_area;    
    % getMagneticStress(rho,Vx,p,bx,by,bz,bnorm,bx0,by0,bz0,b0n,CA,trans_1,trans_2,trans_3)
    [Bstress_x_p(ic_act,jf_act,kc_act), Bstress_y_p(ic_act,jf_act,kc_act), Bstress_z_p(ic_act,jf_act,kc_act), ...
     ~, ~, ~] = getMagneticStress(rho_left(ic_act,jf_act,kc_act),vi_left(ic_act,jf_act,kc_act),p_left(ic_act,jf_act,kc_act),...
                                           bx_left(ic_act,jf_act,kc_act),by_left(ic_act,jf_act,kc_act),bz_left(ic_act,jf_act,kc_act),bnormj(ic_act,jf_act,kc_act),...
                                           fxj(ic_act,jf_act,kc_act),fyj(ic_act,jf_act,kc_act),fzj(ic_act,jf_act,kc_act),fnormj(ic_act,jf_act,kc_act),CA,...
                                           TRANS_j(ic_act,jf_act,kc_act,1),TRANS_j(ic_act,jf_act,kc_act,2),TRANS_j(ic_act,jf_act,kc_act,3));
    
    [~, ~, ~, ...
     Bstress_x_n(ic_act,jf_act,kc_act), Bstress_y_n(ic_act,jf_act,kc_act), Bstress_z_n(ic_act,jf_act,kc_act)] ...
              = getMagneticStress(rho_right(ic_act,jf_act,kc_act),vi_right(ic_act,jf_act,kc_act),p_right(ic_act,jf_act,kc_act),...
                                           bx_right(ic_act,jf_act,kc_act),by_right(ic_act,jf_act,kc_act),bz_right(ic_act,jf_act,kc_act),bnormj(ic_act,jf_act,kc_act),...
                                           fxj(ic_act,jf_act,kc_act),fyj(ic_act,jf_act,kc_act),fzj(ic_act,jf_act,kc_act),fnormj(ic_act,jf_act,kc_act),CA,...
                                           TRANS_j(ic_act,jf_act,kc_act,1),TRANS_j(ic_act,jf_act,kc_act,2),TRANS_j(ic_act,jf_act,kc_act,3));
       
    BstressX_j = Bstress_x_p + Bstress_x_n;
    BstressY_j = Bstress_y_p + Bstress_y_n;
    BstressZ_j = Bstress_z_p + Bstress_z_n;    

    % now add the magnetic HOGS flux
    bsqj(ic_act,jf_act,kc_act) = 0.5*((bx_left(ic_act,jf_act,kc_act)+fxj(ic_act,jf_act,kc_act)).^2+(by_left(ic_act,jf_act,kc_act)+fyj(ic_act,jf_act,kc_act)).^2+(bz_left(ic_act,jf_act,kc_act)+fzj(ic_act,jf_act,kc_act)).^2 + ...
                                      (bx_right(ic_act,jf_act,kc_act)+fxj(ic_act,jf_act,kc_act)).^2+(by_right(ic_act,jf_act,kc_act)+fyj(ic_act,jf_act,kc_act)).^2+(bz_right(ic_act,jf_act,kc_act)+fzj(ic_act,jf_act,kc_act)).^2)/2/CA;
    BstressX_j(ic_act,jf_act,kc_act) = BstressX_j(ic_act,jf_act,kc_act) - bsqj(ic_act,jf_act,kc_act).*(vx_right(ic_act,jf_act,kc_act) - vx_left(ic_act,jf_act,kc_act));
    BstressY_j(ic_act,jf_act,kc_act) = BstressY_j(ic_act,jf_act,kc_act) - bsqj(ic_act,jf_act,kc_act).*(vy_right(ic_act,jf_act,kc_act) - vy_left(ic_act,jf_act,kc_act));
    BstressZ_j(ic_act,jf_act,kc_act) = BstressZ_j(ic_act,jf_act,kc_act) - bsqj(ic_act,jf_act,kc_act).*(vz_right(ic_act,jf_act,kc_act) - vz_left(ic_act,jf_act,kc_act));

    %~~~~~~~~~END OF J-DIRECTION FLUX CALCULATIONS~~~~~~~~~
    
% Step 7. Now do the k-direction flux calculations
    % reconstruct in the k-direction 
    [rho_left, rho_right] = reconstruct_3DV(rho_h,if_act,jf_act,kf_act,PDMB,3,limiter,volume,vol_kface);     
    [vx_left, vx_right] = reconstruct_3DV(vx_h,if_act,jf_act,kf_act,PDMB,3,limiter,volume,vol_kface);    
    [vy_left, vy_right] = reconstruct_3DV(vy_h,if_act,jf_act,kf_act,PDMB,3,limiter,volume,vol_kface);        
    [vz_left, vz_right] = reconstruct_3DV(vz_h,if_act,jf_act,kf_act,PDMB,3,limiter,volume,vol_kface);         
    [p_left, p_right] = reconstruct_3DV(p_h,if_act,jf_act,kf_act,PDMB,3,limiter,volume,vol_kface);   
    
    % transform the velocities to i,j,k coordinate
    vi_left(ic_act,jc_act,kf_act) = vx_left(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,1) + ...
                                    vy_left(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,2) + ...
                                    vz_left(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,3);
    vj_left(ic_act,jc_act,kf_act) = vx_left(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,4) + ...
                                    vy_left(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,5) + ...
                                    vz_left(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,6);
    vk_left(ic_act,jc_act,kf_act) = vx_left(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,7) + ...
                                    vy_left(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,8) + ...
                                    vz_left(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,9);

    vi_right(ic_act,jc_act,kf_act) = vx_right(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,1) + ...
                                     vy_right(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,2) + ...
                                     vz_right(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,3);
    vj_right(ic_act,jc_act,kf_act) = vx_right(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,4) + ...
                                     vy_right(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,5) + ...
                                     vz_right(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,6);
    vk_right(ic_act,jc_act,kf_act) = vx_right(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,7) + ...
                                     vy_right(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,8) + ...
                                     vz_right(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,9);

    % Evaluate the k-interface fluxes                            
    [Frho_p(ic_act,jc_act,kf_act),FrhoVx_p(ic_act,jc_act,kf_act),FrhoVy_p(ic_act,jc_act,kf_act),FrhoVz_p(ic_act,jc_act,kf_act),Feng_p(ic_act,jc_act,kf_act), ...
     ~,~,~,~,~] = getFluidFlux(rho_left(ic_act,jc_act,kf_act),vi_left(ic_act,jc_act,kf_act),vj_left(ic_act,jc_act,kf_act),vk_left(ic_act,jc_act,kf_act),p_left(ic_act,jc_act,kf_act),gamma);
    [~,~,~,~,~, ...
     Frho_n(ic_act,jc_act,kf_act),FrhoVx_n(ic_act,jc_act,kf_act),FrhoVy_n(ic_act,jc_act,kf_act),FrhoVz_n(ic_act,jc_act,kf_act),Feng_n(ic_act,jc_act,kf_act)] ...
     = getFluidFlux(rho_right(ic_act,jc_act,kf_act),vi_right(ic_act,jc_act,kf_act),vj_right(ic_act,jc_act,kf_act),vk_right(ic_act,jc_act,kf_act),p_right(ic_act,jc_act,kf_act),gamma); 
 
    rho_flux_k = Frho_p + Frho_n;
    vx_flux_x = FrhoVx_p + FrhoVx_n;
    vy_flux_x = FrhoVy_p + FrhoVy_n;
    vz_flux_x = FrhoVz_p + FrhoVz_n;
    eng_flux_k = Feng_p + Feng_n;           

    % ~~~~~~~~~~~~~add the HOGS fix, k-direction
    diff_va = 0*Va_eff;
    diff_va(ic_act,jc_act,kf_act) = dif_factor*(Va_eff(ic_act,jc_act,kf_act-1) + Va_eff(ic_act,jc_act,kf_act));
    rho_flux_k(ic_act,jc_act,kf_act) = rho_flux_k(ic_act,jc_act,kf_act) - diff_va(ic_act,jc_act,kf_act).*(rho_right(ic_act,jc_act,kf_act) - rho_left(ic_act,jc_act,kf_act));
    vx_flux_x(ic_act,jc_act,kf_act)  = vx_flux_x(ic_act,jc_act,kf_act)  - diff_va(ic_act,jc_act,kf_act).*(rho_right(ic_act,jc_act,kf_act).*vi_right(ic_act,jc_act,kf_act) - rho_left(ic_act,jc_act,kf_act).*vi_left(ic_act,jc_act,kf_act));
    vy_flux_x(ic_act,jc_act,kf_act)  = vy_flux_x(ic_act,jc_act,kf_act)  - diff_va(ic_act,jc_act,kf_act).*(rho_right(ic_act,jc_act,kf_act).*vj_right(ic_act,jc_act,kf_act) - rho_left(ic_act,jc_act,kf_act).*vj_left(ic_act,jc_act,kf_act));
    vz_flux_x(ic_act,jc_act,kf_act)  = vz_flux_x(ic_act,jc_act,kf_act)  - diff_va(ic_act,jc_act,kf_act).*(rho_right(ic_act,jc_act,kf_act).*vk_right(ic_act,jc_act,kf_act) - rho_left(ic_act,jc_act,kf_act).*vk_left(ic_act,jc_act,kf_act));
    eng_flux_k(ic_act,jc_act,kf_act) = eng_flux_k(ic_act,jc_act,kf_act) - diff_va(ic_act,jc_act,kf_act).*(0.5.*rho_right(ic_act,jc_act,kf_act).*(vi_right(ic_act,jc_act,kf_act).^2+vj_right(ic_act,jc_act,kf_act).^2+vk_right(ic_act,jc_act,kf_act).^2) + p_right(ic_act,jc_act,kf_act)./(gamma-1) ...
                                                                                                        - 0.5.*rho_left(ic_act,jc_act,kf_act).*(vi_left(ic_act,jc_act,kf_act).^2+vj_left(ic_act,jc_act,kf_act).^2+vk_left(ic_act,jc_act,kf_act).^2) - p_left(ic_act,jc_act,kf_act)./(gamma-1));
    
    % transform the momentum fluxes back to x,y,z
    vx_flux_k(ic_act,jc_act,kf_act) = vx_flux_x(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,1) + ...
                                      vy_flux_x(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,4) + ...
                                      vz_flux_x(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,7);
    vy_flux_k(ic_act,jc_act,kf_act) = vx_flux_x(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,2) + ...
                                      vy_flux_x(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,5) + ...
                                      vz_flux_x(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,8);
    vz_flux_k(ic_act,jc_act,kf_act) = vx_flux_x(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,3) + ...
                                      vy_flux_x(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,6) + ...
                                      vz_flux_x(ic_act,jc_act,kf_act).*TRANS_k(ic_act,jc_act,kf_act,9);

   % calculate magnetic stress in the k direction
    [bx_left, bx_right] = reconstruct_3DV(bx_h,if_act,jf_act,kf_act,PDMB,3,limiter,volume,vol_kface);  
    [by_left, by_right] = reconstruct_3DV(by_h,if_act,jf_act,kf_act,PDMB,3,limiter,volume,vol_kface);   
    [bz_left, bz_right] = reconstruct_3DV(bz_h,if_act,jf_act,kf_act,PDMB,3,limiter,volume,vol_kface);   
                               
    bnormk = bk_h./k_face_area;
    
    % getMagneticStress(rho,Vx,p,bx,by,bz,bnorm,bx0,by0,bz0,b0n,CA,trans_1,trans_2,trans_3)
    [Bstress_x_p(ic_act,jc_act,kf_act), Bstress_y_p(ic_act,jc_act,kf_act), Bstress_z_p(ic_act,jc_act,kf_act), ...
     ~, ~, ~] = getMagneticStress(rho_left(ic_act,jc_act,kf_act),vi_left(ic_act,jc_act,kf_act),p_left(ic_act,jc_act,kf_act),...
                                           bx_left(ic_act,jc_act,kf_act),by_left(ic_act,jc_act,kf_act),bz_left(ic_act,jc_act,kf_act),bnormk(ic_act,jc_act,kf_act),...
                                           fxk(ic_act,jc_act,kf_act),fyk(ic_act,jc_act,kf_act),fzk(ic_act,jc_act,kf_act),fnormk(ic_act,jc_act,kf_act),CA,...
                                           TRANS_k(ic_act,jc_act,kf_act,1),TRANS_k(ic_act,jc_act,kf_act,2),TRANS_k(ic_act,jc_act,kf_act,3));
    
    [~, ~, ~, ...
     Bstress_x_n(ic_act,jc_act,kf_act), Bstress_y_n(ic_act,jc_act,kf_act), Bstress_z_n(ic_act,jc_act,kf_act)] ...
              = getMagneticStress(rho_right(ic_act,jc_act,kf_act),vi_right(ic_act,jc_act,kf_act),p_right(ic_act,jc_act,kf_act),...
                                           bx_right(ic_act,jc_act,kf_act),by_right(ic_act,jc_act,kf_act),bz_right(ic_act,jc_act,kf_act),bnormk(ic_act,jc_act,kf_act),...
                                           fxk(ic_act,jc_act,kf_act),fyk(ic_act,jc_act,kf_act),fzk(ic_act,jc_act,kf_act),fnormk(ic_act,jc_act,kf_act),CA,...
                                           TRANS_k(ic_act,jc_act,kf_act,1),TRANS_k(ic_act,jc_act,kf_act,2),TRANS_k(ic_act,jc_act,kf_act,3));
              
    BstressX_k = Bstress_x_p + Bstress_x_n;
    BstressY_k = Bstress_y_p + Bstress_y_n;
    BstressZ_k = Bstress_z_p + Bstress_z_n;    
                                 
    % now add the magnetic HOGS flux
    bsqk(ic_act,jc_act,kf_act) = 0.5*((bx_left(ic_act,jc_act,kf_act)+fxk(ic_act,jc_act,kf_act)).^2+(by_left(ic_act,jc_act,kf_act)+fyk(ic_act,jc_act,kf_act)).^2+(bz_left(ic_act,jc_act,kf_act)+fzk(ic_act,jc_act,kf_act)).^2 + ...
                                      (bx_right(ic_act,jc_act,kf_act)+fxk(ic_act,jc_act,kf_act)).^2+(by_right(ic_act,jc_act,kf_act)+fyk(ic_act,jc_act,kf_act)).^2+(bz_right(ic_act,jc_act,kf_act)+fzk(ic_act,jc_act,kf_act)).^2)/2/CA;
    BstressX_k(ic_act,jc_act,kf_act) = BstressX_k(ic_act,jc_act,kf_act) - bsqk(ic_act,jc_act,kf_act).*(vx_right(ic_act,jc_act,kf_act) - vx_left(ic_act,jc_act,kf_act));
    BstressY_k(ic_act,jc_act,kf_act) = BstressY_k(ic_act,jc_act,kf_act) - bsqk(ic_act,jc_act,kf_act).*(vy_right(ic_act,jc_act,kf_act) - vy_left(ic_act,jc_act,kf_act));
    BstressZ_k(ic_act,jc_act,kf_act) = BstressZ_k(ic_act,jc_act,kf_act) - bsqk(ic_act,jc_act,kf_act).*(vz_right(ic_act,jc_act,kf_act) - vz_left(ic_act,jc_act,kf_act));

    %~~~~~~~~~END OF K-DIRECTION FLUX CALCULATIONS~~~~~~~~~
    
% Step 8. Face-integrals for the flux terms for finite volume update
    % initialize the variables
    drho = rho_h.*0;
    dpxF = rho_h.*0;
    dpyF = rho_h.*0;
    dpzF = rho_h.*0;     
    deng = rho_h.*0;
    dpxB = rho_h.*0;
    dpyB = rho_h.*0;
    dpzB = rho_h.*0;     
    
    % fluid stresses integrated over volume (surface integral): Sum(F*dS)
    % be careful with the index shift in each direction
    drho(ic_act,jc_act,kc_act) = - dt.*(i_face_area(ic_act+1,jc_act,kc_act).*rho_flux_i(ic_act+1,jc_act,kc_act)-i_face_area(ic_act,jc_act,kc_act).*rho_flux_i(ic_act,jc_act,kc_act) +  ...
                                        j_face_area(ic_act,jc_act+1,kc_act).*rho_flux_j(ic_act,jc_act+1,kc_act)-j_face_area(ic_act,jc_act,kc_act).*rho_flux_j(ic_act,jc_act,kc_act) + ...
                                        k_face_area(ic_act,jc_act,kc_act+1).*rho_flux_k(ic_act,jc_act,kc_act+1)-k_face_area(ic_act,jc_act,kc_act).*rho_flux_k(ic_act,jc_act,kc_act))./volume(ic_act,jc_act,kc_act);
    dpxF(ic_act,jc_act,kc_act) = - dt.*(i_face_area(ic_act+1,jc_act,kc_act).*vx_flux_i(ic_act+1,jc_act,kc_act)-i_face_area(ic_act,jc_act,kc_act).*vx_flux_i(ic_act,jc_act,kc_act) +  ...
                                        j_face_area(ic_act,jc_act+1,kc_act).*vx_flux_j(ic_act,jc_act+1,kc_act)-j_face_area(ic_act,jc_act,kc_act).*vx_flux_j(ic_act,jc_act,kc_act) + ...
                                        k_face_area(ic_act,jc_act,kc_act+1).*vx_flux_k(ic_act,jc_act,kc_act+1)-k_face_area(ic_act,jc_act,kc_act).*vx_flux_k(ic_act,jc_act,kc_act))./volume(ic_act,jc_act,kc_act);
    dpyF(ic_act,jc_act,kc_act) = - dt.*(i_face_area(ic_act+1,jc_act,kc_act).*vy_flux_i(ic_act+1,jc_act,kc_act)-i_face_area(ic_act,jc_act,kc_act).*vy_flux_i(ic_act,jc_act,kc_act) +  ...
                                        j_face_area(ic_act,jc_act+1,kc_act).*vy_flux_j(ic_act,jc_act+1,kc_act)-j_face_area(ic_act,jc_act,kc_act).*vy_flux_j(ic_act,jc_act,kc_act) + ...
                                        k_face_area(ic_act,jc_act,kc_act+1).*vy_flux_k(ic_act,jc_act,kc_act+1)-k_face_area(ic_act,jc_act,kc_act).*vy_flux_k(ic_act,jc_act,kc_act))./volume(ic_act,jc_act,kc_act);
    dpzF(ic_act,jc_act,kc_act) = - dt.*(i_face_area(ic_act+1,jc_act,kc_act).*vz_flux_i(ic_act+1,jc_act,kc_act)-i_face_area(ic_act,jc_act,kc_act).*vz_flux_i(ic_act,jc_act,kc_act) +  ...
                                        j_face_area(ic_act,jc_act+1,kc_act).*vz_flux_j(ic_act,jc_act+1,kc_act)-j_face_area(ic_act,jc_act,kc_act).*vz_flux_j(ic_act,jc_act,kc_act) + ...
                                        k_face_area(ic_act,jc_act,kc_act+1).*vz_flux_k(ic_act,jc_act,kc_act+1)-k_face_area(ic_act,jc_act,kc_act).*vz_flux_k(ic_act,jc_act,kc_act))./volume(ic_act,jc_act,kc_act);
    deng(ic_act,jc_act,kc_act) = - dt.*(i_face_area(ic_act+1,jc_act,kc_act).*eng_flux_i(ic_act+1,jc_act,kc_act)-i_face_area(ic_act,jc_act,kc_act).*eng_flux_i(ic_act,jc_act,kc_act) +  ...
                                        j_face_area(ic_act,jc_act+1,kc_act).*eng_flux_j(ic_act,jc_act+1,kc_act)-j_face_area(ic_act,jc_act,kc_act).*eng_flux_j(ic_act,jc_act,kc_act) + ...
                                        k_face_area(ic_act,jc_act,kc_act+1).*eng_flux_k(ic_act,jc_act,kc_act+1)-k_face_area(ic_act,jc_act,kc_act).*eng_flux_k(ic_act,jc_act,kc_act))./volume(ic_act,jc_act,kc_act);

    % magnetic stresses integrated over volume (surface integral)
    dpxB(ic_act,jc_act,kc_act) = - dt.*(i_face_area(ic_act+1,jc_act,kc_act).*BstressX_i(ic_act+1,jc_act,kc_act)-i_face_area(ic_act,jc_act,kc_act).*BstressX_i(ic_act,jc_act,kc_act) +  ...
                                        j_face_area(ic_act,jc_act+1,kc_act).*BstressX_j(ic_act,jc_act+1,kc_act)-j_face_area(ic_act,jc_act,kc_act).*BstressX_j(ic_act,jc_act,kc_act) + ...
                                        k_face_area(ic_act,jc_act,kc_act+1).*BstressX_k(ic_act,jc_act,kc_act+1)-k_face_area(ic_act,jc_act,kc_act).*BstressX_k(ic_act,jc_act,kc_act))./volume(ic_act,jc_act,kc_act);
    dpyB(ic_act,jc_act,kc_act) = - dt.*(i_face_area(ic_act+1,jc_act,kc_act).*BstressY_i(ic_act+1,jc_act,kc_act)-i_face_area(ic_act,jc_act,kc_act).*BstressY_i(ic_act,jc_act,kc_act) +  ...
                                        j_face_area(ic_act,jc_act+1,kc_act).*BstressY_j(ic_act,jc_act+1,kc_act)-j_face_area(ic_act,jc_act,kc_act).*BstressY_j(ic_act,jc_act,kc_act) + ...
                                        k_face_area(ic_act,jc_act,kc_act+1).*BstressY_k(ic_act,jc_act,kc_act+1)-k_face_area(ic_act,jc_act,kc_act).*BstressY_k(ic_act,jc_act,kc_act))./volume(ic_act,jc_act,kc_act);
    dpzB(ic_act,jc_act,kc_act) = - dt.*(i_face_area(ic_act+1,jc_act,kc_act).*BstressZ_i(ic_act+1,jc_act,kc_act)-i_face_area(ic_act,jc_act,kc_act).*BstressZ_i(ic_act,jc_act,kc_act) +  ...
                                        j_face_area(ic_act,jc_act+1,kc_act).*BstressZ_j(ic_act,jc_act+1,kc_act)-j_face_area(ic_act,jc_act,kc_act).*BstressZ_j(ic_act,jc_act,kc_act) + ...
                                        k_face_area(ic_act,jc_act,kc_act+1).*BstressZ_k(ic_act,jc_act,kc_act+1)-k_face_area(ic_act,jc_act,kc_act).*BstressZ_k(ic_act,jc_act,kc_act))./volume(ic_act,jc_act,kc_act);
    
% End of Hydro

end