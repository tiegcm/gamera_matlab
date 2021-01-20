function [Ei,Ej,Ek]=getEfields(rho,vx,vy,vz,bi,bj,bk,BX0_I,BY0_I,BX0_J,BY0_J,BX0_K,BY0_K,...
                              ifaceA_Kedge,jfaceA_Kedge,xnqi_K,ynqi_K,xnqj_K,ynqj_K,DETK,k_edge,vec_k,...
                              jfaceA_Iedge,kfaceA_Iedge,xnqj_I,ynqj_I,xnqk_I,ynqk_I,DETI,i_edge,vec_i,...
                              kfaceA_Jedge,ifaceA_Jedge,xnqk_J,ynqk_J,xnqi_J,ynqi_J,DETJ,j_edge,vec_j,PDMB,limiter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code calculate the electric field
% MAIN ALGORITHM: (Ek as an example)
%     1. get high order interpolated veolcity (v_interp) at k-edges
%        then transform v_interp into the K-edge aligned system use vec_k[]
%     2. get high order reconstructed magnetic field (B_avg) at the k-edge
%        then transform B_avg into the K-edge aligned system use: 
%          bx_cor = DETK.*( ynqj_K.*bi_avg - ynqi_K.*bj_avg) + BX0_K;
%          by_cor = DETK.*(-xnqj_K.*bi_avg + xnqi_K.*bj_avg) + BY0_K;
%     3. compute Ek = v x B + vd*deltaB                          
%     4. convert E-field to E-potential by E*dl
% 
% Ei and Ej are calculate using the same way

global ic_act jc_act kc_act if_act jf_act kf_act CA

% get the electric field at cell edges (k-edge)
% step 1, get high order interpolated veolcities at k-edges and transform
  % to the k-edge coordinates
  vx_avg = center2corner(vx,'ij'); % dimension (if_act,jf_act,kc_act)
  vy_avg = center2corner(vy,'ij'); % dimension (if_act,jf_act,kc_act)
  vz_avg = center2corner(vz,'ij'); % dimension (if_act,jf_act,kc_act)

  % x-component in the k-edge coordinate system
  vx_cor(if_act,jf_act,kc_act) = vx_avg(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,1) + ...
                                 vy_avg(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,2) + ...
                                 vz_avg(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,3);
  % y-component in the k-edge coordinate system
  vy_cor(if_act,jf_act,kc_act) = vx_avg(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,4) + ...
                                 vy_avg(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,5) + ...
                                 vz_avg(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,6);

% step 2, get high order interpolated magnetic flux/field at the k-edge
  % since bi and bj are already face flux, the geometry (face area) is
  % already taken into account in the reconstruction, here we use
  % reconstruct_3D() instead of reconstruct_3DV() without volume invoved
  [bi_left, bi_right] = reconstruct_3DV(bi,if_act,jf_act,kc_act,PDMB, 2,limiter); % dimension (if_act,jf_act,kc_act)
  [bj_left, bj_right] = reconstruct_3DV(bj,if_act,jf_act,kc_act,PDMB, 1,limiter); % dimension (if_act,jf_act,kc_act)

  % b-flux to b-fields
  bi_left = bi_left./ifaceA_Kedge;
  bi_right = bi_right./ifaceA_Kedge;
  bj_left = bj_left./jfaceA_Kedge;
  bj_right = bj_right./jfaceA_Kedge;

  % compute the average/upwind Bi, Bj fields at K-edge. 
  bi_avg = 0.5*(bi_left + bi_right);
  bj_avg = 0.5*(bj_left + bj_right);

  % then transform the edge B field to the K-edge coordinate
  bx_cor = DETK.*( ynqj_K.*bi_avg - ynqi_K.*bj_avg) + BX0_K;
  by_cor = DETK.*(-xnqj_K.*bi_avg + xnqi_K.*bj_avg) + BY0_K;

% step 3. compute the convective part of the Ek
  Ek = -(vx_cor.*by_cor - vy_cor.*bx_cor);

  % then compute the diffusive part for Ek
  rhobar(if_act,jf_act,kc_act) = 0.25*(rho(if_act,jf_act,kc_act) + rho(if_act-1,jf_act,kc_act) + rho(if_act,jf_act-1,kc_act) + rho(if_act-1,jf_act-1,kc_act));
  % This is the diffusive speed for a LLF flux
  dvzz(if_act,jf_act,kc_act) = sqrt((bx_cor(if_act,jf_act,kc_act).^2+by_cor(if_act,jf_act,kc_act).^2)./rhobar(if_act,jf_act,kc_act));
  % limit the dvzz with the Boris speed (if CA>>1 then no limiting)
  dvzz(if_act,jf_act,kc_act) = dvzz(if_act,jf_act,kc_act).*CA./sqrt(dvzz(if_act,jf_act,kc_act).^2+CA.^2) +...
                               sqrt(vx_cor(if_act,jf_act,kc_act).^2 + vy_cor(if_act,jf_act,kc_act).^2);

  % add the convection speed to the diffusive speed, NOTE that the default
  % diffusive speed is the local Alfven speed, in some flow cases if the wave
  % speed is dominated by the convective speed, the Alfven speed may not be
  % large enough to damp out all the numerical oscillations - this could be
  % dvzz(if_act,jf_act,kc_act) = dvzz(if_act,jf_act,kc_act)+sqrt(vx_cor(if_act,jf_act,kc_act).^2 + vy_cor(if_act,jf_act,kc_act).^2);                         
  
  % compute the total Ek with convective part and diffusive part
  Ek(if_act,jf_act,kc_act) = Ek(if_act,jf_act,kc_act) + 0.5*dvzz(if_act,jf_act,kc_act).* (bi_left(if_act,jf_act,kc_act) + ...
                                                                                          bj_right(if_act,jf_act,kc_act)-...
                                                                                          bi_right(if_act,jf_act,kc_act)- ...
                                                                                          bj_left(if_act,jf_act,kc_act));

% step 4. Ek times the edge length - the actually electric potential used in the Faraday's law for updating magnetic fluxes                                                                               
  Ek(if_act,jf_act,kc_act)= Ek(if_act,jf_act,kc_act).*k_edge(if_act,jf_act,kc_act);

  % then clear the temporary variables for Ei and Ej calculations
  clear vx_cor vy_cor bx_cor by_cor;                                                                                
                                                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the electric field at cell edges (i-edge)
% step 1, get high order interpolated veolcities at k-edges and transform
  % to the k-edge coordinates
  vx_avg = center2corner(vx,'jk'); % dimension (if_act,jf_act,kc_act)
  vy_avg = center2corner(vy,'jk'); % dimension (if_act,jf_act,kc_act)
  vz_avg = center2corner(vz,'jk'); % dimension (if_act,jf_act,kc_act)

  % x-component in the i-edge coordinate system
  vx_cor(ic_act,jf_act,kf_act) = vx_avg(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,1) + ...
                                 vy_avg(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,2) + ...
                                 vz_avg(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,3);
  % y-component in the i-edge coordinate system
  vy_cor(ic_act,jf_act,kf_act) = vx_avg(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,4) + ...
                                 vy_avg(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,5) + ...
                                 vz_avg(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,6);

% step 2, get high order interpolated magnetic flux/field at the i-edge    
  [bj_left, bj_right] = reconstruct_3DV(bj,ic_act,jf_act,kf_act,PDMB, 3,limiter); % dimension (if_act,jf_act,kc_act)
  [bk_left, bk_right] = reconstruct_3DV(bk,ic_act,jf_act,kf_act,PDMB, 2,limiter); % dimension (if_act,jf_act,kc_act)

  % b-flux to b-fields
  bj_left = bj_left./jfaceA_Iedge;
  bj_right = bj_right./jfaceA_Iedge;
  bk_left = bk_left./kfaceA_Iedge;
  bk_right = bk_right./kfaceA_Iedge;

  % then compute the B_avg field at I-edge. 
  bj_avg = 0.5*(bj_left + bj_right);
  bk_avg = 0.5*(bk_left + bk_right);

  % then transform the edge B field to the I-edge coordinate
  bx_cor = DETI.*( ynqk_I.*bj_avg - ynqj_I.*bk_avg) + BX0_I;
  by_cor = DETI.*(-xnqk_I.*bj_avg + xnqj_I.*bk_avg) + BY0_I;

% step 3. compute the convective part of the Ei  
  Ei(ic_act,jf_act,kf_act) = -( vx_cor(ic_act,jf_act,kf_act).*by_cor(ic_act,jf_act,kf_act) - vy_cor(ic_act,jf_act,kf_act).*bx_cor(ic_act,jf_act,kf_act));

  % calculate a diffusive speed for eta
  rhobar(ic_act,jf_act,kf_act) = 0.25*(rho(ic_act,jf_act,kf_act) + rho(ic_act,jf_act-1,kf_act) + rho(ic_act,jf_act,kf_act-1) + rho(ic_act,jf_act-1,kf_act-1));
  % This is the local diffusive speed for a LLF flux
  dvzz(ic_act,jf_act,kf_act) = sqrt((bx_cor(ic_act,jf_act,kf_act).^2+by_cor(ic_act,jf_act,kf_act).^2)./rhobar(ic_act,jf_act,kf_act));
  % limit the dvzz with Boris speed by CA (CA>>1 then no limiting)
  dvzz(ic_act,jf_act,kf_act) = dvzz(ic_act,jf_act,kf_act).*CA./sqrt(dvzz(ic_act,jf_act,kf_act).^2+CA.^2) + ...
                               sqrt(vx_cor(ic_act,jf_act,kf_act).^2 + vy_cor(ic_act,jf_act,kf_act).^2);
                           
  % add the convection speed to the diffusive speed, NOTE that the default
  % diffusive speed is the local Alfven speed, in some flow cases if the wave
  % speed is dominated by the convective speed, the Alfven speed may not be
  % large enough to damp out all the numerical oscillations
  % dvzz(ic_act,jf_act,kf_act) = dvzz(ic_act,jf_act,kf_act)+sqrt(vx_cor(ic_act,jf_act,kf_act).^2 + vy_cor(ic_act,jf_act,kf_act).^2);                         
                         
  Ei(ic_act,jf_act,kf_act) = Ei(ic_act,jf_act,kf_act) + 0.5*dvzz(ic_act,jf_act,kf_act).* (bj_left(ic_act,jf_act,kf_act) + ...
                                                                                          bk_right(ic_act,jf_act,kf_act)-...
                                                                                          bj_right(ic_act,jf_act,kf_act)- ...
                                                                                          bk_left(ic_act,jf_act,kf_act));
  
% Step 4. compute the total Ei with convective part and diffusive part                                                                                    
  Ei(ic_act,jf_act,kf_act) = Ei(ic_act,jf_act,kf_act).*i_edge(ic_act,jf_act,kf_act);                                                                                

  clear vx_cor vy_cor bx_cor by_cor;                                                                                

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the electric field at j-edge
% step 1, get high order interpolated veolcities at j-edges and transform
  % to the k-edge coordinates
  vx_avg = center2corner(vx,'ki'); % dimension (if_act,jf_act,kc_act)
  vy_avg = center2corner(vy,'ki'); % dimension (if_act,jf_act,kc_act)
  vz_avg = center2corner(vz,'ki'); % dimension (if_act,jf_act,kc_act)

  % x-component in the j-edge coordinate system
  vx_cor(if_act,jc_act,kf_act) = vx_avg(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,1) + ...
                                 vy_avg(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,2) + ...
                                 vz_avg(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,3);
  % y-component in the j-edge coordinate system
  vy_cor(if_act,jc_act,kf_act) = vx_avg(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,4) + ...
                                 vy_avg(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,5) + ...
                                 vz_avg(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,6);
                           
% step 2, get high order interpolated magnetic flux/field at the j-edge    
  [bk_left, bk_right] = reconstruct_3DV(bk,if_act,jc_act,kf_act,PDMB, 1,limiter); % dimension (if_act,jf_act,kc_act)
  [bi_left, bi_right] = reconstruct_3DV(bi,if_act,jc_act,kf_act,PDMB, 3,limiter); % dimension (if_act,jf_act,kc_act)
                           
  % b-flux to b-fields
  bk_left = bk_left./kfaceA_Jedge;
  bk_right = bk_right./kfaceA_Jedge;
  bi_left = bi_left./ifaceA_Jedge;
  bi_right = bi_right./ifaceA_Jedge;

  % then compute the B_avg field at J-edge. 
  bk_avg = 0.5*(bk_left + bk_right);
  bi_avg = 0.5*(bi_left + bi_right);

  % now transform the B_avg field to the J-edge aligned + background
  bx_cor = DETJ.*( ynqi_J.*bk_avg - ynqk_J.*bi_avg) + BX0_J;
  by_cor = DETJ.*(-xnqi_J.*bk_avg + xnqk_J.*bi_avg) + BY0_J;

% step 3. compute the convective part of the Ej 
  Ej(if_act,jc_act,kf_act) = -( vx_cor(if_act,jc_act,kf_act).*by_cor(if_act,jc_act,kf_act) - ...
                                vy_cor(if_act,jc_act,kf_act).*bx_cor(if_act,jc_act,kf_act) );

  % 7. calculate a diffusive speed for eta
  rhobar(if_act,jc_act,kf_act) = 0.25*(rho(if_act,jc_act,kf_act) + rho(if_act,jc_act,kf_act-1) + rho(if_act-1,jc_act,kf_act) + rho(if_act-1,jc_act,kf_act-1));
  % This is the fastest diffusive speed for a LLF flux
  dvzz(if_act,jc_act,kf_act) = sqrt((bx_cor(if_act,jc_act,kf_act).^2+by_cor(if_act,jc_act,kf_act).^2)./rhobar(if_act,jc_act,kf_act));
  % limit the dvzz with Boris speed by CA (CA>>1 then no limiting)
  dvzz(if_act,jc_act,kf_act) = dvzz(if_act,jc_act,kf_act).*CA./sqrt(dvzz(if_act,jc_act,kf_act).^2+CA.^2) + ...
                               sqrt(vx_cor(if_act,jc_act,kf_act).^2 + vy_cor(if_act,jc_act,kf_act).^2);
                           
  % add the convection speed to the diffusive speed, NOTE that the default
  % diffusive speed is the local Alfven speed, in some flow cases if the wave
  % speed is dominated by the convective speed, the Alfven speed may not be
  % large enough to damp out all the numerical oscillations
  % dvzz(if_act,jc_act,kf_act) = sqrt(vx_cor(if_act,jc_act,kf_act).^2 + vy_cor(if_act,jc_act,kf_act).^2);
     
  % compute the total Ej with convective part and diffusive part
  Ej(if_act,jc_act,kf_act) = Ej(if_act,jc_act,kf_act) + 0.5*dvzz(if_act,jc_act,kf_act).* (bk_left(if_act,jc_act,kf_act) + ...
                                                                                          bi_right(if_act,jc_act,kf_act)-...
                                                                                          bk_right(if_act,jc_act,kf_act)- ...
                                                                                          bi_left(if_act,jc_act,kf_act));
% Step 4. compute the total Ej with convective part and diffusive part                                                                                    
  Ej(if_act,jc_act,kf_act) = Ej(if_act,jc_act,kf_act).*j_edge(if_act,jc_act,kf_act);

  clear vx_cor vy_cor bx_cor by_cor;     
  
  % END of electric field calculations
end
          