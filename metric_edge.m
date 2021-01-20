function [ifaceA_Kedge,jfaceA_Kedge,xnqi_K,ynqi_K,xnqj_K,ynqj_K,DETK,k_edge,vec_k,...
          jfaceA_Iedge,kfaceA_Iedge,xnqj_I,ynqj_I,xnqk_I,ynqk_I,DETI,i_edge,vec_i,...
          kfaceA_Jedge,ifaceA_Jedge,xnqk_J,ynqk_J,xnqi_J,ynqi_J,DETJ,j_edge,vec_j]...
         = metric_edge(x,y,z,...
                       i_face_norm_x,i_face_norm_y,i_face_norm_z,i_face_area,...
                       j_face_norm_x,j_face_norm_y,j_face_norm_z,j_face_area,...
                       k_face_norm_x,k_face_norm_y,k_face_norm_z,k_face_area)

% This code computes the transform matrices for the edge-aligned coordinate
% systems, together with the transform coefficients for the electric field
% calculations
%
% INPUT: x,y,z               - cell corners      
%        i_face_norm_{x,y,z} - face-normal vectors at i-faces           
%        j_face_norm_{x,y,z} - face-normal vectors at j-faces           
%        k_face_norm_{x,y,z} - face-normal vectors at k-faces   
%        {i,j,k}_face_area   - face areas at {i,j,k} interfaces
%
% NOTE: edge-aligned coordinates are only calculated at active edges

global ic_act jc_act kc_act if_act jf_act kf_act PDMB    

disp('Computing Edge-Aligned Coordinates...');

% k-edge aligned coordinate system for Ek
% calculate the coordinate systems at edges for magnetic field evolution
% note that the edge E fields are only calculated at the active regions, so
% the indices are if_act,jf_act,kc_act for k-edge coordinates 
%                 ic_act,jf_act,kf_act for i-edge coordinates
%                 if_act,jc_act,kf_act for j-edge coordinates
% 
% k-edge metrics, dimension (if_act, jf_act, kc_act)
% vec_k(:,:,:,1:9)
% (7-9) vector in the k-edge direction
% (1-3) vector in the direction of sweep
% (4-6) vector in the right handed perp to both k and sweep
% 
% An example transform follows:
% the d0-component in the k-edge coordinate system
% vx_cor(if_act,jf_act,kc_act) = vx_avg(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,1) + ...
%                                vy_avg(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,2) + ...
%                                vz_avg(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,3);
% the d1-component in the k-edge coordinate system
% vy_cor(if_act,jf_act,kc_act) = vx_avg(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,4) + ...
%                                vy_avg(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,5) + ...
%                                vz_avg(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,6);
%
% (dk_edge_x, dk_edge_y, dk_edge_z) is the k-edge vector and 
% (k_edge_x, k_edge_y, k_edge_z) is the k-edge unit vector
dk_edge_x(if_act,jf_act,kc_act) = x(if_act,jf_act,kc_act+1) - x(if_act,jf_act,kc_act);
dk_edge_y(if_act,jf_act,kc_act) = y(if_act,jf_act,kc_act+1) - y(if_act,jf_act,kc_act);
dk_edge_z(if_act,jf_act,kc_act) = z(if_act,jf_act,kc_act+1) - z(if_act,jf_act,kc_act);
k_edge(if_act,jf_act,kc_act) = sqrt(dk_edge_x(if_act,jf_act,kc_act).^2 + dk_edge_y(if_act,jf_act,kc_act).^2 + dk_edge_z(if_act,jf_act,kc_act).^2);
k_edge_x(if_act,jf_act,kc_act) = dk_edge_x(if_act,jf_act,kc_act)./k_edge(if_act,jf_act,kc_act);
k_edge_y(if_act,jf_act,kc_act) = dk_edge_y(if_act,jf_act,kc_act)./k_edge(if_act,jf_act,kc_act);
k_edge_z(if_act,jf_act,kc_act) = dk_edge_z(if_act,jf_act,kc_act)./k_edge(if_act,jf_act,kc_act);

dxd(if_act,jf_act,kc_act) = 0.25.*(x(if_act,jf_act+1,kc_act) + x(if_act,jf_act+1,kc_act+1) - x(if_act,jf_act-1,kc_act) - x(if_act,jf_act-1,kc_act+1));
dyd(if_act,jf_act,kc_act) = 0.25.*(y(if_act,jf_act+1,kc_act) + y(if_act,jf_act+1,kc_act+1) - y(if_act,jf_act-1,kc_act) - y(if_act,jf_act-1,kc_act+1));
dzd(if_act,jf_act,kc_act) = 0.25.*(z(if_act,jf_act+1,kc_act) + z(if_act,jf_act+1,kc_act+1) - z(if_act,jf_act-1,kc_act) - z(if_act,jf_act-1,kc_act+1));

% d{x,y,z}_0 is the first perpendicular unit vector (e1)
[dx_0(if_act,jf_act,kc_act),dy_0(if_act,jf_act,kc_act),dz_0(if_act,jf_act,kc_act)] = ...
    CROSS2(dxd(if_act,jf_act,kc_act),dyd(if_act,jf_act,kc_act),dzd(if_act,jf_act,kc_act),...
           k_edge_x(if_act,jf_act,kc_act),k_edge_y(if_act,jf_act,kc_act),k_edge_z(if_act,jf_act,kc_act));
       
dx_0_norm(if_act,jf_act,kc_act) = sqrt(dx_0(if_act,jf_act,kc_act).^2+dy_0(if_act,jf_act,kc_act).^2+dz_0(if_act,jf_act,kc_act).^2);
dx_0(if_act,jf_act,kc_act) = dx_0(if_act,jf_act,kc_act)./dx_0_norm(if_act,jf_act,kc_act);
dy_0(if_act,jf_act,kc_act) = dy_0(if_act,jf_act,kc_act)./dx_0_norm(if_act,jf_act,kc_act);
dz_0(if_act,jf_act,kc_act) = dz_0(if_act,jf_act,kc_act)./dx_0_norm(if_act,jf_act,kc_act);

% d{x,y,z}_1 is the second perpendicular unit vector (e2), 
% so (e1, e2, dk_edge) forms an orthogonal system at k-edge
[dx_1(if_act,jf_act,kc_act),dy_1(if_act,jf_act,kc_act),dz_1(if_act,jf_act,kc_act)] = ...
    CROSS2(k_edge_x(if_act,jf_act,kc_act),k_edge_y(if_act,jf_act,kc_act),k_edge_z(if_act,jf_act,kc_act),...
           dx_0(if_act,jf_act,kc_act),dy_0(if_act,jf_act,kc_act),dz_0(if_act,jf_act,kc_act));
       
% vec_k matrix is used to transform the velocity and Bupwind vector to the
% k_edge-normal coordinate. NOTE that in the electric field calculation,
% only the d0 (1:3) and d1 (4:6) components are used for vxB
vec_k(if_act,jf_act,kc_act,7) = k_edge_x(if_act,jf_act,kc_act);
vec_k(if_act,jf_act,kc_act,8) = k_edge_y(if_act,jf_act,kc_act);
vec_k(if_act,jf_act,kc_act,9) = k_edge_z(if_act,jf_act,kc_act);

vec_k(if_act,jf_act,kc_act,1) = dx_0(if_act,jf_act,kc_act);
vec_k(if_act,jf_act,kc_act,2) = dy_0(if_act,jf_act,kc_act);
vec_k(if_act,jf_act,kc_act,3) = dz_0(if_act,jf_act,kc_act);

vec_k(if_act,jf_act,kc_act,4) = dx_1(if_act,jf_act,kc_act);
vec_k(if_act,jf_act,kc_act,5) = dy_1(if_act,jf_act,kc_act);
vec_k(if_act,jf_act,kc_act,6) = dz_1(if_act,jf_act,kc_act);

% i-edge metrics, dimension (ic_act, jf_act, kf_act)
% vec_i(:,:,:,1:9)
% (7-9) vector in the i-edge direction
% (1-3) vector in the direction of sweep
% (4-6) vector in the right handed perp to both k and sweep
% 
% (di_edge_x, di_edge_y, di_edge_z) : i-edge vector
% (i_edge_x, i_edge_y, i_edge_z)    : i-edge unit vector
% i_edge                            : i-edge length

di_edge_x(ic_act,jf_act,kf_act) = x(ic_act+1,jf_act,kf_act) - x(ic_act,jf_act,kf_act);
di_edge_y(ic_act,jf_act,kf_act) = y(ic_act+1,jf_act,kf_act) - y(ic_act,jf_act,kf_act);
di_edge_z(ic_act,jf_act,kf_act) = z(ic_act+1,jf_act,kf_act) - z(ic_act,jf_act,kf_act);
i_edge(ic_act,jf_act,kf_act) = sqrt(di_edge_x(ic_act,jf_act,kf_act).^2 + di_edge_y(ic_act,jf_act,kf_act).^2 + di_edge_z(ic_act,jf_act,kf_act).^2);
i_edge_x(ic_act,jf_act,kf_act) = di_edge_x(ic_act,jf_act,kf_act)./i_edge(ic_act,jf_act,kf_act);
i_edge_y(ic_act,jf_act,kf_act) = di_edge_y(ic_act,jf_act,kf_act)./i_edge(ic_act,jf_act,kf_act);
i_edge_z(ic_act,jf_act,kf_act) = di_edge_z(ic_act,jf_act,kf_act)./i_edge(ic_act,jf_act,kf_act);

dxd(ic_act,jf_act,kf_act) = 0.25.*(x(ic_act,jf_act,kf_act+1) + x(ic_act+1,jf_act,kf_act+1) - x(ic_act,jf_act,kf_act-1) - x(ic_act+1,jf_act,kf_act-1));
dyd(ic_act,jf_act,kf_act) = 0.25.*(y(ic_act,jf_act,kf_act+1) + y(ic_act+1,jf_act,kf_act+1) - y(ic_act,jf_act,kf_act-1) - y(ic_act+1,jf_act,kf_act-1));
dzd(ic_act,jf_act,kf_act) = 0.25.*(z(ic_act,jf_act,kf_act+1) + z(ic_act+1,jf_act,kf_act+1) - z(ic_act,jf_act,kf_act-1) - z(ic_act+1,jf_act,kf_act-1));

[dx_0(ic_act,jf_act,kf_act),dy_0(ic_act,jf_act,kf_act),dz_0(ic_act,jf_act,kf_act)] = ...
    CROSS2(dxd(ic_act,jf_act,kf_act),dyd(ic_act,jf_act,kf_act),dzd(ic_act,jf_act,kf_act),...
           i_edge_x(ic_act,jf_act,kf_act),i_edge_y(ic_act,jf_act,kf_act),i_edge_z(ic_act,jf_act,kf_act));

dx_0_norm(ic_act,jf_act,kf_act) = sqrt(dx_0(ic_act,jf_act,kf_act).^2+dy_0(ic_act,jf_act,kf_act).^2+dz_0(ic_act,jf_act,kf_act).^2);
dx_0(ic_act,jf_act,kf_act) = dx_0(ic_act,jf_act,kf_act)./dx_0_norm(ic_act,jf_act,kf_act);
dy_0(ic_act,jf_act,kf_act) = dy_0(ic_act,jf_act,kf_act)./dx_0_norm(ic_act,jf_act,kf_act);
dz_0(ic_act,jf_act,kf_act) = dz_0(ic_act,jf_act,kf_act)./dx_0_norm(ic_act,jf_act,kf_act);

[dx_1(ic_act,jf_act,kf_act),dy_1(ic_act,jf_act,kf_act),dz_1(ic_act,jf_act,kf_act)] = ...
    CROSS2(i_edge_x(ic_act,jf_act,kf_act),i_edge_y(ic_act,jf_act,kf_act),i_edge_z(ic_act,jf_act,kf_act),...
           dx_0(ic_act,jf_act,kf_act),dy_0(ic_act,jf_act,kf_act),dz_0(ic_act,jf_act,kf_act));

vec_i(ic_act,jf_act,kf_act,7) = i_edge_x(ic_act,jf_act,kf_act);
vec_i(ic_act,jf_act,kf_act,8) = i_edge_y(ic_act,jf_act,kf_act);
vec_i(ic_act,jf_act,kf_act,9) = i_edge_z(ic_act,jf_act,kf_act);

vec_i(ic_act,jf_act,kf_act,1) = dx_0(ic_act,jf_act,kf_act);
vec_i(ic_act,jf_act,kf_act,2) = dy_0(ic_act,jf_act,kf_act);
vec_i(ic_act,jf_act,kf_act,3) = dz_0(ic_act,jf_act,kf_act);

vec_i(ic_act,jf_act,kf_act,4) = dx_1(ic_act,jf_act,kf_act);
vec_i(ic_act,jf_act,kf_act,5) = dy_1(ic_act,jf_act,kf_act);
vec_i(ic_act,jf_act,kf_act,6) = dz_1(ic_act,jf_act,kf_act);
     
% j-edge metrics, dimension (if_act, jc_act, kf_act)
% vec_k(:,:,:,1:9)
% (7-9) vector in the j-edge direction
% (1-3) vector in the direction of sweep
% (4-6) vector in the right handed perp to both j and sweep
% 
% (dj_edge_x, dj_edge_y, dj_edge_z) is the j-edge vector and 
% (j_edge_x, j_edge_y, j_edge_z) is the j-edge unit vector

dj_edge_x(if_act,jc_act,kf_act) = x(if_act,jc_act+1,kf_act) - x(if_act,jc_act,kf_act);
dj_edge_y(if_act,jc_act,kf_act) = y(if_act,jc_act+1,kf_act) - y(if_act,jc_act,kf_act);
dj_edge_z(if_act,jc_act,kf_act) = z(if_act,jc_act+1,kf_act) - z(if_act,jc_act,kf_act);
j_edge(if_act,jc_act,kf_act) = sqrt(dj_edge_x(if_act,jc_act,kf_act).^2 + dj_edge_y(if_act,jc_act,kf_act).^2 + dj_edge_z(if_act,jc_act,kf_act).^2);

j_edge_x(if_act,jc_act,kf_act) = dj_edge_x(if_act,jc_act,kf_act)./j_edge(if_act,jc_act,kf_act);
j_edge_y(if_act,jc_act,kf_act) = dj_edge_y(if_act,jc_act,kf_act)./j_edge(if_act,jc_act,kf_act);
j_edge_z(if_act,jc_act,kf_act) = dj_edge_z(if_act,jc_act,kf_act)./j_edge(if_act,jc_act,kf_act);

dxd(if_act,jc_act,kf_act) = 0.25.*(x(if_act+1,jc_act,kf_act) + x(if_act+1,jc_act+1,kf_act) - x(if_act-1,jc_act,kf_act) - x(if_act-1,jc_act+1,kf_act));
dyd(if_act,jc_act,kf_act) = 0.25.*(y(if_act+1,jc_act,kf_act) + y(if_act+1,jc_act+1,kf_act) - y(if_act-1,jc_act,kf_act) - y(if_act-1,jc_act+1,kf_act));
dzd(if_act,jc_act,kf_act) = 0.25.*(z(if_act+1,jc_act,kf_act) + z(if_act+1,jc_act+1,kf_act) - z(if_act-1,jc_act,kf_act) - z(if_act-1,jc_act+1,kf_act));

[dx_0(if_act,jc_act,kf_act),dy_0(if_act,jc_act,kf_act),dz_0(if_act,jc_act,kf_act)] = ...
    CROSS2(dxd(if_act,jc_act,kf_act),dyd(if_act,jc_act,kf_act),dzd(if_act,jc_act,kf_act),...
           j_edge_x(if_act,jc_act,kf_act),j_edge_y(if_act,jc_act,kf_act),j_edge_z(if_act,jc_act,kf_act));
dx_0_norm(if_act,jc_act,kf_act) = sqrt(dx_0(if_act,jc_act,kf_act).^2+dy_0(if_act,jc_act,kf_act).^2+dz_0(if_act,jc_act,kf_act).^2);
dx_0(if_act,jc_act,kf_act) = dx_0(if_act,jc_act,kf_act)./dx_0_norm(if_act,jc_act,kf_act);
dy_0(if_act,jc_act,kf_act) = dy_0(if_act,jc_act,kf_act)./dx_0_norm(if_act,jc_act,kf_act);
dz_0(if_act,jc_act,kf_act) = dz_0(if_act,jc_act,kf_act)./dx_0_norm(if_act,jc_act,kf_act);
[dx_1(if_act,jc_act,kf_act),dy_1(if_act,jc_act,kf_act),dz_1(if_act,jc_act,kf_act)] = ...
    CROSS2(j_edge_x(if_act,jc_act,kf_act),j_edge_y(if_act,jc_act,kf_act),j_edge_z(if_act,jc_act,kf_act),...
           dx_0(if_act,jc_act,kf_act),dy_0(if_act,jc_act,kf_act),dz_0(if_act,jc_act,kf_act));
       
% vec_j matrix is used to transform the velocity and Bupwind vector to the
% k-normal coordinate
vec_j(if_act,jc_act,kf_act,7) = j_edge_x(if_act,jc_act,kf_act);
vec_j(if_act,jc_act,kf_act,8) = j_edge_y(if_act,jc_act,kf_act);
vec_j(if_act,jc_act,kf_act,9) = j_edge_z(if_act,jc_act,kf_act);

vec_j(if_act,jc_act,kf_act,1) = dx_0(if_act,jc_act,kf_act);
vec_j(if_act,jc_act,kf_act,2) = dy_0(if_act,jc_act,kf_act);
vec_j(if_act,jc_act,kf_act,3) = dz_0(if_act,jc_act,kf_act);

vec_j(if_act,jc_act,kf_act,4) = dx_1(if_act,jc_act,kf_act);
vec_j(if_act,jc_act,kf_act,5) = dy_1(if_act,jc_act,kf_act);
vec_j(if_act,jc_act,kf_act,6) = dz_1(if_act,jc_act,kf_act);

%~~~~~~~~~~~~END edge-center edge transform ~~~~~~~~~~~~~~
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the transform coordinates at cell edges (k-edge) for electric field
% calculations, used in getEfields() function

% k-edge metric calculation for transforming magnetic field into the
% (e1,e2) coordinate system for vxB calculation
% high order i-face area. NOTE in very distorted grid '8th' may not work 
[ifaceA_Kedge, ~] = reconstruct_3DV(i_face_area,if_act,jf_act,kc_act,PDMB, 2,'8th'); % dimension (if_act,jf_act,kc_act)
[jfaceA_Kedge, ~] = reconstruct_3DV(j_face_area,if_act,jf_act,kc_act,PDMB, 1,'8th'); % dimension (if_act,jf_act,kc_act)%

% high order estimations for i-face and j-face normal vector. 
% These are the estimated directions for the B_avg (field) at k-edges 
% first interpolate i_face_norm vector in the j-direction
[ifaceNx_Kedge, ~] = reconstruct_3DV(i_face_norm_x.*i_face_area,if_act,jf_act,kc_act,PDMB, 2,'8th'); % dimension (if_act,jf_act,kc_act)
[ifaceNy_Kedge, ~] = reconstruct_3DV(i_face_norm_y.*i_face_area,if_act,jf_act,kc_act,PDMB, 2,'8th'); % dimension (if_act,jf_act,kc_act)
[ifaceNz_Kedge, ~] = reconstruct_3DV(i_face_norm_z.*i_face_area,if_act,jf_act,kc_act,PDMB, 2,'8th'); % dimension (if_act,jf_act,kc_act)
ifaceNx_Kedge = ifaceNx_Kedge./ifaceA_Kedge;
ifaceNy_Kedge = ifaceNy_Kedge./ifaceA_Kedge;
ifaceNz_Kedge = ifaceNz_Kedge./ifaceA_Kedge;
% then interpolate j_face_norm vector in the i-direction
[jfaceNx_Kedge, ~] = reconstruct_3DV(j_face_norm_x.*j_face_area,if_act,jf_act,kc_act,PDMB, 1,'8th'); % dimension (if_act,jf_act,kc_act)
[jfaceNy_Kedge, ~] = reconstruct_3DV(j_face_norm_y.*j_face_area,if_act,jf_act,kc_act,PDMB, 1,'8th'); % dimension (if_act,jf_act,kc_act)
[jfaceNz_Kedge, ~] = reconstruct_3DV(j_face_norm_z.*j_face_area,if_act,jf_act,kc_act,PDMB, 1,'8th'); % dimension (if_act,jf_act,kc_act)
jfaceNx_Kedge = jfaceNx_Kedge./jfaceA_Kedge;
jfaceNy_Kedge = jfaceNy_Kedge./jfaceA_Kedge;
jfaceNz_Kedge = jfaceNz_Kedge./jfaceA_Kedge;

% transform to the K-edge coordinate system using the vec_k matrix
%  B_avg_e1 = DETI.*( ynqk_K.*bi_avg - ynqj_K.*bj_avg) + BX0_K;
%  B_avg_e2 = DETI.*(-xnqk_K.*bi_avg + xnqj_K.*bj_avg) + BY0_K;
% bi direction is ifaceN
% bj direction is jfaceN
% ifaceN vector - x-component in the k-edge coordinate system
dxti_K(if_act,jf_act,kc_act) = ifaceNx_Kedge(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,1) + ...
                               ifaceNy_Kedge(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,2) + ...
                               ifaceNz_Kedge(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,3);
% ifaceN vector - y-component in the k-edge coordinate system
dyti_K(if_act,jf_act,kc_act) = ifaceNx_Kedge(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,4) + ...
                               ifaceNy_Kedge(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,5) + ...
                               ifaceNz_Kedge(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,6);
% normalize to unit vector                         
xnqi_K(if_act,jf_act,kc_act) = dxti_K(if_act,jf_act,kc_act)./sqrt(dxti_K(if_act,jf_act,kc_act).^2 + dyti_K(if_act,jf_act,kc_act).^2);
ynqi_K(if_act,jf_act,kc_act) = dyti_K(if_act,jf_act,kc_act)./sqrt(dxti_K(if_act,jf_act,kc_act).^2 + dyti_K(if_act,jf_act,kc_act).^2);

% jfaceN vector - x-component in the k-edge coordinate system
dxtj_K(if_act,jf_act,kc_act) = jfaceNx_Kedge(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,1) + ...
                               jfaceNy_Kedge(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,2) + ...
                               jfaceNz_Kedge(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,3);
% ifaceN vector - y-component in the k-edge coordinate system
dytj_K(if_act,jf_act,kc_act) = jfaceNx_Kedge(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,4) + ...
                               jfaceNy_Kedge(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,5) + ...
                               jfaceNz_Kedge(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,6);
% normalize to unit vector                         
xnqj_K(if_act,jf_act,kc_act) = dxtj_K(if_act,jf_act,kc_act)./sqrt(dxtj_K(if_act,jf_act,kc_act).^2 + dytj_K(if_act,jf_act,kc_act).^2);
ynqj_K(if_act,jf_act,kc_act) = dytj_K(if_act,jf_act,kc_act)./sqrt(dxtj_K(if_act,jf_act,kc_act).^2 + dytj_K(if_act,jf_act,kc_act).^2);
DETK = 1./(xnqi_K.*ynqj_K - xnqj_K.*ynqi_K);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the transform coordinates at cell edges (i-edge)

% high order j-face,k-face area: j_face_area interpolated in the k-direction
%                                k_face_area interpolated in the j-direction
[jfaceA_Iedge, ~] = reconstruct_3DV(j_face_area,ic_act,jf_act,kf_act,PDMB, 3,'8th'); % dimension (if_act,jf_act,kc_act)
[kfaceA_Iedge, ~] = reconstruct_3DV(k_face_area,ic_act,jf_act,kf_act,PDMB, 2,'8th'); % dimension (if_act,jf_act,kc_act)

% high order j-face normal vector 
% These are the estimated directions for the B_avg (field) at i-edges 
% first interpolate j_face_norm vector in the k-direction
[jfaceNx_Iedge, ~] = reconstruct_3DV(j_face_norm_x.*j_face_area,ic_act,jf_act,kf_act,PDMB, 3,'8th'); % dimension (if_act,jf_act,kc_act)
[jfaceNy_Iedge, ~] = reconstruct_3DV(j_face_norm_y.*j_face_area,ic_act,jf_act,kf_act,PDMB, 3,'8th'); % dimension (if_act,jf_act,kc_act)
[jfaceNz_Iedge, ~] = reconstruct_3DV(j_face_norm_z.*j_face_area,ic_act,jf_act,kf_act,PDMB, 3,'8th'); % dimension (if_act,jf_act,kc_act)
jfaceNx_Iedge = jfaceNx_Iedge./jfaceA_Iedge;
jfaceNy_Iedge = jfaceNy_Iedge./jfaceA_Iedge;
jfaceNz_Iedge = jfaceNz_Iedge./jfaceA_Iedge;
% then interpolate k_face_norm vector in the j-direction
[kfaceNx_Iedge, ~] = reconstruct_3DV(k_face_norm_x.*k_face_area,ic_act,jf_act,kf_act,PDMB, 2,'8th'); % dimension (if_act,jf_act,kc_act)
[kfaceNy_Iedge, ~] = reconstruct_3DV(k_face_norm_y.*k_face_area,ic_act,jf_act,kf_act,PDMB, 2,'8th'); % dimension (if_act,jf_act,kc_act)
[kfaceNz_Iedge, ~] = reconstruct_3DV(k_face_norm_z.*k_face_area,ic_act,jf_act,kf_act,PDMB, 2,'8th'); % dimension (if_act,jf_act,kc_act)
kfaceNx_Iedge = kfaceNx_Iedge./kfaceA_Iedge;
kfaceNy_Iedge = kfaceNy_Iedge./kfaceA_Iedge;
kfaceNz_Iedge = kfaceNz_Iedge./kfaceA_Iedge;
%
% transform coefficients for the (e1,e2) system along I-edges
% the transform should look like this:
%  B_avg_e1 = DETI.*( ynqk_I.*bj_avg - ynqj_I.*bk_avg) + BX0_I;
%  B_avg_e2 = DETI.*(-xnqk_I.*bj_avg + xnqj_I.*bk_avg) + BY0_I;
% jfaceN vector - x-component in the i-edge coordinate system
dxtj_I(ic_act,jf_act,kf_act) = jfaceNx_Iedge(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,1) + ...
                               jfaceNy_Iedge(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,2) + ...
                               jfaceNz_Iedge(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,3);
% jfaceN vector - y-component in the i-edge coordinate system
dytj_I(ic_act,jf_act,kf_act) = jfaceNx_Iedge(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,4) + ...
                               jfaceNy_Iedge(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,5) + ...
                               jfaceNz_Iedge(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,6);
% normalize to unit vector                         
xnqj_I(ic_act,jf_act,kf_act) = dxtj_I(ic_act,jf_act,kf_act)./sqrt(dxtj_I(ic_act,jf_act,kf_act).^2 + dytj_I(ic_act,jf_act,kf_act).^2);
ynqj_I(ic_act,jf_act,kf_act) = dytj_I(ic_act,jf_act,kf_act)./sqrt(dxtj_I(ic_act,jf_act,kf_act).^2 + dytj_I(ic_act,jf_act,kf_act).^2);

% kfaceN vector - x-component in the i-edge coordinate system
dxtk_I(ic_act,jf_act,kf_act) = kfaceNx_Iedge(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,1) + ...
                               kfaceNy_Iedge(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,2) + ...
                               kfaceNz_Iedge(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,3);
% kfaceN vector - y-component in the i-edge coordinate system
dytk_I(ic_act,jf_act,kf_act) = kfaceNx_Iedge(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,4) + ...
                               kfaceNy_Iedge(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,5) + ...
                               kfaceNz_Iedge(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,6);
% normalize to unit vector                         
xnqk_I(ic_act,jf_act,kf_act) = dxtk_I(ic_act,jf_act,kf_act)./sqrt(dxtk_I(ic_act,jf_act,kf_act).^2 + dytk_I(ic_act,jf_act,kf_act).^2);
ynqk_I(ic_act,jf_act,kf_act) = dytk_I(ic_act,jf_act,kf_act)./sqrt(dxtk_I(ic_act,jf_act,kf_act).^2 + dytk_I(ic_act,jf_act,kf_act).^2);
DETI = 1./(xnqj_I.*ynqk_I - xnqk_I.*ynqj_I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the transform coordinates at cell edges (j-edge)

% high order k,i-face area
[kfaceA_Jedge, ~] = reconstruct_3DV(k_face_area,if_act,jc_act,kf_act,PDMB, 1,'8th'); % dimension (if_act,jf_act,kc_act)
[ifaceA_Jedge, ~] = reconstruct_3DV(i_face_area,if_act,jc_act,kf_act,PDMB, 3,'8th'); % dimension (if_act,jf_act,kc_act)

% high order k,i-face normal vector
% These are the estimated directions for the B_avg (field) at j-edges 
% first interpolate k_face_norm vector in the i-direction
[kfaceNx_Jedge, ~] = reconstruct_3DV(k_face_norm_x.*k_face_area,if_act,jc_act,kf_act,PDMB, 1,'8th'); % dimension (if_act,jf_act,kc_act)
[kfaceNy_Jedge, ~] = reconstruct_3DV(k_face_norm_y.*k_face_area,if_act,jc_act,kf_act,PDMB, 1,'8th'); % dimension (if_act,jf_act,kc_act)
[kfaceNz_Jedge, ~] = reconstruct_3DV(k_face_norm_z.*k_face_area,if_act,jc_act,kf_act,PDMB, 1,'8th'); % dimension (if_act,jf_act,kc_act)
kfaceNx_Jedge = kfaceNx_Jedge./kfaceA_Jedge;
kfaceNy_Jedge = kfaceNy_Jedge./kfaceA_Jedge;
kfaceNz_Jedge = kfaceNz_Jedge./kfaceA_Jedge;
% then interpolate i_face_norm vector in the k-direction
[ifaceNx_Jedge, ~] = reconstruct_3DV(i_face_norm_x.*i_face_area,if_act,jc_act,kf_act,PDMB, 3,'8th'); % dimension (if_act,jf_act,kc_act)
[ifaceNy_Jedge, ~] = reconstruct_3DV(i_face_norm_y.*i_face_area,if_act,jc_act,kf_act,PDMB, 3,'8th'); % dimension (if_act,jf_act,kc_act)
[ifaceNz_Jedge, ~] = reconstruct_3DV(i_face_norm_z.*i_face_area,if_act,jc_act,kf_act,PDMB, 3,'8th'); % dimension (if_act,jf_act,kc_act)
ifaceNx_Jedge = ifaceNx_Jedge./ifaceA_Jedge;
ifaceNy_Jedge = ifaceNy_Jedge./ifaceA_Jedge;
ifaceNz_Jedge = ifaceNz_Jedge./ifaceA_Jedge;

% transform coefficients for the (e1,e2) system along J-edges
% the transform should look like this:
%  B_avg_e1 = DETJ.*( ynqk_J.*bk_avg - ynqj_J.*bi_avg) + BX0_J;
%  B_avg_e2 = DETJ.*(-xnqk_J.*bk_avg + xnqj_J.*bi_avg) + BY0_J;
% bk direction is kfaceN
% bi direction is ifaceN
% kfaceN vector - x-component in the k-edge coordinate system
dxtk_J(if_act,jc_act,kf_act) = kfaceNx_Jedge(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,1) + ...
                               kfaceNy_Jedge(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,2) + ...
                               kfaceNz_Jedge(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,3);
% kfaceN vector - y-component in the k-edge coordinate system
dytk_J(if_act,jc_act,kf_act) = kfaceNx_Jedge(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,4) + ...
                               kfaceNy_Jedge(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,5) + ...
                               kfaceNz_Jedge(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,6);
% normalize to unit vector                         
xnqk_J(if_act,jc_act,kf_act) = dxtk_J(if_act,jc_act,kf_act)./sqrt(dxtk_J(if_act,jc_act,kf_act).^2 + dytk_J(if_act,jc_act,kf_act).^2);
ynqk_J(if_act,jc_act,kf_act) = dytk_J(if_act,jc_act,kf_act)./sqrt(dxtk_J(if_act,jc_act,kf_act).^2 + dytk_J(if_act,jc_act,kf_act).^2);

% ifaceN vector - x-component in the j-edge coordinate system
dxti_J(if_act,jc_act,kf_act) = ifaceNx_Jedge(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,1) + ...
                               ifaceNy_Jedge(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,2) + ...
                               ifaceNz_Jedge(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,3);
% ifaceN vector - y-component in the k-edge coordinate system
dyti_J(if_act,jc_act,kf_act) = ifaceNx_Jedge(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,4) + ...
                               ifaceNy_Jedge(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,5) + ...
                               ifaceNz_Jedge(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,6);
% normalize to unit vector                         
xnqi_J(if_act,jc_act,kf_act) = dxti_J(if_act,jc_act,kf_act)./sqrt(dxti_J(if_act,jc_act,kf_act).^2 + dyti_J(if_act,jc_act,kf_act).^2);
ynqi_J(if_act,jc_act,kf_act) = dyti_J(if_act,jc_act,kf_act)./sqrt(dxti_J(if_act,jc_act,kf_act).^2 + dyti_J(if_act,jc_act,kf_act).^2);
DETJ = 1./(xnqk_J.*ynqi_J - xnqi_J.*ynqk_J);

% END of edge metric
end