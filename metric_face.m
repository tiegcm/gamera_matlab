function [TRANS_i,TRANS_j,TRANS_k,...
          i_face_norm_x,i_face_norm_y,i_face_norm_z,...
          j_face_norm_x,j_face_norm_y,j_face_norm_z,...
          k_face_norm_x,k_face_norm_y,k_face_norm_z]=metric_face(x,y,z)

% This code computes the transform matrices for the face-normal coordinate
% system, together with the unit vectors for the face normal directions at
% face centers
%
% INPUT: x,y,z - cell corners      
      
global I J K Ip1 Jp1 Kp1

disp('Computing Face-Normal Coordinates');

% average i-unit vector and j-unit vector at i-j face
di_face_x(I,J,Kp1) = 0.5*(x(I+1,J+1,Kp1) - x(I,J+1,Kp1)) + 0.5*(x(I+1,J,Kp1) - x(I,J,Kp1));
di_face_y(I,J,Kp1) = 0.5*(y(I+1,J+1,Kp1) - y(I,J+1,Kp1)) + 0.5*(y(I+1,J,Kp1) - y(I,J,Kp1));
di_face_z(I,J,Kp1) = 0.5*(z(I+1,J+1,Kp1) - z(I,J+1,Kp1)) + 0.5*(z(I+1,J,Kp1) - z(I,J,Kp1));

dj_face_x(I,J,Kp1) = 0.5*(x(I+1,J+1,Kp1) - x(I+1,J,Kp1)) + 0.5*(x(I,J+1,Kp1) - x(I,J,Kp1));
dj_face_y(I,J,Kp1) = 0.5*(y(I+1,J+1,Kp1) - y(I+1,J,Kp1)) + 0.5*(y(I,J+1,Kp1) - y(I,J,Kp1));
dj_face_z(I,J,Kp1) = 0.5*(z(I+1,J+1,Kp1) - z(I+1,J,Kp1)) + 0.5*(z(I,J+1,Kp1) - z(I,J,Kp1));

di_face(I,J,Kp1) = sqrt(di_face_x(I,J,Kp1).^2+di_face_y(I,J,Kp1).^2+di_face_z(I,J,Kp1).^2);
i_face_x(I,J,Kp1) = di_face_x(I,J,Kp1)./di_face(I,J,Kp1);
i_face_y(I,J,Kp1) = di_face_y(I,J,Kp1)./di_face(I,J,Kp1);
i_face_z(I,J,Kp1) = di_face_z(I,J,Kp1)./di_face(I,J,Kp1);

dj_face(I,J,Kp1) = sqrt(dj_face_x(I,J,Kp1).^2+dj_face_y(I,J,Kp1).^2+dj_face_z(I,J,Kp1).^2);
j_face_x(I,J,Kp1) = dj_face_x(I,J,Kp1)./dj_face(I,J,Kp1);
j_face_y(I,J,Kp1) = dj_face_y(I,J,Kp1)./dj_face(I,J,Kp1);
j_face_z(I,J,Kp1) = dj_face_z(I,J,Kp1)./dj_face(I,J,Kp1);

% k-face normal direction - i_face cross j_face. 
kfx(I,J,Kp1) =  di_face_y(I,J,Kp1).* dj_face_z(I,J,Kp1) - di_face_z(I,J,Kp1).* dj_face_y(I,J,Kp1);
kfy(I,J,Kp1) =-(di_face_x(I,J,Kp1).* dj_face_z(I,J,Kp1) - di_face_z(I,J,Kp1).* dj_face_x(I,J,Kp1));
kfz(I,J,Kp1) =  di_face_x(I,J,Kp1).* dj_face_y(I,J,Kp1) - di_face_y(I,J,Kp1).* dj_face_x(I,J,Kp1);

% unit k-face normal direction
k_face_norm_x(I,J,Kp1) =  kfx(I,J,Kp1)./sqrt(kfx(I,J,Kp1).^2 + kfy(I,J,Kp1).^2 + kfz(I,J,Kp1).^2);
k_face_norm_y(I,J,Kp1) =  kfy(I,J,Kp1)./sqrt(kfx(I,J,Kp1).^2 + kfy(I,J,Kp1).^2 + kfz(I,J,Kp1).^2);
k_face_norm_z(I,J,Kp1) =  kfz(I,J,Kp1)./sqrt(kfx(I,J,Kp1).^2 + kfy(I,J,Kp1).^2 + kfz(I,J,Kp1).^2);

% j-like face direction, which is k cross i, forms an orthogonal system
[j_face_x(I,J,Kp1),j_face_y(I,J,Kp1),j_face_z(I,J,Kp1)] = ...
    CROSS2(k_face_norm_x(I,J,Kp1),k_face_norm_y(I,J,Kp1),k_face_norm_z(I,J,Kp1), ...
           i_face_x(I,J,Kp1),i_face_y(I,J,Kp1),i_face_z(I,J,Kp1));

% define the TRANS_k matrix as used in LFM and gamera
% to transpose to the i,j,k coordinates, use
% vk(:,:,:) = TRANS_k(:,:,:,1).*vx(:,:,:) +  TRANS_k(:,:,:,2).*vy(:,:,:)+ TRANS_k(:,:,:,3).*vz(:,:,:)
% vi(:,:,:) = TRANS_k(:,:,:,4).*vx(:,:,:) +  TRANS_k(:,:,:,5).*vy(:,:,:)+ TRANS_k(:,:,:,6).*vz(:,:,:)
% vj(:,:,:) = TRANS_k(:,:,:,7).*vx(:,:,:) +  TRANS_k(:,:,:,8).*vy(:,:,:)+ TRANS_k(:,:,:,9).*vz(:,:,:)
% since TRANS_k is a unity matrix, the inverse is just TRANS_k', so
% vx(:,:,:) = TRANS_k(:,:,:,1).*vk(:,:,:) +  TRANS_k(:,:,:,4).*vi(:,:,:)+ TRANS_k(:,:,:,7).*vj(:,:,:)
% vy(:,:,:) = TRANS_k(:,:,:,2).*vk(:,:,:) +  TRANS_k(:,:,:,5).*vi(:,:,:)+ TRANS_k(:,:,:,8).*vj(:,:,:)
% vz(:,:,:) = TRANS_k(:,:,:,3).*vk(:,:,:) +  TRANS_k(:,:,:,6).*vi(:,:,:)+ TRANS_k(:,:,:,9).*vj(:,:,:)
%
% Note that the actual direction k,j,i doesn't matter in the 1-D flux function,
% i.e., for the flux function, as long as the normal velocity is correctly used 
% Thus in the code, all the normal velocity conponent is called "vi", then
% the next one is called "vj" and the last one is called "vk"
TRANS_k(I,J,Kp1,1) = k_face_norm_x(I,J,Kp1);
TRANS_k(I,J,Kp1,2) = k_face_norm_y(I,J,Kp1);
TRANS_k(I,J,Kp1,3) = k_face_norm_z(I,J,Kp1);
TRANS_k(I,J,Kp1,4) = i_face_x(I,J,Kp1);
TRANS_k(I,J,Kp1,5) = i_face_y(I,J,Kp1);
TRANS_k(I,J,Kp1,6) = i_face_z(I,J,Kp1);
TRANS_k(I,J,Kp1,7) = j_face_x(I,J,Kp1);
TRANS_k(I,J,Kp1,8) = j_face_y(I,J,Kp1);
TRANS_k(I,J,Kp1,9) = j_face_z(I,J,Kp1);

% average j-unit vector and k-unit vector at j-k face
dj_face_x(Ip1,J,K) = 0.5*(x(Ip1,J+1,K+1) - x(Ip1,J,K+1)) + 0.5*(x(Ip1,J+1,K) - x(Ip1,J,K));
dj_face_y(Ip1,J,K) = 0.5*(y(Ip1,J+1,K+1) - y(Ip1,J,K+1)) + 0.5*(y(Ip1,J+1,K) - y(Ip1,J,K));
dj_face_z(Ip1,J,K) = 0.5*(z(Ip1,J+1,K+1) - z(Ip1,J,K+1)) + 0.5*(z(Ip1,J+1,K) - z(Ip1,J,K));

dk_face_x(Ip1,J,K) = 0.5*(x(Ip1,J+1,K+1) - x(Ip1,J+1,K)) + 0.5*(x(Ip1,J,K+1) - x(Ip1,J,K));
dk_face_y(Ip1,J,K) = 0.5*(y(Ip1,J+1,K+1) - y(Ip1,J+1,K)) + 0.5*(y(Ip1,J,K+1) - y(Ip1,J,K));
dk_face_z(Ip1,J,K) = 0.5*(z(Ip1,J+1,K+1) - z(Ip1,J+1,K)) + 0.5*(z(Ip1,J,K+1) - z(Ip1,J,K));

dj_face(Ip1,J,K) = sqrt(dj_face_x(Ip1,J,K).^2 + dj_face_y(Ip1,J,K).^2 + dj_face_z(Ip1,J,K).^2);
j_face_x(Ip1,J,K) = dj_face_x(Ip1,J,K)./dj_face(Ip1,J,K);
j_face_y(Ip1,J,K) = dj_face_y(Ip1,J,K)./dj_face(Ip1,J,K);
j_face_z(Ip1,J,K) = dj_face_z(Ip1,J,K)./dj_face(Ip1,J,K);

% i-face normal direction - i_face cross j_face
ifx(Ip1,J,K) =  dj_face_y(Ip1,J,K).* dk_face_z(Ip1,J,K) - dj_face_z(Ip1,J,K).* dk_face_y(Ip1,J,K);
ify(Ip1,J,K) =-(dj_face_x(Ip1,J,K).* dk_face_z(Ip1,J,K) - dj_face_z(Ip1,J,K).* dk_face_x(Ip1,J,K));
ifz(Ip1,J,K) =  dj_face_x(Ip1,J,K).* dk_face_y(Ip1,J,K) - dj_face_y(Ip1,J,K).* dk_face_x(Ip1,J,K);

% normalized i-face normal direction
i_face_norm_x(Ip1,J,K) =  ifx(Ip1,J,K)./sqrt(ifx(Ip1,J,K).^2 + ify(Ip1,J,K).^2 + ifz(Ip1,J,K).^2);
i_face_norm_y(Ip1,J,K) =  ify(Ip1,J,K)./sqrt(ifx(Ip1,J,K).^2 + ify(Ip1,J,K).^2 + ifz(Ip1,J,K).^2);
i_face_norm_z(Ip1,J,K) =  ifz(Ip1,J,K)./sqrt(ifx(Ip1,J,K).^2 + ify(Ip1,J,K).^2 + ifz(Ip1,J,K).^2);

% the thir direction that forms an orthogonal system
[k_face_x(Ip1,J,K),k_face_y(Ip1,J,K),k_face_z(Ip1,J,K)] = ...
    CROSS2(i_face_norm_x(Ip1,J,K),i_face_norm_y(Ip1,J,K),i_face_norm_z(Ip1,J,K), ...
    j_face_x(Ip1,J,K),j_face_y(Ip1,J,K),j_face_z(Ip1,J,K));

TRANS_i(Ip1,J,K,1) = i_face_norm_x(Ip1,J,K);
TRANS_i(Ip1,J,K,2) = i_face_norm_y(Ip1,J,K);
TRANS_i(Ip1,J,K,3) = i_face_norm_z(Ip1,J,K);
TRANS_i(Ip1,J,K,4) = j_face_x(Ip1,J,K);
TRANS_i(Ip1,J,K,5) = j_face_y(Ip1,J,K);
TRANS_i(Ip1,J,K,6) = j_face_z(Ip1,J,K);
TRANS_i(Ip1,J,K,7) = k_face_x(Ip1,J,K);
TRANS_i(Ip1,J,K,8) = k_face_y(Ip1,J,K);
TRANS_i(Ip1,J,K,9) = k_face_z(Ip1,J,K);

% calculate (k,i) face area, dk_face_center cross di_face_center
% average k-unit vector and i-unit vector at k-i face
dk_face_x(I,Jp1,K) = 0.5*(x(I+1,Jp1,K+1) - x(I+1,Jp1,K)) + 0.5*(x(I,Jp1,K+1) - x(I,Jp1,K));
dk_face_y(I,Jp1,K) = 0.5*(y(I+1,Jp1,K+1) - y(I+1,Jp1,K)) + 0.5*(y(I,Jp1,K+1) - y(I,Jp1,K));
dk_face_z(I,Jp1,K) = 0.5*(z(I+1,Jp1,K+1) - z(I+1,Jp1,K)) + 0.5*(z(I,Jp1,K+1) - z(I,Jp1,K));

di_face_x(I,Jp1,K) = 0.5*(x(I+1,Jp1,K+1) - x(I,Jp1,K+1)) + 0.5*(x(I+1,Jp1,K) - x(I,Jp1,K));
di_face_y(I,Jp1,K) = 0.5*(y(I+1,Jp1,K+1) - y(I,Jp1,K+1)) + 0.5*(y(I+1,Jp1,K) - y(I,Jp1,K));
di_face_z(I,Jp1,K) = 0.5*(z(I+1,Jp1,K+1) - z(I,Jp1,K+1)) + 0.5*(z(I+1,Jp1,K) - z(I,Jp1,K));

dk_face(I,Jp1,K) = sqrt(dk_face_x(I,Jp1,K).^2+dk_face_y(I,Jp1,K).^2+dk_face_z(I,Jp1,K).^2);
k_face_x(I,Jp1,K) = dk_face_x(I,Jp1,K)./dk_face(I,Jp1,K);
k_face_y(I,Jp1,K) = dk_face_y(I,Jp1,K)./dk_face(I,Jp1,K);
k_face_z(I,Jp1,K) = dk_face_z(I,Jp1,K)./dk_face(I,Jp1,K);

jfx(I,Jp1,K) =  dk_face_y(I,Jp1,K).* di_face_z(I,Jp1,K) - dk_face_z(I,Jp1,K).* di_face_y(I,Jp1,K);
jfy(I,Jp1,K) =-(dk_face_x(I,Jp1,K).* di_face_z(I,Jp1,K) - dk_face_z(I,Jp1,K).* di_face_x(I,Jp1,K));
jfz(I,Jp1,K) =  dk_face_x(I,Jp1,K).* di_face_y(I,Jp1,K) - dk_face_y(I,Jp1,K).* di_face_x(I,Jp1,K);

% j-face normal direction - k_face cross i_face
j_face_norm_x(I,Jp1,K) = jfx(I,Jp1,K)./sqrt(jfx(I,Jp1,K).^2 + jfy(I,Jp1,K).^2 + jfz(I,Jp1,K).^2);
j_face_norm_y(I,Jp1,K) = jfy(I,Jp1,K)./sqrt(jfx(I,Jp1,K).^2 + jfy(I,Jp1,K).^2 + jfz(I,Jp1,K).^2);
j_face_norm_z(I,Jp1,K) = jfz(I,Jp1,K)./sqrt(jfx(I,Jp1,K).^2 + jfy(I,Jp1,K).^2 + jfz(I,Jp1,K).^2);

[i_face_x(I,Jp1,K), i_face_y(I,Jp1,K),i_face_z(I,Jp1,K),] = ...
    CROSS2(j_face_norm_x(I,Jp1,K),j_face_norm_y(I,Jp1,K),j_face_norm_z(I,Jp1,K),...
           k_face_x(I,Jp1,K),k_face_y(I,Jp1,K),k_face_z(I,Jp1,K));

TRANS_j(I,Jp1,K,1) = j_face_norm_x(I,Jp1,K);
TRANS_j(I,Jp1,K,2) = j_face_norm_y(I,Jp1,K);
TRANS_j(I,Jp1,K,3) = j_face_norm_z(I,Jp1,K);
TRANS_j(I,Jp1,K,4) = k_face_x(I,Jp1,K);
TRANS_j(I,Jp1,K,5) = k_face_y(I,Jp1,K);
TRANS_j(I,Jp1,K,6) = k_face_z(I,Jp1,K);
TRANS_j(I,Jp1,K,7) = i_face_x(I,Jp1,K);
TRANS_j(I,Jp1,K,8) = i_face_y(I,Jp1,K);
TRANS_j(I,Jp1,K,9) = i_face_z(I,Jp1,K);

end