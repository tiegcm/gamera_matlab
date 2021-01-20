function [dpxB_G,dpyB_G,dpzB_G]=background_force(fsqi,fnxi,fnyi,fnzi,TRANS_i,i_face_area,...
                                                 fsqj,fnxj,fnyj,fnzj,TRANS_j,j_face_area,...
                                                 fsqk,fnxk,fnyk,fnzk,TRANS_k,k_face_area,volume)

% This subroutine calculates the Lorentz force from the background magnetic field, 
% For force-free background fields, the Lorentz force computed by this
% algorithm should be essentially zero O(dx^12) since 12-th order Gaussian
% quadrature is used to compute all the face integrals 
%
% INPUT: face-integrated terms from background_field()
%        TRANS_i, TRANS_j, TRANS_k 
%        {i,j,k}_fae_area
%        volume
                                             
global I J K Ip1 Jp1 Kp1
global nx_total ny_total nz_total

disp('Initialize Background Lorentz Force...');

% Initialize Gaussian integrals based Lorentz terms for B0
% Here dpxB_G is the x-component of the J0xB0 term integrated over volume(I,J,K)    
%      dpyB_G is the y-component of the J0xB0 term integrated over volume(I,J,K)    
%      dpzB_G is the z-component of the J0xB0 term integrated over volume(I,J,K)    
dpxB_G = zeros(nx_total-1,ny_total-1,nz_total-1);
dpyB_G = zeros(nx_total-1,ny_total-1,nz_total-1);
dpzB_G = zeros(nx_total-1,ny_total-1,nz_total-1);

I_stress_X(Ip1,J,K) = TRANS_i(Ip1,J,K,1).*0.5.*fsqi(Ip1,J,K) - fnxi(Ip1,J,K);
I_stress_Y(Ip1,J,K) = TRANS_i(Ip1,J,K,2).*0.5.*fsqi(Ip1,J,K) - fnyi(Ip1,J,K);
I_stress_Z(Ip1,J,K) = TRANS_i(Ip1,J,K,3).*0.5.*fsqi(Ip1,J,K) - fnzi(Ip1,J,K);

J_stress_X(I,Jp1,K) = TRANS_j(I,Jp1,K,1).*0.5.*fsqj(I,Jp1,K) - fnxj(I,Jp1,K);
J_stress_Y(I,Jp1,K) = TRANS_j(I,Jp1,K,2).*0.5.*fsqj(I,Jp1,K) - fnyj(I,Jp1,K);
J_stress_Z(I,Jp1,K) = TRANS_j(I,Jp1,K,3).*0.5.*fsqj(I,Jp1,K) - fnzj(I,Jp1,K);

K_stress_X(I,J,Kp1) = TRANS_k(I,J,Kp1,1).*0.5.*fsqk(I,J,Kp1) - fnxk(I,J,Kp1);
K_stress_Y(I,J,Kp1) = TRANS_k(I,J,Kp1,2).*0.5.*fsqk(I,J,Kp1) - fnyk(I,J,Kp1);
K_stress_Z(I,J,Kp1) = TRANS_k(I,J,Kp1,3).*0.5.*fsqk(I,J,Kp1) - fnzk(I,J,Kp1);

dpxB_G(I,J,K) = (I_stress_X(I+1,J,K).*i_face_area(I+1,J,K) - I_stress_X(I,J,K).*i_face_area(I,J,K) + ...
                 J_stress_X(I,J+1,K).*j_face_area(I,J+1,K) - J_stress_X(I,J,K).*j_face_area(I,J,K) + ...
                 K_stress_X(I,J,K+1).*k_face_area(I,J,K+1) - K_stress_X(I,J,K).*k_face_area(I,J,K) )./volume(I,J,K);
dpyB_G(I,J,K) = (I_stress_Y(I+1,J,K).*i_face_area(I+1,J,K) - I_stress_Y(I,J,K).*i_face_area(I,J,K) + ...
                 J_stress_Y(I,J+1,K).*j_face_area(I,J+1,K) - J_stress_Y(I,J,K).*j_face_area(I,J,K) + ...
                 K_stress_Y(I,J,K+1).*k_face_area(I,J,K+1) - K_stress_Y(I,J,K).*k_face_area(I,J,K) )./volume(I,J,K);
dpzB_G(I,J,K) = (I_stress_Z(I+1,J,K).*i_face_area(I+1,J,K) - I_stress_Z(I,J,K).*i_face_area(I,J,K) + ...
                 J_stress_Z(I,J+1,K).*j_face_area(I,J+1,K) - J_stress_Z(I,J,K).*j_face_area(I,J,K) + ...
                 K_stress_Z(I,J,K+1).*k_face_area(I,J,K+1) - K_stress_Z(I,J,K).*k_face_area(I,J,K) )./volume(I,J,K);
