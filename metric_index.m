function [I,J,K,...
          Ip1,Jp1,Kp1,...
          ic_act,jc_act,kc_act,...
          if_act,jf_act,kf_act,...
          ic_lb,jc_lb,kc_lb,...
          ic_rb,jc_rb,kc_rb,...
          if_lb,jf_lb,kf_lb,...
          if_rb,jf_rb,kf_rb] = metric_index(nx,ny,nz,nx_total,ny_total,nz_total,NO2)

% This code computes the indices for evolving the MHD equation and
% enforcing boundary conditions
% INPUT: nx,ny,nz                   - # of active cells
%        nx_total,ny_total,nz_total - # of total cells
%        NO2                        - # of ghost cells
      
% I,J,K are the indices for all cell centers
I = 1:nx_total-1;
J = 1:ny_total-1;
K = 1:nz_total-1;

% Ip1,Jp1,Kp1 are the indices for all cell corners
Ip1 = 1:nx_total;
Jp1 = 1:ny_total;
Kp1 = 1:nz_total;

% index of active cell centers
ic_act = NO2+1:NO2+nx;
jc_act = NO2+1:NO2+ny;
kc_act = NO2+1:NO2+nz;

% index of active face centers
if_act = NO2+1:NO2+nx+1;
jf_act = NO2+1:NO2+ny+1;
kf_act = NO2+1:NO2+nz+1;

% index of left-boudary for cell centers
ic_lb = 1:NO2;
jc_lb = 1:NO2;
kc_lb = 1:NO2;

% index of right-boundary for cell centers
ic_rb = nx+NO2+1: nx+NO2+NO2;
jc_rb = ny+NO2+1: ny+NO2+NO2;
kc_rb = nz+NO2+1: nz+NO2+NO2;

% index of left-boundary for face centers
if_lb = 1:NO2;
jf_lb = 1:NO2;
kf_lb = 1:NO2;

% index of right-boundary for face centers
if_rb = nx+NO2+1+1: nx+NO2+NO2+1;
jf_rb = ny+NO2+1+1: ny+NO2+NO2+1;
kf_rb = nz+NO2+1+1: nz+NO2+NO2+1;