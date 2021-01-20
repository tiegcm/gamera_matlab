function vx_avg = center2corner(vx,direction)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the averaged values at cell corners, from cell centers.
% Note that the dimensions of arrays in this function varies depending on 
% the location of the values, which is a bit tricky.
%
% MAIN ALGORITHM (in the x-y plane):
%   STEP 1: interpolate cell-centered values in the x-direction to cell
%           faces - dimensions from (I,J,Kc) to (if_act,J,Kc)
%   STEP 2: interpolate the cell-faced values in the y-direction to cell
%           cornor- dimensions from (if_act,J,Kc) to (if_act,if_act,Kc)
% 
% The y-z and z-x plane calculations simply follow the same algorithm
%
% INPUT: velocity, plane of sweep

global ic_act jc_act kc_act if_act jf_act kf_act I J K

if (strcmp(direction,'ij')==1)
    % step 1 & 2: interpolate velocity from cell center to corner
    %         first interpolate in the x-direction, so the index for
    %         x-dimension changes from ic_act to if_act, for the y-direction,
    %         interpolation is done for the whole y-domain, i.e., J
    %         then interpolate in the y-direction, so the index for y-dimension
    %         changes from J to jf_act
    
    % interpolate vx from cell center to cell corner, vx_interp dimension is (if_act, jf_act, kc_act)
    vx_interp_x(if_act,J,kc_act) = (-3*vx(if_act-4,J,kc_act)+29*vx(if_act-3,J,kc_act)-139*vx(if_act-2,J,kc_act)+533*vx(if_act-1,J,kc_act)+...
        533*vx(if_act,J,kc_act)-139*vx(if_act+1,J,kc_act)+29*vx(if_act+2,J,kc_act)-3*vx(if_act+3,J,kc_act))/840;
    vx_avg(if_act,jf_act,kc_act) = (-3*vx_interp_x(if_act,jf_act-4,kc_act)+29*vx_interp_x(if_act,jf_act-3,kc_act)-139*vx_interp_x(if_act,jf_act-2,kc_act)+533*vx_interp_x(if_act,jf_act-1,kc_act)+...
        533*vx_interp_x(if_act,jf_act,kc_act)-139*vx_interp_x(if_act,jf_act+1,kc_act)+29*vx_interp_x(if_act,jf_act+2,kc_act)-3*vx_interp_x(if_act,jf_act+3,kc_act))/840;
    
elseif (strcmp(direction,'jk')==1) %
    % step 1 & 2: interpolate velocity from cell center to corner
    %         first interpolate in the y-direction, so the index for
    %         y-dimension changes from jc_act to jf_act, for the y-direction,
    %         interpolation is done for the whole z-domain, i.e., K
    %         then interpolate in the z-direction, so the index for z-dimension
    %         changes from K to kf_act
    
    % interpolate vx from cell center to cell corner, vx_interp dimension is (ic_act, jf_act, kf_act)
    vx_interp_x(ic_act,jf_act,K) = (-3*vx(ic_act,jf_act-4,K)+29*vx(ic_act,jf_act-3,K)-139*vx(ic_act,jf_act-2,K)+533*vx(ic_act,jf_act-1,K)+...
        533*vx(ic_act,jf_act,K)-139*vx(ic_act,jf_act+1,K)+29*vx(ic_act,jf_act+2,K)-3*vx(ic_act,jf_act+3,K))/840;
    vx_avg(ic_act,jf_act,kf_act) = (-3*vx_interp_x(ic_act,jf_act,kf_act-4)+29*vx_interp_x(ic_act,jf_act,kf_act-3)-139*vx_interp_x(ic_act,jf_act,kf_act-2)+533*vx_interp_x(ic_act,jf_act,kf_act-1)+...
        533*vx_interp_x(ic_act,jf_act,kf_act)-139*vx_interp_x(ic_act,jf_act,kf_act+1)+29*vx_interp_x(ic_act,jf_act,kf_act+2)-3*vx_interp_x(ic_act,jf_act,kf_act+3))/840;
    
elseif (strcmp(direction,'ki')==1) 
    % step 1 & 2: interpolate velocity from cell center to corner
    %         first interpolate in the z-direction, so the index for
    %         z-dimension changes from kc_act to kf_act, for the x-direction,
    %         interpolation is done for the whole x-domain, i.e., J
    %         then interpolate in the x-direction, so the index for x-dimension
    %         changes from J to jf_act
    
    % interpolate vx from cell center to cell corner, vx_interp dimension is (if_act, jc_act, kf_act)
    vx_interp_x(I,jc_act,kf_act) = (-3*vx(I,jc_act,kf_act-4)+29*vx(I,jc_act,kf_act-3)-139*vx(I,jc_act,kf_act-2)+533*vx(I,jc_act,kf_act-1)+...
        533*vx(I,jc_act,kf_act)-139*vx(I,jc_act,kf_act+1)+29*vx(I,jc_act,kf_act+2)-3*vx(I,jc_act,kf_act+3))/840;
    vx_avg(if_act,jc_act,kf_act) = (-3*vx_interp_x(if_act-4,jc_act,kf_act)+29*vx_interp_x(if_act-3,jc_act,kf_act)-139*vx_interp_x(if_act-2,jc_act,kf_act)+533*vx_interp_x(if_act-1,jc_act,kf_act)+...
        533*vx_interp_x(if_act,jc_act,kf_act)-139*vx_interp_x(if_act+1,jc_act,kf_act)+29*vx_interp_x(if_act+2,jc_act,kf_act)-3*vx_interp_x(if_act+3,jc_act,kf_act))/840;
end

end


