function [rho_left, rho_right] = reconstruct_3DV(rho_h,if_act,jf_act,kf_act,PDMB, direction,type,volume,volume_face)

% This code computes the interface L/Rs through high-order reconstruction.
% The reconstruction applies to conserved variables, e.g., mass to
% incorporate the influence of geometry (especially important for
% non-orthogonal, rapidly changing geometries)
%
% MAIN ALGORITHMS:
%   Step 1. Interpolate conserved quantities along a stencil
%   Step 2. Splitting the interpolated interface state into L and R
%
% NOTE: 1. the DEFAULT limiter is the PDM limiter, including both centered and
%          upwind reconstruction. This subroutine should be easily adapted to
%          include other reconstruction method such as TVD, WENO, etc.
%       2. if no "volume" information is given as input, the subroutine
%          sets "volume = 1" and only reconstructs rho_h (density).

if (nargin<=6) % set the default limiter, other options include '8th' (no limiting) and 'PDMU' (upwind)
    type = 'PDM';
end

if (nargin<=7)             % if no "volume" is provided as input, then set volume=1 
    volume = rho_h*0+1.0;
    volume_face = rho_h*0+1.0;
    rho_h0 = rho_h;
    rho_h = rho_h.*volume; 
else                       % if "volume" is provided as input, then compute the conserved quantities (mass)
    rho_h0 = rho_h;        % this is the DENSITY used in the nonlinear switchers
    rho_h = rho_h.*volume; % this is the MASS used for high-order interpolation
end

if(direction==1) % reconstruction in the I-DIRECTION (first dimension)
    %8-th order reconstruction
    switch type
        
        case '8th' % 8-th order centered reconstruction w/o limiting
            
            % interpolate/reconstruct the conserved quantity (rho*volume)
            rho_left(if_act,jf_act,kf_act) = (-3*rho_h(if_act-4,jf_act,kf_act)+29*rho_h(if_act-3,jf_act,kf_act)-139*rho_h(if_act-2,jf_act,kf_act)+533*rho_h(if_act-1,jf_act,kf_act)+...
                                             533*rho_h(if_act,jf_act,kf_act)-139*rho_h(if_act+1,jf_act,kf_act)+29*rho_h(if_act+2,jf_act,kf_act)-3*rho_h(if_act+3,jf_act,kf_act))/840;
            % divide the interface volume to get interface rho
            rho_left(if_act,jf_act,kf_act) = rho_left(if_act,jf_act,kf_act)./volume_face(if_act,jf_act,kf_act); % since no limiting, left- and right-states are identical
            rho_right(if_act,jf_act,kf_act) = rho_left(if_act,jf_act,kf_act); % since no limiting, left- and right-states are identical

        case 'PDM' % 8-th order centered reconstruction w/ PDM limiting
            
            % interpolate/reconstruct the conserved quantity (rho*volume),
            % default 8th-order, centered
            rho_interp(if_act,jf_act,kf_act) = (-3*rho_h(if_act-4,jf_act,kf_act)+29*rho_h(if_act-3,jf_act,kf_act)-139*rho_h(if_act-2,jf_act,kf_act)+533*rho_h(if_act-1,jf_act,kf_act)+...
                                               533*rho_h(if_act,jf_act,kf_act)-139*rho_h(if_act+1,jf_act,kf_act)+29*rho_h(if_act+2,jf_act,kf_act)-3*rho_h(if_act+3,jf_act,kf_act))/840;
            % divide the interface volume to get interface rho
            rho_interp(if_act,jf_act,kf_act) = rho_interp(if_act,jf_act,kf_act)./volume_face(if_act,jf_act,kf_act);
            % limiting on rho (NOT rho*volume!) to get L/R use the PDM limiter
            [rho_left(if_act,jf_act,kf_act),rho_right(if_act,jf_act,kf_act)]= ...
                PDM2(rho_h0(if_act-2,jf_act,kf_act),rho_h0(if_act-1,jf_act,kf_act),rho_h0(if_act,jf_act,kf_act),rho_h0(if_act+1,jf_act,kf_act),rho_interp(if_act,jf_act,kf_act),PDMB);

        case 'PDMU' % 7-th order upwind reconstruction w/ PDM limiting & non-clipping (optional)
            
            % upwind interpolation for left-state (stencil shifted towards left)
            rho_interp(if_act,jf_act,kf_act) = -1/140*rho_h(if_act-4,jf_act,kf_act)+5/84*rho_h(if_act-3,jf_act,kf_act)-101/420*rho_h(if_act-2,jf_act,kf_act)+319/420*rho_h(if_act-1,jf_act,kf_act)+...
                                              107/210*rho_h(if_act,jf_act,kf_act)-19/210*rho_h(if_act+1,jf_act,kf_act)+1/105*rho_h(if_act+2,jf_act,kf_act);
            % divide the interface volume to get interface rho
            rho_interp(if_act,jf_act,kf_act) = rho_interp(if_act,jf_act,kf_act)./volume_face(if_act,jf_act,kf_act);
            % limiting on rho (NOT rho*volume!) to get L use the PDMU limiter
            rho_left(if_act,jf_act,kf_act) = PDMU7(rho_h0(if_act-3,jf_act,kf_act),rho_h0(if_act-2,jf_act,kf_act),rho_h0(if_act-1,jf_act,kf_act),rho_h0(if_act,jf_act,kf_act),rho_h0(if_act+1,jf_act,kf_act),rho_interp(if_act,jf_act,kf_act),PDMB);

            % upwind interpolation for right-state (stencil shifted towards right)
            rho_interp(if_act,jf_act,kf_act) = -1/140*rho_h(if_act+3,jf_act,kf_act)+5/84*rho_h(if_act+2,jf_act,kf_act)-101/420*rho_h(if_act+1,jf_act,kf_act)+319/420*rho_h(if_act,jf_act,kf_act)+...
                                              107/210*rho_h(if_act-1,jf_act,kf_act)-19/210*rho_h(if_act-2,jf_act,kf_act)+1/105*rho_h(if_act-3,jf_act,kf_act);
            % divide the interface volume to get interface rho
            rho_interp(if_act,jf_act,kf_act) = rho_interp(if_act,jf_act,kf_act)./volume_face(if_act,jf_act,kf_act);
            % limiting on rho (NOT rho*volume!) to get L use the PDMU limiter
            % NOTE that in order to use the same limiter subroutine for the
            % right state, the indices goes in to the PDMU7() function is
            % reversed.
            rho_right(if_act,jf_act,kf_act)= PDMU7(rho_h0(if_act+2,jf_act,kf_act),rho_h0(if_act+1,jf_act,kf_act),rho_h0(if_act,jf_act,kf_act),rho_h0(if_act-1,jf_act,kf_act),rho_h0(if_act-2,jf_act,kf_act),rho_interp(if_act,jf_act,kf_act),PDMB);
            
    end
    
elseif(direction==2)

    switch type % reconstruction in the J-DIRECTION (second dimension)
        
        case '8th' % 8-th order centered reconstruction w/o limiting
            
            rho_left(if_act,jf_act,kf_act) = (-3*rho_h(if_act,jf_act-4,kf_act)+29*rho_h(if_act,jf_act-3,kf_act)-139*rho_h(if_act,jf_act-2,kf_act)+533*rho_h(if_act,jf_act-1,kf_act)+...
                                             533*rho_h(if_act,jf_act,kf_act)-139*rho_h(if_act,jf_act+1,kf_act)+29*rho_h(if_act,jf_act+2,kf_act)-3*rho_h(if_act,jf_act+3,kf_act))/840;
            rho_left(if_act,jf_act,kf_act) = rho_left(if_act,jf_act,kf_act)./volume_face(if_act,jf_act,kf_act);
            rho_right(if_act,jf_act,kf_act) = rho_left(if_act,jf_act,kf_act);
        
        case 'PDM' % 8-th order centered reconstruction w/ PDM limiting
            
            % interpolation
            rho_interp(if_act,jf_act,kf_act) = (-3*rho_h(if_act,jf_act-4,kf_act)+29*rho_h(if_act,jf_act-3,kf_act)-139*rho_h(if_act,jf_act-2,kf_act)+533*rho_h(if_act,jf_act-1,kf_act)+...
                                               533*rho_h(if_act,jf_act,kf_act)-139*rho_h(if_act,jf_act+1,kf_act)+29*rho_h(if_act,jf_act+2,kf_act)-3*rho_h(if_act,jf_act+3,kf_act))/840;
            rho_interp(if_act,jf_act,kf_act) = rho_interp(if_act,jf_act,kf_act)./volume_face(if_act,jf_act,kf_act);
            % limiting
            [rho_left(if_act,jf_act,kf_act),rho_right(if_act,jf_act,kf_act)]= ...
                PDM2(rho_h0(if_act,jf_act-2,kf_act),rho_h0(if_act,jf_act-1,kf_act),rho_h0(if_act,jf_act,kf_act),...
                     rho_h0(if_act,jf_act+1,kf_act),rho_interp(if_act,jf_act,kf_act),PDMB);
       
        case 'PDMU' % 7-th order upwind reconstruction w/ PDM limiting & non-clipping (optional)
            
            % upwind interpolation and limiting for left-state
            rho_interp(if_act,jf_act,kf_act) = -1/140*rho_h(if_act,jf_act-4,kf_act)+5/84*rho_h(if_act,jf_act-3,kf_act)-101/420*rho_h(if_act,jf_act-2,kf_act)+319/420*rho_h(if_act,jf_act-1,kf_act)+...
                                              107/210*rho_h(if_act,jf_act,kf_act)-19/210*rho_h(if_act,jf_act+1,kf_act)+1/105*rho_h(if_act,jf_act+2,kf_act);
            rho_interp(if_act,jf_act,kf_act) = rho_interp(if_act,jf_act,kf_act)./volume_face(if_act,jf_act,kf_act);
            rho_left(if_act,jf_act,kf_act) = PDMU7(rho_h0(if_act,jf_act-3,kf_act),rho_h0(if_act,jf_act-2,kf_act),rho_h0(if_act,jf_act-1,kf_act),rho_h0(if_act,jf_act,kf_act),rho_h0(if_act,jf_act+1,kf_act),rho_interp(if_act,jf_act,kf_act),PDMB);

            % upwind interpolation and limiting for right-state
            rho_interp(if_act,jf_act,kf_act) = -1/140*rho_h(if_act,jf_act+3,kf_act)+5/84*rho_h(if_act,jf_act+2,kf_act)-101/420*rho_h(if_act,jf_act+1,kf_act)+319/420*rho_h(if_act,jf_act,kf_act)+...
                                              107/210*rho_h(if_act,jf_act-1,kf_act)-19/210*rho_h(if_act,jf_act-2,kf_act)+1/105*rho_h(if_act,jf_act-3,kf_act);
            rho_interp(if_act,jf_act,kf_act) = rho_interp(if_act,jf_act,kf_act)./volume_face(if_act,jf_act,kf_act);
            rho_right(if_act,jf_act,kf_act)= PDMU7(rho_h0(if_act,jf_act+2,kf_act),rho_h0(if_act,jf_act+1,kf_act),rho_h0(if_act,jf_act,kf_act),rho_h0(if_act,jf_act-1,kf_act),rho_h0(if_act,jf_act-2,kf_act),rho_interp(if_act,jf_act,kf_act),PDMB);
 
    end
    
elseif(direction==3) % reconstruction in the K-DIRECTION (second dimension)
    
    switch type
        
        case '8th' % 8-th order centered reconstruction w/o limiting
            
            rho_left(if_act,jf_act,kf_act) = (-3*rho_h(if_act,jf_act,kf_act-4)+29*rho_h(if_act,jf_act,kf_act-3)-139*rho_h(if_act,jf_act,kf_act-2)+533*rho_h(if_act,jf_act,kf_act-1)+...
                                             533*rho_h(if_act,jf_act,kf_act)-139*rho_h(if_act,jf_act,kf_act+1)+29*rho_h(if_act,jf_act,kf_act+2)-3*rho_h(if_act,jf_act,kf_act+3))/840;
            rho_left(if_act,jf_act,kf_act) = rho_left(if_act,jf_act,kf_act)./volume_face(if_act,jf_act,kf_act);
            rho_right(if_act,jf_act,kf_act) = rho_left(if_act,jf_act,kf_act);

        case 'PDM'  % 8-th order centered reconstruction w/ PDM limiting
            
            % interpolation
            rho_interp(if_act,jf_act,kf_act) = (-3*rho_h(if_act,jf_act,kf_act-4)+29*rho_h(if_act,jf_act,kf_act-3)-139*rho_h(if_act,jf_act,kf_act-2)+533*rho_h(if_act,jf_act,kf_act-1)+...
                                               533*rho_h(if_act,jf_act,kf_act)-139*rho_h(if_act,jf_act,kf_act+1)+29*rho_h(if_act,jf_act,kf_act+2)-3*rho_h(if_act,jf_act,kf_act+3))/840;
            rho_interp(if_act,jf_act,kf_act) = rho_interp(if_act,jf_act,kf_act)./volume_face(if_act,jf_act,kf_act);

            % limiting
            [rho_left(if_act,jf_act,kf_act),rho_right(if_act,jf_act,kf_act)]= ...
                PDM2(rho_h0(if_act,jf_act,kf_act-2),rho_h0(if_act,jf_act,kf_act-1),rho_h0(if_act,jf_act,kf_act),...
                     rho_h0(if_act,jf_act,kf_act+1),rho_interp(if_act,jf_act,kf_act),PDMB);
                 
        case 'PDMU' % 7-th order upwind reconstruction w/ PDM limiting & non-clipping (optional)
            
            % upwind interpolation and limiting for left-state
            rho_interp(if_act,jf_act,kf_act) = -1/140*rho_h(if_act,jf_act,kf_act-4)+5/84*rho_h(if_act,jf_act,kf_act-3)-101/420*rho_h(if_act,jf_act,kf_act-2)+319/420*rho_h(if_act,jf_act,kf_act-1)+...
                                              107/210*rho_h(if_act,jf_act,kf_act)-19/210*rho_h(if_act,jf_act,kf_act+1)+1/105*rho_h(if_act,jf_act,kf_act+2);
            rho_interp(if_act,jf_act,kf_act) = rho_interp(if_act,jf_act,kf_act)./volume_face(if_act,jf_act,kf_act);
            rho_left(if_act,jf_act,kf_act) = PDMU7(rho_h0(if_act,jf_act,kf_act-3),rho_h0(if_act,jf_act,kf_act-2),rho_h0(if_act,jf_act,kf_act-1),rho_h0(if_act,jf_act,kf_act),rho_h0(if_act,jf_act,kf_act+1),rho_interp(if_act,jf_act,kf_act),PDMB);

            % upwind interpolation and limiting for right-state
            rho_interp(if_act,jf_act,kf_act) = -1/140*rho_h(if_act,jf_act,kf_act+3)+5/84*rho_h(if_act,jf_act,kf_act+2)-101/420*rho_h(if_act,jf_act,kf_act+1)+319/420*rho_h(if_act,jf_act,kf_act)+...
                                              107/210*rho_h(if_act,jf_act,kf_act-1)-19/210*rho_h(if_act,jf_act,kf_act-2)+1/105*rho_h(if_act,jf_act,kf_act-3);
            rho_interp(if_act,jf_act,kf_act) = rho_interp(if_act,jf_act,kf_act)./volume_face(if_act,jf_act,kf_act);
            rho_right(if_act,jf_act,kf_act)= PDMU7(rho_h0(if_act,jf_act,kf_act+2),rho_h0(if_act,jf_act,kf_act+1),rho_h0(if_act,jf_act,kf_act),rho_h0(if_act,jf_act,kf_act-1),rho_h0(if_act,jf_act,kf_act-2),rho_interp(if_act,jf_act,kf_act),PDMB);
  
    end

end
    
% END of reconstruction
end