function [Bstress_x_p, Bstress_y_p, Bstress_z_p, ...
          Bstress_x_n, Bstress_y_n, Bstress_z_n] = getMagneticStress(rho,Vx,p,bx,by,bz,bnorm,bx0,by0,bz0,b0n,CA,trans_1,trans_2,trans_3)

    % This code calculates the magnetic stresses on a cell face with the
    % normal direction x. The bx, by, bz are the purturbed B fields thats
    % tracked in the solver; the bx0, by0, bz0 are the background B fields
    % thats not changing with time. Note that for highly non-linear
    % background B fields, the bx0, by0, bz0 must be the high-order (usually
    % 8th) face integrated average fields rather than point values

    % note that the b0^2 terms are removed, this is actually curl b x (b+b0), 
    % assuming curl b0 x b0 is zero, which is calculated in Initialize() and 
    % added to the total stress in Hydro() after the magnetic stress calculation
    
    % INPUT:
    % rho, p: density, pressure
    % Vx: normal velocity 
    % bx,by,bz: magnetic fields
    % bnorm: interface magnetic flux - will be the same for both left and right stress calculations, which is very important
    % bx0,by0,bz0: background field components at interface
    % b0n: magnetic flux of background fields at interface
    % CA: Speed of light
    % trans_1, trans_2, trans_3: the x, y, z component of the face-normal vector
    
    % NOTE: 
    % this flux function computes the stress terms in the based (x,y,z)
    % coordinate system directly, so no transform is needed after the
    % calculation, which is different from the fluid stress calculations.
    
    va2 = ((bx+bx0).^2+(by+by0).^2+(bz+bz0).^2)./rho;
    valf = va2.*CA.^2./(va2+CA.^2); % effective Alfven speed (if CA>>1, than valf = va2)
    ptot = 2*p./rho + valf; % total pressure used in the Maxwellian distribution
    lamda = 1./(ptot); % width of the Maxwellian distribution
    bsq = (bx.^2+by.^2+bz.^2)+2.*(bx.*bx0+by.*by0+bz.*bz0);
    Vx0_p = 0.5*erfc(-sqrt(lamda).*Vx); % zeroth-order positive velocity moment
    Vx0_n = 0.5*erfc(+sqrt(lamda).*Vx); % zeroth-order negative velocity moment
    
    % stress in the positive x-direction (not the i,j,k normal direction!!)
    % since the magnetic stress is not a function of velocity, only
    % zeroth-moment integrals are needed in calculating the stresses
    Bstress_x_p = (0.5*trans_1.*bsq - (bx+bx0).*bnorm - bx.*b0n).*Vx0_p;
    Bstress_y_p = (0.5*trans_2.*bsq - (by+by0).*bnorm - by.*b0n).*Vx0_p;
    Bstress_z_p = (0.5*trans_3.*bsq - (bz+bz0).*bnorm - bz.*b0n).*Vx0_p;
    
    % stress in the negative x-direction
    Bstress_x_n = (0.5*trans_1.*bsq - (bx+bx0).*bnorm - bx.*b0n).*Vx0_n;
    Bstress_y_n = (0.5*trans_2.*bsq - (by+by0).*bnorm - by.*b0n).*Vx0_n;
    Bstress_z_n = (0.5*trans_3.*bsq - (bz+bz0).*bnorm - bz.*b0n).*Vx0_n;
    
end
