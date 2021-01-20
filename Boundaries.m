function [rho,p,vx,vy,vz,bi,bj,bk]= Boundaries(rho,p,vx,vy,vz,bi,bj,bk,NO,...
                                              x_bound,y_bound,z_bound,...
                                              ic_act,jc_act,kc_act,...
                                              if_act,jf_act,kf_act,...
                                              ic_lb,ic_rb,jc_lb,jc_rb,kc_lb,kc_rb,...
                                              if_lb,if_rb,jf_lb,jf_rb,kf_lb,kf_rb)
           
% This code specify simple boundary conditions, currently implemented
% options include: 
%                'periodic'  - periodic b.c.
%                'symmetric' - symmetric b.c. (normal gradient = 0)
%                'extrap'    - zeroth-order extrapolation b.c. 
%
% NOTE: the code only specify the primary variables include the plasma
%       variables (rho, vx,vy,vz,p) and magnetic fluxes (bi,bj,bk).
%       Magnetic fields (bx,by,bz) are derived from (bi,bj,bk) after B.Cs
%       are specified. In some simulation it could be a good idea to
%       specify b.c. for the magnetic fields (bx, by, bz) as well due to
%       the complexity in the boundary conditions.

NO2 = NO/2;
[nx,ny,nz]=size(rho);
nx = nx - NO;
ny = ny - NO;
nz = nz - NO;

switch x_bound
    case 'periodic'
        % --periodic boundary - x direction
        rho(ic_rb,:,:) = rho(ic_act(1):ic_act(NO2),:,:);
        p(ic_rb,:,:)   = p(ic_act(1):ic_act(NO2),:,:);
        vx(ic_rb,:,:)  = vx(ic_act(1):ic_act(NO2),:,:);
        vy(ic_rb,:,:)  = vy(ic_act(1):ic_act(NO2),:,:);
        vz(ic_rb,:,:)  = vz(ic_act(1):ic_act(NO2),:,:);
        bi(if_rb,:,:)  = bi(if_act(2):if_act(NO2+1),:,:); % <--- 'bi' has dimension (Ip1,J,K), so the I boundary use i-face active index
        bj(ic_rb,:,:)  = bj(ic_act(1):ic_act(NO2),:,:);   % <--- 'bj' has dimension (I,Jp1,K), so the I boundary use i-center active index
        bk(ic_rb,:,:)  = bk(ic_act(1):ic_act(NO2),:,:);   % <--- 'bk' has dimension (I,J,Kp1), so the I boundary use i-center active index
        % --periodic boundary - x direction
        rho(ic_lb,:,:) = rho(ic_act(end-NO2+1):ic_act(end),:,:);
        vx(ic_lb,:,:)  = vx(ic_act(end-NO2+1):ic_act(end),:,:);
        vy(ic_lb,:,:)  = vy(ic_act(end-NO2+1):ic_act(end),:,:);
        vz(ic_lb,:,:)  = vz(ic_act(end-NO2+1):ic_act(end),:,:);
        p(ic_lb,:,:)   =  p(ic_act(end-NO2+1):ic_act(end),:,:);
        bi(if_lb,:,:)  = bi(if_act(end-NO2):if_act(end-1),:,:);
        bj(ic_lb,:,:)  = bj(ic_act(end-NO2+1):ic_act(end),:,:);
        bk(ic_lb,:,:)  = bk(ic_act(end-NO2+1):ic_act(end),:,:);
        
    case 'symmetric'
        % -symmetric boundary - x direction
        rho(ic_lb,:,:) =rho(NO+1-ic_lb,:,:);
        vx(ic_lb,:,:)  = vx(NO+1-ic_lb,:,:);
        vy(ic_lb,:,:)  = vy(NO+1-ic_lb,:,:);
        vz(ic_lb,:,:)  = vz(NO+1-ic_lb,:,:);
        p(ic_lb,:,:)   =  p(NO+1-ic_lb,:,:);
        bi(if_lb,:,:)  = bi(NO+1-if_lb,:,:);
        bj(ic_lb,:,:)  = bj(NO+1-ic_lb,:,:);
        bk(ic_lb,:,:)  = bk(NO+1-ic_lb,:,:);
        % -symmetric boundary - x direction
        rho(ic_rb,:,:) = rho(2*nx+NO+1-ic_rb,:,:);
        p(ic_rb,:,:) = p(2*nx+NO+1-ic_rb,:,:);
        vx(ic_rb,:,:) = vx(2*nx+NO+1-ic_rb,:,:);
        vy(ic_rb,:,:) = vy(2*nx+NO+1-ic_rb,:,:);
        vz(ic_rb,:,:) = vz(2*nx+NO+1-ic_rb,:,:);
        bi(if_rb,:,:) = bi(2*nx+NO+1-if_rb+2,:,:); % <--- 'bi' has dimension (Ip1,J,K), so the I boundary use i-face active index
        bj(ic_rb,:,:) = bj(2*nx+NO+1-ic_rb,:,:);   % <--- 'bj' has dimension (I,Jp1,K), so the I boundary use i-center active index
        bk(ic_rb,:,:) = bk(2*nx+NO+1-ic_rb,:,:);   % <--- 'bk' has dimension (I,J,Kp1), so the I boundary use i-center active index

    case 'reflect'
        % -reflect boundary - x direction
        rho(ic_lb,:,:) =rho(NO+1-ic_lb,:,:);
        vx(ic_lb,:,:)  = -vx(NO+1-ic_lb,:,:);
        vy(ic_lb,:,:)  = -vy(NO+1-ic_lb,:,:);
        vz(ic_lb,:,:)  = -vz(NO+1-ic_lb,:,:);
        p(ic_lb,:,:)   =  p(NO+1-ic_lb,:,:);
        bi(if_lb,:,:)  = bi(NO+1-if_lb,:,:);
        bj(ic_lb,:,:)  = bj(NO+1-ic_lb,:,:);
        bk(ic_lb,:,:)  = bk(NO+1-ic_lb,:,:);
        % -symmetric boundary - x direction
        rho(ic_rb,:,:) = rho(2*nx+NO+1-ic_rb,:,:);
        p(ic_rb,:,:) = p(2*nx+NO+1-ic_rb,:,:);
        vx(ic_rb,:,:) = -vx(2*nx+NO+1-ic_rb,:,:);
        vy(ic_rb,:,:) = -vy(2*nx+NO+1-ic_rb,:,:);
        vz(ic_rb,:,:) = -vz(2*nx+NO+1-ic_rb,:,:);
        bi(if_rb,:,:) = bi(2*nx+NO+1-if_rb+2,:,:); % <--- 'bi' has dimension (Ip1,J,K), so the I boundary use i-face active index
        bj(ic_rb,:,:) = bj(2*nx+NO+1-ic_rb,:,:);   % <--- 'bj' has dimension (I,Jp1,K), so the I boundary use i-center active index
        bk(ic_rb,:,:) = bk(2*nx+NO+1-ic_rb,:,:);   % <--- 'bk' has dimension (I,J,Kp1), so the I boundary use i-center active index
        
        
    case 'extrap'
        
        % zeroth-order extrapolation in x
        ic_b1 = [ic_act(1) ic_act(1) ic_act(1) kc_act(1)];
        if_b1 = [if_act(1) if_act(1) if_act(1) if_act(1)];
        rho(ic_lb,:,:) = rho(ic_b1,:,:);
        vx(ic_lb,:,:)  = vx(ic_b1,:,:);
        vy(ic_lb,:,:)  = vy(ic_b1,:,:);
        vz(ic_lb,:,:)  = vz(ic_b1,:,:);
        p(ic_lb,:,:)   =  p(ic_b1,:,:);
        bi(if_lb,:,:)  = bi(if_b1,:,:);
        bj(ic_lb,:,:)  = bj(ic_b1,:,:);
        bk(ic_lb,:,:)  = bk(ic_b1,:,:);
        % --zeroth order extrapolation boundary - x direction
        ic_b2 = [ic_act(end) ic_act(end) ic_act(end) ic_act(end)];
        if_b2 = [if_act(end) if_act(end) if_act(end) if_act(end)];
        rho(ic_rb,:,:) = rho(ic_b2,:,:);
        vx(ic_rb,:,:)  = vx(ic_b2,:,:);
        vy(ic_rb,:,:)  = vy(ic_b2,:,:);
        vz(ic_rb,:,:)  = vz(ic_b2,:,:);
        p(ic_rb,:,:)   =  p(ic_b2,:,:);
        bi(if_rb,:,:)  = bi(if_b2,:,:);
        bj(ic_rb,:,:)  = bj(ic_b2,:,:);
        bk(ic_rb,:,:)  = bk(ic_b2,:,:);
        
end

switch y_bound
    case 'periodic'
        % % --periodic boundary - y direction;
        rho(:,jc_lb,:) =rho(:,jc_act(end-NO2+1):jc_act(end),:);
        vx(:,jc_lb,:)  = vx(:,jc_act(end-NO2+1):jc_act(end),:);
        vy(:,jc_lb,:)  = vy(:,jc_act(end-NO2+1):jc_act(end),:);
        vz(:,jc_lb,:)  = vz(:,jc_act(end-NO2+1):jc_act(end),:);
        p(:,jc_lb,:)   =  p(:,jc_act(end-NO2+1):jc_act(end),:);
        bi(:,jc_lb,:)  = bi(:,jc_act(end-NO2+1):jc_act(end),:);   % <--- 'bi' has dimension (Ip1,J,K), so the J boundary use j-center active index
        bj(:,jf_lb,:)  = bj(:,jf_act(end-NO2):jf_act(end-1),:);   % <--- 'bj' has dimension (I,Jp1,K), so the J boundary use j-face active index
        bk(:,jc_lb,:)  = bk(:,jc_act(end-NO2+1):jc_act(end),:);   % <--- 'bk' has dimension (I,J,Kp1), so the J boundary use j-center active index
        
        % --periodic boundary - y direction
        rho(:,jc_rb,:) = rho(:,jc_act(1):jc_act(NO2),:);
        vx(:,jc_rb,:)  = vx(:,jc_act(1):jc_act(NO2),:);
        vy(:,jc_rb,:)  = vy(:,jc_act(1):jc_act(NO2),:);
        vz(:,jc_rb,:)  = vz(:,jc_act(1):jc_act(NO2),:);
        p(:,jc_rb,:)   =  p(:,jc_act(1):jc_act(NO2),:);
        bi(:,jc_rb,:)  = bi(:,jc_act(1):jc_act(NO2),:);
        bj(:,jf_rb,:)  = bj(:,jf_act(2):jf_act(NO2+1),:);
        bk(:,jc_rb,:)  = bk(:,jc_act(1):jc_act(NO2),:);
        
    case 'symmetric'
        % --symmetric boundary - y direction;
        rho(:,jc_lb,:) =rho(:,NO+1-jc_lb,:);
        vx(:,jc_lb,:)  = vx(:,NO+1-jc_lb,:);
        vy(:,jc_lb,:)  = vy(:,NO+1-jc_lb,:);
        vz(:,jc_lb,:)  = -vz(:,NO+1-jc_lb,:);
        p(:,jc_lb,:)   =  p(:,NO+1-jc_lb,:);
        bi(:,jc_lb,:)  = bi(:,NO+1-jc_lb,:);   % <--- 'bi' has dimension (Ip1,J,K), so the J boundary use j-center active index
        bj(:,jf_lb,:)  = bj(:,NO+1-jf_lb,:);   % <--- 'bj' has dimension (I,Jp1,K), so the J boundary use j-face active index
        bk(:,jc_lb,:)  = bk(:,NO+1-jc_lb,:);   % <--- 'bk' has dimension (I,J,Kp1), so the J boundary use j-center active index
        
        % --symmetric boundary - y direction
        rho(:,jc_rb,:) =rho(:,2*ny+NO+1-jc_rb,:);
        vx(:,jc_rb,:)  = vx(:,2*ny+NO+1-jc_rb,:);
        vy(:,jc_rb,:)  = vy(:,2*ny+NO+1-jc_rb,:);
        vz(:,jc_rb,:)  = vz(:,2*ny+NO+1-jc_rb,:);
        p(:,jc_rb,:)   =  p(:,2*ny+NO+1-jc_rb,:);
        bi(:,jc_rb,:)  = bi(:,2*ny+NO+1-jc_rb,:);
        bj(:,jf_rb,:)  = bj(:,2*ny+NO+1-jf_rb+2,:);
        bk(:,jc_rb,:)  = bk(:,2*ny+NO+1-jc_rb,:);

    case 'reflect'
        % --symmetric boundary - y direction;
        rho(:,jc_lb,:) =rho(:,NO+1-jc_lb,:);
        vx(:,jc_lb,:)  = vx(:,NO+1-jc_lb,:);
        vy(:,jc_lb,:)  = -vy(:,NO+1-jc_lb,:);
        vz(:,jc_lb,:)  = vz(:,NO+1-jc_lb,:);
        p(:,jc_lb,:)   =  p(:,NO+1-jc_lb,:);
        bi(:,jc_lb,:)  = bi(:,NO+1-jc_lb,:);   % <--- 'bi' has dimension (Ip1,J,K), so the J boundary use j-center active index
        bj(:,jf_lb,:)  = bj(:,NO+1-jf_lb,:);   % <--- 'bj' has dimension (I,Jp1,K), so the J boundary use j-face active index
        bk(:,jc_lb,:)  = bk(:,NO+1-jc_lb,:);   % <--- 'bk' has dimension (I,J,Kp1), so the J boundary use j-center active index
        
        % --symmetric boundary - y direction
        rho(:,jc_rb,:) =rho(:,2*ny+NO+1-jc_rb,:);
        vx(:,jc_rb,:)  = vx(:,2*ny+NO+1-jc_rb,:);
        vy(:,jc_rb,:)  = -vy(:,2*ny+NO+1-jc_rb,:);
        vz(:,jc_rb,:)  = vz(:,2*ny+NO+1-jc_rb,:);
        p(:,jc_rb,:)   =  p(:,2*ny+NO+1-jc_rb,:);
        bi(:,jc_rb,:)  = bi(:,2*ny+NO+1-jc_rb,:);
        bj(:,jf_rb,:)  = bj(:,2*ny+NO+1-jf_rb+2,:);
        bk(:,jc_rb,:)  = bk(:,2*ny+NO+1-jc_rb,:);
        
        
    case 'extrap'
        jc_b1 = [jc_act(1) jc_act(1) jc_act(1) jc_act(1)];
        jf_b1 = [jf_act(1) jf_act(1) jf_act(1) jf_act(1)];
        
        rho(:,jc_lb,:) =rho(:,jc_b1,:);
        vx(:,jc_lb,:)  = vx(:,jc_b1,:);
        vy(:,jc_lb,:)  = vy(:,jc_b1,:);
        vz(:,jc_lb,:)  = vz(:,jc_b1,:);
        p(:,jc_lb,:)   =  p(:,jc_b1,:);
        bi(:,jc_lb,:)  = bi(:,jc_b1,:);
        bj(:,jf_lb,:)  = bj(:,jf_b1,:);
        bk(:,jc_lb,:)  = bk(:,jc_b1,:);
        
        % --zeroth order extrapolation boundary - y direction
        jc_b2 = [jc_act(end) jc_act(end) jc_act(end) jc_act(end)];
        jf_b2 = [jf_act(end) jf_act(end) jf_act(end) jf_act(end)];
        
        rho(:,jc_rb,:) =rho(:,jc_b2,:);
        vx(:,jc_rb,:)  = vx(:,jc_b2,:);
        vy(:,jc_rb,:)  = vy(:,jc_b2,:);
        vz(:,jc_rb,:)  = vz(:,jc_b2,:);
        p(:,jc_rb,:)   =  p(:,jc_b2,:);
        bi(:,jc_rb,:)  = bi(:,jc_b2,:);
        bj(:,jf_rb,:)  = bj(:,jf_b2,:);
        bk(:,jc_rb,:)  = bk(:,jc_b2,:);
end

switch z_bound
    case 'symmetric'
        % --routflow boundary - z direction;
        rho(:,:,kc_lb) =rho(:,:,NO+1-kc_lb);
        vx(:,:,kc_lb)  = vx(:,:,NO+1-kc_lb);
        vy(:,:,kc_lb)  = vy(:,:,NO+1-kc_lb);
        vz(:,:,kc_lb)  = vz(:,:,NO+1-kc_lb);
        p(:,:,kc_lb)   =  p(:,:,NO+1-kc_lb);
        bi(:,:,kc_lb)  = bi(:,:,NO+1-kc_lb);   % <--- 'bi' has dimension (Ip1,J,K), so the J boundary use j-center active index
        bj(:,:,kc_lb)  = bj(:,:,NO+1-kc_lb);   % <--- 'bj' has dimension (I,Jp1,K), so the J boundary use j-face active index
        bk(:,:,kf_lb)  = bk(:,:,NO+1-kf_lb);   % <--- 'bk' has dimension (I,J,Kp1), so the J boundary use j-center active index
        
        % --routflow boundary - z direction
        rho(:,:,kc_rb) =rho(:,:,2*nz+NO+1-kc_rb);
        vx(:,:,kc_rb)  = vx(:,:,2*nz+NO+1-kc_rb);
        vy(:,:,kc_rb)  = vy(:,:,2*nz+NO+1-kc_rb);
        vz(:,:,kc_rb)  = vz(:,:,2*nz+NO+1-kc_rb);
        p(:,:,kc_rb)   =  p(:,:,2*nz+NO+1-kc_rb);
        bi(:,:,kc_rb)  = bi(:,:,2*nz+NO+1-kc_rb);
        bj(:,:,kc_rb)  = bj(:,:,2*nz+NO+1-kc_rb);
        bk(:,:,kf_rb)  = bk(:,:,2*nz+NO+1-kf_rb+2);
        
    case 'periodic'
        % periodic boundary - z direction
        rho(:,:,kc_lb) =rho(:,:,kc_act(end-NO2+1):kc_act(end));
        vx(:,:,kc_lb)  = vx(:,:,kc_act(end-NO2+1):kc_act(end));
        vy(:,:,kc_lb)  = vy(:,:,kc_act(end-NO2+1):kc_act(end));
        vz(:,:,kc_lb)  = vz(:,:,kc_act(end-NO2+1):kc_act(end));
        p(:,:,kc_lb)   =  p(:,:,kc_act(end-NO2+1):kc_act(end));
        bi(:,:,kc_lb)  = bi(:,:,kc_act(end-NO2+1):kc_act(end));   % <--- 'bi' has dimension (Ip1,J,K), so the J boundary use j-center active index
        bj(:,:,kc_lb)  = bj(:,:,kc_act(end-NO2+1):kc_act(end));   % <--- 'bj' has dimension (I,Jp1,K), so the J boundary use j-center active index
        bk(:,:,kf_lb)  = bk(:,:,kf_act(end-NO2+1):kf_act(end));   % <--- 'bk' has dimension (I,J,Kp1), so the J boundary use j-face active index
        
        % periodic boundary - z direction
        rho(:,:,kc_rb) =rho(:,:,kc_act(1):kc_act(NO2));
        vx(:,:,kc_rb)  = vx(:,:,kc_act(1):kc_act(NO2));
        vy(:,:,kc_rb)  = vy(:,:,kc_act(1):kc_act(NO2));
        vz(:,:,kc_rb)  = vz(:,:,kc_act(1):kc_act(NO2));
        p(:,:,kc_rb)   =  p(:,:,kc_act(1):kc_act(NO2));
        bi(:,:,kc_rb)  = bi(:,:,kc_act(1):kc_act(NO2));
        bj(:,:,kc_rb)  = bj(:,:,kc_act(1):kc_act(NO2));
        bk(:,:,kf_rb)  = bk(:,:,kf_act(1):kf_act(NO2));
        
    case 'extrap'
        % % --zeroth order extraplation boundary - z direction
        kc_b1 = [kc_act(1) kc_act(1) kc_act(1) kc_act(1)];
        kf_b1 = [kf_act(1) kf_act(1) kf_act(1) kf_act(1)];
        
        rho(:,:,kc_lb) =rho(:,:,kc_b1);
        vx(:,:,kc_lb)  = vx(:,:,kc_b1);
        vy(:,:,kc_lb)  = vy(:,:,kc_b1);
        vz(:,:,kc_lb)  = vz(:,:,kc_b1);
        p(:,:,kc_lb)   =  p(:,:,kc_b1);
        bi(:,:,kc_lb)  = bi(:,:,kc_b1);
        bj(:,:,kc_lb)  = bj(:,:,kc_b1);
        bk(:,:,kf_lb)  = bk(:,:,kf_b1);
        
        kc_b2 = [kc_act(end) kc_act(end) kc_act(end) kc_act(end)];
        kf_b2 = [kf_act(end) kf_act(end) kf_act(end) kf_act(end)];
        % --zeroth order extrapolation boundary - z direction
        rho(:,:,kc_rb) =rho(:,:,kc_b2);
        vx(:,:,kc_rb)  = vx(:,:,kc_b2);
        vy(:,:,kc_rb)  = vy(:,:,kc_b2);
        vz(:,:,kc_rb)  = vz(:,:,kc_b2);
        p(:,:,kc_rb)   =  p(:,:,kc_b2);
        bi(:,:,kc_rb)  = bi(:,:,kc_b2);   % <--- 'bi' has dimension (Ip1,J,K), so the K boundary use k-center active index
        bj(:,:,kc_rb)  = bj(:,:,kc_b2);   % <--- 'bj' has dimension (I,Jp1,K), so the K boundary use k-center active index
        bk(:,:,kf_rb)  = bk(:,:,kf_b2);   % <--- 'bk' has dimension (I,J,Kp1), so the K boundary use k-face active index
        
end

end