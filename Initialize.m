function [rho,vx,vy,vz,p,bi,bj,bk,...
          BX0,BY0,BZ0,...
          x_bound,y_bound,z_bound]=Initialize(x,y,z,xc,yc,zc,xi,xj,xk,i_face_area,j_face_area,k_face_area,prob)

% This subroutine initialize the plasma variables and magnetic fluxes, 
% together with the specification of boundary conditions
% the function form of background field is also specified in the subroutine
%
% INPUT: grids, face_areas, etc., prob
      
global I J K Ip1 Jp1 Kp1 
 
disp('Initializing Plasma Variables...');

% Define premitive plasma variables at cell centers
rho= zeros(size(xc));
vx = zeros(size(xc));
vy = zeros(size(xc));
vz = zeros(size(xc));
p = zeros(size(xc));
% Define magnetic fluxes at cell interfaces
bi = zeros(size(xi));
bj = zeros(size(xj));
bk = zeros(size(xk));

Vec_init = 1; % default method using vector potential for initialization
              % if Vec_init ~= 1, then the initialization should specify
              % the magnetic fluxes bi, bj, bk, directly (not recommended)

switch prob
    
    case 'OT2D'
        % Background field functions
        BX0 = @(x,y,z) 0;
        BY0 = @(x,y,z) 0;
        BZ0 = @(x,y,z) 0;
        
        % boundary condition options
        x_bound = 'periodic';
        y_bound = 'periodic';
        z_bound = 'extrap';
        
        % Vector potential for initial bi, bj, bk
        Ax = @(x,y,z) 0;
        Ay = @(x,y,z) 0;
        Az = @(x,y,z) (cos(4*pi*x)/pi/4 + cos(2*pi*y)/pi/2)/sqrt(pi*4);
        
        % plasma variables for 2-D Orzag-Tang vortex
        rho = rho+25/(pi*36);
        p   = p+5/(pi*12);
        vx  = -sin(2*pi*yc);
        vy  =  sin(2*pi*xc);
        vz(:) = 0;
        
    case 'BW2D'
        % Background field functions
        BX0 = @(x,y,z) 1/sqrt(2);
        BY0 = @(x,y,z) 1/sqrt(2);
        BZ0 = @(x,y,z) 0.0;
        
        % boundary condition options
        x_bound = 'periodic';
        y_bound = 'periodic';
        z_bound = 'extrap';
        
        % Vector potential for initial bi, bj, bk
        Ax = @(x,y,z) (x+y+z).*0;
        Ay = @(x,y,z) (x+y+z).*0;
        Az = @(x,y,z) (x+y+z).*0;
        
        % plasma variables for 2-D MHD Blast Wave
        rho(:)=1;
        p(:)=0.1;
        rc = sqrt((xc-0.5).^2+(yc-0.5).^2);
        p(rc<0.1) = 10.0;
        vx(:) = 0;
        vy(:) = 0;
        vz(:) = 0;
        
    case 'LOOP2D'
        % Background field functions
        BX0 = @(x,y,z) 0.0;
        BY0 = @(x,y,z) 0.0;
        BZ0 = @(x,y,z) 0.0;
        
        % boundary condition options
        x_bound = 'periodic';
        y_bound = 'periodic';
        z_bound = 'extrap';
        
        % Vector potential for initial bi, bj, bk
        Ax = @(x,y,z) 0.0;
        Ay = @(x,y,z) 0.0;
        Az = @(x,y,z) max(1e-3*(0.2-sqrt((x-0.5).^2+(y-0.5).^2)),0);

        % plasma variables for 2-D field-loop advection
        rho(:)=1;
        p(:)=1;
        vx(:) = sqrt(2)/2;
        vy(:) = sqrt(2)/2;
        vz(:) = 0.0;

    case 'KH2D'
        % Background field functions
        BX0 = @(x,y,z) 0.0;
        BY0 = @(x,y,z) 0.0;
        BZ0 = @(x,y,z) 0.0;
        
        % boundary condition options
        x_bound = 'periodic';
        y_bound = 'periodic';
        z_bound = 'extrap';
        
        % Vector potential for initial bi, bj, bk
        Ax = @(x,y,z) 0.0;
        Ay = @(x,y,z) 0.0;
        Az = @(x,y,z) 0.5*y;

        % plasma variables for 2-D Kelvin-Helmholtz instability
        rho(:) = 1;
        rho(abs(yc-0.5)<0.25)=2;
        p(:)=2.5;
        vx(:)=0.5;
        vx(abs(yc-0.5)<0.25)=-0.5;
        vy(:)=0;
        vx=vx+0.01*rand(length(I),length(J));
        vy=vy+0.01*rand(length(I),length(J));
        vz(:) = 0.0;
    
    case 'RT2D'
        % Background field functions
        BX0 = @(x,y,z) 0.0;
        BY0 = @(x,y,z) 0.0;
        BZ0 = @(x,y,z) 0.0;
        
        % boundary condition options
        x_bound = 'periodic';
        y_bound = 'reflect';
        z_bound = 'extrap';
        
        % Vector potential for initial bi, bj, bk
        Ax = @(x,y,z) 0.0;
        Ay = @(x,y,z) 0.0;
        Az = @(x,y,z) 0.0;

        % constants
        global G0x G0y G0z; %gravity acceleration
        P0 = 2.5;
        
        % plasma variables for 2-D Kelvin-Helmholtz instability
        rho = yc*0+1;
        rho(yc>0.0)=2;
        p=P0 + G0y.*rho.*yc;
        vx=yc*0;
        vy = 0.01.*(1+cos(4*pi.*yc)).*(1+cos(4*pi.*xc))/4;
        vz = yc*0;
        
    otherwise
        error('PROB not identified, exiting..');
end

% initialize bi, bj and bk using vector potential defined by function
% handles Ax, Ay and Az specified in the initialization
% integrate Ax along i edge

disp('Initializing Magnetic Fluxes...');

if(Vec_init==1)
    
    % Integrate A dot dl along i-edges
    for i=I(1):I(end)
        for j=Jp1(1):Jp1(end)
            for k=Kp1(1):Kp1(end)
                LAi(i,j,k) = (x(i+1,j,k)-x(i,j,k)).*GaussianLineIntegral(Ax,x(i+1,j,k),y(i+1,j,k),z(i+1,j,k),x(i,j,k),y(i,j,k),z(i,j,k)) + ...
                             (y(i+1,j,k)-y(i,j,k)).*GaussianLineIntegral(Ay,x(i+1,j,k),y(i+1,j,k),z(i+1,j,k),x(i,j,k),y(i,j,k),z(i,j,k)) + ...
                             (z(i+1,j,k)-z(i,j,k)).*GaussianLineIntegral(Az,x(i+1,j,k),y(i+1,j,k),z(i+1,j,k),x(i,j,k),y(i,j,k),z(i,j,k));
            end
        end
    end
    % integrate A dot dl along j-edges
    for i=Ip1(1):Ip1(end)
        for j=J(1):J(end)
            for k=Kp1(1):Kp1(end)
                LAj(i,j,k) = (y(i,j+1,k)-y(i,j,k)).*GaussianLineIntegral(Ay,x(i,j+1,k),y(i,j+1,k),z(i,j+1,k),x(i,j,k),y(i,j,k),z(i,j,k)) + ...
                             (x(i,j+1,k)-x(i,j,k)).*GaussianLineIntegral(Ax,x(i,j+1,k),y(i,j+1,k),z(i,j+1,k),x(i,j,k),y(i,j,k),z(i,j,k)) + ...
                             (z(i,j+1,k)-z(i,j,k)).*GaussianLineIntegral(Az,x(i,j+1,k),y(i,j+1,k),z(i,j+1,k),x(i,j,k),y(i,j,k),z(i,j,k));
            end
        end
    end
    % integrate A dot dl along k-edges
    for i=Ip1(1):Ip1(end)
        for j=Jp1(1):Jp1(end)
            for k=K(1):K(end)
                LAk(i,j,k) = (z(i,j,k+1)-z(i,j,k)).*GaussianLineIntegral(Az,x(i,j,k+1),y(i,j,k+1),z(i,j,k+1),x(i,j,k),y(i,j,k),z(i,j,k)) + ...
                             (x(i,j,k+1)-x(i,j,k)).*GaussianLineIntegral(Ax,x(i,j,k+1),y(i,j,k+1),z(i,j,k+1),x(i,j,k),y(i,j,k),z(i,j,k)) + ...
                             (y(i,j,k+1)-y(i,j,k)).*GaussianLineIntegral(Ay,x(i,j,k+1),y(i,j,k+1),z(i,j,k+1),x(i,j,k),y(i,j,k),z(i,j,k)) ;
            end
        end
    end
    
    % use Stokes theorm for bi, bj, bk - face-integrated fluxes -
    % note that in this code bi, bj, bk are magnetic fluxes throuth i, j, k faces, respectively
    % rather than magneitic field vecotrs at the face center. 
    i=Ip1;j=J;k=K;
    bi(i,j,k) = LAj(i,j,k) - LAj(i,j,k+1) + LAk(i,j+1,k) - LAk(i,j,k);
    
    i=I;j=Jp1;k=K;
    bj(i,j,k) = -( LAi(i,j,k) - LAi(i,j,k+1) + LAk(i+1,j,k) - LAk(i,j,k) );
    
    i=I;j=J;k=Kp1;
    bk(i,j,k) = LAi(i,j,k) - LAi(i,j+1,k) + LAj(i+1,j,k) - LAj(i,j,k);
end

end
