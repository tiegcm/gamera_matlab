%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                The Matlab version of GAMERA MHD Solver
%            Authors B. Zhang, J. Lyon, K. Sorathia, Aug. 2018
% 
%  MAIN FEATURES:
%     physics:
%        1. Single-fluid, ideal MHD equations
%        2. Plasma energy equation (semi-conservative)
%        3. Semi-relativistic (Boris/Alfven) correction
%     numerics:
%        1. Finite-volume method
%        2. Non-orthogonal Curvilinear grids
%        3. Constraint transport for magnetic flux
%        4. second-order predictor/corrector in time
%        5. seventh-/eighth-order spatial reconstruction
%        6. PDM limiter with non-clipping option
%        7. Gas-kinetic flux functions
%        8. Background magnetic field splitting
%
%  REMARKS:
%        1. Go GAMERA! 
%        2. Working is better than consistancy
%        3. It is what it is (nothing to hide)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

% key simulation parameters
global NO NO2 CA PDMB 

NO=8;           % order of reconstruction - not actually used in the reconstruction module
NO2 = NO/2;     % # of ghost cells on each side of a stencil, default 4
CFL = 0.3;      % the Courant number, related to the PDMB value, see Kain [1987] JCP for detail
gamma= 1.4;     % ratio of specific heats
PDMB = 4;       % PDMB value for the PDM limiter, the smaller PDMB the more diffusion, default 4
CA = 1e10;      % upper limit of the speed of light, in normalized units
limiter = 'PDM';% the choice of limiting algorithm: 'PDM'  = 8th-order centered scheme
                %                                   'PDMU' = 7th-order upwind scheme

% Default problems: 1. OT2D   - Orszag-Tang vortex
%                   2. BW2D   - spherical blast wave
%                   3. LOOP2D - field-loop advection
%                   4. KH2D   - Kelvin-Helmholtz
% prob specifications are in the Initialize() function
% all the problem speifications are from the Athana MHD code test suite by
% Stone et al. [2008] (https://www.astro.princeton.edu/~jstone/Athena/tests/)
prob = 'RT2D'; 
T_stop = 10.0; % maximum sim time

global nx ny nz nx_total ny_total nz_total
global G0x G0y G0z; % gravitational acceleration
G0x = 0.0;
G0y = -0.1;
G0z = 0.0;

% generate a grid with (nx,ny,nz) points. For 2-D simulations nz = 1
nx = 50; % # of active cell centers in the first dimension (i-direction)
ny = 150; % # of active cell centers in the second dimension (j-direction)
nz = 1;   % # of active cell centers in the third dimension (k-direction)
v0 = 0.0; % level or distortion, v0=0 (no distortion); v0=0.2 (extreme distortion)

% generate a distorted Cartesian grid as the default grid
[x,y,z]=Generate_Distorted_Cartesian(nx,ny,nz,v0,NO2);
[nx_total,ny_total,nz_total]=size(x);

% scale the grid to [-0.25, 0.25]x[-0.75, 0.75]
x = (x-0.5)*2*0.25;
y = (y-0.5)*2*0.75;

% calculate the indices for active and the whole domain, all the indices
% are global variables accessible from any subroutines
global I J K Ip1 Jp1 Kp1 ic_act jc_act kc_act if_act jf_act kf_act
global ic_lb jc_lb kc_lb ic_rb jc_rb kc_rb if_lb jf_lb kf_lb if_rb jf_rb kf_rb

% metric calculation 1: grid index
[I,J,K,...                % I,J,K are the indices for cell centers
 Ip1,Jp1,Kp1,...          % Ip1,Jp1,Kp1 are the indices for cell corners
 ic_act,jc_act,kc_act,... % index of active cell centers
 if_act,jf_act,kf_act,... % index of active face centers
 ic_lb,jc_lb,kc_lb,...    % index of left for cell centers and faces
 ic_rb,jc_rb,kc_rb,...
 if_lb,jf_lb,kf_lb,...    % index of right for cell centers and faces
 if_rb,jf_rb,kf_rb] = metric_index(nx,ny,nz,nx_total,ny_total,nz_total,NO2);

% metric calculation 2: grid locations and lengths
[xc,yc,zc,...             % cell centers
 xi,yi,zi,...             % i-face centers 
 xj,yj,zj,...             % j-face centers 
 xk,yk,zk,...             % k-face centers
 dx,dy,dz,...             % cell lengths
 di,dj,dk,... 
 i_face_area,j_face_area,k_face_area,volume]=metric_grid(x,y,z);
% 
% metric calculation 3: face-normal transforms
% TRANS_{i,j,k} are the transform matrix at cell face, used in the flux
% calculations transforming the (x,y,z) component into fae normal
% coordinates at i,j,k interfaces, respectively.
% (i_face_norm_x,i_face_norm_y,i_face_norm_x): i-face normal vector
% (j_face_norm_x,j_face_norm_y,j_face_norm_x): j-face normal vector
% (k_face_norm_x,k_face_norm_y,k_face_norm_x): k-face normal vector
[TRANS_i,TRANS_j,TRANS_k,...
 i_face_norm_x,i_face_norm_y,i_face_norm_z,...
 j_face_norm_x,j_face_norm_y,j_face_norm_z,...
 k_face_norm_x,k_face_norm_y,k_face_norm_z]=metric_face(x,y,z);

% metric calculation 4: edge-aligned transforms
% vec_{i,j,k} are transform matrix for edge-aligned coordinate transforms
% (ifaceA_Kedge,jfaceA_Kedge): i- and j-face at k-edges for magnetic field estimation at cell edges
% xnqi_K,ynqi_K,xnqj_K,ynqj_K,DETK: transform elements for B^{e1} and B^{e2} at k-edges
% xnqi_J,ynqi_J,xnqj_J,ynqj_J,DETJ: transform elements for B^{e1} and B^{e2} at j-edges
% xnqi_I,ynqi_I,xnqj_I,ynqj_I,DETI: transform elements for B^{e1} and B^{e2} at i-edges
% {i,j,k}_edeg: corresponding lengths of cell edges
[ifaceA_Kedge,jfaceA_Kedge,xnqi_K,ynqi_K,xnqj_K,ynqj_K,DETK,k_edge,vec_k,...
 jfaceA_Iedge,kfaceA_Iedge,xnqj_I,ynqj_I,xnqk_I,ynqk_I,DETI,i_edge,vec_i,...
 kfaceA_Jedge,ifaceA_Jedge,xnqk_J,ynqk_J,xnqi_J,ynqi_J,DETJ,j_edge,vec_j]...
         = metric_edge(x,y,z,i_face_norm_x,i_face_norm_y,i_face_norm_z,i_face_area,...
                             j_face_norm_x,j_face_norm_y,j_face_norm_z,j_face_area,...
                             k_face_norm_x,k_face_norm_y,k_face_norm_z,k_face_area);

% intialization, which specifies:
% rho,vx,vy,vz,p - cell cenered plasma variables
% bi,bj,bk - face magnetic fluxes using vector potential (Ax,Ay,Az)
% BX0,BY0,BZ0 - background magnetic field functions
% {x,y,z}_bound - boundary condition choices: 'periodic'  - periodic
%                                             'symmetric' - zero derivative
%                                             'extrap'    - zeroth-order extrapolation
[rho,vx,vy,vz,p,bi,bj,bk,...
 BX0,BY0,BZ0,...
 x_bound,y_bound,z_bound]=Initialize(x,y,z,xc,yc,zc,xi,xj,xk,i_face_area,j_face_area,k_face_area,prob);

% impose boundary conditions for plasma variables and magnetic fluxes
% x_bound, y_bound and z_bound are specified in the 'Initialize' subroutine
[rho,p,vx,vy,vz,bi,bj,bk]= Boundaries(rho,p,vx,vy,vz,bi,bj,bk,NO,...
                                      x_bound,y_bound,z_bound,...
                                      ic_act,jc_act,kc_act,...
                                      if_act,jf_act,kf_act,...
                                      ic_lb,ic_rb,jc_lb,jc_rb,kc_lb,kc_rb,...
                                      if_lb,if_rb,jf_lb,jf_rb,kf_lb,kf_rb);

% compute bx, by, bz from bi, bj, bk (Eqn. 28)
bx(I,J,K) = ( bi(I+1,J,K).*xi(I+1,J,K)-bi(I,J,K).*xi(I,J,K) + bj(I,J+1,K).*xj(I,J+1,K)-bj(I,J,K).*xj(I,J,K) + bk(I,J,K+1).*xk(I,J,K+1)-bk(I,J,K).*xk(I,J,K) )./volume(I,J,K);
by(I,J,K) = ( bi(I+1,J,K).*yi(I+1,J,K)-bi(I,J,K).*yi(I,J,K) + bj(I,J+1,K).*yj(I,J+1,K)-bj(I,J,K).*yj(I,J,K) + bk(I,J,K+1).*yk(I,J,K+1)-bk(I,J,K).*yk(I,J,K) )./volume(I,J,K);
bz(I,J,K) = ( bi(I+1,J,K).*zi(I+1,J,K)-bi(I,J,K).*zi(I,J,K) + bj(I,J+1,K).*zj(I,J+1,K)-bj(I,J,K).*zj(I,J,K) + bk(I,J,K+1).*zk(I,J,K+1)-bk(I,J,K).*zk(I,J,K) )./volume(I,J,K);

% compute total magnetic field. BX0, BY0, BZ0 are specifield in the 'Initialize' subroutine
% bx_total, by_total and bz_total are not directlly used in the MHD solver, they
% are only used in the time step calculation (Alfven speed)
bx_total(I,J,K) = bx(I,J,K) + BX0(xc,yc,zc);
by_total(I,J,K) = by(I,J,K) + BY0(xc,yc,zc);
bz_total(I,J,K) = bz(I,J,K) + BZ0(xc,yc,zc);

% compute Background fields using 12-th orderGaussian quadrature
% use the i-interfaces as an example:
% f{x,y,z}i : the {x-,y-,z-}component of B0 integrated on i-interfaces 
% fnormi    : the (B0 dot n_i) normal component of B0 integrated on i-interfaces 
% fsqi      : the B0^2 component integrated on i-interfaces 
% fn{x,y,z}i: the (B0 dot n_i)*B0{x,y,z} integrated on i-interfaces 
% ALONG i-edges, BX0_I and BY0_I are the two perpendicular components in
% the edge-aligned system for vxB calculations
[fxi,fyi,fzi,fnormi,fsqi,fnxi,fnyi,fnzi,...
 fxj,fyj,fzj,fnormj,fsqj,fnxj,fnyj,fnzj,...
 fxk,fyk,fzk,fnormk,fsqk,fnxk,fnyk,fnzk,...
 BX0_I,BY0_I,BX0_J,BY0_J,BX0_K,BY0_K]=background_field(x,y,z,BX0,BY0,BZ0,vec_i,vec_j,vec_k);
      
% compute the Lorentz force from the background magnetic field
[dpxB_G,dpyB_G,dpzB_G]=background_force(fsqi,fnxi,fnyi,fnzi,TRANS_i,i_face_area,...
                                        fsqj,fnxj,fnyj,fnzj,TRANS_j,j_face_area,...
                                        fsqk,fnxk,fnyk,fnzk,TRANS_k,k_face_area,volume);

disp('O.K. to here, Bill!'); % ask John...

% Compute conserved hydrodynamic variables
[rhovx,rhovy,rhovz,eng] = getConservedVariables(rho,vx,vy,vz,p,gamma);

% compute the time step
dt = getDT(rho,vx,vy,vz,bx_total,by_total,bz_total,p,gamma,di,dj,dk,CFL);

% Adam-Bashforth predictor step
% Save the initial states for the first Adam-Bashforth time stepping
% plasma variables
rho_p = rho;
vx_p = vx;
vy_p = vy;
vz_p = vz;
p_p= p;
% magnetic variables, here we do predictor for both flux and fields
bx_p = bx;
by_p = by;
bz_p = bz;
bi_p = bi;
bj_p = bj;
bk_p = bk;

dt0 = dt;

Time = 0; % simulation time
RealT = 0;% real time used
N_max = 1000000; % max simulation steps

% Main Loop
for n=1:N_max
    
    Tstart=tic;
    
    % compute the time step
    dt = getDT(rho,vx,vy,vz,bx_total,by_total,bz_total,p,gamma,di,dj,dk,CFL);
    
    Time=Time+dt; % advance the simulation time
    
    % compute conserved variables from primitive variables
    [rhovx,rhovy,rhovz,eng] = getConservedVariables(rho,vx,vy,vz,p,gamma);
    
    % Step 1: get the half time step values at t = N+1/2 - linear extrapolation in time
    % plasma variables
    rho_h = rho + dt./dt0./2.*(rho-rho_p);
    vx_h = vx + dt./dt0./2.*(vx-vx_p);
    vy_h = vy + dt./dt0./2.*(vy-vy_p);
    vz_h = vz + dt./dt0./2.*(vz-vz_p);
    p_h = p + dt./dt0./2.*(p-p_p);  
    % magnetic fluxes and fields
    bx_h = bx + dt./dt0./2.*(bx-bx_p); 
    by_h = by + dt./dt0./2.*(by-by_p);
    bz_h = bz + dt./dt0./2.*(bz-bz_p);
    bi_h = bi + dt./dt0./2.*(bi-bi_p); 
    bj_h = bj + dt./dt0./2.*(bj-bj_p);
    bk_h = bk + dt./dt0./2.*(bk-bk_p);    

    % save the base step at t=N for updates
    rho0 = rho;
    rhovx0 = rhovx;
    rhovy0 = rhovy;
    rhovz0 = rhovz;
    eng0 = eng; 
    bx0 = bx;
    by0 = by;
    bz0 = bz;    
    vx0 = vx;
    vy0 = vy;
    vz0 = vz;

    % save the state variables (include dt) at t=N for the next A-B predicotr        
    rho_p = rho;
    vx_p = vx;
    vy_p = vy;
    vz_p = vz;   
    p_p= p;
    bx_p = bx;
    by_p = by;
    bz_p = bz;
    bi_p = bi;
    bj_p = bj;
    bk_p = bk;    
    dt0 = dt;   
    
    % calculat the electric fields using predictor values
    [Ei,Ej,Ek]=getEfields(rho_h,vx_h,vy_h,vz_h,bi_h,bj_h,bk_h,BX0_I,BY0_I,BX0_J,BY0_J,BX0_K,BY0_K,...
                         ifaceA_Kedge,jfaceA_Kedge,xnqi_K,ynqi_K,xnqj_K,ynqj_K,DETK,k_edge,vec_k,...
                         jfaceA_Iedge,kfaceA_Iedge,xnqj_I,ynqj_I,xnqk_I,ynqk_I,DETI,i_edge,vec_i,...
                         kfaceA_Jedge,ifaceA_Jedge,xnqk_J,ynqk_J,xnqi_J,ynqi_J,DETJ,j_edge,vec_j,PDMB,limiter);
    
    % compute effective Alfven speeds (if CA>>1, va_eff == Valfvn)
    Valfvn = sqrt(bx_total.^2+by_total.^2+bz_total.^2)./sqrt(rho); % Alfven speed VA
    Va_eff = Valfvn.*CA./sqrt(Valfvn.^2+CA.^2); % Boris-corrected (limited) VA
                  
    % calculate the mass flux (drho), fluid stress (dpxF,dpyF,dpzF), energy flux (deng) 
    % and magnetic stress (dpxB,dpyB,dpzB) using predictor values
    [drho,dpxF,dpyF,dpzF,deng,dpxB,dpyB,dpzB]=Hydro(rho_h,vx_h,vy_h,vz_h,p_h,bi_h,bj_h,bk_h,bx_h,by_h,bz_h,TRANS_i,TRANS_j,TRANS_k,...
                                                    fxi,fyi,fzi,fnormi,fxj,fyj,fzj,fnormj,fxk,fyk,fzk,fnormk,...
                                                    i_face_area,j_face_area,k_face_area,...
                                                    CA,Va_eff,volume,PDMB,limiter,gamma,n,dt);

    % add the lorentz force terms from the background field B0 - not the
    % deltas are -dt*volume_integral./volume, only need active cells
    dpxB(ic_act,jc_act,kc_act) = dpxB(ic_act,jc_act,kc_act) -dt.*dpxB_G(ic_act,jc_act,kc_act)*0;
    dpyB(ic_act,jc_act,kc_act) = dpyB(ic_act,jc_act,kc_act) -dt.*dpyB_G(ic_act,jc_act,kc_act)*0;    
    dpzB(ic_act,jc_act,kc_act) = dpzB(ic_act,jc_act,kc_act) -dt.*dpzB_G(ic_act,jc_act,kc_act)*0;
    
    % update of face-center magnetic fluxes (bi,bj,bk) through Faraday's law
    bi(if_act,jc_act,kc_act) = bi(if_act,jc_act,kc_act) - dt.*( (Ek(if_act,jc_act+1,kc_act)-Ek(if_act,jc_act,kc_act)) ...
                                                               -(Ej(if_act,jc_act,kc_act+1)-Ej(if_act,jc_act,kc_act)) );
    bj(ic_act,jf_act,kc_act) = bj(ic_act,jf_act,kc_act) - dt.*( (Ei(ic_act,jf_act,kc_act+1)-Ei(ic_act,jf_act,kc_act)) ...
                                                               -(Ek(ic_act+1,jf_act,kc_act)-Ek(ic_act,jf_act,kc_act)) );
    bk(ic_act,jc_act,kf_act) = bk(ic_act,jc_act,kf_act) - dt.*( (Ej(ic_act+1,jc_act,kf_act)-Ej(ic_act,jc_act,kf_act)) ...
                                                               -(Ei(ic_act,jc_act+1,kf_act)-Ei(ic_act,jc_act,kf_act)) );    
                                                            
    % compute bx, by, bz from bi, bj, bk (Eqn. 28) - now B fields are solved
    bx(ic_act,jc_act,kc_act) = ( bi(ic_act+1,jc_act,kc_act).*xi(ic_act+1,jc_act,kc_act)-bi(ic_act,jc_act,kc_act).*xi(ic_act,jc_act,kc_act) + bj(ic_act,jc_act+1,kc_act).*xj(ic_act,jc_act+1,kc_act)-bj(ic_act,jc_act,kc_act).*xj(ic_act,jc_act,kc_act) + bk(ic_act,jc_act,kc_act+1).*xk(ic_act,jc_act,kc_act+1)-bk(ic_act,jc_act,kc_act).*xk(ic_act,jc_act,kc_act) )./volume(ic_act,jc_act,kc_act);
    by(ic_act,jc_act,kc_act) = ( bi(ic_act+1,jc_act,kc_act).*yi(ic_act+1,jc_act,kc_act)-bi(ic_act,jc_act,kc_act).*yi(ic_act,jc_act,kc_act) + bj(ic_act,jc_act+1,kc_act).*yj(ic_act,jc_act+1,kc_act)-bj(ic_act,jc_act,kc_act).*yj(ic_act,jc_act,kc_act) + bk(ic_act,jc_act,kc_act+1).*yk(ic_act,jc_act,kc_act+1)-bk(ic_act,jc_act,kc_act).*yk(ic_act,jc_act,kc_act) )./volume(ic_act,jc_act,kc_act);
    bz(ic_act,jc_act,kc_act) = ( bi(ic_act+1,jc_act,kc_act).*zi(ic_act+1,jc_act,kc_act)-bi(ic_act,jc_act,kc_act).*zi(ic_act,jc_act,kc_act) + bj(ic_act,jc_act+1,kc_act).*zj(ic_act,jc_act+1,kc_act)-bj(ic_act,jc_act,kc_act).*zj(ic_act,jc_act,kc_act) + bk(ic_act,jc_act,kc_act+1).*zk(ic_act,jc_act,kc_act+1)-bk(ic_act,jc_act,kc_act).*zk(ic_act,jc_act,kc_act) )./volume(ic_act,jc_act,kc_act);
    
    %initialize intermediate variables for implementing the Boris correction
    if(n==1)
        bx1 = bx.*0;
        by1 = by.*0;
        bz1 = bz.*0;
        bx2 = bx.*0;
        by2 = by.*0;
        bz2 = bz.*0;
    end
    
    % get mid-time step total B field, using volume-averaged B at cell
    % this is total B(n+1/2), f{x,y,z}k are the Gaussian integrated fluxes
    % at k-interfaces, here we use a simple face average for cell-center field estimations 
    bx1(ic_act,jc_act,kc_act) = 0.5*(bx(ic_act,jc_act,kc_act)+bx0(ic_act,jc_act,kc_act)) + 0.5.*(fxk(ic_act,jc_act,kc_act)+fxk(ic_act,jc_act,kc_act+1));
    by1(ic_act,jc_act,kc_act) = 0.5*(by(ic_act,jc_act,kc_act)+by0(ic_act,jc_act,kc_act)) + 0.5.*(fyk(ic_act,jc_act,kc_act)+fyk(ic_act,jc_act,kc_act+1));
    bz1(ic_act,jc_act,kc_act) = 0.5*(bz(ic_act,jc_act,kc_act)+bz0(ic_act,jc_act,kc_act)) + 0.5.*(fzk(ic_act,jc_act,kc_act)+fzk(ic_act,jc_act,kc_act+1));
    % this is total B(n+1)
    bx2(ic_act,jc_act,kc_act) = bx(ic_act,jc_act,kc_act) + 0.5.*(fxk(ic_act,jc_act,kc_act)+fxk(ic_act,jc_act,kc_act+1));
    by2(ic_act,jc_act,kc_act) = by(ic_act,jc_act,kc_act) + 0.5.*(fyk(ic_act,jc_act,kc_act)+fyk(ic_act,jc_act,kc_act+1));
    bz2(ic_act,jc_act,kc_act) = bz(ic_act,jc_act,kc_act) + 0.5.*(fzk(ic_act,jc_act,kc_act)+fzk(ic_act,jc_act,kc_act+1));
    
    % update the bulk hydrodynamic equation (without magnetic stress) to get density and pressure                                
    rho   = rho0   + drho; 
    rhovx = rhovx0 + dpxF;
    rhovy = rhovy0 + dpyF;
    rhovz = rhovz0 + dpzF;
    eng   = eng0   + deng;
    
    % solve for the thermal pressure 
    vx = rhovx./rho;
    vy = rhovy./rho;
    vz = rhovz./rho;
    p = (eng - 0.5.*rho.*(vx.^2+vy.^2+vz.^2)).*(gamma-1); 
        
    % now density and pressure are solved
    % next apply the magnetic stress WITH Boris correction to get velocity
    % for ideal MHD simulations WITHOUT Boris correction, CA >> 1, here we
    % use the default CA=1e10, which means the Boris correction is OFF
    % NOTE: Boric correction is off, the operator splitting step is just
    %              rhovx = rhovx + dpxB
    %              rhovy = rhovy + dpxB
    %              rhovz = rhovz + dpxB
    
    % first estimate the Alfven constant for perp momentum
    b1 = sqrt(bx1.^2+by1.^2+bz1.^2);
    b2 = sqrt(bx2.^2+by2.^2+bz2.^2);
    alfn_ratio = b1.^2./rho./CA^2;  % use the updated rho at T(n+1) for the alfven correction 
    perp_ratio = 1./(1+alfn_ratio);
    bdx = dpxB./rho.*perp_ratio;
    bdy = dpyB./rho.*perp_ratio;
    bdz = dpzB./rho.*perp_ratio;
    % then update the momentum at t(N+1) with Alfven correction 
    dv_alf = alfn_ratio.*drho;
    vx_tmp = alfn_ratio.*(dpxF - drho.*vx0);
    vy_tmp = alfn_ratio.*(dpyF - drho.*vy0);
    vz_tmp = alfn_ratio.*(dpzF - drho.*vz0);
    bdotv = bx1.*vx_tmp + by1.*vy_tmp + bz1.*vz_tmp;
    bdotv = bdotv./(b1.^2+1e-10);
    rhovx = rhovx0 + perp_ratio.*( dpxF + dpxB + dv_alf.*vx0 + bx1.*bdotv); % if CA>>1, then rhovx = rhovx0 + dpxF + dpxB
    rhovy = rhovy0 + perp_ratio.*( dpyF + dpyB + dv_alf.*vy0 + by1.*bdotv); %                rhovy = rhovy0 + dpyF + dpyB            
    rhovz = rhovz0 + perp_ratio.*( dpzF + dpzB + dv_alf.*vz0 + bz1.*bdotv); %                rhovz = rhovz0 + dpzF + dpzB 
                                                                  
    % get velocities with magnetic stress - now bulk vx, vy, vz are solved,
    vx = rhovx./rho; % NOTE: rho here is already updated at t(N+1)
    vy = rhovy./rho;
    vz = rhovz./rho;
  
    % now add gravitational acceleration - FIXME: need flux method for
    % conservation constrains
    vx = vx + dt.*G0x;
    vy = vy + dt.*G0y;
    vz = vz + dt.*G0z;
   
    % apply boundary conditions to plasma and magnetic flux
    [rho,p,vx,vy,vz,bi,bj,bk] = Boundaries(rho,p,vx,vy,vz,bi,bj,bk,NO,...
                                           x_bound,y_bound,z_bound,...
                                           ic_act,jc_act,kc_act,...
                                           if_act,jf_act,kf_act,...
                                           ic_lb,ic_rb,jc_lb,jc_rb,kc_lb,kc_rb,...
                                           if_lb,if_rb,jf_lb,jf_rb,kf_lb,kf_rb);
                                          
    % compute bx, by, bz from bi, bj, bk (Eqn. 28)
    bx(I,J,K) = ( bi(I+1,J,K).*xi(I+1,J,K)-bi(I,J,K).*xi(I,J,K) + bj(I,J+1,K).*xj(I,J+1,K)-bj(I,J,K).*xj(I,J,K) + bk(I,J,K+1).*xk(I,J,K+1)-bk(I,J,K).*xk(I,J,K) )./volume(I,J,K);
    by(I,J,K) = ( bi(I+1,J,K).*yi(I+1,J,K)-bi(I,J,K).*yi(I,J,K) + bj(I,J+1,K).*yj(I,J+1,K)-bj(I,J,K).*yj(I,J,K) + bk(I,J,K+1).*yk(I,J,K+1)-bk(I,J,K).*yk(I,J,K) )./volume(I,J,K);
    bz(I,J,K) = ( bi(I+1,J,K).*zi(I+1,J,K)-bi(I,J,K).*zi(I,J,K) + bj(I,J+1,K).*zj(I,J+1,K)-bj(I,J,K).*zj(I,J,K) + bk(I,J,K+1).*zk(I,J,K+1)-bk(I,J,K).*zk(I,J,K) )./volume(I,J,K);
    
    % adding the background field to cell center B fields, bx_total,
    % by_total and bz_total are not directlly used in the MHD solver, they
    % are only used in the time step calculation (Alfven speed)
    bx_total(I,J,K) = bx(I,J,K) + BX0(xc,yc,zc);
    by_total(I,J,K) = by(I,J,K) + BY0(xc,yc,zc);
    bz_total(I,J,K) = bz(I,J,K) + BZ0(xc,yc,zc);
    
    RealDT = toc(Tstart);   % Real time for each step (diagnostics only)
    RealT = RealT + RealDT; % Total real time used    (diagnostics only)
    
    % diagnostic plots
    if(mod(n,20)==0)
        
        % Check divergence of B
        divB = bi(ic_act+1,jc_act,kc_act) - bi(ic_act,jc_act,kc_act) + ...
               bj(ic_act,jc_act+1,kc_act) - bj(ic_act,jc_act,kc_act) + ...
               bk(ic_act,jc_act,kc_act+1) - bk(ic_act,jc_act,kc_act);

        % Magnetic energy
        b2 = (bx_total.^2+by_total.^2+bz_total.^2);
        
        % Mach and Beta   
        Vsound = sqrt(gamma.*p./rho);               % sound speed
        Vfluid = sqrt(vx.^2+vy.^2+vz.^2);           % fluid speed
        P_beta = 2./gamma.*Vsound.^2./Va_eff.^2;    % plasma beta
        Mach = Vfluid./Vsound;                      % Mach number
        
        % plot density and pressure
        figure(1);
        subplot(1,2,1)
        pcolor(xc(ic_act,jc_act,kc_act(1)),yc(ic_act,jc_act,kc_act(1)),squeeze(rho(ic_act,jc_act,kc_act(1))));
        shading flat;xlabel('x'),ylabel('y');
        title('rho')
        
        subplot(1,2,2)
        Temper = p(ic_act,jc_act,kc_act(1))./rho(ic_act,jc_act,kc_act(1));
        pcolor(xc(ic_act,jc_act,kc_act(1)),yc(ic_act,jc_act,kc_act(1)),squeeze(Temper));
        shading flat;xlabel('x'),ylabel('y');
        title('P')
        
        % output information
        disp(['Loop ',num2str(n),', Sim Time = ',num2str(Time), ', Real Time = ',num2str(RealT),...
            ', Sim DT = ', num2str(dt), ', Real DT = ',num2str(RealDT)]);
        disp(['           ','Max(divB) = ',num2str(max(abs(divB(:)))),', Max(Mach) = ',num2str(max(Mach(:))),', Min(rho) = ',num2str(min(rho(:))),...
              ', Min(p) = ',num2str(min(p(:))),', Min(beta) = ',num2str(min(P_beta(:)))]);
    
        pause(0.05);
        
    end
    
    if(Time > T_stop)
        break;
    end

end
